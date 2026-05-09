#!/usr/bin/env python3
"""
Full pairwise anchor validation — BLAST+clasp every anchor match against
whole-genome databases and classify hits as:

  A  — best hit overlaps expected anchor region (confirmed)
  B1 — best hit is elsewhere but same query segment overlaps pipeline coords
  B2 — best hit is elsewhere and different query segment

Reports per-pair, per-direction statistics and a global summary.

Located in utils/util_code/. Default paths assume this location relative to
the project root (../../) and the working directory (../../template/).

Usage:
  python utils/util_code/validate_anchor_scores.py                      # defaults
  python utils/util_code/validate_anchor_scores.py --cores 16
  python utils/util_code/validate_anchor_scores.py --work-dir template/ --root-dir ./
  python utils/util_code/validate_anchor_scores.py --keep-temp          # keep BLAST/clasp files
"""

import pickle
import sys
import subprocess
import argparse
import shutil
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict
from statistics import median

# clasp columns — same as pipeline
CLASP_COLUMNS_C = ['7', '8', '9', '10', '12']
CLASP_COLUMNS_BIG_C = ['1', '2']

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ROOT_DIR = str(SCRIPT_DIR / '..' / '..')
DEFAULT_WORK_DIR = str(SCRIPT_DIR / '..' / '..' / 'template')


# Utility functions

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
            'N': 'N', 'n': 'n',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D',
            'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
            'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
            'd': 'h', 'h': 'd'}
    return ''.join(comp.get(c, 'N') for c in reversed(seq))


def coords_overlap(s1, e1, s2, e2):
    """Half-open interval overlap check: [s1, e1) and [s2, e2)."""
    return s1 < e2 and s2 < e1


def overlap_length(s1, e1, s2, e2):
    """Length of overlap between [s1, e1) and [s2, e2). 0 if none."""
    return max(0, min(e1, e2) - max(s1, s2))


def abs_pos_from_hit(hit, seqids, seqlen):
    """Convert clasp hit subject coords to absolute genome position [start, end)."""
    chromo = hit['subject_name']
    if chromo not in seqids:
        return None, None
    pos = seqids.index(chromo)
    offset = seqlen[pos - 1] if pos > 0 else 0
    # BLAST 1-based inclusive -> 0-based half-open
    abs_start = hit['subject_start'] + offset - 1
    abs_end = hit['subject_end'] + offset
    return abs_start, abs_end


def write_multi_fasta(filepath, seqs):
    """Write multi-sequence FASTA. seqs: list of (header, sequence)."""
    with open(filepath, 'w') as f:
        for header, seq in seqs:
            f.write(f'>{header}\n')
            for i in range(0, len(seq), 80):
                f.write(str(seq[i:i + 80]) + '\n')


def run_blast(db, query, outfile, word_size=4, evalue='1e-3'):
    """Run blastn with -strand plus. Returns True on success."""
    cmd = ['blastn', '-db', db, '-query', query,
           '-strand', 'plus', '-outfmt', '6',
           '-evalue', evalue, '-word_size', str(word_size),
           '-out', outfile]
    return subprocess.run(cmd, capture_output=True, text=True).returncode == 0


def run_clasp(infile, outfile, l_param=2, e_param=0.1):
    """Run clasp.x with -C 1 2 for multi-query chaining. Returns True on success."""
    cmd = (['clasp.x', '-m', '-i', infile, '-c'] + CLASP_COLUMNS_C +
           ['-C'] + CLASP_COLUMNS_BIG_C +
           ['-l', str(l_param), '-e', str(e_param), '-o', outfile])
    return subprocess.run(cmd, capture_output=True, text=True).returncode == 0


def parse_clasp_multi_query(filepath):
    """
    Parse clasp output from a multi-query run.

    Returns {anchor_idx: [hit_dicts]} grouped by query anchor index.
    Hit dicts contain: query_name, subject_name, query_start, query_end,
    subject_start, subject_end, score.
    """
    results = defaultdict(list)
    path = Path(filepath)
    if not path.exists() or path.stat().st_size == 0:
        return results
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 8:
                continue
            qname = parts[1]
            if not qname.startswith('anchor_'):
                continue
            try:
                aidx = int(qname.split('_', 1)[1])
            except ValueError:
                continue
            results[aidx].append({
                'query_name': qname,
                'subject_name': parts[2],
                'query_start': int(parts[3]),
                'query_end': int(parts[4]),
                'subject_start': int(parts[5]),
                'subject_end': int(parts[6]),
                'score': float(parts[7]),
            })
    return results


def ensure_blastdb(org, root_dir):
    """Check for whole-genome blastdb; create from FASTA if missing. Returns True if available."""
    db_prefix = Path(root_dir) / 'utils' / 'blastdbs' / org
    for ext in ('.nin', '.ndb', '.nsq'):
        if Path(str(db_prefix) + ext).exists():
            return True
    fasta = Path(root_dir) / 'utils' / 'genomes' / f'{org}.fasta'
    if not fasta.exists():
        return False
    Path(db_prefix).parent.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        ['makeblastdb', '-in', str(fasta), '-out', str(db_prefix), '-dbtype', 'nucl'],
        capture_output=True, text=True).returncode == 0


def load_orgs(work_dir):
    """Load organism list from work_dir/orgs."""
    with open(Path(work_dir) / 'orgs') as f:
        return [line.strip() for line in f if line.strip()]


# Phase 1: Preload data

def preload_data(orgs, root_dir):
    """
    Load aligned maps for all organisms, extract anchor sequences for
    every position involved in a pairwise match, store small metadata
    (seqids, seqlen), and free genome strings.

    Returns (anchor_seqs, anchor_info, small_metas, aligned_maps).
    """
    anchor_dir = Path(root_dir) / 'anchors'

    # Step 1: Load all aligned maps
    aligned_maps = {}
    for org in orgs:
        p = anchor_dir / 'aligned' / org
        if not p.exists():
            print(f"  SKIP {org}: no aligned map", file=sys.stderr)
            continue
        with open(p, 'rb') as f:
            aligned_maps[org] = pickle.load(f)

    # Step 2: Identify all genome positions needing sequence extraction.
    # This includes own anchors with matches AND target positions from other orgs'
    # matches (since we need the partner anchor's sequence for direction 2).
    needed = defaultdict(set)
    for org_src, amap in aligned_maps.items():
        for i, data in amap.items():
            if 'matches' not in data:
                continue
            for org_tgt, tgt_data in data['matches'].items():
                if 'matches' in tgt_data:
                    for j in tgt_data['matches']:
                        needed[org_src].add(i)
                        needed[org_tgt].add(j)
                if 'dups_matches' in tgt_data:
                    syn = tgt_data['dups_matches'].get('syntenic', set())
                    for j in tgt_data['dups_matches']:
                        if j != 'syntenic' and j in syn:
                            needed[org_src].add(i)
                            needed[org_tgt].add(j)

    # Step 3: Load genomes, extract sequences, free genome strings
    anchor_seqs = {}
    anchor_info = {}
    small_metas = {}
    ready = set()

    for org in orgs:
        if org not in needed:
            continue
        meta_path = Path(root_dir) / 'utils' / 'metadata_genomes' / org
        if not meta_path.exists():
            print(f"  SKIP {org}: no metadata_genomes", file=sys.stderr)
            continue
        if not ensure_blastdb(org, root_dir):
            print(f"  SKIP {org}: no blastdb and cannot create", file=sys.stderr)
            continue

        with open(meta_path, 'rb') as f:
            seqids, seqlen, genome_seq = pickle.load(f)
        small_metas[org] = (list(seqids), list(seqlen))

        anchor_seqs[org] = {}
        anchor_info[org] = {}
        amap = aligned_maps.get(org, {})

        for pos in needed[org]:
            if pos not in amap:
                continue
            ad = amap[pos]
            alen = ad['end'] - ad['start']
            if alen <= 0:
                continue
            seq = str(genome_seq[pos:pos + alen])
            if seq:
                anchor_seqs[org][pos] = seq
                anchor_info[org][pos] = {
                    'len': alen,
                    'chromo': ad.get('chromosome', '?'),
                }

        del genome_seq
        ready.add(org)
        print(f"  {org}: {len(anchor_seqs[org])} anchor seqs extracted", file=sys.stderr)

    aligned_maps = {k: v for k, v in aligned_maps.items() if k in ready}
    return anchor_seqs, anchor_info, small_metas, aligned_maps


# Phase 2: Generate pair tasks

def _extract_matches(amap, org_src, org_tgt):
    """
    Extract all match tuples from org_src's aligned_map targeting org_tgt.

    Returns [(i, j, score, coords_src, coords_tgt, is_reverse, is_dup), ...]
    where i is in org_src, j is in org_tgt.
    coords_src = pipeline hit coords relative to anchor i (0-based inclusive).
    coords_tgt = pipeline hit coords relative to anchor j (0-based inclusive).
    """
    out = []
    own_key = f'hit coordinates in (own) {org_src} candidate'
    tgt_key = f'hit coordinates in {org_tgt} candidate'

    for i, data in amap.items():
        if 'matches' not in data:
            continue
        if org_tgt not in data['matches']:
            continue
        tgt_data = data['matches'][org_tgt]

        # Regular matches
        if 'matches' in tgt_data:
            for j, md in tgt_data['matches'].items():
                out.append((i, j,
                            md.get('match score', 0),
                            md.get(own_key),
                            md.get(tgt_key),
                            md.get('match is on other strand in other genome', False),
                            False))
        # Syntenic dup matches
        if 'dups_matches' in tgt_data:
            syn = tgt_data['dups_matches'].get('syntenic', set())
            for j, md in tgt_data['dups_matches'].items():
                if j == 'syntenic':
                    continue
                if j not in syn:
                    continue
                out.append((i, j,
                            md.get('match score', 0),
                            md.get(own_key),
                            md.get(tgt_key),
                            md.get('match is on other strand in other genome', False),
                            True))
    return out


def generate_pair_tasks(aligned_maps, anchor_seqs, anchor_info, small_metas,
                        root_dir, tmp_dir):
    """
    Build one task dict per unique organism pair containing all their matches.

    For pair (org_a, org_b) with org_a < org_b alphabetically, matches are
    collected from both aligned maps and deduplicated by (pos_a, pos_b).
    Each match tuple is normalized to (pos_in_org_a, pos_in_org_b, ...).
    """
    all_orgs = sorted(aligned_maps.keys())
    tasks = []

    for idx_a, org_a in enumerate(all_orgs):
        for org_b in all_orgs[idx_a + 1:]:
            # Matches from org_a->org_b: (i_in_a, j_in_b, score, coords_a, coords_b, ...)
            matches_ab = _extract_matches(aligned_maps[org_a], org_a, org_b)

            # Matches from org_b->org_a: (i_in_b, j_in_a, score, coords_b, coords_a, ...)
            # Normalize to (pos_in_a, pos_in_b, score, coords_a, coords_b, ...)
            matches_ba_raw = _extract_matches(aligned_maps[org_b], org_b, org_a)
            matches_ba_norm = [(j, i, sc, ca, cb, rev, dup)
                               for i, j, sc, cb, ca, rev, dup in matches_ba_raw]

            # Deduplicate by (pos_a, pos_b) — same match seen from both sides
            seen = {}
            for m in matches_ab + matches_ba_norm:
                key = (m[0], m[1])
                if key not in seen:
                    seen[key] = m
            matches = list(seen.values())

            if not matches:
                continue

            # Filter to matches where we have both anchor sequences
            filtered = []
            for pos_a, pos_b, sc, ca, cb, rev, dup in matches:
                if (org_a in anchor_seqs and pos_a in anchor_seqs[org_a] and
                        org_b in anchor_seqs and pos_b in anchor_seqs[org_b]):
                    filtered.append((pos_a, pos_b, sc, ca, cb, rev, dup))

            if not filtered:
                continue

            # Collect only the sequences needed for this pair
            seqs_a = {}
            seqs_b = {}
            lens_a = {}
            lens_b = {}
            for pos_a, pos_b, *_ in filtered:
                seqs_a[pos_a] = anchor_seqs[org_a][pos_a]
                seqs_b[pos_b] = anchor_seqs[org_b][pos_b]
                lens_a[pos_a] = anchor_info[org_a][pos_a]['len']
                lens_b[pos_b] = anchor_info[org_b][pos_b]['len']

            tasks.append({
                'org_a': org_a,
                'org_b': org_b,
                'matches': filtered,
                'seqs_a': seqs_a,
                'seqs_b': seqs_b,
                'lens_a': lens_a,
                'lens_b': lens_b,
                'meta_a': small_metas[org_a],
                'meta_b': small_metas[org_b],
                'root_dir': root_dir,
                'tmp_dir': tmp_dir,
            })

    return tasks


# Phase 3: Process pair (runs in worker)

def categorize_hits(all_hits, target_pos, target_len, own_coords, seqids, seqlen):
    """
    Classify BLAST+clasp hits for one anchor in one direction.

    Checks whether the best overall hit overlaps the expected anchor region
    [target_pos, target_pos + target_len) in the target genome.

    Returns (category, best_A_score, best_B_score, overlap_len_if_B1)
    where category is 'A', 'B1', 'B2', or 'no_hit'.
    """
    if not all_hits:
        return 'no_hit', None, None, None

    target_end = target_pos + target_len

    # Parse pipeline own coords (0-based inclusive) -> half-open
    have_own = isinstance(own_coords, (list, tuple)) and len(own_coords) >= 2
    if have_own:
        own_s = int(own_coords[0])
        own_e = int(own_coords[1]) + 1  # inclusive -> exclusive

    cat_A_hits = []
    cat_B_hits = []

    for hit in all_hits:
        abs_s, abs_e = abs_pos_from_hit(hit, seqids, seqlen)
        if abs_s is None:
            continue
        if coords_overlap(abs_s, abs_e, target_pos, target_end):
            cat_A_hits.append(hit)
        else:
            cat_B_hits.append(hit)

    best_A = max((h['score'] for h in cat_A_hits), default=None)
    best_B = max((h['score'] for h in cat_B_hits), default=None)

    # Direction confirmed if best overall hit is in Category A
    if best_A is not None and (best_B is None or best_A >= best_B):
        return 'A', best_A, best_B, None

    # Best hit is in B — classify as B1 or B2
    if best_B is not None:
        best_B_hit = max(cat_B_hits, key=lambda h: h['score'])
        if have_own:
            # BLAST query coords: 1-based inclusive -> 0-based half-open
            q_s = best_B_hit['query_start'] - 1
            q_e = best_B_hit['query_end']
            if coords_overlap(q_s, q_e, own_s, own_e):
                ol = overlap_length(q_s, q_e, own_s, own_e)
                return 'B1', best_A, best_B, ol
            else:
                return 'B2', best_A, best_B, None
        else:
            return 'B2', best_A, best_B, None

    return 'no_hit', None, None, None


def process_pair(task):
    """
    Full validation pipeline for one organism pair.

    1. Write 4 multi-FASTA files (fwd/rev × org_a/org_b anchors)
    2. Run 4 BLAST+clasp (2 directions × 2 orientations)
    3. Parse clasp output grouped by query anchor
    4. Categorize each match in both directions as A/B1/B2
    """
    org_a = task['org_a']
    org_b = task['org_b']
    matches = task['matches']
    seqs_a = task['seqs_a']
    seqs_b = task['seqs_b']
    lens_a = task['lens_a']
    lens_b = task['lens_b']
    meta_a = task['meta_a']  # (seqids, seqlen)
    meta_b = task['meta_b']
    root_dir = task['root_dir']
    tmp_dir = task['tmp_dir']

    pair_dir = Path(tmp_dir) / f'{org_a}___{org_b}'
    pair_dir.mkdir(parents=True, exist_ok=True)

    db_a = f'{root_dir}/utils/blastdbs/{org_a}'
    db_b = f'{root_dir}/utils/blastdbs/{org_b}'

    # Step 1: Write multi-FASTA files with deduplicated anchor indices
    indices_a = sorted(set(m[0] for m in matches))
    indices_b = sorted(set(m[1] for m in matches))

    fwd_a_fa = str(pair_dir / 'fwd_a.fa')
    rev_a_fa = str(pair_dir / 'rev_a.fa')
    fwd_b_fa = str(pair_dir / 'fwd_b.fa')
    rev_b_fa = str(pair_dir / 'rev_b.fa')

    write_multi_fasta(fwd_a_fa,
                      [(f'anchor_{i}', seqs_a[i]) for i in indices_a])
    write_multi_fasta(rev_a_fa,
                      [(f'anchor_{i}', reverse_complement(seqs_a[i])) for i in indices_a])
    write_multi_fasta(fwd_b_fa,
                      [(f'anchor_{j}', seqs_b[j]) for j in indices_b])
    write_multi_fasta(rev_b_fa,
                      [(f'anchor_{j}', reverse_complement(seqs_b[j])) for j in indices_b])

    # Step 2: Run 4 BLAST+clasp (2 directions × 2 orientations)
    runs = [
        (fwd_a_fa, db_b, 'dir1_fwd'),  # org_a anchors fwd -> org_b genome
        (rev_a_fa, db_b, 'dir1_rev'),  # org_a anchors rev -> org_b genome
        (fwd_b_fa, db_a, 'dir2_fwd'),  # org_b anchors fwd -> org_a genome
        (rev_b_fa, db_a, 'dir2_rev'),  # org_b anchors rev -> org_a genome
    ]

    parsed = {}
    for query_fa, db, tag in runs:
        blast_out = str(pair_dir / f'{tag}_blast')
        clasp_out = str(pair_dir / f'{tag}_clasp')
        if run_blast(db, query_fa, blast_out) and run_clasp(blast_out, clasp_out):
            parsed[tag] = parse_clasp_multi_query(clasp_out)
        else:
            parsed[tag] = defaultdict(list)

    # Step 3: Combine fwd+rev hits per anchor per direction
    def combined_hits(tag_fwd, tag_rev, idx):
        return parsed[tag_fwd].get(idx, []) + parsed[tag_rev].get(idx, [])

    # Step 4: Categorize each match in both directions
    dir1_results = {'A': 0, 'B1': 0, 'B2': 0, 'no_hit': 0, 'B1_overlaps': []}
    dir2_results = {'A': 0, 'B1': 0, 'B2': 0, 'no_hit': 0, 'B1_overlaps': []}
    both_confirmed = 0
    details = []

    for pos_a, pos_b, score, coords_a, coords_b, is_reverse, is_dup in matches:
        # Direction 1: anchor at pos_a -> org_b genome
        # Expected landing zone: [pos_b, pos_b + len_b)
        # Pipeline own coords for this direction: coords_a (in anchor at pos_a)
        hits_d1 = combined_hits('dir1_fwd', 'dir1_rev', pos_a)
        cat1, bestA1, bestB1, ol1 = categorize_hits(
            hits_d1, pos_b, lens_b[pos_b], coords_a,
            meta_b[0], meta_b[1])
        dir1_results[cat1] += 1
        if cat1 == 'B1' and ol1 is not None:
            dir1_results['B1_overlaps'].append(ol1)

        # Direction 2: anchor at pos_b -> org_a genome
        # Expected landing zone: [pos_a, pos_a + len_a)
        # Pipeline own coords for this direction: coords_b (in anchor at pos_b)
        hits_d2 = combined_hits('dir2_fwd', 'dir2_rev', pos_b)
        cat2, bestA2, bestB2, ol2 = categorize_hits(
            hits_d2, pos_a, lens_a[pos_a], coords_b,
            meta_a[0], meta_a[1])
        dir2_results[cat2] += 1
        if cat2 == 'B1' and ol2 is not None:
            dir2_results['B1_overlaps'].append(ol2)

        if cat1 == 'A' and cat2 == 'A':
            both_confirmed += 1

        # Record details for non-confirmed matches
        if cat1 != 'A' or cat2 != 'A':
            details.append({
                'pos_a': pos_a, 'pos_b': pos_b,
                'cat1': cat1, 'bestA1': bestA1, 'bestB1': bestB1, 'ol1': ol1,
                'cat2': cat2, 'bestA2': bestA2, 'bestB2': bestB2, 'ol2': ol2,
                'pipeline_score': score, 'is_dup': is_dup,
            })

    return {
        'org_a': org_a, 'org_b': org_b,
        'n_matches': len(matches),
        'dir1': dir1_results,
        'dir2': dir2_results,
        'both_confirmed': both_confirmed,
        'details': details,
    }


# Phase 4: Reporting

def _fmt_score(s):
    return f'{s:.1f}' if s is not None else '-'


def _fmt_dir_detail(cat, best_other, ol):
    if cat == 'A':
        return 'A'
    if cat == 'B1':
        return f'B1(best_other={_fmt_score(best_other)}, overlap={ol}bp)'
    if cat == 'B2':
        return f'B2(best_other={_fmt_score(best_other)})'
    return 'no_hit'


def print_pair_report(r):
    """Print per-pair per-direction breakdown."""
    n = r['n_matches']
    d1 = r['dir1']
    d2 = r['dir2']

    print(f"PAIR {r['org_a']} vs {r['org_b']}: {n} matches")

    # Direction 1
    print(f"  Dir1 ({r['org_a']} anchors \u2192 {r['org_b']} genome):")
    print(f"    best_in_A: {d1['A']}/{n}    not_in_A: {n - d1['A']}/{n}")
    if n - d1['A'] > 0:
        ols = d1['B1_overlaps']
        ol_str = (f" [overlap_lengths: {','.join(str(x) for x in ols)}bp]"
                  if ols else '')
        print(f"    of non-A: B1 (overlaps query coords): {d1['B1']}{ol_str}")
        print(f"              B2 (different query region): {d1['B2']}")
        if d1['no_hit'] > 0:
            print(f"              no_hit: {d1['no_hit']}")

    # Direction 2
    print(f"  Dir2 ({r['org_b']} anchors \u2192 {r['org_a']} genome):")
    print(f"    best_in_A: {d2['A']}/{n}    not_in_A: {n - d2['A']}/{n}")
    if n - d2['A'] > 0:
        ols = d2['B1_overlaps']
        ol_str = (f" [overlap_lengths: {','.join(str(x) for x in ols)}bp]"
                  if ols else '')
        print(f"    of non-A: B1 (overlaps query coords): {d2['B1']}{ol_str}")
        print(f"              B2 (different query region): {d2['B2']}")
        if d2['no_hit'] > 0:
            print(f"              no_hit: {d2['no_hit']}")

    print(f"  Both directions confirmed: {r['both_confirmed']}/{n}")

    # Non-confirmed match details
    if r['details']:
        print(f"  DETAILS (non-confirmed matches):")
        for d in r['details'][:20]:
            d1_str = _fmt_dir_detail(d['cat1'], d['bestB1'], d['ol1'])
            d2_str = _fmt_dir_detail(d['cat2'], d['bestB2'], d['ol2'])
            dup_tag = ' [dup]' if d.get('is_dup') else ''
            print(f"    i={d['pos_a']} j={d['pos_b']}: "
                  f"dir1={d1_str} dir2={d2_str} "
                  f"pipeline={_fmt_score(d['pipeline_score'])}{dup_tag}")
        if len(r['details']) > 20:
            print(f"    ... and {len(r['details']) - 20} more")
    print()


def print_global_summary(results):
    """Aggregate and print global statistics across all pairs."""
    total = sum(r['n_matches'] for r in results)
    both = sum(r['both_confirmed'] for r in results)
    # Inclusion-exclusion: d1_only = dir1_A - both, d2_only = dir2_A - both
    d1_only = sum(r['dir1']['A'] - r['both_confirmed'] for r in results)
    d2_only = sum(r['dir2']['A'] - r['both_confirmed'] for r in results)
    neither = total - both - d1_only - d2_only

    # Aggregate non-confirmed direction breakdown
    all_b1 = sum(r['dir1']['B1'] + r['dir2']['B1'] for r in results)
    all_b2 = sum(r['dir1']['B2'] + r['dir2']['B2'] for r in results)
    all_no = sum(r['dir1']['no_hit'] + r['dir2']['no_hit'] for r in results)
    all_b1_ols = []
    for r in results:
        all_b1_ols.extend(r['dir1']['B1_overlaps'])
        all_b1_ols.extend(r['dir2']['B1_overlaps'])

    print("=" * 70)
    print("=== GLOBAL SUMMARY ===")
    print("=" * 70)
    print(f"Total matches checked:     {total}")
    if total > 0:
        print(f"Both directions confirmed: {both} ({100 * both / total:.1f}%)")
        print(f"Dir1 only:                 {d1_only}")
        print(f"Dir2 only:                 {d2_only}")
        print(f"Neither:                   {neither}")
        print()
        print("Non-confirmed directions breakdown:")
        if all_b1_ols:
            print(f"  B1 (same query region, higher score elsewhere):  {all_b1}"
                  f"  overlap_lengths: min={min(all_b1_ols)}"
                  f" median={median(all_b1_ols):.0f}"
                  f" max={max(all_b1_ols)}")
        else:
            print(f"  B1 (same query region, higher score elsewhere):  {all_b1}")
        print(f"  B2 (different query region):                     {all_b2}")
        print(f"  no_hit:                                          {all_no}")
    print("=" * 70)


# Main

def main():
    parser = argparse.ArgumentParser(
        description='Full pairwise anchor validation with A/B1/B2 categorization')
    parser.add_argument('--work-dir', default=DEFAULT_WORK_DIR,
                        help=f'Working directory with orgs file (default: {DEFAULT_WORK_DIR})')
    parser.add_argument('--root-dir', default=DEFAULT_ROOT_DIR,
                        help=f'Project root with anchors/, utils/ (default: {DEFAULT_ROOT_DIR})')
    parser.add_argument('--cores', type=int, default=32,
                        help='Parallel workers (default: 32)')
    parser.add_argument('--keep-temp', action='store_true',
                        help='Keep BLAST/clasp temp files')
    args = parser.parse_args()

    work_dir = str(Path(args.work_dir).resolve())
    root_dir = str(Path(args.root_dir).resolve())

    tmp_dir = str(Path(root_dir) / 'to_del_tmp' / 'validate')
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    orgs = load_orgs(work_dir)

    print(f"=== Full Pairwise Anchor Validation ===")
    print(f"work_dir:  {work_dir}")
    print(f"root_dir:  {root_dir}")
    print(f"cores:     {args.cores}")
    print(f"organisms: {len(orgs)}")
    print()

    # Phase 1: Preload
    print("Phase 1: Preloading data...", file=sys.stderr)
    anchor_seqs, anchor_info, small_metas, aligned_maps = preload_data(orgs, root_dir)
    n_seqs = sum(len(v) for v in anchor_seqs.values())
    print(f"  {len(aligned_maps)} organisms ready, "
          f"{n_seqs} total anchor seqs", file=sys.stderr)

    # Phase 2: Generate pair tasks
    print("Phase 2: Generating pair tasks...", file=sys.stderr)
    tasks = generate_pair_tasks(aligned_maps, anchor_seqs, anchor_info,
                                small_metas, root_dir, tmp_dir)
    total_matches = sum(len(t['matches']) for t in tasks)
    print(f"  {len(tasks)} pairs, {total_matches} total matches", file=sys.stderr)

    if not tasks:
        print("No pairs with matches found. Nothing to validate.")
        return

    # Free preloaded data no longer needed by main process
    del aligned_maps, anchor_seqs, anchor_info

    # Phase 3: Parallel execution
    print(f"Phase 3: Processing {len(tasks)} pairs on {args.cores} cores...",
          file=sys.stderr)

    if args.cores == 1:
        results = [process_pair(t) for t in tasks]
    else:
        with Pool(args.cores) as pool:
            results = pool.map(process_pair, tasks)

    # Phase 4: Reporting
    print()
    print("=" * 70)
    print("=== PER-PAIR RESULTS ===")
    print("=" * 70)
    print()

    for r in sorted(results, key=lambda x: (x['org_a'], x['org_b'])):
        print_pair_report(r)

    print_global_summary(results)

    # Cleanup
    if not args.keep_temp:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        print("Temp files cleaned up.", file=sys.stderr)
    else:
        print(f"Temp files kept at: {tmp_dir}", file=sys.stderr)

    print("\n=== Validation complete ===")


if __name__ == '__main__':
    main()
