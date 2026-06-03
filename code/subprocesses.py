#!/usr/bin/env python3
"""
Wrapper functions for external bioinformatics tools.
Reads default parameters from pipeline_config.yaml

This module provides:
- blast(): Run BLAST with configurable parameters
- clasp(): Run clasp with configurable parameters
- get_config(): Get the loaded configuration dictionary
- get_score_threshold(): Get score threshold for filtering
- get_tolerance(): Get tolerance for multiple matches
- get_min_anchor_size(): Get minimum anchor/candidate size
- get_overlap_threshold(): Get overlap threshold for dups integration
- get_gap_size(): Get gap size for dups shrinking
- get_split_size(): Get FASTA splitter part sizes
"""

import resource
import platform
from subprocess import run, TimeoutExpired
from os.path import isfile
from os import getcwd
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import yaml
import pathlib

# Load pipeline configuration
def load_pipeline_config():
    """Load pipeline configuration from YAML file."""
    try:
        # Try to find config relative to this file
        code_dir = pathlib.Path(__file__).parent.resolve()
        root_dir = code_dir.parent
        config_file = root_dir / 'template' / 'pipeline_config.yaml'

        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        # Return defaults if config not found (maintains backwards compatibility)
        return {
            'blast': {'word_size': 11, 'evalue': 1e-3, 'strand': 'plus', 'outfmt': 6},
            'clasp': {'l': 0.5, 'e': 0, 'columns_c': [7, 8, 9, 10, 12], 'columns_C': [1, 2]},
            'scoring': {'score_threshold': 40},
            'matching': {'tolerance': 10000},
            'splitting': {'self_blast_part_size': 100000, 'pairwise_part_size': 10000000},
            'regions': {'min_anchor_size': 300},
            'dups_integration': {'overlap_threshold': 0.3, 'gap_size': 100},
        }

# Load config once at module import
_CONFIG = load_pipeline_config()

def get_config():
    """Return the loaded configuration dictionary."""
    return _CONFIG

def get_score_threshold():
    """Get score threshold for BLAST/clasp filtering (default: 40)."""
    return _CONFIG.get('scoring', {}).get('score_threshold', 40)

def get_tolerance():
    """Get tolerance for multiple match distance in bp (default: 10000)."""
    return _CONFIG.get('matching', {}).get('tolerance', 10000)

def get_min_anchor_size():
    """Get minimum anchor/candidate size in bp (default: 300)."""
    return _CONFIG.get('regions', {}).get('min_anchor_size', 300)

def get_overlap_threshold():
    """Get overlap threshold for dups integration (default: 0.3 = 30%)."""
    return _CONFIG.get('dups_integration', {}).get('overlap_threshold', 0.3)

def get_gap_size():
    """Get gap size in bp for dups shrinking (default: 100)."""
    return _CONFIG.get('dups_integration', {}).get('gap_size', 100)

def get_split_size(split_type='self_blast'):
    """Get FASTA splitter part size.

    Args:
        split_type: 'self_blast' (default: 100000) or 'pairwise' (default: 10000000)
    """
    splitting = _CONFIG.get('splitting', {})
    if split_type == 'pairwise':
        return splitting.get('pairwise_part_size', 10000000)
    return splitting.get('self_blast_part_size', 100000)

def get_mem_limit_bytes():
    mem = _CONFIG.get('memory', {})
    max_mem_mb = mem.get('max_mem_mb', None)
    if max_mem_mb is None:
        return None
    cores = max(1, mem.get('cores', 1))
    return int((max_mem_mb * 1024 * 1024) / cores)



def get_timeout_seconds(mode='self'):
    key = 'timeout_bcamm_minutes' if mode == 'pairwise' else 'timeout_self_minutes'
    mins = _CONFIG.get('memory', {}).get(key)
    if mins is None or mins <= 0:
        return None
    return mins * 60

def _run_with_limit(cmd, mem_limit_bytes, timeout_sec=None):
    kwargs = {}
    if timeout_sec is not None:
        kwargs['timeout'] = timeout_sec
    if mem_limit_bytes and platform.system() == 'Linux':
        limit = int(mem_limit_bytes)
        def _set_limit():
            resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
        kwargs['preexec_fn'] = _set_limit
    try:
        return run(cmd, **kwargs).returncode
    except TimeoutExpired:
        print(f'[TIMEOUT] killed after {timeout_sec}s: {" ".join(str(c) for c in cmd[:3])}...')
        return -99

def split_fasta_for_retry(fasta_path, sub_chunk_size, retry_dir):
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    n = len(records)
    if n <= 1:
        return None
    sub_chunk_size = max(1, sub_chunk_size)
    n_chunks = (n + sub_chunk_size - 1) // sub_chunk_size
    width = len(str(n_chunks))
    base = os.path.basename(fasta_path)
    os.makedirs(retry_dir, exist_ok=True)
    chunks = []
    for i in range(0, n, sub_chunk_size):
        chunk_records = records[i:i + sub_chunk_size]
        chunk_path = os.path.join(retry_dir, f"{base}.sub-{len(chunks)+1:0{width}d}")
        SeqIO.write(chunk_records, chunk_path, 'fasta')
        chunks.append(chunk_path)
    return chunks

def split_single_sequence(fwd_fasta_path, rev_fasta_path, retry_dir):
    part_size = _CONFIG.get('memory', {}).get('single_seq_part_size', 300)
    gap = 100

    fwd_records = list(SeqIO.parse(fwd_fasta_path, 'fasta'))
    if len(fwd_records) != 1:
        return None
    rec = fwd_records[0]
    seq_str = str(rec.seq)
    seq_len = len(seq_str)
    if seq_len <= part_size:
        return None

    hdr_parts = rec.id.split('$')
    abs_start = int(hdr_parts[1])
    abs_end = int(hdr_parts[-1])
    if abs_end != abs_start + seq_len:
        raise RuntimeError(
            f'Header abs_end ({abs_end}) != abs_start + seq_len ({abs_start} + {seq_len} = {abs_start + seq_len})'
            f' in {fwd_fasta_path} — header/sequence length mismatch')
    chromo = hdr_parts[2]
    chromo_rel_start = int(hdr_parts[3])

    fwd_new = []
    rev_new = []
    sub_windows = []
    pos = 0
    while pos < seq_len:
        end = pos + part_size
        if end > seq_len:
            end = seq_len
        if end < seq_len and (seq_len - end) < (part_size + gap):
            end = seq_len
        sub_seq = seq_str[pos:end]
        new_abs_start = abs_start + pos
        new_abs_end = abs_start + end
        new_rel = chromo_rel_start + pos
        hdr = f"kmer${new_abs_start}${chromo}${new_rel}${new_abs_end}"
        fwd_new.append(SeqRecord(Seq(sub_seq), id=hdr, description=''))
        rev_new.append(SeqRecord(Seq(sub_seq).reverse_complement(), id=hdr, description=''))
        sub_windows.append((new_abs_start, new_abs_end))
        if end >= seq_len:
            break
        pos = end + gap

    os.makedirs(retry_dir, exist_ok=True)
    fwd_base = os.path.basename(fwd_fasta_path)
    rev_base = os.path.basename(rev_fasta_path)
    fwd_split = os.path.join(retry_dir, 'fwd_' + fwd_base)
    rev_split = os.path.join(retry_dir, 'rev_' + rev_base)
    SeqIO.write(fwd_new, fwd_split, 'fasta')
    SeqIO.write(rev_new, rev_split, 'fasta')
    # marker: parent -> sub-windows (read by parse_clasp_out to fix orig_idx)
    marker = os.path.join(retry_dir, 'split_windows.tsv')
    with open(marker, 'a') as f:
        subs = '\t'.join(f'{s},{e}' for s, e in sub_windows)
        f.write(f'{abs_start},{abs_end}\t{subs}\n')
    print(f'[OOM-SPLIT] Split single seq ({seq_len}bp) into {len(fwd_new)} parts of ~{part_size}bp, gap={gap}bp')
    return fwd_split, rev_split

def split_sequences_to_windows(fasta_path, is_reverse=False):
    part_size = _CONFIG.get('memory', {}).get('single_seq_part_size', 300)
    gap = 100
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    if not records:
        return None

    new_records = []
    any_split = False
    split_entries = []
    for rec in records:
        seq_str = str(rec.seq)
        seq_len = len(seq_str)
        if seq_len <= part_size:
            new_records.append(rec)
            continue
        any_split = True
        # if reverse: RC back to forward for correct header computation, then RC each window back
        if is_reverse:
            seq_str = str(Seq(seq_str).reverse_complement())
        hdr_parts = rec.id.split('$')
        abs_start = int(hdr_parts[1])
        abs_end = int(hdr_parts[-1])
        if abs_end != abs_start + seq_len:
            raise RuntimeError(
                f'Header abs_end ({abs_end}) != abs_start + seq_len ({abs_start} + {seq_len} = {abs_start + seq_len})'
                f' in {fasta_path} — header/sequence length mismatch')
        chromo = hdr_parts[2]
        chromo_rel_start = int(hdr_parts[3])
        sub_windows = []
        pos = 0
        while pos < seq_len:
            end = pos + part_size
            if end > seq_len:
                end = seq_len
            if end < seq_len and (seq_len - end) < (part_size + gap):
                end = seq_len
            sub_seq = seq_str[pos:end]
            new_abs_start = abs_start + pos
            new_abs_end = abs_start + end
            new_rel = chromo_rel_start + pos
            hdr = f"kmer${new_abs_start}${chromo}${new_rel}${new_abs_end}"
            if is_reverse:
                new_records.append(SeqRecord(Seq(sub_seq).reverse_complement(), id=hdr, description=''))
            else:
                new_records.append(SeqRecord(Seq(sub_seq), id=hdr, description=''))
            sub_windows.append((new_abs_start, new_abs_end))
            if end >= seq_len:
                break
            pos = end + gap
        split_entries.append((abs_start, abs_end, sub_windows))

    if not any_split:
        return None
    SeqIO.write(new_records, fasta_path, 'fasta')
    # marker: parent -> sub-windows (read by parse_clasp_out to fix orig_idx)
    marker = fasta_path + '.split_windows.tsv'
    with open(marker, 'w') as f:
        for p_start, p_end, subs in split_entries:
            sub_str = '\t'.join(f'{s},{e}' for s, e in subs)
            f.write(f'{p_start},{p_end}\t{sub_str}\n')
    print(f'[OOM-SPLIT] Split {len(records)} sequences into {len(new_records)} windows of ~{part_size}bp, gap={gap}bp')
    return len(new_records)

def _log_split_failure(blast_outfile, message):
    # blast_outfile = {work_dir}/blast_out_XXX/{org_or_tuple}/{file} -> 3 dirs up = work_dir
    work_dir = os.path.dirname(os.path.dirname(os.path.dirname(blast_outfile)))
    split_log_dir = os.path.join(work_dir, 'SPLIT_LOG')
    os.makedirs(split_log_dir, exist_ok=True)
    parts = os.path.normpath(blast_outfile).split(os.sep)
    job_name = '__'.join(parts[-3:])
    log_file = os.path.join(split_log_dir, f'{job_name}.log')
    with open(log_file, 'a') as f:
        f.write(message + '\n')

def _write_clasp_fail(clasp_outfile, msg):
    # write # comment line so awk ($1 != "#") skips it; preserves partial output if any
    os.makedirs(os.path.dirname(clasp_outfile), exist_ok=True)
    with open(clasp_outfile, 'a') as f:
        f.write(f'# {msg}\n')

def run_blast_clasp_with_retry(db, query_fasta, blast_outfile, clasp_outfile,
                                word_size, mem_limit_bytes, clasp_mode='self',
                                rev_fasta_path=None, orientation='forward'):
    timeout_sec = get_timeout_seconds('pairwise' if clasp_mode == 'pairwise' else 'self')
    is_rev = (orientation == 'reverse')

    # no limits and no timeout -> original simple behavior (with binary error check)
    if not mem_limit_bytes and timeout_sec is None:
        rc = blast(db, query_fasta, blast_outfile, word_size)
        if rc != 0:
            raise RuntimeError(f'blast failed (rc={rc}) on {blast_outfile} — check blastn binary and target database')
        rc = clasp(blast_outfile, clasp_outfile, mode=clasp_mode)
        if rc != 0:
            raise RuntimeError(f'clasp failed (rc={rc}) on {clasp_outfile} — check clasp.x binary and conda environment')
        return 0

    # first attempt with limit/timeout
    rc = blast(db, query_fasta, blast_outfile, word_size, mem_limit_bytes, timeout_sec)
    blast_ok = (rc == 0)
    if blast_ok:
        rc = clasp(blast_outfile, clasp_outfile, mode=clasp_mode,
                   mem_limit_bytes=mem_limit_bytes, timeout_sec=timeout_sec)
    if rc == 0:
        return 0

    # rc > 0 = application error from the tool itself.
    # blast rc=4 (official "out of memory") AND rc=3 (generic "Error in BLAST engine")
    # can both be RLIMIT_AS / memory induced: under a tight RLIMIT_AS, blastn's engine
    # catches the allocation failure internally and exits rc=3 ("Unknown error code -1")
    # instead of the clean rc=4 path. So route blast rc=3 through the splitter too — if it
    # keeps failing it ends in a len*2 placeholder (a single window) or, for a systematic
    # cause (corrupt DB mid-search, disk I/O), a flood of placeholders caught downstream by
    # the ">=5 error scores = systematic failure" heuristic. Genuinely-fatal blast errors
    # have their own codes and STILL raise: rc=1 (query/options), rc=2 (bad DB), rc=6
    # (output). A missing binary is FileNotFoundError upstream; a wrong-arch binary is SIGILL
    # (rc<0, handled in the retry path below). Note: this only runs in the limit/timeout path
    # — a true no-limits run keeps rc=3 fatal via the simple path above.
    if rc > 0:
        which = 'clasp' if blast_ok else 'blast'
        if (which == 'clasp' and rc != 255) or (which == 'blast' and rc not in (3, 4)):
            raise RuntimeError(
                f'{which} application error (rc={rc}) on {blast_outfile}'
                f' — check binary/conda environment (not OOM/timeout, splitting cannot help)')

    # retryable failure: timeout (rc=-99), known OOM codes (blast rc=4, clasp rc=255),
    # or signal kill (rc<0) — could be system OOM, RLIMIT_AS, or other signal
    which = 'clasp' if blast_ok else 'blast'
    if rc == -99:
        label = 'TIMEOUT'
    elif rc > 0:
        label = f'{which}-rc={rc}'
    else:
        label = f'signal={-rc}'
    print(f'[RETRY {label}] {which} failed rc={rc} on {blast_outfile}')

    # --- splitting retry path (only OOM/timeout reach here, never application errors) ---
    # all retry files go under blast_outfile dir (per-job, no race conditions)
    blast_base = os.path.basename(blast_outfile)
    retry_dir = os.path.join(os.path.dirname(blast_outfile), 'retry', blast_base)
    os.makedirs(retry_dir, exist_ok=True)

    n_seqs = sum(1 for _ in SeqIO.parse(query_fasta, 'fasta'))

    # single sequence
    if n_seqs == 1:
        # pairwise: window split would create headers not in candidates dict
        if clasp_mode == 'pairwise':
            msg = f'[{label}-FAIL] Single seq pairwise, cannot split: {query_fasta}'
            print(msg)
            _log_split_failure(blast_outfile, msg)
            _write_clasp_fail(clasp_outfile, msg)
            return 0
        # split_single_sequence must always read forward content first
        if is_rev and rev_fasta_path:
            result = split_single_sequence(rev_fasta_path, query_fasta, retry_dir)
        else:
            result = split_single_sequence(query_fasta, rev_fasta_path, retry_dir)
        if result is None:
            msg = f'[{label}-FAIL] Single seq too short to split: {query_fasta}'
            print(msg)
            _log_split_failure(blast_outfile, msg)
            _write_clasp_fail(clasp_outfile, msg)
            return 0
        fwd_split, rev_split = result
        query_split = rev_split if is_rev else fwd_split
        rc = blast(db, query_split, blast_outfile, word_size, mem_limit_bytes, timeout_sec)
        if rc == 0:
            rc = clasp(blast_outfile, clasp_outfile, mode=clasp_mode,
                       mem_limit_bytes=mem_limit_bytes, timeout_sec=timeout_sec)
        if rc != 0:
            msg = f'[{label}-FAIL] Still fails (rc={rc}) after single-seq split: {query_fasta}'
            print(msg)
            _log_split_failure(blast_outfile, msg)
            _write_clasp_fail(clasp_outfile, msg)
        return 0

    # multi-sequence -> Round 1: sub-chunks
    if clasp_mode == 'pairwise':
        sub_chunk_size = 1
    else:
        sub_chunk_size = max(1, n_seqs // _CONFIG.get('memory', {}).get('oom_resplit_divisor', 100))
    sub_chunks = split_fasta_for_retry(query_fasta, sub_chunk_size, retry_dir)
    if sub_chunks is None:
        raise RuntimeError(
            f'split_fasta_for_retry returned None for {query_fasta} with {n_seqs} sequences'
            f' — file may be corrupt or unreadable')

    print(f'[{label}-RETRY] Split {query_fasta} ({n_seqs} seqs) into {len(sub_chunks)} sub-chunks')

    clasp_sub_outputs = []
    width = len(str(len(sub_chunks)))
    for idx, sub_fasta in enumerate(sub_chunks):
        sub_blast = os.path.join(retry_dir, f'blast.sub-{idx+1:0{width}d}')
        sub_clasp = os.path.join(retry_dir, f'clasp.sub-{idx+1:0{width}d}')

        rc = blast(db, sub_fasta, sub_blast, word_size, mem_limit_bytes, timeout_sec)
        if rc == 0:
            rc = clasp(sub_blast, sub_clasp, mode=clasp_mode,
                       mem_limit_bytes=mem_limit_bytes, timeout_sec=timeout_sec)

        if rc != 0:
            if clasp_mode == 'pairwise':
                # pairwise: 300bp window splitting creates headers not in candidates dict
                msg = f'[{label}-FAIL] Sub-chunk {idx+1}/{len(sub_chunks)} of {query_fasta} failed (rc={rc}), pairwise cannot window-split'
                print(msg)
                _log_split_failure(blast_outfile, msg)
                with open(sub_clasp, 'w') as f:
                    f.write(f'# {msg}\n')
            else:
                # Round 2: split sub-chunk sequences into 300bp windows
                print(f'[{label}-RETRY] Sub-chunk {idx+1}/{len(sub_chunks)} failed (rc={rc}), splitting to windows...')
                n_windows = split_sequences_to_windows(sub_fasta, is_reverse=is_rev)
                if n_windows is not None:
                    print(f'[{label}-RETRY] Sub-chunk {idx+1} -> {n_windows} windows')
                    rc = blast(db, sub_fasta, sub_blast, word_size, mem_limit_bytes, timeout_sec)
                    if rc == 0:
                        rc = clasp(sub_blast, sub_clasp, mode=clasp_mode,
                                   mem_limit_bytes=mem_limit_bytes, timeout_sec=timeout_sec)

                if rc != 0:
                    msg = f'[{label}-FAIL] Sub-chunk {idx+1}/{len(sub_chunks)} of {query_fasta} failed after window split (rc={rc})'
                    print(msg)
                    _log_split_failure(blast_outfile, msg)
                    with open(sub_clasp, 'w') as f:
                        f.write(f'# {msg}\n')

        clasp_sub_outputs.append(sub_clasp)

    # concatenate all sub-clasp outputs (# lines filtered by downstream awk)
    with open(clasp_outfile, 'w') as outf:
        for sub_out in clasp_sub_outputs:
            if os.path.exists(sub_out):
                with open(sub_out) as inf:
                    outf.write(inf.read())

    return 0

def _read_round2_markers(retry_dir):
    markers = {}
    if not os.path.isdir(retry_dir):
        return markers
    for marker_file in glob.glob(os.path.join(retry_dir, '*split_windows*')):
        with open(marker_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                p = parts[0].split(',')
                p_start, p_end = int(p[0]), int(p[1])
                subs = []
                for s in parts[1:]:
                    ss, se = s.split(',')
                    subs.append((int(ss), int(se)))
                markers[p_start] = (p_end, subs)
    return markers

def _filter_clasp_entries(clasp_file, abs_start, abs_end):
    if not os.path.exists(clasp_file):
        return
    kept = []
    removed = 0
    with open(clasp_file) as f:
        for line in f:
            if line.startswith('#'):
                kept.append(line)
                continue
            parts = line.split()
            if len(parts) >= 2:
                hdr = parts[1].split('$')
                try:
                    if int(hdr[1]) == abs_start and int(hdr[-1]) == abs_end:
                        removed += 1
                        continue
                except (IndexError, ValueError):
                    pass
            kept.append(line)
    with open(clasp_file, 'w') as f:
        f.writelines(kept)
    if removed:
        print(f'[ROUND2-COORD] Removed {removed} parent entries ({abs_start}-{abs_end}) from {os.path.basename(clasp_file)}')

def coordinate_round2_markers(work_dir, org, split_file, db, word_size, mem_limit_bytes):
    timeout_sec = get_timeout_seconds('self')

    markers = {}
    for orient in ['forward', 'reverse']:
        retry_dir = os.path.join(work_dir, f'blast_out_{orient}', org, 'retry', split_file)
        markers[orient] = _read_round2_markers(retry_dir)

    fwd_starts = set(markers['forward'].keys())
    rev_starts = set(markers['reverse'].keys())
    if fwd_starts == rev_starts:
        return

    for orient_with, orient_without in [('forward', 'reverse'), ('reverse', 'forward')]:
        one_sided = set(markers[orient_with].keys()) - set(markers[orient_without].keys())
        if not one_sided:
            continue

        is_rev = (orient_without == 'reverse')
        split_out = os.path.join(work_dir, f'split_out_{orient_without}', org, split_file)
        clasp_outfile = os.path.join(work_dir, f'clasp_out_{orient_without}', org, split_file)
        retry_dir = os.path.join(work_dir, f'blast_out_{orient_without}', org, 'retry', split_file)
        os.makedirs(retry_dir, exist_ok=True)

        # read all records once
        all_records = {int(rec.id.split('$')[1]): rec for rec in SeqIO.parse(split_out, 'fasta')}

        for p_start in sorted(one_sided):
            p_end, subs = markers[orient_with][p_start]
            if p_start not in all_records:
                print(f'[ROUND2-COORD] WARNING: seq abs_start={p_start} not found in {split_out}')
                continue

            print(f'[ROUND2-COORD] {orient_with} has Round 2 marker for {p_start}-{p_end}, '
                  f'window-splitting {orient_without}')

            # write single parent seq to temp, then window-split in place
            temp_fasta = os.path.join(retry_dir, f'coord_{p_start}.fasta')
            SeqIO.write([all_records[p_start]], temp_fasta, 'fasta')
            n_windows = split_sequences_to_windows(temp_fasta, is_reverse=is_rev)

            if n_windows is None:
                print(f'[ROUND2-COORD] Seq {p_start}-{p_end} too short for window split, skipping')
                continue

            # blast + clasp on the windowed file
            temp_blast = os.path.join(retry_dir, f'coord_{p_start}.blast')
            temp_clasp = os.path.join(retry_dir, f'coord_{p_start}.clasp')
            rc = blast(db, temp_fasta, temp_blast, word_size, mem_limit_bytes, timeout_sec)
            if rc == 0:
                rc = clasp(temp_blast, temp_clasp, mode='self',
                           mem_limit_bytes=mem_limit_bytes, timeout_sec=timeout_sec)
            if rc != 0:
                msg = f'[ROUND2-COORD] blast/clasp failed (rc={rc}) for coord windows of {p_start}'
                print(msg)
                _write_clasp_fail(temp_clasp, msg)

            # remove parent entries from main clasp output, append window entries
            _filter_clasp_entries(clasp_outfile, p_start, p_end)
            with open(clasp_outfile, 'a') as outf:
                if os.path.exists(temp_clasp):
                    with open(temp_clasp) as inf:
                        outf.write(inf.read())

            print(f'[ROUND2-COORD] Updated {orient_without} clasp for seq {p_start}-{p_end} ({n_windows} windows)')

def makeblastdb(org):
    path = getcwd()
    if not (isfile(path + "/blastdbs/{}.nin".format(org)) and isfile(path + "/blastdbs/{}.nhr".format(org)) and isfile(path + "/blastdbs/{}.nsq".format(org))):
        cmd = [
            'makeblastdb',
            '-in', f'genomes/{org}.fasta',
            '-out', f'blastdbs/{org}',
            '-dbtype', 'nucl'
        ]
        run(cmd)

def blast(db, query, outfile, word_size=None, mem_limit_bytes=None, timeout_sec=None):
    if word_size is None:
        word_size = _CONFIG['blast']['word_size']
    evalue = _CONFIG['blast']['evalue']
    strand = _CONFIG['blast']['strand']
    outfmt = _CONFIG['blast']['outfmt']
    cmd = [
        'blastn',
        '-db', db,
        '-query', query,
        '-strand', strand,
        '-outfmt', str(outfmt),
        '-evalue', str(evalue),
        '-word_size', str(word_size),
        '-out', outfile
    ]
    return _run_with_limit(cmd, mem_limit_bytes, timeout_sec)

def clasp(infile, outfile, l=None, e=None, mode='self', mem_limit_bytes=None, timeout_sec=None):
    if l is None:
        if mode == 'pairwise':
            l = _CONFIG['clasp'].get('l_pairwise', 2)
        else:
            l = _CONFIG['clasp']['l']
    if e is None:
        if mode == 'pairwise':
            e = _CONFIG['clasp'].get('e_pairwise', 0.1)
        else:
            e = _CONFIG['clasp']['e']
    columns_c = [str(c) for c in _CONFIG['clasp']['columns_c']]
    columns_C = [str(c) for c in _CONFIG['clasp']['columns_C']]
    cmd = ['clasp.x', '-m', '-i', infile, '-c'] + columns_c + ['-C'] + columns_C + ['-l', str(l), '-e', str(e), '-o', outfile]
    return _run_with_limit(cmd, mem_limit_bytes, timeout_sec)
