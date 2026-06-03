#! /usr/bin/env python3

import os
import sys
import argparse
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_genome(genomes_dir,org):
    seqs = SeqIO.parse(f'{genomes_dir}/{org}.fasta', "fasta")
    s = {}
    for seq in seqs:
        s[seq.id] = seq.seq
    return(s)

def translate_seq(seq,frame,orientation):
    # handle reverse orientation
    if orientation == 'reverse':
        seq = seq.reverse_complement()

    # translate based on frame
    if frame is None:
        # no frame - try all 3 and pick longest
        proteins = []
        for start_pos in [0,1,2]:
            prot = seq[start_pos:].translate(to_stop=False)
            proteins.append(prot)
        return max(proteins, key=len)
    else:
        # use specified frame
        frame_abs = abs(frame)
        start_pos = (frame_abs - 1) % 3
        return seq[start_pos:].translate(to_stop=False)

def deduplicate_overlapping_hits(hits, min_reciprocal_overlap=0.5, max_size_ratio=2.0):
    """
    Deduplicate overlapping hits that likely represent the same genomic element.

    Strategy:
    - Group hits by chromosome
    - Cluster overlapping hits based on reciprocal overlap and size ratio
    - Keep best representative (longest sequence) per cluster
    - Use union of coordinates for genomic range

    Args:
        hits: List of hit dictionaries with keys: chromo, start, end, seq, gene_id, etc.
        min_reciprocal_overlap: Minimum reciprocal overlap to consider merging (default 0.5)
        max_size_ratio: Maximum size ratio to consider merging (default 2.0)

    Returns:
        Deduplicated list of hits
    """
    from collections import defaultdict

    # Group hits by chromosome
    by_chrom = defaultdict(list)
    for hit in hits:
        by_chrom[hit['chromo']].append(hit)

    deduplicated = []

    for chromo, chrom_hits in by_chrom.items():
        # Sort by start position
        chrom_hits.sort(key=lambda h: h['start'])

        # Track which hits have been merged
        merged = [False] * len(chrom_hits)

        for i in range(len(chrom_hits)):
            if merged[i]:
                continue

            # Start a new cluster with hit i
            cluster = [i]
            max_end = chrom_hits[i]['end']

            # Find all hits that should be merged with hit i
            for j in range(i + 1, len(chrom_hits)):
                if merged[j]:
                    continue

                # Early exit: hits sorted by start, so no more overlaps possible
                if chrom_hits[j]['start'] > max_end:
                    break

                # Check if any hit in cluster overlaps with hit j
                should_merge = False
                for cluster_idx in cluster:
                    h1 = chrom_hits[cluster_idx]
                    h2 = chrom_hits[j]

                    # Calculate overlap
                    overlap_start = max(h1['start'], h2['start'])
                    overlap_end = min(h1['end'], h2['end'])
                    overlap_length = max(0, overlap_end - overlap_start)

                    if overlap_length > 0:
                        # Calculate reciprocal overlap
                        len1 = h1['end'] - h1['start']
                        len2 = h2['end'] - h2['start']
                        recip_overlap = overlap_length / min(len1, len2)

                        # Calculate size ratio
                        size_ratio = max(len1, len2) / min(len1, len2)

                        # Merge if reciprocal overlap is high AND size ratio is reasonable
                        if recip_overlap >= min_reciprocal_overlap and size_ratio <= max_size_ratio:
                            should_merge = True
                            break

                if should_merge:
                    cluster.append(j)
                    merged[j] = True
                    max_end = max(max_end, chrom_hits[j]['end'])

            # Process cluster: pick best representative
            cluster_hits = [chrom_hits[idx] for idx in cluster]

            # Score by sequence length (longer = better/more complete)
            # For hits without sequences, use genomic length
            def score_hit(h):
                if h['seq'] is not None and len(h['seq']) > 0:
                    return len(h['seq'])
                else:
                    return h['end'] - h['start']

            best_hit = max(cluster_hits, key=score_hit)

            # Create merged hit with:
            # - Union of coordinates (earliest start, latest end)
            # - Best hit's sequence
            # - Best hit's gene_id
            merged_hit = {
                'org': best_hit['org'],
                'chromo': best_hit['chromo'],
                'start': min(h['start'] for h in cluster_hits),
                'end': max(h['end'] for h in cluster_hits),
                'orientation': best_hit['orientation'],
                'gene_id': best_hit['gene_id'],
                'seq': best_hit['seq'],
                'frame': best_hit.get('frame')
            }

            deduplicated.append(merged_hit)
            merged[i] = True

    return deduplicated

def process_org(org,input_dir,genomes_dir,output_dir,mode,deduplicate_overlaps=False):

    temp_coords_file = f'{input_dir}/temp_coords_{org}'
    if not os.path.isfile(temp_coords_file):
        raise FileNotFoundError(
            f"ERROR: temp_coords file not found: {temp_coords_file}\n"
            f"  This file should have been created by run_blast_on_regions.py for {org}.\n"
            f"  Expected format: org\\tchromo\\tstart\\tend\\torientation\\tgene_id\\t[sequence]"
        )

    # parse temp_coords
    # Format: org chromo start end orientation gene_id [sequence]
    # sequence column can be empty for blastn (will extract from genome)
    hits = []
    with open(temp_coords_file,'r') as f:
        for line in f:
            if not line.strip():
                continue
            # Use rstrip('\n') not strip() — strip() eats the trailing tab that
            # marks an empty sequence column in blastn output, turning 7 cols -> 6.
            cols = line.rstrip('\n').split('\t')

            if len(cols) < 6:
                raise ValueError(f"temp_coords line has {len(cols)} columns, expected at least 6 "
                                 f"(org chromo start end orientation gene_id [sequence]). "
                                 f"Line: {line.rstrip()}")

            hit = {
                'org': cols[0],
                'chromo': cols[1],
                'start': int(cols[2]),
                'end': int(cols[3]),
                'orientation': cols[4],
                'gene_id': cols[5]
            }
            # sequence column is optional: present for tblastn, empty/absent for blastn
            seq_val = cols[6] if len(cols) > 6 else ''
            if seq_val:
                # sequence provided - use it directly (tblastn case)
                hit['seq'] = seq_val
                hit['frame'] = None
            else:
                # no sequence - will extract from genome (blastn case)
                hit['seq'] = None
                hit['frame'] = None
            hits.append(hit)

    if len(hits) == 0:
        print(f"WARNING: {org}: temp_coords file has 0 parseable hits: {temp_coords_file}\n"
              f"  Expected format: org\\tchromo\\tstart\\tend\\torientation\\tgene_id\\t[sequence]")
        return(0)

    # Deduplicate overlapping hits if requested
    if deduplicate_overlaps:
        original_count = len(hits)
        hits = deduplicate_overlapping_hits(hits)
        if len(hits) < original_count:
            print(f"  {org}: Deduplicated {original_count} hits to {len(hits)} ({original_count - len(hits)} merged)")

    # load genome lazily (skip if all hits have pre-provided sequences)
    genome = None

    # extract sequences
    sequences = {}
    n_empty_seq = 0
    n_translation_failed = 0
    for hit in hits:
        chromo = hit['chromo']
        start = hit['start']
        end = hit['end']
        gene_id = hit['gene_id']
        orientation = hit['orientation']
        frame = hit.get('frame')
        provided_seq = hit.get('seq')

        # if sequence provided in temp_coords and protein mode, use it directly
        if provided_seq is not None and mode == 'protein':
            seq_obj = Seq(provided_seq)
            record = SeqRecord(seq_obj, id=gene_id, description="")
            sequences[gene_id] = (record,chromo,start,end,orientation)
            continue

        # otherwise extract from genome
        if genome is None:
            genome = get_genome(genomes_dir, org)
        if chromo not in genome:
            raise ValueError(
                f"ERROR: Chromosome '{chromo}' not found in genome for {org}.\n"
                f"  Hit: {gene_id} at {chromo}:{start}-{end}\n"
                f"  Genome chromosomes: {list(genome.keys())[:10]}{'...' if len(genome) > 10 else ''}\n"
                f"  temp_coords file: {temp_coords_file}"
            )

        nucl_seq = genome[chromo][start-1:end]  # start is 1-based (from BLAST/CLASP); Python slicing is 0-based
        if len(nucl_seq) == 0:
            n_empty_seq += 1
            print(f"WARNING: {org}/{gene_id}: 0-length nucleotide sequence at {chromo}:{start}-{end} "
                  f"(chromosome length: {len(genome[chromo])}). This should not happen — check coordinates.")
            continue

        if mode == 'protein':
            # no provided sequence - translate from genome
            try:
                prot_seq = translate_seq(nucl_seq,frame,orientation)
            except Exception as e:
                n_translation_failed += 1
                print(f"WARNING: {org}/{gene_id}: Translation failed at {chromo}:{start}-{end} "
                      f"(seq length: {len(nucl_seq)}, orientation: {orientation}, frame: {frame}): {e}")
                continue
            record = SeqRecord(prot_seq, id=gene_id, description="")
            sequences[gene_id] = (record,chromo,start,end,orientation)
        else:
            # nucleotide mode - always use full genomic region
            if orientation == 'reverse':
                nucl_seq = nucl_seq.reverse_complement()
            # gene_id already contains org (e.g. hit-1_Org_Query from run_blast_on_regions.py)
            # Use it directly — same as protein mode — so GFF ID= and FASTA header match.
            record = SeqRecord(nucl_seq, id=gene_id, description="")
            sequences[gene_id] = (record,chromo,start,end,orientation)

    if len(sequences) == 0:
        print(f"WARNING: {org}: 0 sequences extracted from {len(hits)} hits"
              f" (empty_seq: {n_empty_seq}, translation_failed: {n_translation_failed}).\n"
              f"  Expected: temp_coords entries with valid genome coordinates or pre-provided sequences.")
        return(0)

    # write gff
    os.makedirs(f'{output_dir}/gff', exist_ok=True)
    with open(f'{output_dir}/gff/{org}.gff','w') as f:
        for gene_id,(record,chromo,start,end,orientation) in sorted(sequences.items()):
            # Write proper GFF3 format: seqid source type start end score strand phase attributes
            # Use CDS as type so run_synthology_prot.py will parse it
            strand_symbol = '+' if orientation == 'forward' else '-'
            f.write(f'{chromo}\t.\tCDS\t{start}\t{end}\t.\t{strand_symbol}\t0\tID={gene_id}\n')

    # write sequences
    extension = 'faa' if mode == 'protein' else 'fna'
    seq_dir = 'proteins' if mode == 'protein' else 'seqs'
    os.makedirs(f'{output_dir}/{seq_dir}', exist_ok=True)
    records = [record for gene_id,(record,chromo,start,end,orientation) in sorted(sequences.items())]
    SeqIO.write(records, f'{output_dir}/{seq_dir}/{org}.{extension}', 'fasta')

    return(len(sequences))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert BLAST results to synthology format',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir', type=str, help="directory with temp_coords files")
    parser.add_argument('--genomes_dir', type=str, nargs='?', default='../../utils/genomes', help="directory with genome FASTA files")
    parser.add_argument('--output_dir', type=str, nargs='?', default='annotations_for_synthology_from_regions', help="output directory")
    parser.add_argument('--mode', type=str, nargs='?', default='protein', choices=['protein','nucleotide'], help="protein or nucleotide mode")
    parser.add_argument('--deduplicate_overlaps', action='store_true', help="merge overlapping hits that likely represent the same element (useful for multi-copy gene families)")
    parser.add_argument('--cores', type=int, default=1, help="number of parallel processes (default: 1)")

    args = parser.parse_args()

    input_dir = args.input_dir
    genomes_dir = args.genomes_dir
    output_dir = args.output_dir
    mode = args.mode
    deduplicate_overlaps = args.deduplicate_overlaps

    path = os.getcwd()

    orgs = []
    if os.path.isfile('orgs'):
        with open('orgs') as f:
            for line in f:
                org = line.strip()
                if os.path.isfile(f'{input_dir}/temp_coords_{org}'):
                    orgs.append(org)
    else:
        # fallback: list all temp_coords files
        if os.path.isdir(input_dir):
            for filename in os.listdir(input_dir):
                if filename.startswith('temp_coords_'):
                    org = filename.replace('temp_coords_','')
                    orgs.append(org)

    print(f'{len(orgs)} organisms with temp_coords files')

    if deduplicate_overlaps:
        print("Deduplication of overlapping hits is ENABLED")

    # Pre-create output directories
    os.makedirs(f'{output_dir}/gff', exist_ok=True)
    seq_dir = 'proteins' if mode == 'protein' else 'seqs'
    os.makedirs(f'{output_dir}/{seq_dir}', exist_ok=True)

    total_sequences = 0
    if args.cores > 1:
        with mp.Pool(processes=args.cores) as pool:
            results = pool.starmap(process_org, [
                (org, input_dir, genomes_dir, output_dir, mode, deduplicate_overlaps)
                for org in orgs
            ])
        for org, n in zip(orgs, results):
            print(f'{org}: {n} sequences extracted')
            total_sequences += n
    else:
        for org in orgs:
            n = process_org(org, input_dir, genomes_dir, output_dir, mode, deduplicate_overlaps)
            print(f'{org}: {n} sequences extracted')
            total_sequences += n

    if total_sequences == 0:
        print(f"\nERROR: 0 total sequences extracted across all {len(orgs)} organisms. "
              f"Something is fundamentally wrong with the input data or genome files.")
        sys.exit(1)

    print(f'\nDone. Total: {total_sequences} sequences from {len(orgs)} organisms.')
