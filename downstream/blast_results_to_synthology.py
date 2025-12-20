#! /usr/bin/env python3

from subprocess import run
import os
import argparse
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

            # Find all hits that should be merged with hit i
            for j in range(i + 1, len(chrom_hits)):
                if merged[j]:
                    continue

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
        return(0)

    # parse temp_coords
    # Format: org chromo start end orientation gene_id [sequence]
    # sequence column can be empty for blastn (will extract from genome)
    hits = []
    with open(temp_coords_file,'r') as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split('\t')

            # temp_coords must have exactly 7 columns
            if len(cols) != 7:
                raise ValueError(f"temp_coords line has {len(cols)} columns, expected 7 (org chromo start end orientation gene_id sequence). Line: {line.strip()}")

            hit = {
                'org': cols[0],
                'chromo': cols[1],
                'start': int(cols[2]),
                'end': int(cols[3]),
                'orientation': cols[4],
                'gene_id': cols[5]
            }
            # check if sequence is provided (column 6, can be empty for blastn)
            if len(cols[6]) > 0:
                # sequence provided - use it directly (tblastn case)
                hit['seq'] = cols[6]
                hit['frame'] = None
            else:
                # no sequence - will extract from genome (blastn case)
                hit['seq'] = None
                hit['frame'] = None
            hits.append(hit)

    if len(hits) == 0:
        return(0)

    # Deduplicate overlapping hits if requested
    if deduplicate_overlaps:
        original_count = len(hits)
        hits = deduplicate_overlapping_hits(hits)
        if len(hits) < original_count:
            print(f"  {org}: Deduplicated {original_count} hits to {len(hits)} ({original_count - len(hits)} merged)")

    # load genome
    genome = get_genome(genomes_dir,org)

    # extract sequences
    sequences = {}
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
            sequences[gene_id] = (record,chromo,start,end)
            continue

        # otherwise extract from genome
        if chromo not in genome:
            continue

        nucl_seq = genome[chromo][start:end]
        if len(nucl_seq) == 0:
            continue

        if mode == 'protein':
            # no provided sequence - translate from genome
            try:
                prot_seq = translate_seq(nucl_seq,frame,orientation)
            except:
                continue
            record = SeqRecord(prot_seq, id=gene_id, description="")
            sequences[gene_id] = (record,chromo,start,end)
        else:
            # nucleotide mode - always use full genomic region
            if orientation == 'reverse':
                nucl_seq = nucl_seq.reverse_complement()
            record = SeqRecord(nucl_seq, id=gene_id, description="")
            sequences[gene_id] = (record,chromo,start,end)

    if len(sequences) == 0:
        return(0)

    # write gff
    run(f'mkdir -p {output_dir}/gff',shell=True)
    with open(f'{output_dir}/gff/{org}.gff','w') as f:
        for gene_id,(record,chromo,start,end) in sorted(sequences.items()):
            # Write proper GFF3 format: seqid source type start end score strand phase attributes
            # Use CDS as type so run_synthology_prot.py will parse it
            f.write(f'{chromo}\t.\tCDS\t{start}\t{end}\t.\t+\t.\tID={gene_id}\n')

    # write sequences
    extension = 'faa' if mode == 'protein' else 'fna'
    seq_dir = 'proteins' if mode == 'protein' else 'nucleotides'
    run(f'mkdir -p {output_dir}/{seq_dir}',shell=True)
    records = [record for gene_id,(record,chromo,start,end) in sorted(sequences.items())]
    SeqIO.write(records, f'{output_dir}/{seq_dir}/{org}.{extension}', 'fasta')

    return(len(sequences))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert BLAST results to synthology format',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_dir', type=str, help="directory with temp_coords files")
    parser.add_argument('--genomes_dir', type=str, nargs='?', default='../../utils/genomes', help="directory with genome FASTA files")
    parser.add_argument('--output_dir', type=str, nargs='?', default='annotations_for_synthology_from_regions', help="output directory")
    parser.add_argument('--mode', type=str, nargs='?', default='protein', choices=['protein','nucleotide'], help="protein or nucleotide mode")
    parser.add_argument('--deduplicate_overlaps', action='store_true', help="merge overlapping hits that likely represent the same element (useful for multi-copy gene families)")

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

    for org in orgs:
        n = process_org(org,input_dir,genomes_dir,output_dir,mode,deduplicate_overlaps)
        print(f'{org}: {n} sequences extracted')

    print('Done')
