#!/usr/bin/env python3
"""
annotation_hits_to_synthology.py

Filter annotation BLAST hits by bitscore and convert to synthology format.

This script:
1. Filters hits by bitscore threshold (default: 30)
2. Removes overlapping genes (deduplication)
3. Copies passing hits from annotation files to synthology format
4. Writes GFF3 + protein FASTA compatible with run_synthology_prot.py

Usage:
    python3 annotation_hits_to_synthology.py \\
        --blast_hits blast_results \\
        --annotations_dir annotations_raw \\
        --output_dir annotations_for_synthology_from_annotation \\
        --min_bitscore 30 \\
        --deduplicate_overlaps

Output:
    - annotations_for_synthology_from_annotation/gff/{org}.gff
    - annotations_for_synthology_from_annotation/proteins/{org}.faa

Next step:
    Use with run_synthology_prot.py:
    python3 run_synthology_prot.py --annotation_dir annotations_for_synthology_from_annotation
"""

import os
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_gff_attributes(attr_string):
    """
    Parse GFF3 attributes field.

    Args:
        attr_string: Semicolon-separated key=value pairs

    Returns:
        Dictionary of attributes
    """
    attrs = {}
    if not attr_string or attr_string == '.':
        return attrs

    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key] = value
    return attrs


def load_gff_genes(gff_file):
    """
    Load gene coordinates from GFF3 file.

    Args:
        gff_file: Path to GFF3 file

    Returns:
        Dictionary mapping gene ID to coordinate dict (seqid, start, end, strand)
    """
    genes = {}

    if not os.path.isfile(gff_file):
        return genes

    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = cols

            # Only process CDS/gene features
            if feature_type not in ['CDS', 'gene', 'mRNA', 'transcript']:
                continue

            attrs = parse_gff_attributes(attributes)
            gene_id = attrs.get('ID', '')

            if gene_id:
                genes[gene_id] = {
                    'seqid': seqid,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand
                }

    return genes


def load_proteins(faa_file):
    """
    Load protein sequences from FASTA file.

    Args:
        faa_file: Path to protein FASTA file

    Returns:
        Dictionary mapping gene ID to protein sequence
    """
    proteins = {}

    if not os.path.isfile(faa_file):
        return proteins

    for record in SeqIO.parse(faa_file, 'fasta'):
        proteins[record.id] = record.seq

    return proteins


def parse_region_seqid(seqid, org):
    """
    Parse region-based seqid to extract chromosome, coordinates, and direction.

    Region seqids are created by extract_syntenic_fastas.py with format:
        {org}_{chromo}_{start}_{end}_{direction}
    e.g.: Capsicum_annuum_Chr07_67354775_67647425_forward

    Args:
        seqid: Sequence identifier from GFF
        org: Organism name (prefix to strip)

    Returns:
        Tuple (chr_name, region_start, region_end, direction) or None if not a region seqid
    """
    prefix = f'{org}_'
    if not seqid.startswith(prefix):
        return None

    remainder = seqid[len(prefix):]
    parts = remainder.split('_')

    # Need at least: chromo, start, end, direction
    if len(parts) < 4:
        return None

    direction = parts[-1]
    if direction not in ('forward', 'reverse'):
        return None

    try:
        region_start = int(parts[-3])
        region_end = int(parts[-2])
    except ValueError:
        return None

    chr_name = '_'.join(parts[:-3])
    return (chr_name, region_start, region_end, direction)


def deduplicate_overlapping_hits(hits, min_reciprocal_overlap=0.5, max_size_ratio=2.0):
    """
    Deduplicate overlapping hits that likely represent the same genomic element.

    Reuses logic from blast_results_to_synthology.py (lines 36-142).

    Args:
        hits: List of hit dictionaries with keys: seqid, start, end, seq, gene_id, etc.
        min_reciprocal_overlap: Minimum reciprocal overlap to consider merging
        max_size_ratio: Maximum size ratio to consider merging

    Returns:
        Deduplicated list of hits
    """
    # Group hits by chromosome
    by_chrom = defaultdict(list)
    for hit in hits:
        by_chrom[hit['seqid']].append(hit)

    deduplicated = []

    for seqid, chrom_hits in by_chrom.items():
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

            # Process cluster: pick best representative (longest protein or highest bitscore)
            cluster_hits = [chrom_hits[idx] for idx in cluster]

            def score_hit(h):
                # Prefer longer proteins, then higher bitscore
                seq_len = len(h['seq']) if h['seq'] is not None else 0
                return (seq_len, h.get('bitscore', 0))

            best_hit = max(cluster_hits, key=score_hit)

            # Create merged hit with union of coordinates
            merged_hit = {
                'org': best_hit['org'],
                'seqid': best_hit['seqid'],
                'start': min(h['start'] for h in cluster_hits),
                'end': max(h['end'] for h in cluster_hits),
                'strand': best_hit['strand'],
                'gene_id': best_hit['gene_id'],
                'seq': best_hit['seq'],
                'query_id': best_hit['query_id'],
                'bitscore': best_hit.get('bitscore', 0)
            }

            deduplicated.append(merged_hit)
            merged[i] = True

    return deduplicated


def process_organism(org, blast_hits_dir, annotations_dir, output_dir, min_bitscore, deduplicate_overlaps):
    """
    Process one organism: filter hits, format for synthology.

    Args:
        org: Organism identifier
        blast_hits_dir: Directory with BLAST hit files
        annotations_dir: Directory with annotation files
        output_dir: Output directory for synthology format
        min_bitscore: Minimum bitscore threshold
        deduplicate_overlaps: Whether to merge overlapping hits

    Returns:
        Number of sequences written
    """
    hits_file = f'{blast_hits_dir}/hits_{org}.tsv'

    if not os.path.isfile(hits_file):
        raise FileNotFoundError(
            f"ERROR: BLAST hits file not found: {hits_file}\n"
            f"  Expected output from blast_query_on_annotations.py for {org}.\n"
            f"  Expected format: TSV with header 'org\\tquery_id\\tgene_id\\tbitscore'"
        )

    # Load GFF and proteins from annotation results
    gff_file = f'{annotations_dir}/gff_raw/{org}_predictions.gff'
    faa_file = f'{annotations_dir}/proteins_raw/{org}_predictions.faa'

    genes = load_gff_genes(gff_file)
    proteins = load_proteins(faa_file)

    # Convert region-relative coordinates to genome-absolute coordinates.
    # Annotation GFF uses region seqids (from gene predictor running on extracted FASTAs).
    # Output GFF should use original chromosome names and genome-level coordinates.
    for gene_id in list(genes.keys()):
        gene = genes[gene_id]
        region_info = parse_region_seqid(gene['seqid'], org)
        if region_info:
            chr_name, region_start, region_end, direction = region_info
            if direction == 'forward':
                # genome_pos_1based = region_start_0based + local_pos_1based
                gene['seqid'] = chr_name
                gene['start'] = region_start + gene['start']
                gene['end'] = region_start + gene['end']
            else:  # reverse
                orig_start, orig_end = gene['start'], gene['end']
                gene['seqid'] = chr_name
                gene['start'] = region_end - orig_end + 1
                gene['end'] = region_end - orig_start + 1
                gene['strand'] = '-' if gene['strand'] == '+' else '+'

    if len(genes) == 0 or len(proteins) == 0:
        # GFF files can be empty if the predictor found nothing for this organism's regions.
        # But if a hits file exists, BLAST found matches, so the annotation must have had content.
        print(f"WARNING: {org}: hits file exists but annotations are empty "
              f"(genes: {len(genes)}, proteins: {len(proteins)}).\n"
              f"  GFF: {gff_file}\n  FAA: {faa_file}\n"
              f"  This is unexpected — if BLAST found hits, the annotation files should have content.")
        return 0

    # Parse BLAST hits
    hits = []
    n_filtered_bitscore = 0
    with open(hits_file) as f:
        header = f.readline()  # Skip header
        for line in f:
            if not line.strip():
                continue

            cols = line.strip().split('\t')
            if len(cols) < 4:
                raise ValueError(
                    f"ERROR: Malformed line in hits file {hits_file}: {len(cols)} columns, expected 4.\n"
                    f"  Expected format: org\\tquery_id\\tgene_id\\tbitscore\n"
                    f"  This file is produced by blast_query_on_annotations.py and should always have 4 columns.\n"
                    f"  Line: {line.strip()}"
                )

            # Format: org query_id gene_id bitscore
            org_id = cols[0]
            query_id = cols[1]
            gene_id = cols[2]
            bitscore = float(cols[3])

            # Filter by bitscore
            if bitscore < min_bitscore:
                n_filtered_bitscore += 1
                continue

            # Get gene coordinates
            if gene_id not in genes:
                raise ValueError(
                    f"ERROR: gene_id '{gene_id}' from BLAST hits not found in GFF for {org}.\n"
                    f"  Hits file: {hits_file}\n"
                    f"  GFF file: {gff_file}\n"
                    f"  GFF has {len(genes)} genes. Sample IDs: {list(genes.keys())[:5]}\n"
                    f"  The gene IDs in the hits file must match those in the GFF from annotate_syn_regions.py."
                )

            gene_coords = genes[gene_id]

            # Get protein sequence
            if gene_id not in proteins:
                raise ValueError(
                    f"ERROR: gene_id '{gene_id}' found in GFF but not in protein FASTA for {org}.\n"
                    f"  GFF file: {gff_file}\n"
                    f"  FAA file: {faa_file}\n"
                    f"  Proteins has {len(proteins)} entries. Sample IDs: {list(proteins.keys())[:5]}\n"
                    f"  The protein FASTA should contain all genes from the GFF."
                )

            protein_seq = proteins[gene_id]

            hits.append({
                'org': org,
                'query_id': query_id,
                'gene_id': gene_id,
                'bitscore': bitscore,
                'seqid': gene_coords['seqid'],
                'start': gene_coords['start'],
                'end': gene_coords['end'],
                'strand': gene_coords['strand'],
                'seq': protein_seq
            })

    if len(hits) == 0:
        print(f"WARNING: {org}: 0 hits passed filters "
              f"(filtered by bitscore<{min_bitscore}: {n_filtered_bitscore}).")
        return 0

    # Deduplicate if requested
    if deduplicate_overlaps:
        original_count = len(hits)
        hits = deduplicate_overlapping_hits(hits)
        if len(hits) < original_count:
            print(f"  {org}: Deduplicated {original_count} hits to {len(hits)} ({original_count - len(hits)} merged)")

    # Renumber hits for synthology format
    for idx, hit in enumerate(hits, start=1):
        # Create new gene ID: hit-{N}_{org}_{query_protein}
        hit['new_gene_id'] = f"hit-{idx}_{org}_{hit['query_id']}"

    # Write GFF file
    os.makedirs(f'{output_dir}/gff', exist_ok=True)
    with open(f'{output_dir}/gff/{org}.gff', 'w') as f:
        for hit in sorted(hits, key=lambda h: (h['seqid'], h['start'])):
            f.write(f"{hit['seqid']}\t.\tCDS\t{hit['start']}\t{hit['end']}\t.\t{hit['strand']}\t0\tID={hit['new_gene_id']}\n")

    # Write protein FASTA file
    os.makedirs(f'{output_dir}/proteins', exist_ok=True)
    records = []
    for hit in sorted(hits, key=lambda h: h['new_gene_id']):
        record = SeqRecord(hit['seq'], id=hit['new_gene_id'], description='')
        records.append(record)

    SeqIO.write(records, f'{output_dir}/proteins/{org}.faa', 'fasta')

    return len(hits)


def main():
    parser = argparse.ArgumentParser(
        description='Filter annotation BLAST hits and convert to synthology format',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--blast_hits', type=str, required=True,
                       help="Directory with BLAST hit files (from blast_query_on_annotations.py)")
    parser.add_argument('--annotations_dir', type=str, default='annotations_raw',
                       help="Directory with annotation files (from annotate_syn_regions.py)")
    parser.add_argument('--output_dir', type=str, default='annotations_for_synthology_from_annotation',
                       help="Output directory for synthology format")
    parser.add_argument('--min_bitscore', type=float, default=30,
                       help="Minimum bitscore threshold to keep hits")
    parser.add_argument('--deduplicate_overlaps', action='store_true',
                       help="Merge overlapping hits (useful for multi-copy gene families)")
    parser.add_argument('--cores', type=int, default=1,
                       help="Number of parallel processes (default: 1)")

    args = parser.parse_args()

    # Validate inputs
    if not os.path.isdir(args.blast_hits):
        print(f"ERROR: BLAST hits directory not found: {args.blast_hits}")
        sys.exit(1)

    if not os.path.isdir(args.annotations_dir):
        print(f"ERROR: Annotations directory not found: {args.annotations_dir}")
        sys.exit(1)

    # Get list of organisms
    orgs = []
    if os.path.isfile('orgs'):
        with open('orgs') as f:
            for line in f:
                org = line.strip()
                if org and os.path.isfile(f'{args.blast_hits}/hits_{org}.tsv'):
                    orgs.append(org)
    else:
        # Fallback: list all hits_*.tsv files
        if os.path.isdir(args.blast_hits):
            for filename in os.listdir(args.blast_hits):
                if filename.startswith('hits_') and filename.endswith('.tsv'):
                    org = filename.replace('hits_', '').replace('.tsv', '')
                    orgs.append(org)

    if len(orgs) == 0:
        print("ERROR: No organisms found with BLAST hit files")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Annotation Hits to Synthology Conversion")
    print(f"{'='*60}")
    print(f"BLAST hits dir:    {args.blast_hits}")
    print(f"Annotations dir:   {args.annotations_dir}")
    print(f"Output directory:  {args.output_dir}")
    print(f"Organisms:         {len(orgs)}")
    print(f"Min bitscore:      {args.min_bitscore}")
    print(f"Deduplication:     {'enabled' if args.deduplicate_overlaps else 'disabled'}")
    print(f"{'='*60}\n")

    # Pre-create output directories
    os.makedirs(f'{args.output_dir}/gff', exist_ok=True)
    os.makedirs(f'{args.output_dir}/proteins', exist_ok=True)

    # Process each organism
    total_sequences = 0
    if args.cores > 1:
        with mp.Pool(processes=args.cores) as pool:
            results = pool.starmap(process_organism, [
                (org, args.blast_hits, args.annotations_dir, args.output_dir,
                 args.min_bitscore, args.deduplicate_overlaps)
                for org in orgs
            ])
        for org, n in zip(orgs, results):
            if n > 0:
                print(f"{org}: {n} sequences written")
            total_sequences += n
    else:
        for org in orgs:
            n = process_organism(
                org, args.blast_hits, args.annotations_dir, args.output_dir,
                args.min_bitscore, args.deduplicate_overlaps
            )
            if n > 0:
                print(f"{org}: {n} sequences written")
            total_sequences += n

    if total_sequences == 0:
        print(f"\nERROR: 0 total sequences across all {len(orgs)} organisms. "
              f"No output was produced. Check warnings above.")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Conversion Complete")
    print(f"{'='*60}")
    print(f"Total sequences:   {total_sequences}")
    print(f"GFF files:         {args.output_dir}/gff/")
    print(f"Protein files:     {args.output_dir}/proteins/")
    print(f"\nNext step:")
    print(f"  python3 run_synthology_prot.py --annotation_dir {args.output_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
