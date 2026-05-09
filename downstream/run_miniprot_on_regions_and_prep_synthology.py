#!/usr/bin/env python3
"""
run_miniprot_on_regions_and_prep_synthology.py

Run miniprot on extracted genomic regions and produce synthology-format output.

This script combines alignment and post-processing in one step:
1. Runs miniprot protein-to-genome alignment on extracted region FASTAs (with --trans)
2. Converts coordinates from region-relative to genome-absolute
3. Parses multi-exon genes and translated proteins from miniprot GFF3 output
4. Filters by identity and frameshift thresholds
5. Optionally deduplicates overlapping hits
6. Writes synthology-format GFF + protein FASTA

Miniprot's --trans flag provides translated protein sequences directly via ##STA
comment lines in GFF output, eliminating the need for genome extraction/re-translation.

Usage:
    python3 run_miniprot_on_regions_and_prep_synthology.py query.faa \\
        --fastas_dir fastas \\
        --output_dir annotations_for_synthology_from_miniprot \\
        --cores 8 \\
        --min_identity 0.3 \\
        --max_frameshifts 3

Output:
    - {output_dir}/gff/{org}.gff       : Synthology-format GFF files
    - {output_dir}/proteins/{org}.faa  : Protein FASTA files
    - {output_dir}/miniprot_gff/{org}_hits.gff : Intermediate converted GFF (for debugging)

Next step:
    python3 run_synthology_prot.py --annotation_dir {output_dir}
"""

import os
import sys
import argparse
import multiprocessing as mp
from subprocess import run as subprocess_run, PIPE, CalledProcessError
from os.path import isfile, isdir
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Coordinate conversion (region-relative -> genome-absolute)
# ---------------------------------------------------------------------------

def parse_region_seqid(seqid, org):
    """
    Parse region-based seqid to extract chromosome, coordinates, and direction.

    Region seqids are created by extract_syntenic_fastas.py with format:
        {org}_{chromo}_{start}_{end}_{direction}
    e.g.: Capsicum_annuum_Chr07_67354775_67647425_forward

    Args:
        seqid: Sequence identifier from FASTA header
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


def convert_gff_line(line, org):
    """
    Convert a miniprot GFF3 line from region-relative to genome-absolute coordinates.

    Region FASTAs have 0-based start coordinates (Python slicing).
    Miniprot GFF features use 1-based coordinates relative to the region.
    Conversion: genome_pos_1based = region_start_0based + local_pos_1based

    Args:
        line: A GFF3 line (may be feature or comment)
        org: Organism name

    Returns:
        Converted line string, or None to drop the line (e.g. ##sequence-region)
    """
    if not line.strip():
        return line

    if line.startswith('#'):
        # Drop ##sequence-region lines — they reference region seqids and
        # their coordinates are meaningless after conversion
        if line.startswith('##sequence-region'):
            return None
        return line

    cols = line.rstrip('\n').split('\t')
    if len(cols) != 9:
        return line  # Not a valid GFF feature line, pass through

    seqid = cols[0]
    region_info = parse_region_seqid(seqid, org)
    if region_info is None:
        return line  # Not a region seqid, pass through unchanged

    chr_name, region_start, region_end, direction = region_info
    local_start = int(cols[3])
    local_end = int(cols[4])
    strand = cols[6]

    if direction == 'forward':
        genome_start = region_start + local_start
        genome_end = region_start + local_end
        genome_strand = strand
    else:  # reverse
        genome_start = region_end - local_end + 1
        genome_end = region_end - local_start + 1
        genome_strand = '-' if strand == '+' else '+'

    cols[0] = chr_name
    cols[3] = str(genome_start)
    cols[4] = str(genome_end)
    cols[6] = genome_strand

    return '\t'.join(cols) + '\n'


# ---------------------------------------------------------------------------
# Miniprot execution
# ---------------------------------------------------------------------------

def check_miniprot_installation(miniprot_path):
    """
    Verify miniprot is installed and accessible.

    Args:
        miniprot_path: Path to miniprot binary

    Returns:
        True if installed, False otherwise
    """
    try:
        result = subprocess_run([miniprot_path, '--version'],
                                capture_output=True, text=True, timeout=5)
        version_info = result.stdout + result.stderr
        print(f"Found miniprot: {version_info.strip()}")
        return True
    except (FileNotFoundError, CalledProcessError, Exception) as e:
        print(f"ERROR: miniprot not found at '{miniprot_path}'")
        print(f"Details: {e}")
        print("\nInstallation instructions:")
        print("  1. From source: git clone https://github.com/lh3/miniprot && cd miniprot && make")
        print("  2. Via conda: conda install -c bioconda miniprot")
        print("  3. Specify path: --miniprot_path /path/to/miniprot")
        return False


# ---------------------------------------------------------------------------
# GFF3 parsing
# ---------------------------------------------------------------------------

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


def parse_miniprot_gff_lines(gff_lines):
    """
    Parse miniprot GFF3 lines (in memory) and group features by parent mRNA.

    Miniprot outputs (with --trans flag):
    - ##STA comment lines with translated protein sequences
    - mRNA features (parent alignments)
    - CDS features (individual exons, linked by Parent attribute)
    - stop_codon features (if present)

    Each ##STA line appears immediately before its corresponding mRNA feature.

    Args:
        gff_lines: List of GFF3 line strings

    Returns:
        Dictionary mapping mRNA ID to gene data dict with 'mRNA', 'CDS',
        and 'protein' keys
    """
    genes = defaultdict(lambda: {
        'mRNA': None,
        'CDS': [],
        'protein': None
    })

    pending_protein = None  # ##STA protein waiting for its mRNA

    for line in gff_lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('##STA\t'):
            # Miniprot --trans output: ##STA\t<translated_protein_sequence>
            pending_protein = line.split('\t', 1)[1].strip()
            continue

        if line.startswith('#'):
            continue

        cols = line.split('\t')
        if len(cols) != 9:
            continue

        seqid, source, feature_type, start, end, score, strand, phase, attributes = cols
        start, end = int(start), int(end)

        attrs = parse_gff_attributes(attributes)

        if feature_type == 'mRNA':
            mrna_id = attrs.get('ID', '')
            if mrna_id:
                genes[mrna_id]['mRNA'] = {
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'attributes': attrs
                }
                # Associate pending ##STA protein with this mRNA
                if pending_protein is not None:
                    genes[mrna_id]['protein'] = pending_protein
                    pending_protein = None

        elif feature_type == 'CDS':
            parent_id = attrs.get('Parent', '')
            if parent_id:
                genes[parent_id]['CDS'].append({
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'phase': phase,
                    'attributes': attrs
                })

    # Filter genes with CDS features
    result = {}
    for mrna_id, gene_data in genes.items():
        if len(gene_data['CDS']) > 0:
            result[mrna_id] = gene_data

    return result


def extract_quality_metrics(gene_data):
    """
    Extract quality metrics from miniprot alignment.

    Args:
        gene_data: Gene dictionary from parse_miniprot_gff_lines

    Returns:
        Dictionary with identity, frameshift, stop_codon counts
    """
    if gene_data['mRNA']:
        attrs = gene_data['mRNA']['attributes']
        identity = float(attrs.get('Identity', '0'))
        frameshift = int(attrs.get('Frameshift', '0'))
        stop_codons = int(attrs.get('StopCodon', '0'))
        rank = int(attrs.get('Rank', '0'))
        target = attrs.get('Target', '')
        return {
            'identity': identity,
            'frameshift': frameshift,
            'stop_codons': stop_codons,
            'rank': rank,
            'target': target
        }
    else:
        return {
            'identity': 0.0,
            'frameshift': 0,
            'stop_codons': 0,
            'rank': 0,
            'target': ''
        }


def merge_cds_to_single_range(cds_list):
    """
    Merge multi-exon CDS features into single genomic range.

    Strategy: Use earliest start and latest end coordinates.
    This captures the full gene span including introns.

    Args:
        cds_list: List of CDS feature dictionaries

    Returns:
        Dictionary with merged coordinates: seqid, start, end, strand
    """
    if not cds_list:
        return None

    seqid = cds_list[0]['seqid']
    strand = cds_list[0]['strand']

    min_start = min(cds['start'] for cds in cds_list)
    max_end = max(cds['end'] for cds in cds_list)

    return {
        'seqid': seqid,
        'start': min_start,
        'end': max_end,
        'strand': strand
    }


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def deduplicate_overlapping_hits(hits, min_reciprocal_overlap=0.5, max_size_ratio=2.0):
    """
    Deduplicate overlapping hits that likely represent the same genomic element.

    Args:
        hits: List of hit dictionaries with keys: seqid, start, end, seq, gene_id, etc.
        min_reciprocal_overlap: Minimum reciprocal overlap to consider merging
        max_size_ratio: Maximum size ratio to consider merging

    Returns:
        Deduplicated list of hits
    """
    by_chrom = defaultdict(list)
    for hit in hits:
        by_chrom[hit['seqid']].append(hit)

    deduplicated = []

    for seqid, chrom_hits in by_chrom.items():
        chrom_hits.sort(key=lambda h: h['start'])

        merged = [False] * len(chrom_hits)

        for i in range(len(chrom_hits)):
            if merged[i]:
                continue

            cluster = [i]
            max_end = chrom_hits[i]['end']

            for j in range(i + 1, len(chrom_hits)):
                if merged[j]:
                    continue

                if chrom_hits[j]['start'] > max_end:
                    break

                should_merge = False
                for cluster_idx in cluster:
                    h1 = chrom_hits[cluster_idx]
                    h2 = chrom_hits[j]

                    overlap_start = max(h1['start'], h2['start'])
                    overlap_end = min(h1['end'], h2['end'])
                    overlap_length = max(0, overlap_end - overlap_start)

                    if overlap_length > 0:
                        len1 = h1['end'] - h1['start']
                        len2 = h2['end'] - h2['start']
                        recip_overlap = overlap_length / min(len1, len2)
                        size_ratio = max(len1, len2) / min(len1, len2)

                        if recip_overlap >= min_reciprocal_overlap and size_ratio <= max_size_ratio:
                            should_merge = True
                            break

                if should_merge:
                    cluster.append(j)
                    merged[j] = True
                    max_end = max(max_end, chrom_hits[j]['end'])

            cluster_hits = [chrom_hits[idx] for idx in cluster]

            def score_hit(h):
                if h['seq'] is not None and len(h['seq']) > 0:
                    return len(h['seq'])
                else:
                    return h['end'] - h['start']

            best_hit = max(cluster_hits, key=score_hit)

            merged_hit = {
                'seqid': best_hit['seqid'],
                'start': min(h['start'] for h in cluster_hits),
                'end': max(h['end'] for h in cluster_hits),
                'strand': best_hit['strand'],
                'gene_id': best_hit['gene_id'],
                'seq': best_hit['seq'],
                'identity': best_hit['identity'],
                'frameshift': best_hit['frameshift']
            }

            deduplicated.append(merged_hit)
            merged[i] = True

    return deduplicated


# ---------------------------------------------------------------------------
# Per-organism processing (miniprot + synthology in one step)
# ---------------------------------------------------------------------------

def process_org(org, query, fastas_dir, output_dir,
                miniprot_path, threads, miniprot_params,
                min_identity, max_frameshifts, deduplicate_overlaps):
    """
    Run miniprot (with --trans) and convert to synthology format for one organism.

    Uses miniprot's native --trans flag to obtain translated protein sequences
    directly via ##STA comment lines, avoiding genome re-extraction.

    Args:
        org: Organism identifier
        query: Path to protein query FASTA file
        fastas_dir: Directory containing extracted region FASTAs
        output_dir: Base output directory
        miniprot_path: Path to miniprot binary
        threads: Number of threads for miniprot
        miniprot_params: Additional miniprot parameters (string)
        min_identity: Minimum alignment identity (0-1)
        max_frameshifts: Maximum frameshifts allowed (None = unlimited)
        deduplicate_overlaps: Whether to merge overlapping hits

    Returns:
        Number of sequences written to synthology output
    """
    # --- 1. Check input ---
    genome_fasta = f'{fastas_dir}/{org}/all.fasta'
    if not isfile(genome_fasta):
        print(f"WARNING: {org} - region FASTA not found at {genome_fasta}")
        return 0

    # --- 2. Run miniprot ---
    cmd_parts = [miniprot_path, '-t', str(threads), '--gff', '--trans']
    if miniprot_params:
        cmd_parts.extend(miniprot_params.split())
    cmd_parts.extend([genome_fasta, query])

    try:
        result = subprocess_run(cmd_parts, stdout=PIPE, stderr=PIPE,
                                text=True, timeout=3600)
    except Exception as e:
        print(f"ERROR: {org} - miniprot failed: {e}")
        return 0

    # --- 3. Convert coordinates to genome-absolute ---
    converted_lines = []
    for line in result.stdout.splitlines(keepends=True):
        converted = convert_gff_line(line, org)
        if converted is not None:
            converted_lines.append(converted)

    # Write intermediate GFF for debugging
    os.makedirs(f'{output_dir}/miniprot_gff', exist_ok=True)
    with open(f'{output_dir}/miniprot_gff/{org}_hits.gff', 'w') as f:
        f.writelines(converted_lines)

    # --- 4. Parse GFF from memory ---
    genes = parse_miniprot_gff_lines(converted_lines)

    if len(genes) == 0:
        print(f"WARNING: {org}: 0 mRNA/CDS features from miniprot.")
        return 0

    # --- 5. Filter and collect hits ---
    hits = []
    n_filtered_identity = 0
    n_filtered_frameshift = 0
    n_no_protein = 0

    for gene_idx, (mrna_id, gene_data) in enumerate(genes.items(), start=1):
        metrics = extract_quality_metrics(gene_data)

        # Apply filters
        if metrics['identity'] < min_identity:
            n_filtered_identity += 1
            continue
        if max_frameshifts is not None and metrics['frameshift'] > max_frameshifts:
            n_filtered_frameshift += 1
            continue

        # Merge CDS to single range for GFF output coordinates
        merged_coords = merge_cds_to_single_range(gene_data['CDS'])
        if not merged_coords:
            raise ValueError(
                f"ERROR: CDS merge failed for {mrna_id} in {org} despite having "
                f"{len(gene_data['CDS'])} CDS features.\n"
                f"  CDS features: {gene_data['CDS']}"
            )

        genome_seqid = merged_coords['seqid']
        genome_start = merged_coords['start']
        genome_end = merged_coords['end']
        genome_strand = merged_coords['strand']

        # Get protein from miniprot --trans output (##STA line)
        protein_str = gene_data.get('protein')
        if not protein_str:
            n_no_protein += 1
            print(f"WARNING: {org}/{mrna_id}: no ##STA protein from miniprot at "
                  f"{genome_seqid}:{genome_start}-{genome_end} ({genome_strand}). "
                  f"This should not happen with --trans flag.")
            continue

        # Strip trailing stop codon '*' if present
        protein_clean = protein_str.rstrip('*')
        if not protein_clean:
            n_no_protein += 1
            print(f"WARNING: {org}/{mrna_id}: ##STA protein is empty after stripping stop codons at "
                  f"{genome_seqid}:{genome_start}-{genome_end} ({genome_strand}).")
            continue
        protein_seq = Seq(protein_clean)

        # Generate gene ID from Target attribute
        target_info = metrics.get('target', '')
        query_name = target_info.split()[0] if target_info else 'unknown'
        gene_id = f'hit-{gene_idx}_{org}_{query_name}'

        hits.append({
            'seqid': genome_seqid,
            'start': genome_start,
            'end': genome_end,
            'strand': genome_strand,
            'gene_id': gene_id,
            'seq': protein_seq,
            'identity': metrics['identity'],
            'frameshift': metrics['frameshift']
        })

    if len(hits) == 0:
        print(f"WARNING: {org}: 0 hits passed filters from {len(genes)} parsed genes "
              f"(identity: {n_filtered_identity}, frameshift: {n_filtered_frameshift}, "
              f"no_protein: {n_no_protein})")
        return 0

    # --- 6. Deduplicate ---
    if deduplicate_overlaps:
        original_count = len(hits)
        hits = deduplicate_overlapping_hits(hits)
        if len(hits) < original_count:
            print(f"  {org}: Deduplicated {original_count} hits to {len(hits)} "
                  f"({original_count - len(hits)} merged)")

    # --- 7. Write synthology output ---
    os.makedirs(f'{output_dir}/gff', exist_ok=True)
    with open(f'{output_dir}/gff/{org}.gff', 'w') as f:
        for hit in sorted(hits, key=lambda h: (h['seqid'], h['start'])):
            f.write(f"{hit['seqid']}\t.\tCDS\t{hit['start']}\t{hit['end']}\t.\t"
                    f"{hit['strand']}\t0\tID={hit['gene_id']}\n")

    os.makedirs(f'{output_dir}/proteins', exist_ok=True)
    records = [SeqRecord(hit['seq'], id=hit['gene_id'], description='')
               for hit in sorted(hits, key=lambda h: h['gene_id'])]
    SeqIO.write(records, f'{output_dir}/proteins/{org}.faa', 'fasta')

    return len(hits)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Run miniprot on regions and prepare synthology format output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Miniprot arguments
    parser.add_argument('query', type=str,
                       help="Protein query FASTA file")
    parser.add_argument('--fastas_dir', type=str, default='fastas',
                       help="Directory with organism FASTAs (each in {org}/all.fasta)")
    parser.add_argument('--output_dir', type=str,
                       default='annotations_for_synthology_from_miniprot',
                       help="Output directory for synthology format files")
    parser.add_argument('--cores', type=int, default=1,
                       help="Number of parallel organism processes")
    parser.add_argument('--threads_per_job', type=int, default=1,
                       help="Threads per miniprot job")
    parser.add_argument('--miniprot_path', type=str, default='miniprot',
                       help="Path to miniprot binary")
    parser.add_argument('--miniprot_params', type=str, default='-I -j 1',
                       help="Additional miniprot parameters (e.g., '-G 50000 -I -j 1')")

    # Filtering arguments
    parser.add_argument('--min_identity', type=float, default=0.3,
                       help="Minimum alignment identity (0-1) to keep hits")
    parser.add_argument('--max_frameshifts', type=int, default=3,
                       help="Maximum frameshifts allowed (-1 for unlimited)")
    parser.add_argument('--deduplicate_overlaps', action='store_true',
                       help="Merge overlapping hits at the same locus (e.g., redundant splice variants)")

    args = parser.parse_args()

    # Handle -1 = unlimited
    max_frameshifts = None if args.max_frameshifts < 0 else args.max_frameshifts

    # Validate inputs
    if not isfile(args.query):
        print(f"ERROR: Query file not found: {args.query}")
        sys.exit(1)

    if not isdir(args.fastas_dir):
        print(f"ERROR: FASTAs directory not found: {args.fastas_dir}")
        sys.exit(1)

    # Check miniprot installation
    if not check_miniprot_installation(args.miniprot_path):
        sys.exit(1)

    # Create output directories
    os.makedirs(f'{args.output_dir}/miniprot_gff', exist_ok=True)
    os.makedirs(f'{args.output_dir}/gff', exist_ok=True)
    os.makedirs(f'{args.output_dir}/proteins', exist_ok=True)

    # Get list of organisms
    orgs = []
    if isfile('orgs'):
        with open('orgs') as f:
            for line in f:
                org = line.strip()
                if org and isfile(f'{args.fastas_dir}/{org}/all.fasta'):
                    orgs.append(org)
    else:
        # Fallback: list all subdirectories with all.fasta
        if isdir(args.fastas_dir):
            for org in os.listdir(args.fastas_dir):
                org_dir = f'{args.fastas_dir}/{org}'
                if isdir(org_dir) and isfile(f'{org_dir}/all.fasta'):
                    orgs.append(org)

    if len(orgs) == 0:
        print("ERROR: No organisms found with FASTA files")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Miniprot Alignment + Synthology Preparation")
    print(f"{'='*60}")
    print(f"Query:             {args.query}")
    print(f"FASTAs directory:  {args.fastas_dir}")
    print(f"Output directory:  {args.output_dir}")
    print(f"Organisms:         {len(orgs)}")
    print(f"Parallel jobs:     {args.cores}")
    print(f"Threads per job:   {args.threads_per_job}")
    print(f"Miniprot params:   {args.miniprot_params}")
    print(f"Min identity:      {args.min_identity}")
    print(f"Max frameshifts:   {max_frameshifts if max_frameshifts is not None else 'unlimited'}")
    print(f"Deduplication:     {'enabled' if args.deduplicate_overlaps else 'disabled'}")
    print(f"{'='*60}\n")

    if len(orgs) <= 10:
        print(f"Organisms: {', '.join(orgs)}\n")

    # Process all organisms
    total_sequences = 0

    if args.cores > 1:
        with mp.Pool(processes=args.cores) as pool:
            results = pool.starmap(process_org, [
                (org, args.query, args.fastas_dir, args.output_dir,
                 args.miniprot_path, args.threads_per_job, args.miniprot_params,
                 args.min_identity, max_frameshifts, args.deduplicate_overlaps)
                for org in orgs
            ])
        for org, n in zip(orgs, results):
            if n > 0:
                print(f"{org}: {n} sequences")
            total_sequences += n
    else:
        for org in orgs:
            n = process_org(
                org, args.query, args.fastas_dir, args.output_dir,
                args.miniprot_path, args.threads_per_job, args.miniprot_params,
                args.min_identity, max_frameshifts, args.deduplicate_overlaps
            )
            if n > 0:
                print(f"{org}: {n} sequences")
            total_sequences += n

    if total_sequences == 0:
        print(f"\nERROR: 0 total sequences across all {len(orgs)} organisms. "
              f"No output was produced. Check warnings above.")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Complete")
    print(f"{'='*60}")
    print(f"Total sequences:   {total_sequences}")
    print(f"GFF files:         {args.output_dir}/gff/")
    print(f"Protein files:     {args.output_dir}/proteins/")
    print(f"Debug GFF:         {args.output_dir}/miniprot_gff/")
    print(f"\nNext step:")
    print(f"  python3 run_synthology_prot.py --annotation_dir {args.output_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
