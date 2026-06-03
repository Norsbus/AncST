#!/usr/bin/env python3

"""
Filter Overlapping Sequences

Removes overlapping genes from GFF/FAA input files, keeping only the longest sequence.
This is useful for cleaning up annotations with multiple isoforms or alternative transcripts
before running synthology.

Usage:
    filter_overlapping_sequences.py <input_dir> <output_dir>

Input directory structure:
    input_dir/
    ├── gff/
    │   ├── species1.gff
    │   ├── species2.gff
    │   └── ...
    └── proteins/
        ├── species1.faa
        ├── species2.faa
        └── ...

Output directory structure (same as input):
    output_dir/
    ├── gff/
    │   ├── species1.gff
    │   ├── species2.gff
    │   └── ...
    └── proteins/
        ├── species1.faa
        ├── species2.faa
        └── ...

The script will:
1. Parse each species' GFF and FAA files
2. Merge CDS exons to get full coding regions
3. Remove shorter genes that overlap with longer ones
4. Write filtered GFF (only kept genes) and FAA (only kept sequences) to output
"""

import argparse
import os
import sys
from collections import defaultdict
from Bio import SeqIO


def parse_gff_faa(gff_path, faa_path, species):
    """
    Parse GFF and FAA files for a species.

    Args:
        gff_path: Path to GFF file
        faa_path: Path to FAA file
        species: Species name

    Returns:
        tuple: (parsed_genes, sequences, filtered_genes)
            - parsed_genes: dict {seqid: [(start, end, strand, gene_id)]}
            - sequences: dict {gene_id: sequence_string}
            - filtered_genes: list of (gene_id, reason, details) for filtered genes
    """
    print(f"  Parsing {species}...")

    # Load sequences from FAA
    sequences = {}
    gene_ids = set()
    for record in SeqIO.parse(faa_path, "fasta"):
        gene_id = record.id.split()[0]
        sequences[gene_id] = str(record.seq)
        gene_ids.add(gene_id)

    print(f"    Found {len(sequences)} sequences in FAA")

    # Parse GFF and match to sequences
    gene_intervals = defaultdict(lambda: defaultdict(list))

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = cols

            # Only parse CDS features
            if ftype != "CDS":
                continue

            start, end = int(start), int(end)

            # Extract gene ID from attributes
            attr_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in attrs.split(";") if "=" in kv}

            candidates = set()
            for key, value in attr_dict.items():
                candidates.add(value)
                if ":" in value:
                    candidates.add(value.split(":", 1)[1])

            matched_id = None
            for candidate in candidates:
                if candidate in gene_ids:
                    matched_id = candidate
                    break

            if not matched_id:
                continue

            gene_intervals[seqid][matched_id].append((start, end, strand))

    # Merge CDS exons to get full coding regions
    parsed = {}
    for seqid in gene_intervals:
        for gene_id in gene_intervals[seqid]:
            intervals = gene_intervals[seqid][gene_id]
            # Always merge CDS exons to get full coding region span
            min_start = min(s for s, e, strand in intervals)
            max_end = max(e for s, e, strand in intervals)
            strand = intervals[0][2]

            if seqid not in parsed:
                parsed[seqid] = []
            parsed[seqid].append((min_start, max_end, strand, gene_id))

    total_before = sum(len(genes) for genes in parsed.values())
    print(f"    Parsed {total_before} genes from GFF")

    # Filter overlapping genes (keep longest)
    filtered_genes = []
    for seqid in parsed:
        coords = parsed[seqid]
        coords.sort()

        kept = []
        for start, end, strand, pid in coords:
            overlap = False
            for i, (ks, ke, kstrand, kid) in enumerate(kept):
                if not (end < ks or start > ke):
                    overlap = True
                    new_len = end - start
                    old_len = ke - ks
                    if new_len > old_len:
                        # Current gene is longer - filter out the kept gene
                        filtered_genes.append((
                            kid,
                            'filtered_out_due_to_overlap_with_larger_seq',
                            {
                                'chromosome': seqid,
                                'position': f"{ks}-{ke}",
                                'length': old_len,
                                'overlaps_with': pid,
                                'overlapping_position': f"{start}-{end}",
                                'kept_length': new_len
                            }
                        ))
                        kept[i] = (start, end, strand, pid)
                    else:
                        # Current gene is shorter - filter it out
                        filtered_genes.append((
                            pid,
                            'filtered_out_due_to_overlap_with_larger_seq',
                            {
                                'chromosome': seqid,
                                'position': f"{start}-{end}",
                                'length': new_len,
                                'overlaps_with': kid,
                                'overlapping_position': f"{ks}-{ke}",
                                'kept_length': old_len
                            }
                        ))
                    break

            if not overlap:
                kept.append((start, end, strand, pid))

        parsed[seqid] = kept

    total_after = sum(len(genes) for genes in parsed.values())
    print(f"    After overlap filtering: {total_after} genes ({total_before - total_after} filtered)")

    return parsed, sequences, filtered_genes


def write_filtered_gff(original_gff_path, output_gff_path, kept_gene_ids):
    """
    Write filtered GFF file containing only kept genes.

    Args:
        original_gff_path: Path to original GFF file
        output_gff_path: Path to output GFF file
        kept_gene_ids: Set of gene IDs to keep
    """
    with open(original_gff_path) as infile, open(output_gff_path, 'w') as outfile:
        for line in infile:
            # Keep header lines
            if line.startswith("#"):
                outfile.write(line)
                continue

            if not line.strip():
                continue

            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = cols

            # Only filter CDS features
            if ftype != "CDS":
                outfile.write(line)
                continue

            # Extract gene ID from attributes
            attr_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in attrs.split(";") if "=" in kv}

            # Check if this gene should be kept
            keep = False
            for value in attr_dict.values():
                if value in kept_gene_ids:
                    keep = True
                    break
                if ":" in value and value.split(":", 1)[1] in kept_gene_ids:
                    keep = True
                    break

            if keep:
                outfile.write(line)


def write_filtered_faa(original_faa_path, output_faa_path, kept_gene_ids):
    """
    Write filtered FAA file containing only kept genes.

    Args:
        original_faa_path: Path to original FAA file
        output_faa_path: Path to output FAA file
        kept_gene_ids: Set of gene IDs to keep
    """
    kept_records = []
    for record in SeqIO.parse(original_faa_path, "fasta"):
        gene_id = record.id.split()[0]
        if gene_id in kept_gene_ids:
            kept_records.append(record)

    SeqIO.write(kept_records, output_faa_path, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description='Filter overlapping genes from GFF/FAA files, keeping longest',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  filter_overlapping_sequences.py input_annotations/ filtered_annotations/

Input directory must contain:
  input_dir/gff/        - GFF files (*.gff)
  input_dir/proteins/   - FAA files (*.faa)

Output directory will be created with same structure.

Filtering logic:
  - Merges CDS exons to get full coding region
  - When genes overlap, keeps the longer one
  - Writes filtered GFF (only kept CDS features) and FAA (only kept sequences)
  - Logs all filtered genes to output_dir/filtered_genes_log.txt
        """
    )

    parser.add_argument('input_dir', help='Input directory with gff/ and proteins/ subdirectories')
    parser.add_argument('output_dir', help='Output directory for filtered files')

    args = parser.parse_args()

    # Validate input directory
    input_gff_dir = os.path.join(args.input_dir, 'gff')
    input_faa_dir = os.path.join(args.input_dir, 'proteins')

    if not os.path.isdir(input_gff_dir):
        print(f"ERROR: GFF directory not found: {input_gff_dir}")
        sys.exit(1)

    if not os.path.isdir(input_faa_dir):
        print(f"ERROR: Proteins directory not found: {input_faa_dir}")
        sys.exit(1)

    # Create output directories
    output_gff_dir = os.path.join(args.output_dir, 'gff')
    output_faa_dir = os.path.join(args.output_dir, 'proteins')

    os.makedirs(output_gff_dir, exist_ok=True)
    os.makedirs(output_faa_dir, exist_ok=True)

    print("\n" + "="*80)
    print("FILTER OVERLAPPING SEQUENCES")
    print("="*80)
    print(f"\nInput:  {args.input_dir}")
    print(f"Output: {args.output_dir}\n")

    # Process each species
    all_filtered = []
    species_stats = []

    gff_files = [f for f in os.listdir(input_gff_dir) if f.endswith('.gff')]

    if not gff_files:
        print(f"ERROR: No GFF files (*.gff) found in {input_gff_dir}")
        sys.exit(1)

    print(f"Found {len(gff_files)} species to process\n")

    for gff_file in sorted(gff_files):
        species = gff_file.replace('.gff', '')

        gff_path = os.path.join(input_gff_dir, gff_file)
        faa_file = species + '.faa'
        faa_path = os.path.join(input_faa_dir, faa_file)

        if not os.path.exists(faa_path):
            print(f"  WARNING: Skipping {species} - FAA file not found: {faa_file}")
            continue

        # Parse and filter
        parsed_genes, sequences, filtered = parse_gff_faa(gff_path, faa_path, species)

        # Get kept gene IDs
        kept_gene_ids = set()
        for seqid, genes in parsed_genes.items():
            for start, end, strand, gene_id in genes:
                kept_gene_ids.add(gene_id)

        # Write filtered files
        output_gff_path = os.path.join(output_gff_dir, gff_file)
        output_faa_path = os.path.join(output_faa_dir, faa_file)

        write_filtered_gff(gff_path, output_gff_path, kept_gene_ids)
        write_filtered_faa(faa_path, output_faa_path, kept_gene_ids)

        # Track statistics
        total_genes = len(kept_gene_ids) + len(filtered)
        species_stats.append((species, total_genes, len(kept_gene_ids), len(filtered)))

        # Add species to filtered genes
        for gene_id, reason, details in filtered:
            all_filtered.append((species, gene_id, reason, details))

    # Write filtered genes log
    log_path = os.path.join(args.output_dir, 'filtered_genes_log.txt')
    with open(log_path, 'w') as f:
        f.write("# Genes filtered due to overlaps\n")
        f.write(f"# Total filtered: {len(all_filtered)}\n\n")

        for species, gene_id, reason, details in all_filtered:
            f.write(f"Species: {species}\n")
            f.write(f"Gene: {gene_id}\n")
            f.write(f"  Chromosome: {details['chromosome']}\n")
            f.write(f"  Position: {details['position']}\n")
            f.write(f"  Length: {details['length']}\n")
            f.write(f"  Overlaps with: {details['overlaps_with']}\n")
            f.write(f"  Overlapping position: {details['overlapping_position']}\n")
            f.write(f"  Kept gene length: {details['kept_length']}\n")
            f.write(f"  Reason: {reason}\n\n")

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\n{'Species':<30} {'Input':<10} {'Kept':<10} {'Filtered':<10} {'%Kept':<10}")
    print("-" * 80)

    total_input = 0
    total_kept = 0
    total_filtered = 0

    for species, input_count, kept_count, filtered_count in species_stats:
        pct_kept = (kept_count / input_count * 100) if input_count > 0 else 0
        print(f"{species:<30} {input_count:<10} {kept_count:<10} {filtered_count:<10} {pct_kept:>6.1f}%")
        total_input += input_count
        total_kept += kept_count
        total_filtered += filtered_count

    print("-" * 80)
    total_pct = (total_kept / total_input * 100) if total_input > 0 else 0
    print(f"{'TOTAL':<30} {total_input:<10} {total_kept:<10} {total_filtered:<10} {total_pct:>6.1f}%")

    print(f"\nFiltered genes logged to: {log_path}")
    print(f"\nOutput written to:")
    print(f"  GFF: {output_gff_dir}/")
    print(f"  FAA: {output_faa_dir}/")
    print("\n" + "="*80 + "\n")


if __name__ == '__main__':
    main()
