#!/usr/bin/env python3
"""
blast_query_on_annotations.py

Find which annotated proteins match user-provided query proteins using blastp.

This script:
1. Runs blastp: query proteins (user-provided) -> annotated proteins (subject)
2. Optionally runs CLASP to chain multi-domain hits
3. Parses BLAST/CLASP output to extract hits with coordinates and scores
4. Maps hits back to gene IDs from GFF annotations

Usage:
    python3 blast_query_on_annotations.py \\
        --query_proteins query.faa \\
        --annotations_dir annotations_raw \\
        --output_dir blast_results \\
        --cores 8 \\
        --blast_evalue 1e-3 \\
        --use_clasp y

Output:
    - blast_results/hits_{org}.tsv : BLAST hits with gene coordinates

Next step:
    Use with annotation_hits_to_synthology.py to filter and format for synthology
"""

import os
import sys
import argparse
import multiprocessing as mp
from subprocess import run
from os.path import isfile


def blastp(query, subject, outfile, evalue):
    """
    Run blastp: query proteins vs subject proteins.

    Args:
        query: Query protein FASTA (user-provided proteins of interest)
        subject: Subject protein FASTA (annotated proteins)
        outfile: Output file path
        evalue: E-value threshold

    Returns:
        0 if successful
    """
    cmd = f'blastp -query {query} -subject {subject} ' \
          f'-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" ' \
          f'-evalue {evalue} -out {outfile}'
    run(cmd, shell=True)
    return 0


def clasp(infile, outfile, clasp_path, l, e):
    """
    Run CLASP to chain BLAST hits.

    Args:
        infile: BLAST output file
        outfile: CLASP output file
        clasp_path: Path to clasp.x executable
        l: CLASP -l parameter (gap length penalty)
        e: CLASP -e parameter (gap extension penalty)

    Returns:
        0 if successful
    """
    cmd = [clasp_path, '-m', '-f', '-O',
           '-i', infile,
           '-c', '7', '8', '9', '10', '12',  # qstart qend sstart send bitscore
           '-C', '1', '2',  # qseqid sseqid
           '-l', str(l),
           '-e', str(e),
           '-o', outfile]
    run(cmd)
    return 0


def parse_clasp_output(clasp_file, org):
    """
    Parse CLASP output to extract chained hits.

    Args:
        clasp_file: Path to CLASP output file
        org: Organism identifier

    Returns:
        List of hit dictionaries with keys: query_id, gene_id, bitscore, qstart, qend, sstart, send
    """
    hits = []

    if not isfile(clasp_file):
        raise FileNotFoundError(
            f"ERROR: CLASP output file not found: {clasp_file}\n"
            f"  CLASP was expected to produce output for {org}.\n"
            f"  This file should have been created by the CLASP chaining step."
        )

    with open(clasp_file) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        if line[0] == '#':
            # CLASP header/comment
            i += 1
            continue

        cols = line.strip().split('\t')

        # Only C and F lines expected in CLASP output body
        if cols[0] not in ['C', 'F']:
            raise ValueError(f"Unexpected CLASP line type '{cols[0]}', expected 'C' or 'F'. Line: {line.strip()}")

        if cols[0] == 'C':
            # C-line: chained hit
            # Format: C qseqid sseqid bitscore num_fragments qstart qend sstart send
            if len(cols) != 8:
                raise ValueError(f"CLASP C-line has {len(cols)} columns, expected 8. Line: {line.strip()}")

            query_id = cols[1]
            gene_id = cols[2]  # This is the gene ID from annotation (e.g., "g1.t1")
            bitscore = float(cols[3])
            # Note: We don't need qstart/qend/sstart/send for this pipeline
            # Script 3 only needs gene_id and bitscore, then looks up coords from GFF

            hits.append({
                'org': org,
                'query_id': query_id,
                'gene_id': gene_id,
                'bitscore': bitscore
            })

        i += 1

    return hits


def parse_blast_output(blast_file, org):
    """
    Parse raw BLAST output (no CLASP).

    Args:
        blast_file: Path to BLAST output file
        org: Organism identifier

    Returns:
        List of hit dictionaries
    """
    hits = []

    if not isfile(blast_file):
        raise FileNotFoundError(
            f"ERROR: BLAST output file not found: {blast_file}\n"
            f"  blastp was expected to produce output for {org}.\n"
            f"  Expected format: tabular -outfmt 6 with 12 columns (qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore)"
        )

    with open(blast_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line[0] == '#':
                continue

            cols = line.strip().split('\t')
            if len(cols) != 12:
                raise ValueError(
                    f"ERROR: Malformed BLAST output line in {blast_file}: {len(cols)} columns, expected 12.\n"
                    f"  Expected format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\n"
                    f"  Line: {line.strip()}"
                )

            # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
            query_id = cols[0]
            gene_id = cols[1]  # Gene ID from annotation
            bitscore = float(cols[11])
            # Note: We don't need qstart/qend/sstart/send for this pipeline
            # Script 3 only needs gene_id and bitscore, then looks up coords from GFF

            hits.append({
                'org': org,
                'query_id': query_id,
                'gene_id': gene_id,
                'bitscore': bitscore
            })

    return hits


def exec_blast_clasp(org, query, annotations_dir, output_dir, blast_evalue, use_clasp, clasp_path, clasp_l, clasp_e):
    """
    Execute blastp (and optionally CLASP) for one organism.

    Args:
        org: Organism identifier
        query: Query protein FASTA file
        annotations_dir: Directory with annotation results
        output_dir: Output directory for BLAST results
        blast_evalue: BLAST e-value threshold
        use_clasp: Whether to use CLASP chaining
        clasp_path: Path to clasp.x
        clasp_l: CLASP -l parameter
        clasp_e: CLASP -e parameter

    Returns:
        Number of hits found
    """
    # Input: annotated proteins from Script 1
    subject_faa = f'{annotations_dir}/proteins_raw/{org}_predictions.faa'

    if not isfile(subject_faa):
        raise FileNotFoundError(
            f"ERROR: Annotated protein file not found: {subject_faa}\n"
            f"  Expected output from annotate_syn_regions.py for {org}.\n"
            f"  The file was present when the organism list was built but is now missing."
        )

    # Run blastp
    blast_out = f'{output_dir}/blast_out/{org}'
    blastp(query, subject_faa, blast_out, blast_evalue)

    # Parse results
    if use_clasp:
        # Run CLASP
        if isfile(blast_out) and os.path.getsize(blast_out) > 0:
            clasp_out = f'{output_dir}/clasp_out/{org}'
            clasp(blast_out, clasp_out, clasp_path, clasp_l, clasp_e)

            # Parse CLASP output
            hits = parse_clasp_output(clasp_out, org)
        else:
            hits = []
    else:
        # Parse raw BLAST output
        hits = parse_blast_output(blast_out, org)

    # Write hits to file
    if len(hits) > 0:
        with open(f'{output_dir}/hits_{org}.tsv', 'w') as f:
            # Header
            f.write('org\tquery_id\tgene_id\tbitscore\n')
            for hit in hits:
                f.write(f"{hit['org']}\t{hit['query_id']}\t{hit['gene_id']}\t{hit['bitscore']}\n")

    return len(hits)


def main():
    parser = argparse.ArgumentParser(
        description='BLAST query proteins against annotated proteins',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--query_proteins', type=str, required=True,
                       help="Query protein FASTA file (proteins of interest)")
    parser.add_argument('--annotations_dir', type=str, default='annotations_raw',
                       help="Directory with annotation results (from annotate_syn_regions.py)")
    parser.add_argument('--output_dir', type=str, default='blast_results',
                       help="Output directory for BLAST results")
    parser.add_argument('--cores', type=int, default=1,
                       help="Number of parallel processes")
    parser.add_argument('--blast_evalue', type=float, default=1e-3,
                       help="BLAST e-value threshold")
    parser.add_argument('--use_clasp', type=str, default='y', choices=['y', 'n'],
                       help="Use CLASP for chaining hits (y or n)")
    parser.add_argument('--clasp_path', type=str, default='../clasp.x',
                       help="Path to clasp.x executable")
    parser.add_argument('--clasp_l', type=float, default=0.5,
                       help="CLASP -l parameter (gap length penalty)")
    parser.add_argument('--clasp_e', type=float, default=0,
                       help="CLASP -e parameter (gap extension penalty)")

    args = parser.parse_args()

    # Validate inputs
    if not isfile(args.query_proteins):
        print(f"ERROR: Query protein file not found: {args.query_proteins}")
        sys.exit(1)

    if not os.path.isdir(args.annotations_dir):
        print(f"ERROR: Annotations directory not found: {args.annotations_dir}")
        sys.exit(1)

    if not os.path.isdir(f'{args.annotations_dir}/proteins_raw'):
        print(f"ERROR: No proteins_raw subdirectory in {args.annotations_dir}")
        print(f"Expected output from annotate_syn_regions.py")
        sys.exit(1)

    use_clasp = args.use_clasp == 'y'

    # Create output directories
    os.makedirs(f'{args.output_dir}/blast_out', exist_ok=True)
    if use_clasp:
        os.makedirs(f'{args.output_dir}/clasp_out', exist_ok=True)

    # Get list of organisms
    orgs = []
    if isfile('orgs'):
        with open('orgs') as f:
            for line in f:
                org = line.strip()
                if org and isfile(f'{args.annotations_dir}/proteins_raw/{org}_predictions.faa'):
                    orgs.append(org)
    else:
        # Fallback: list all _predictions.faa files
        protein_dir = f'{args.annotations_dir}/proteins_raw'
        if os.path.isdir(protein_dir):
            for filename in os.listdir(protein_dir):
                if filename.endswith('_predictions.faa'):
                    org = filename.replace('_predictions.faa', '')
                    orgs.append(org)

    if len(orgs) == 0:
        print("ERROR: No organisms found with annotation files")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"BLAST Query Proteins Against Annotations")
    print(f"{'='*60}")
    print(f"Query proteins:    {args.query_proteins}")
    print(f"Annotations dir:   {args.annotations_dir}")
    print(f"Output directory:  {args.output_dir}")
    print(f"Organisms:         {len(orgs)}")
    print(f"Cores:             {args.cores}")
    print(f"BLAST e-value:     {args.blast_evalue}")
    print(f"Use CLASP:         {'yes' if use_clasp else 'no'}")
    print(f"{'='*60}\n")

    # Process organisms
    total_hits = 0

    if args.cores > 1:
        with mp.Pool(processes=args.cores) as pool:
            results = pool.starmap(
                exec_blast_clasp,
                [(org, args.query_proteins, args.annotations_dir, args.output_dir,
                  args.blast_evalue, use_clasp, args.clasp_path, args.clasp_l, args.clasp_e)
                 for org in orgs]
            )
        for org, n in zip(orgs, results):
            if n > 0:
                print(f"{org}: {n} hits found")
            total_hits += n
    else:
        for org in orgs:
            n = exec_blast_clasp(
                org, args.query_proteins, args.annotations_dir, args.output_dir,
                args.blast_evalue, use_clasp, args.clasp_path, args.clasp_l, args.clasp_e
            )
            if n > 0:
                print(f"{org}: {n} hits found")
            total_hits += n

    print(f"\n{'='*60}")
    print(f"BLAST Search Complete")
    print(f"{'='*60}")
    print(f"Total hits:        {total_hits}")
    print(f"Hit files:         {args.output_dir}/hits_*.tsv")
    print(f"\nNext step:")
    print(f"  python3 annotation_hits_to_synthology.py \\")
    print(f"      --blast_hits {args.output_dir} \\")
    print(f"      --annotations_dir {args.annotations_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
