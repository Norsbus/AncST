#!/usr/bin/env python3
"""
annotate_syn_regions.py

Run ab initio gene prediction on extracted syntenic regions using Augustus or Helixer.

This script:
1. Runs gene prediction on genomic regions (DNA sequences)
2. Parses GFF3 output to extract CDS features
3. Groups CDS by parent gene/mRNA (handles multi-exon genes)
4. Concatenates ALL CDS exons (no introns), translates to protein
5. Selects longest protein isoform per gene

Usage:
    python3 annotate_syn_regions.py \\
        --predictor augustus \\
        --fastas_dir fastas \\
        --output_dir annotations_raw \\
        --cores 8 \\
        --augustus_species fly

    python3 annotate_syn_regions.py \\
        --predictor helixer \\
        --fastas_dir fastas \\
        --output_dir annotations_raw \\
        --cores 8 \\
        --helixer_model invertebrate

Output:
    - annotations_raw/gff_raw/{org}_predictions.gff : Raw predictor GFF3
    - annotations_raw/proteins_raw/{org}_predictions.faa : Longest protein per gene

Next step:
    Use with blast_query_on_annotations.py to find which genes match query proteins
"""

import os
import sys
import argparse
import multiprocessing as mp
from collections import defaultdict
from subprocess import run, PIPE
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Helixer uses S50 (50-byte fixed ASCII) for seqids in its HDF5 files.
# Our FASTA seqids can exceed 50 chars, causing truncation and KeyError downstream.
HELIXER_SEQID_LIMIT = 50


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


def run_augustus(fasta_file, species, output_gff):
    """
    Run Augustus gene prediction on FASTA file.

    Args:
        fasta_file: Input FASTA (genomic DNA)
        species: Augustus species model (e.g., 'fly', 'human', 'arabidopsis')
        output_gff: Output GFF3 file path

    Returns:
        0 if successful, non-zero otherwise
    """
    cmd = [
        'augustus',
        '--species=' + species,
        '--gff3=on',
        '--strand=both',
        '--noInFrameStop=true',
        fasta_file
    ]

    try:
        with open(output_gff, 'w') as out:
            result = run(cmd, stdout=out, stderr=PIPE, check=True)
        return 0
    except Exception as e:
        print(f"ERROR running Augustus: {e}", file=sys.stderr)
        return 1


def run_helixer(fasta_file, model, output_gff):
    """
    Run Helixer gene prediction on FASTA file.

    Args:
        fasta_file: Input FASTA (genomic DNA)
        model: Helixer model (fungi, land_plant, vertebrate, invertebrate, auto)
        output_gff: Output GFF3 file path

    Returns:
        0 if successful, non-zero otherwise
    """
    cmd = [
        'Helixer.py',
        '--lineage', model,
        '--fasta-path', fasta_file,
        '--gff-output-path', output_gff
    ]

    try:
        result = run(cmd, check=True, capture_output=True)
        return 0
    except Exception as e:
        print(f"ERROR running Helixer: {e}", file=sys.stderr)
        return 1


def parse_gene_predictions_gff(gff_file):
    """
    Parse gene prediction GFF3 output and group CDS by parent gene/mRNA.

    Works with both Augustus and Helixer output formats.

    Args:
        gff_file: Path to prediction GFF3 file

    Returns:
        Dictionary mapping gene ID to list of CDS features
        Each CDS feature is a dict with: seqid, start, end, strand, phase
    """
    genes = defaultdict(lambda: {
        'gene': None,
        'mRNAs': defaultdict(list)  # mRNA_id -> list of CDS
    })

    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            cols = line.split('\t')
            if len(cols) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = cols
            start, end = int(start), int(end)

            attrs = parse_gff_attributes(attributes)

            if feature_type == 'gene':
                gene_id = attrs.get('ID', '')
                if gene_id:
                    genes[gene_id]['gene'] = {
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand
                    }

            elif feature_type in ['mRNA', 'transcript']:
                # mRNA/transcript - link to parent gene
                gene_id = attrs.get('Parent', '')
                mrna_id = attrs.get('ID', '')
                if gene_id and mrna_id:
                    # Store mRNA metadata (currently unused but available)
                    pass

            elif feature_type == 'CDS':
                # CDS exon - link to parent mRNA or gene
                parent_id = attrs.get('Parent', '')
                if parent_id:
                    # Determine gene ID (parent might be mRNA or gene)
                    # For Augustus: Parent is usually transcript ID
                    # For Helixer: Parent might be gene ID directly
                    # We'll use Parent as the grouping key
                    genes[parent_id]['mRNAs'][parent_id].append({
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'phase': phase
                    })

    # Flatten: one representative transcript per gene
    # Strategy: for each gene, pick the transcript with most CDS (or longest total length)
    result = {}

    for gene_id, gene_data in genes.items():
        # If this is a transcript ID (has CDS), use it directly
        if len(gene_data['mRNAs'].get(gene_id, [])) > 0:
            result[gene_id] = gene_data['mRNAs'][gene_id]
        # Otherwise, pick best transcript
        elif len(gene_data['mRNAs']) > 0:
            # Pick transcript with most exons (or longest total)
            best_mrna_id = max(
                gene_data['mRNAs'].keys(),
                key=lambda m: sum(cds['end'] - cds['start'] for cds in gene_data['mRNAs'][m])
            )
            result[gene_id] = gene_data['mRNAs'][best_mrna_id]

    return result


def extract_and_translate_cds(genome, cds_list):
    """
    Extract CDS sequences from genome and translate to protein.

    Concatenates ALL exons (no introns), handles strand orientation,
    applies GFF3 phase, trims trailing partial codons, translates.

    Args:
        genome: Dictionary of chromosome sequences (seqid -> Seq)
        cds_list: List of CDS feature dictionaries

    Returns:
        Tuple of (protein_seq, trim_info) where:
        - protein_seq: Protein Seq object, or None if extraction/translation fails
        - trim_info: dict with trim diagnostics if any trimming occurred, else None
          Keys: 'phase', 'trailing', 'seqid', 'min_start', 'max_end', 'seq_len',
                'dist_from_start', 'dist_from_end'
    """
    if not cds_list:
        return None, None

    # All CDS should be on same chromosome and strand
    seqid = cds_list[0]['seqid']
    strand = cds_list[0]['strand']

    if seqid not in genome:
        raise ValueError(
            f"ERROR: seqid '{seqid}' from prediction GFF not found in input FASTA.\n"
            f"  FASTA seqids: {list(genome.keys())[:10]}{'...' if len(genome) > 10 else ''}\n"
            f"  The predictor was run on this FASTA, so all seqids must match."
        )

    # Sort CDS by genomic position
    sorted_cds = sorted(cds_list, key=lambda x: x['start'])

    # Extract exon sequences and concatenate
    exon_seqs = []
    for cds in sorted_cds:
        # GFF3 uses 1-based inclusive coordinates
        # Python uses 0-based exclusive: genome[start-1:end]
        exon_start = cds['start'] - 1
        exon_end = cds['end']
        exon_seq = genome[seqid][exon_start:exon_end]
        exon_seqs.append(exon_seq)

    # Concatenate exons (spliced, no introns)
    concatenated = Seq(''.join(str(s) for s in exon_seqs))

    if len(concatenated) == 0:
        print(f"WARNING: 0-length CDS sequence after exon concatenation for {seqid}. "
              f"CDS count: {len(cds_list)}, positions: {[(c['start'],c['end']) for c in cds_list]}")
        return None, None

    # Handle strand orientation BEFORE translation
    if strand == '-':
        concatenated = concatenated.reverse_complement()
        # For minus strand, the first exon in coding order is the last in genomic order
        first_coding_cds = sorted_cds[-1]
    else:
        first_coding_cds = sorted_cds[0]

    # Apply GFF3 phase: skip N bases at start of first coding exon to reach
    # the first complete codon. Phase 0 = already on boundary, 1 = skip 1, 2 = skip 2.
    # This handles partial genes at region boundaries and exon-boundary frame offsets.
    phase_str = first_coding_cds.get('phase', '0')
    phase = int(phase_str) if phase_str in ('0', '1', '2') else 0
    if phase > 0:
        concatenated = concatenated[phase:]

    # Trim trailing partial codon (can happen with truncated genes at region edges)
    remainder = len(concatenated) % 3
    if remainder != 0:
        concatenated = concatenated[:len(concatenated) - remainder]

    # Build trim diagnostics when any trimming occurred.
    # Boundary determination is done by the caller (process_organism) which can see
    # whether this gene is the first/last predicted gene on its sequence.
    trim_info = None
    if phase > 0 or remainder > 0:
        seq_len = len(genome[seqid])
        min_start = min(cds['start'] for cds in sorted_cds)
        max_end = max(cds['end'] for cds in sorted_cds)
        trim_info = {
            'phase': phase,
            'trailing': remainder,
            'seqid': seqid,
            'min_start': min_start,
            'max_end': max_end,
            'seq_len': seq_len,
            'dist_from_start': min_start - 1,        # bp from sequence start
            'dist_from_end': seq_len - max_end,       # bp from sequence end
        }

    if len(concatenated) == 0:
        print(f"WARNING: 0-length CDS after phase/trim for {seqid}. "
              f"Original length before trim: {len(concatenated) + remainder + phase}")
        return None, trim_info

    # Translate to protein
    try:
        protein = concatenated.translate(to_stop=False)
        return protein, trim_info
    except Exception as e:
        print(f"WARNING: Translation failed for {seqid} "
              f"(seq length: {len(concatenated)}, strand: {strand}): {e}")
        return None, trim_info


def load_genome(fasta_file):
    """
    Load genome sequences into dictionary.

    Args:
        fasta_file: Path to genome FASTA file

    Returns:
        Dictionary mapping chromosome ID to sequence
    """
    if not os.path.isfile(fasta_file):
        return {}

    genome = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        genome[record.id] = record.seq

    return genome


def get_augustus_species(org, species_arg):
    """
    Determine Augustus species model for organism.

    Args:
        org: Organism identifier
        species_arg: Either a species name (same for all) or path to mapping file

    Returns:
        Species model name
    """
    # Check if species_arg is a file (organism -> species mapping)
    if os.path.isfile(species_arg):
        with open(species_arg) as f:
            for line in f:
                if line.strip():
                    cols = line.strip().split('\t')
                    if len(cols) >= 2 and cols[0] == org:
                        return cols[1]
        # Not found in mapping file, use default
        return 'generic'
    else:
        # Use same species for all organisms
        return species_arg


def shorten_fasta_for_helixer(fasta_file, output_dir, org):
    """
    Create a temporary FASTA with shortened seqids if any exceed Helixer's S50 limit.

    Args:
        fasta_file: Path to original FASTA file
        output_dir: Directory to write temporary FASTA into
        org: Organism name (used to make temp file unique for multiprocessing)

    Returns:
        (fasta_path, mapping) where:
        - If no seqids exceed 50 chars: (original fasta_file, None) — no-op
        - If shortening needed: (tmp_fasta_path, {short_id: original_id})
    """
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    needs_shortening = any(len(rec.id) > HELIXER_SEQID_LIMIT for rec in records)

    if not needs_shortening:
        return fasta_file, None

    long_ids = [rec.id for rec in records if len(rec.id) > HELIXER_SEQID_LIMIT]
    print(f"  {org}: Helixer S50 workaround: {len(long_ids)}/{len(records)} seqids exceed "
          f"{HELIXER_SEQID_LIMIT} chars, using temporary short IDs")

    mapping = {}  # short_id -> original_id
    tmp_fasta = os.path.join(output_dir, f'_helixer_tmp_{org}.fasta')

    shortened_records = []
    for i, rec in enumerate(records):
        short_id = f's{i}'
        mapping[short_id] = rec.id
        shortened_records.append(SeqRecord(rec.seq, id=short_id, description=''))

    SeqIO.write(shortened_records, tmp_fasta, 'fasta')
    return tmp_fasta, mapping


def restore_gff_seqids(gff_file, mapping):
    """
    Replace shortened seqids in GFF column 1 back to original seqids.

    Args:
        gff_file: Path to GFF file (modified in-place)
        mapping: Dict of {short_id: original_id}
    """
    reverse_map = {short: orig for short, orig in mapping.items()}

    with open(gff_file, 'r') as f:
        lines = f.readlines()

    with open(gff_file, 'w') as f:
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                f.write(line)
                continue
            cols = line.split('\t')
            if len(cols) >= 1 and cols[0] in reverse_map:
                cols[0] = reverse_map[cols[0]]
                f.write('\t'.join(cols))
            else:
                f.write(line)


def process_organism(org, fastas_dir, output_dir, predictor, augustus_species, helixer_model):
    """
    Process one organism: run gene prediction, extract proteins.

    Args:
        org: Organism identifier
        fastas_dir: Directory with extracted FASTA files
        output_dir: Output directory
        predictor: 'augustus' or 'helixer'
        augustus_species: Species model for Augustus
        helixer_model: Model for Helixer

    Returns:
        Number of proteins extracted
    """
    fasta_file = f'{fastas_dir}/{org}/all.fasta'

    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(
            f"ERROR: FASTA file not found: {fasta_file}\n"
            f"  {org} was in the organism list but its FASTA is missing."
        )

    # Create output directories
    os.makedirs(f'{output_dir}/gff_raw', exist_ok=True)
    os.makedirs(f'{output_dir}/proteins_raw', exist_ok=True)

    output_gff = f'{output_dir}/gff_raw/{org}_predictions.gff'
    output_faa = f'{output_dir}/proteins_raw/{org}_predictions.faa'

    # Run gene prediction
    if predictor == 'augustus':
        species = get_augustus_species(org, augustus_species)
        ret = run_augustus(fasta_file, species, output_gff)
        seqid_mapping = None
    else:  # helixer
        # Helixer uses S50 (50-byte) seqids internally — shorten if needed
        helixer_fasta, seqid_mapping = shorten_fasta_for_helixer(fasta_file, output_dir, org)
        ret = run_helixer(helixer_fasta, helixer_model, output_gff)
        # Restore original seqids in GFF and clean up temp FASTA
        if seqid_mapping and ret == 0:
            restore_gff_seqids(output_gff, seqid_mapping)
        if seqid_mapping and helixer_fasta != fasta_file:
            try:
                os.remove(helixer_fasta)
            except OSError:
                pass

    if ret != 0:
        print(f"WARNING: {org} - gene prediction failed")
        return 0

    # Parse GFF3 to get genes with CDS
    genes = parse_gene_predictions_gff(output_gff)

    if len(genes) == 0:
        print(f"WARNING: {org}: {predictor} produced GFF with 0 genes containing CDS features.\n"
              f"  GFF file: {output_gff}\n"
              f"  Expected: GFF3 with gene/mRNA/CDS features from {predictor}.")
        return 0

    # Load genome for sequence extraction
    genome = load_genome(fasta_file)
    if len(genome) == 0:
        raise ValueError(
            f"ERROR: {org}: FASTA file is empty or unparseable: {fasta_file}\n"
            f"  This FASTA was just used as input for {predictor}, so it must be valid."
        )

    # Identify the first and last predicted gene per sequence (seqid).
    # These are the genes most likely to be truncated at region boundaries.
    gene_extents = {}  # gene_id -> (seqid, min_start, max_end)
    for gene_id, cds_list in genes.items():
        if cds_list:
            sid = cds_list[0]['seqid']
            gmin = min(c['start'] for c in cds_list)
            gmax = max(c['end'] for c in cds_list)
            gene_extents[gene_id] = (sid, gmin, gmax)

    first_gene_per_seq = {}  # seqid -> gene_id with smallest start
    last_gene_per_seq = {}   # seqid -> gene_id with largest end
    for gene_id, (sid, gmin, gmax) in gene_extents.items():
        if sid not in first_gene_per_seq or gmin < gene_extents[first_gene_per_seq[sid]][1]:
            first_gene_per_seq[sid] = gene_id
        if sid not in last_gene_per_seq or gmax > gene_extents[last_gene_per_seq[sid]][2]:
            last_gene_per_seq[sid] = gene_id
    boundary_gene_ids = set(first_gene_per_seq.values()) | set(last_gene_per_seq.values())

    # Extract and translate proteins
    proteins = []
    n_translation_failed = 0
    n_boundary_trimmed = 0
    n_interior_trimmed = 0
    interior_trim_examples = []
    for gene_id, cds_list in genes.items():
        protein_seq, trim_info = extract_and_translate_cds(genome, cds_list)

        # Track trimming diagnostics
        if trim_info is not None:
            is_boundary = gene_id in boundary_gene_ids
            if is_boundary:
                n_boundary_trimmed += 1
            else:
                n_interior_trimmed += 1
                if len(interior_trim_examples) < 5:
                    interior_trim_examples.append(
                        f"    {gene_id}: phase={trim_info['phase']}, "
                        f"trailing={trim_info['trailing']}bp trimmed, "
                        f"CDS span={trim_info['min_start']}-{trim_info['max_end']}, "
                        f"seq_len={trim_info['seq_len']}, "
                        f"dist_from_edges={trim_info['dist_from_start']}bp/"
                        f"{trim_info['dist_from_end']}bp"
                    )

        if protein_seq is None or len(protein_seq) == 0:
            n_translation_failed += 1
            continue

        # Create SeqRecord
        record = SeqRecord(
            protein_seq,
            id=gene_id,
            description=f'organism={org} predictor={predictor}'
        )
        proteins.append(record)

    # Report trimming diagnostics per organism
    if n_boundary_trimmed > 0:
        print(f"  {org}: {n_boundary_trimmed} genes trimmed (first/last gene on their sequence, "
              f"expected for partial genes at region boundaries)")
    if n_interior_trimmed > 0:
        print(f"  *** WARNING *** {org}: {n_interior_trimmed} genes with non-multiple-of-3 CDS "
              f"that are NOT the first/last gene on their sequence!")
        print(f"  These are interior genes — partial codons here may indicate a GFF parsing "
              f"bug or coordinate error. Examples:")
        for example in interior_trim_examples:
            print(example)
        if n_interior_trimmed > 5:
            print(f"    ... and {n_interior_trimmed - 5} more")

    if len(proteins) == 0:
        print(f"WARNING: {org}: 0 proteins from {len(genes)} predicted genes "
              f"(translation_failed: {n_translation_failed}).\n"
              f"  All translations failed. Check if genome sequences match GFF coordinates.")
        return 0

    if n_translation_failed > 0:
        print(f"  {org}: {len(proteins)} proteins translated, {n_translation_failed} failed")

    # Write protein FASTA
    SeqIO.write(proteins, output_faa, 'fasta')

    return len(proteins)


def main():
    parser = argparse.ArgumentParser(
        description='Run ab initio gene prediction on syntenic regions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--predictor', type=str, required=True,
                       choices=['augustus', 'helixer'],
                       help="Gene prediction tool (augustus or helixer)")
    parser.add_argument('--fastas_dir', type=str, default='fastas',
                       help="Directory with extracted FASTA files ({org}/all.fasta)")
    parser.add_argument('--output_dir', type=str, default='annotations_raw',
                       help="Output directory for predictions")
    parser.add_argument('--cores', type=int, default=1,
                       help="Number of parallel processes")
    parser.add_argument('--augustus_species', type=str, default='generic',
                       help="Augustus species model or path to mapping file (org<tab>species)")
    parser.add_argument('--helixer_model', type=str, default='invertebrate',
                       choices=['fungi', 'land_plant', 'vertebrate', 'invertebrate', 'auto'],
                       help="Helixer model")

    args = parser.parse_args()

    # Validate inputs
    if not os.path.isdir(args.fastas_dir):
        print(f"ERROR: FASTA directory not found: {args.fastas_dir}")
        sys.exit(1)

    # Get list of organisms
    orgs = []
    if os.path.isfile('orgs'):
        with open('orgs') as f:
            for line in f:
                org = line.strip()
                if org and os.path.isfile(f'{args.fastas_dir}/{org}/all.fasta'):
                    orgs.append(org)
    else:
        # Fallback: list all subdirectories with all.fasta
        if os.path.isdir(args.fastas_dir):
            for org in os.listdir(args.fastas_dir):
                if os.path.isdir(f'{args.fastas_dir}/{org}') and os.path.isfile(f'{args.fastas_dir}/{org}/all.fasta'):
                    orgs.append(org)

    if len(orgs) == 0:
        print("ERROR: No organisms found with FASTA files")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Ab Initio Gene Prediction on Syntenic Regions")
    print(f"{'='*60}")
    print(f"Predictor:         {args.predictor}")
    print(f"FASTA directory:   {args.fastas_dir}")
    print(f"Output directory:  {args.output_dir}")
    print(f"Organisms:         {len(orgs)}")
    print(f"Cores:             {args.cores}")
    if args.predictor == 'augustus':
        print(f"Augustus species:  {args.augustus_species}")
    else:
        print(f"Helixer model:     {args.helixer_model}")
    print(f"{'='*60}\n")

    # Process organisms
    total_proteins = 0

    if args.cores > 1:
        with mp.Pool(processes=args.cores) as pool:
            results = pool.starmap(
                process_organism,
                [(org, args.fastas_dir, args.output_dir, args.predictor,
                  args.augustus_species, args.helixer_model) for org in orgs]
            )
        for org, n in zip(orgs, results):
            if n > 0:
                print(f"{org}: {n} proteins extracted")
            total_proteins += n
    else:
        for org in orgs:
            n = process_organism(
                org, args.fastas_dir, args.output_dir, args.predictor,
                args.augustus_species, args.helixer_model
            )
            if n > 0:
                print(f"{org}: {n} proteins extracted")
            total_proteins += n

    if total_proteins == 0:
        print(f"\nERROR: 0 total proteins across all {len(orgs)} organisms. "
              f"Gene prediction produced no usable output.")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Gene Prediction Complete")
    print(f"{'='*60}")
    print(f"Total proteins:    {total_proteins}")
    print(f"GFF files:         {args.output_dir}/gff_raw/")
    print(f"Protein files:     {args.output_dir}/proteins_raw/")
    print(f"\nNext step:")
    print(f"  python3 blast_query_on_annotations.py \\")
    print(f"      --query_proteins <query.faa> \\")
    print(f"      --annotations_dir {args.output_dir}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
