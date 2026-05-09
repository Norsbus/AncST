#!/usr/bin/env python3
"""
coords_to_synthology.py

Convert a TSV of genomic coordinates to GFF3 + FASTA files suitable for
downstream/run_synthology_nucl.py and downstream/run_synthology_prot.py.

Input TSV format (tab-separated, with header row):
    Genome AC       Chr/Contig AC      Strand   Start    end
    GCA_000166975.1 AABZ01000048.1     minus    19809    18596

Start/end may be in descending order for minus-strand rows; the script
normalises them so that gff_start <= gff_end as required by GFF3.

Output directory layout:
    {output_dir}/gff/{Genome AC}.gff           (always)
    {output_dir}/proteins/{Genome AC}.faa      (--mode prot)
    {output_dir}/seqs/{Genome AC}.fna          (--mode nucl)

Gene IDs use the scheme: {Chr_AC}_{gff_start}_{gff_end}_{strand_char}
    e.g.: AABZ01000048.1_18596_19809_m

For prot mode, FASTA headers equal the GFF ID= value (matched directly by
run_synthology_prot.py).  For nucl mode, FASTA headers are prefixed with
the genome accession ({genome_ac}_{gene_id}) because run_synthology_nucl.py
prepends the species name when looking up records.

Usage:
    python utils/util_code/coords_to_synthology.py tsv_file --mode prot
    python utils/util_code/coords_to_synthology.py tsv_file --mode nucl \\
        --genomes-dir utils/genomes --output-dir annotations_for_synthology

Options:
    tsv_file              Input TSV (header + 5-col rows)
    --mode {nucl,prot}    Write nucleotide (.fna) or translated protein (.faa)
    --genomes-dir DIR     Path to genome FASTA files (default: <script>/../../genomes)
    --output-dir DIR      Output root directory (default: annotations_for_synthology)
    --genetic-code INT    NCBI translation table (default: 1, standard code)
    --source STR          GFF3 source column value (default: .)
"""

import argparse
import sys
from collections import defaultdict, namedtuple
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_GENOMES_DIR = str(SCRIPT_DIR / '..' / 'genomes')
DEFAULT_OUTPUT_DIR = 'annotations_for_synthology'

Row = namedtuple('Row', ['genome_ac', 'chr_ac', 'strand', 'raw_start', 'raw_end'])


def parse_tsv(tsv_file):
    """Parse input TSV, return list of Row namedtuples (header skipped)."""
    rows = []
    with open(tsv_file, 'r') as f:
        f.readline()  # skip header
        for lineno, line in enumerate(f, start=2):
            line = line.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 5:
                print(f"WARNING: line {lineno}: expected 5 columns, got {len(parts)} — skipping")
                continue
            genome_ac = parts[0].strip()
            chr_ac = parts[1].strip()
            strand = parts[2].strip()
            try:
                raw_start = int(parts[3].strip())
                raw_end = int(parts[4].strip())
            except ValueError:
                print(f"WARNING: line {lineno}: invalid coordinates "
                      f"'{parts[3].strip()}'/'{parts[4].strip()}' — skipping")
                continue
            if strand not in ('plus', 'minus'):
                print(f"WARNING: line {lineno}: unrecognised strand '{strand}' "
                      f"(expected 'plus' or 'minus') — skipping")
                continue
            rows.append(Row(genome_ac, chr_ac, strand, raw_start, raw_end))
    return rows


def group_by_genome(rows):
    """Group rows by genome accession, preserving insertion order."""
    groups = defaultdict(list)
    for row in rows:
        groups[row.genome_ac].append(row)
    return groups


def make_gene_id(chr_ac, gff_start, gff_end, strand):
    """Build stable unique gene ID: {chr_ac}_{gff_start}_{gff_end}_{strand_char}."""
    strand_char = 'm' if strand == 'minus' else 'p'
    return f"{chr_ac}_{gff_start}_{gff_end}_{strand_char}"


def process_genomes(groups, genomes_dir, output_dir, mode, genetic_code, source):
    """
    For each genome group: load FASTA, extract + orient sequences, write GFF + FASTA.

    Returns (n_genomes_written, n_genes_written, n_rows_skipped).
    """
    gff_dir = output_dir / 'gff'
    if mode == 'prot':
        seq_dir = output_dir / 'proteins'
        seq_ext = '.faa'
    else:
        seq_dir = output_dir / 'seqs'
        seq_ext = '.fna'

    gff_dir.mkdir(parents=True, exist_ok=True)
    seq_dir.mkdir(parents=True, exist_ok=True)

    n_genomes = 0
    n_genes = 0
    n_skipped = 0

    for genome_ac, rows in groups.items():
        fasta_path = Path(genomes_dir) / f"{genome_ac}.fasta"
        if not fasta_path.exists():
            print(f"WARNING: genome FASTA not found: {fasta_path} "
                  f"— skipping {len(rows)} row(s) for {genome_ac}")
            n_skipped += len(rows)
            continue

        # Load one genome at a time to keep memory manageable
        genome_dict = SeqIO.to_dict(SeqIO.parse(str(fasta_path), 'fasta'))

        gff_lines = []
        seq_records = []

        for row in rows:
            if row.chr_ac not in genome_dict:
                print(f"WARNING: contig '{row.chr_ac}' not in {fasta_path.name} — skipping")
                n_skipped += 1
                continue

            # GFF3 requires start <= end, both 1-based
            gff_start = min(row.raw_start, row.raw_end)
            gff_end = max(row.raw_start, row.raw_end)
            strand_gff = '+' if row.strand == 'plus' else '-'
            gene_id = make_gene_id(row.chr_ac, gff_start, gff_end, row.strand)

            # Extract sequence — BioPython uses 0-based half-open intervals
            raw_seq = genome_dict[row.chr_ac].seq[gff_start - 1 : gff_end]
            if row.strand == 'minus':
                raw_seq = raw_seq.reverse_complement()

            # GFF3 line: exactly 9 tab-separated columns
            gff_line = '\t'.join([
                row.chr_ac,       # col 1: seqid
                source,           # col 2: source
                'CDS',            # col 3: type
                str(gff_start),   # col 4: start (1-based)
                str(gff_end),     # col 5: end (1-based)
                '.',              # col 6: score
                strand_gff,       # col 7: strand
                '0',              # col 8: phase
                f"ID={gene_id}",  # col 9: attributes
            ])
            gff_lines.append(gff_line)

            if mode == 'prot':
                # Translate; warn on partial codons and internal stops
                if len(raw_seq) % 3 != 0:
                    print(f"WARNING: {gene_id}: length {len(raw_seq)} not divisible by 3 "
                          f"(partial codon — translating anyway)")
                protein = raw_seq.translate(table=genetic_code, to_stop=False)
                prot_str = str(protein)
                if len(prot_str) > 1 and '*' in prot_str[:-1]:
                    print(f"WARNING: {gene_id}: internal stop codon(s) present")
                # prot mode: FASTA header = bare gene_id (matched directly by prot parser)
                seq_records.append(SeqRecord(Seq(prot_str), id=gene_id, description=''))
            else:
                # nucl mode: FASTA header = {genome_ac}_{gene_id}
                # run_synthology_nucl.py prepends species when matching, so header must
                # equal f"{species}_{ID_value_from_GFF}"
                fasta_header = f"{genome_ac}_{gene_id}"
                seq_records.append(SeqRecord(raw_seq, id=fasta_header, description=''))

            n_genes += 1

        if not seq_records:
            print(f"WARNING: no genes written for {genome_ac} — skipping output files")
            continue

        # Write GFF (with version pragma as comment — parsers skip # lines)
        gff_path = gff_dir / f"{genome_ac}.gff"
        with open(str(gff_path), 'w') as fh:
            fh.write("##gff-version 3\n")
            for line in gff_lines:
                fh.write(line + '\n')

        # Write sequence file
        seq_path = seq_dir / f"{genome_ac}{seq_ext}"
        SeqIO.write(seq_records, str(seq_path), 'fasta')

        n_genomes += 1

    return n_genomes, n_genes, n_skipped


def main():
    parser = argparse.ArgumentParser(
        description='Convert genomic coordinate TSV to GFF3 + FASTA for run_synthology',
    )
    parser.add_argument('tsv_file',
                        help='Input TSV file (header row + 5-column data rows)')
    parser.add_argument('--mode', choices=['nucl', 'prot'], required=True,
                        help='Write nucleotide sequences (.fna) or translated proteins (.faa)')
    parser.add_argument('--genomes-dir', default=DEFAULT_GENOMES_DIR,
                        help=f'Directory containing {{genome_ac}}.fasta files '
                             f'(default: {DEFAULT_GENOMES_DIR})')
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR,
                        help=f'Output root directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--genetic-code', type=int, default=1,
                        help='NCBI translation table number, prot mode only (default: 1)')
    parser.add_argument('--source', default='.',
                        help='GFF3 source column value (default: .)')
    args = parser.parse_args()

    tsv_path = Path(args.tsv_file)
    if not tsv_path.exists():
        print(f"ERROR: input file not found: {tsv_path}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output_dir)

    print(f"Input TSV:    {tsv_path}")
    print(f"Mode:         {args.mode}")
    print(f"Genomes dir:  {args.genomes_dir}")
    print(f"Output dir:   {output_dir}")
    if args.mode == 'prot':
        print(f"Genetic code: {args.genetic_code}")

    rows = parse_tsv(str(tsv_path))
    print(f"Parsed {len(rows)} data rows from TSV")

    groups = group_by_genome(rows)
    print(f"Found {len(groups)} unique genome accessions")
    print()

    n_genomes, n_genes, n_skipped = process_genomes(
        groups, args.genomes_dir, output_dir,
        args.mode, args.genetic_code, args.source,
    )

    print()
    print("Done.")
    print(f"  Genomes processed: {n_genomes}")
    print(f"  Genes written:     {n_genes}")
    print(f"  Rows skipped:      {n_skipped}")

    if n_genes == 0:
        print("ERROR: zero genes written — check genome paths and TSV content", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
