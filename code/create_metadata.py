#!/usr/bin/env python3
"""
Create genome metadata pickles for AncST pipeline.
Creates two pickle files:
- small_meta/{org} - (seqids, seqlen) for chromosome info
- metadata_genomes/{org} - (seqids, seqlen, full_sequence) for detailed info

Called from Snakemake Phase 0 preprocessing.
"""

from Bio import SeqIO
import pickle
from sys import argv
import os
from pathlib import Path


def make_small_metadata(org, root_dir):
    """
    Create small metadata pickle with chromosome IDs and cumulative lengths.

    Args:
        org: Organism identifier
        root_dir: Root directory containing utils/

    Returns:
        Path to created small_meta file
    """
    output_path = Path(root_dir) / "utils" / "small_meta" / org

    if output_path.exists():
        print(f"Small metadata already exists: {output_path}")
        return output_path

    genome_path = Path(root_dir) / "utils" / "genomes" / f"{org}.fasta"

    seqs = SeqIO.parse(str(genome_path), "fasta")
    seqids = []
    seqlen = []
    tot = 0

    for seq in seqs:
        seqids.append(seq.id)
        tot += len(seq.seq)
        seqlen.append(tot)

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'wb') as f:
        pickle.dump((seqids, seqlen), f)

    print(f"Created small metadata: {output_path}")
    print(f"  Chromosomes: {len(seqids)}, Total length: {tot:,} bp")

    return output_path


def make_metadata(org, root_dir):
    """
    Create full metadata pickle with chromosome IDs, lengths, and full sequence.

    Args:
        org: Organism identifier
        root_dir: Root directory containing utils/

    Returns:
        Path to created metadata file
    """
    output_path = Path(root_dir) / "utils" / "metadata_genomes" / org

    if output_path.exists():
        print(f"Metadata already exists: {output_path}")
        return output_path

    genome_path = Path(root_dir) / "utils" / "genomes" / f"{org}.fasta"

    seqs = SeqIO.parse(str(genome_path), "fasta")
    s = ''
    seqids = []
    seqlen = []
    tot = 0

    for seq in seqs:
        s += seq.seq
        seqids.append(seq.id)
        tot += len(seq.seq)
        seqlen.append(tot)

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'wb') as f:
        pickle.dump((seqids, seqlen, s), f)

    print(f"Created full metadata: {output_path}")
    print(f"  Chromosomes: {len(seqids)}, Total length: {tot:,} bp")
    print(f"  Full sequence stored: {len(s):,} characters")

    return output_path


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Create genome metadata pickles for AncST pipeline')
    parser.add_argument('organism', help='Organism identifier (e.g., GCF_000001215.4)')
    parser.add_argument('root_dir', help='Root directory containing utils/ subdirectory')
    parser.add_argument('--task', choices=['small_meta', 'metadata', 'both'], default='both',
                        help='Which metadata to create: small_meta, metadata, or both (default)')

    args = parser.parse_args()

    org = args.organism
    root_dir = args.root_dir
    task = args.task

    print(f"Creating metadata for {org}")
    print(f"Root directory: {root_dir}")
    print(f"Task: {task}")

    # Create requested metadata files
    if task in ['small_meta', 'both']:
        make_small_metadata(org, root_dir)

    if task in ['metadata', 'both']:
        make_metadata(org, root_dir)

    print(f"Metadata creation complete for {org}")
