#! /usr/bin/env python3

from sys import argv
from subprocess import run
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_genome(genomes_dir,org):
    seqs = SeqIO.parse(f'{genomes_dir}/{org}.fasta', "fasta")
    s = {}
    for seq in seqs:
        s[seq.id] = seq.seq
    return(s)

def extract_fastas(org,regions,genomes_dir,output_dir):

    genome = get_genome(genomes_dir,org)

    records = []
    for chromo,start,end in regions:
        if chromo not in genome:
            continue
        seq = genome[chromo][start:end]
        record = SeqRecord(seq, id=f"{org}_{chromo}_{start}_{end}_forward", description="")
        records.append(record)

    run(f'mkdir -p {output_dir}/{org}',shell=True)
    SeqIO.write(records, f'{output_dir}/{org}/all.fasta', 'fasta')

    return(len(records))

if __name__ == "__main__":

    # args: syntenic_regions_file genomes_dir output_dir
    if len(argv) < 2:
        print("Usage: python extract_syntenic_fastas.py syntenic_regions_succinct [genomes_dir] [output_dir]")
        exit(1)

    regions_file = argv[1]
    genomes_dir = argv[2] if len(argv) > 2 else '../../utils/genomes'
    output_dir = argv[3] if len(argv) > 3 else 'fastas'

    path = os.getcwd()

    # parse syntenic regions
    regions_by_org = {}
    with open(regions_file,'r') as f:
        for line in f:
            if line.startswith('org') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 4:
                continue
            org,chromo,start,end = cols[0],cols[1],int(cols[2]),int(cols[3])
            if org not in regions_by_org:
                regions_by_org[org] = []
            regions_by_org[org].append((chromo,start,end))

    # extract fastas
    for org in regions_by_org:
        regions = regions_by_org[org]
        n = extract_fastas(org,regions,genomes_dir,output_dir)
        print(f'{org}: {n} regions extracted')
