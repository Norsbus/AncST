#! /usr/bin/env python3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from sys import argv
from subprocess import run
from pprint import pprint

def write_scaffolds(genome,scaffolds):

    mapping = {}
    
    new_genome = {}
    with open(scaffolds) as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if '>' in line:
                new_genome[line.strip()] = []
                last = line.strip()
                place = 0
            elif len(line.split()) > 0:
                new_genome[last].append((line.strip()[:-1].strip(),line.strip()[-1].strip()))
                mapping[line.strip()[:-1].strip()] = (last,place,line.strip()[-1].strip())
                place += 1

    seq_lists = {}

    for new_scaffold,contigs in new_genome.items():
        seq_lists[new_scaffold] = ['' for i in range(len(contigs))]

    done = {}

    for seq in SeqIO.parse(genome,'fasta'):
        if seq.id in mapping:
            scaff,place,ori = mapping[seq.id]
            if ori == '+' or ori == '0':
                seq_lists[scaff][place] = str(seq.seq)
            elif ori == '-' or ori == '1':
                seq_lists[scaff][place] = str(seq.reverse_complement().seq)
            done[seq.id] = 1

    records = []
    join_str = "N"*100
    for new_scaffold,contigs in seq_lists.items():
        new_seq = join_str.join(contigs)
        records.append(SeqRecord(Seq(new_seq),id=f"{new_scaffold[1:]}",description=""))
    c = 0
    for seq in SeqIO.parse(genome,'fasta'):
        if seq.id not in done:
            c += 1
            records.append(seq)
    SeqIO.write(records, f'fastas/{outname}', "fasta")

    return(0)

if __name__ == "__main__":
    genome = argv[1].strip()
    scaffolds = argv[2].strip()
    outname = argv[3]
    run('mkdir -p fastas',shell=True)
    write_scaffolds(genome,scaffolds)
