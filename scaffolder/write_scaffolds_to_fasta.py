#! /usr/bin/env python3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from sys import argv
from subprocess import run

def write_scaffolds(genome,scaffolds):

    new_genome = {}
    with open(scaffolds) as f:
        for line in f:
            if len(line.strip()) == 0:
                continue
            if '>scaffold' in line:
                new_genome[line.strip()] = []
                last = line.strip()
            else:
                new_genome[last].append((line.strip()[:-1],line.strip()[-1]))
    
    records = []
    for new_scaffold,contigs in new_genome.items():
        print(f'writing {new_scaffold} with contigs:')
        print(contigs)
        new_seq = ''
        for contig,ori in contigs:
            contig = contig.strip()
            ori = ori.strip()
            print(f'looking for contig {contig} with ori {ori}')
            for seq in SeqIO.parse(genome,'fasta'):
                if seq.id == contig:
                    print(f'found contig {contig} with ori {ori}')
                    if ori == '+':
                        new_seq += str(seq.seq)
                    else:
                        new_seq += str(seq.reverse_complement().seq)
                    break
            new_seq += 100*'N'
        records.append(SeqRecord(Seq(new_seq[:-100]),id=f"{new_scaffold[1:]}",description=""))
    SeqIO.write(records, f'scaffolds_fastas/{genome.split("/")[-1]}2{scaffolds.split("/")[-1]}.fasta', "fasta")

    return(0)

if __name__ == "__main__":
    genome = argv[1].strip()
    scaffolds = argv[2].strip()
    run('mkdir -p scaffolds_fastas',shell=True)
    write_scaffolds(genome,scaffolds)
