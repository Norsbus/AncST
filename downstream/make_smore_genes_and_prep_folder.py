#! /usr/bin/env python3

#NT_033779.5	Drosophilidae_0	1518061	1517735	-	94	766	NA	AATGGCCGAGCTACTTTTTCCACAACTGGAGCGTCCTTTGCCCTCGCTGCCATCGCTGCACTACACCCTGTTTGCCTACAGGGAGGAGTTGCGACGCCGCGATGCCCCGTTCATGAAGATGTCCACCATCAAGCTGCATCTCACGGACAACCTCATCCTGCAGACGATCAAGAACATCCGGCAGTATGACACCATCGAGATTATGAATCTCAACCAGGAGATAAACTTCAAGCGGCGGCTGACCAAACAAATGAGGAAGGTGCGCAAACTGGAGAAGCTCGGCTTGCACATCGATCCCCGGAAACTCAACGAGGACGGCAAATGGCA	NA	NA	NA

import re,pickle
from Bio import SeqIO
import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
from subprocess import run

def get_genome(genomes_dir,org):
    
    seqs = SeqIO.parse(f'{genomes_dir}/{org}.fasta', "fasta")
    s = {}
    for seq in seqs:
        s[seq.id] = seq.seq

    return(s)

if __name__ == '__main__':

    genomes_dir = sys.argv[1]

    org_mapping,chr_mapping = get_mapping()

    coords = {}
    orgs = set()
    with open('../coords','r') as f:
        for line in f:
            line = line.split()
            org = line[0]
            orgs.add(org)
            chromo = line[1]
            if org not in coords:
                coords[org] = []
            start = int(line[2])
            end = int(line[3])
            if line[4] == 'forward':
                strand = '+'
            else:
                strand = '-'
            coords[org].append((chromo,start,end,strand))

    run('rm -rf smore_genes && mkdir -p smore_genes',shell=True)
    for org in orgs:
        org_orig = org
        org = org_mapping[org]
        count = 0
        out = open(f'smore_genes/{org}.bed','w+')
        s = get_genome(genomes_dir,org_orig)
        for coord in coords[org_orig]:
            count += 1
            chromo,start,end,strand = coord
            if strand == 'reverse':
                seq = s[chromo][start:end].reverse_complement()
                strand = '-'
            else:
                seq = s[chromo][start:end]
                strand = '+'
            chromo = chr_mapping[org][chromo]
            closest_upstream = ('NA',1e15)
            closest_downstream = ('NA',1e15)
            with open(f'smore_anchors/{org}.bed','r') as f:
                for line in f:
                    a_chromo = line.split()[0]
                    a_strand = line.split()[4]
                    if chromo != a_chromo or strand != a_strand:
                        continue
                    ident = line.split()[1].split('_')[-1]
                    a_start = int(line.split()[2])
                    a_length = int(line.split()[3])
                    a_end = a_start + a_length
                    if a_end < start:
                        if start - a_end < closest_upstream[1]:
                            closest_upstream = (ident,start-a_end)
                    elif a_start > end:
                        if a_start - end < closest_downstream[1]:
                            closest_downstream = (ident,a_start-end)

            if closest_upstream[0] == 'NA':
                print('no upstream ')
            if closest_downstream[0] == 'NA':
                print('no downstream')
            out.write('\t'.join([str(x) for x in [chromo,org+f'_{count}',start,end,strand,closest_upstream[0],closest_downstream[0],'NA',seq,'NA','NA','NA']])+'\n')
        out.close()

    run(f'rm -rf prep && mkdir -p prep',shell=True)
    run(f'mv smore_genes genes && mv genes prep',shell=True)
    run(f'cp -r smore_anchors temp && gzip temp/* && mv temp prep',shell=True)
    run(f'mkdir -p out',shell=True)
