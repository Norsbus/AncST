#! /usr/bin/env python3

import pickle
from bisect import bisect_left
from Bio import SeqIO
from subprocess import run
import sys
sys.path.append('../../utils/')
from get_mapping import get_mapping


if __name__ == '__main__':

    org_mapping,chr_mapping = get_mapping()

    coords = {}
    seqid = 1

    margin = int(sys.argv[1])
    genomes = sys.argv[2]

    to_write = []

    with open('../coords') as f:
        for line in f:
            line = line.strip().split()
            org,chromo,start,end,orientation,hit_name = line
            if org not in coords:
                coords[org] = {}
            if chromo not in coords[org]:
                coords[org][chromo] = []
            genome = SeqIO.parse(f'{genomes}/{org}.fasta','fasta')
            for seq in genome:
                if seq.id == chromo:
                    if orientation == 'forward':
                        seq = seq[max(int(start)-margin,0):min(int(end)+margin,len(seq))]
                        seq.id = f'{seqid} {org} {chromo} {start} {end}'
                        seq.description = ''
                        to_write.append(seq)
                    else:
                        seq = seq[max(int(start)-margin,0):min(int(end)+margin,len(seq))].reverse_complement()
                        seq.id = f'{seqid} {org} {chromo} {start} {end}'
                        seq.description = ''
                        to_write.append(seq)
                    coords[org][chromo].append((max(int(start)-margin,0),min(int(end)+margin,len(seq)),orientation,seqid))
                    break
            seqid += 1

    SeqIO.write(to_write,'to_align','fasta')

    lines = set()
    with open('dialign.homology') as f:
        for line in f:
            #input('next')
            name1,chromo1,start1,end1,hit_start1,hit_end1,name2,chromo2,start2,end2,hit_start2,hit_end2,score,length,flipped = line.strip().split()
            org1 = org_mapping[name1.split('ele')[0].split('chr')[0]]
            org2 = org_mapping[name2.split('ele')[0].split('chr')[0]]
            chromo1 = chr_mapping[org1][chromo1]
            chromo2 = chr_mapping[org2][chromo2]
            start1 = int(start1)
            start2 = int(start2)
            end1 = int(end1)
            end2 = int(end2)
            length = int(length)
            hit_start1 = int(hit_start1)
            hit_start2 = int(hit_start2)
            hit_end1 = int(hit_end1)
            hit_end2 = int(hit_end2)
            real_start1 = start1 + hit_start1
            real_end1 = real_start1 + length
            real_start2 = start2 + hit_start2
            real_end2 = real_start2 + length

            if org1 not in coords or chromo1 not in coords[org1] or org2 not in coords or chromo2 not in coords[org2]:
                continue
            for coord in coords[org1][chromo1]:
                start_coord1 = int(coord[0])
                end_coord1 = int(coord[1])
                if (real_start1 > start_coord1 and real_start1 < end_coord1) and (real_end1 > start_coord1 and real_end1 < end_coord1):
                    pass
                else:
                    continue
                ori1 = coord[2]
                seqid1 = coord[3]
                for coord2 in coords[org2][chromo2]:
                    start_coord2 = int(coord2[0])
                    end_coord2 = int(coord2[1])
                    if (real_start2 > start_coord2 and real_start2 < end_coord2) and (real_end2 > start_coord2 and real_end2 < end_coord2):
                        pass
                    else:
                        continue
                    ori2 = coord2[2]
                    seqid2 = coord2[3]
                    if ori1 == ori2 and flipped == 'forward':
                        if ori1 == 'reverse':
                            start1 = end_coord1 - real_end1
                            start2 = end_coord2 - real_end2
                        else:
                            start1 = real_start1 - start_coord1
                            start2 = real_start2 - start_coord2
                    elif ori1 != ori2 and flipped == 'reverse':
                        if ori1 == 'reverse':
                            start1 = end_coord1 - real_end1
                            start2 = real_start2 - start_coord2
                        else:
                            start1 = real_start1 - start_coord1
                            start2 = end_coord2 - real_end2
                    else:
                        continue
                    lines.add(f'{seqid1}\t{seqid2}\t{start1}\t{start2}\t{length}\t{score}\n')

                
    with open('to_align.anc','w') as f:
        for line in lines:
            f.write(line)
