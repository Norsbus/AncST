#! /usr/bin/env python3

from bisect import bisect_left
import sys
sys.path.append('../utils/')
from get_mapping import get_mapping

def write():

    out = open('pairwise_alignments_table_only_relevant','w')
    with open('pairwise_alignments_table') as f:
        for line in f:
            # 1	Anopheles_gambiae_AgamP4_dna_toplevel	2L	1289026	1289622	GCF_013141755.1	seq3	43110301	43110897	reverse	527.0	NA	NA
            id,org1,chromo1,start1,end1,org2,chromo2,start2,end2,ori_alignment,score,alt_score1,alt_score2 = line.strip().split()
            start1 = int(start1)
            end1 = int(end1)
            start2 = int(start2)
            end2 = int(end2)
            if org1 not in coords or chromo1 not in coords[org1] or org2 not in coords or chromo2 not in coords[org2]:
                continue
            valid1 = 0
            valid2 = 0
            for s,e in coords[org1][chromo1]:
                if abs(s-start1) < margin or abs(e-end1) < margin:
                    valid1 = 1
            for s,e in coords[org2][chromo2]:
                if abs(s-start2) < margin or abs(e-end2) < margin:
                    valid2 = 1
            if valid1 == 1 and valid2 == 1:
                out.write(line)

    return 0

if __name__ == "__main__":
    
    margin = int(sys.argv[1])
   
    coords = {}
    with open('../utils/coords_both_genomic_and_regional_nono') as f:
        for line in f:
            org,chromo,start,end,ori,hit_name = line.strip().split()
            if org not in coords:
                coords[org] = {}
            if chromo not in coords[org]:
                coords[org][chromo] = []
            coords[org][chromo].append((int(start),int(end)))

    write()
