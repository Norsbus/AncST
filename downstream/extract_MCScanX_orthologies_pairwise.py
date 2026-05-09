#! /usr/bin/env python3

import pickle
import sys
sys.path.append('../../utils/')
from get_mapping import get_mapping
from subprocess import run
from pprint import pprint

def merge(c1,c2):
    c = c1 + c2
    return(list(set(c)))


if __name__ == "__main__":

    margin = int(sys.argv[1])

    org_mapping,chr_mapping = get_mapping()

    alloc = {}

    skip = {}
    lines = 0
    tot_eles = 0
    no_customs = 0

    with open(f'MCScanX_plus_coords_only_relevant_{margin}.collinearity','r') as f:
        for line in f:
            if 'Alignment' in line:
                if lines != 0 and no_customs < tot_eles / 2:#valid == 0:
                    for l in lines:
                        skip[l] = 1
                lines = []
                #valid = 0
                no_customs = 0
                tot_eles = 0
            else:
                tot_eles += 1
                if line.count('CUSTOM') != 2:
                    #valid = 1
                    no_customs += 1
                else:
                    lines.append(line)
    if no_customs < tot_eles / 2:#valid == 0:
        for l in lines:
            skip[l] = 1

    with open(f'MCScanX_plus_coords_only_relevant_{margin}.collinearity','r') as f:
        for line in f:
            if 'CUSTOM' in line and line not in skip:
                line = line.split(':')[1].split()
                first = line[0]
                org1 = first.split('chr')[0]
                if org1 not in alloc:
                    alloc[org1] = {}
                if first not in alloc[org1]:
                    alloc[org1][first] = {}
                second = line[1]
                org2 = second.split('chr')[0]
                if org2 not in alloc:
                    alloc[org2] = {}
                if second not in alloc[org2]:
                    alloc[org2][second] = {}
                if org1 not in alloc[org2][second]:
                    alloc[org2][second][org1] = set()
                if org2 not in alloc[org1][first]:
                    alloc[org1][first][org2] = set()
                alloc[org1][first][org2].add(second)
                alloc[org2][second][org1].add(first)
    
    run(f'rm -rf pairwise_orthologies_{margin} && mkdir -p pairwise_orthologies_{margin}',shell=True) 
    for org,bib in alloc.items():
        org = org_mapping[org]
        with open(f'pairwise_orthologies_{margin}/{org}','w') as f:
            for ele,bib2 in bib.items():
                f.write(f'------------------------ element: {ele} ------------------------\n')
                f.write(f'------------------------ aligns with ------------------------\n')
                for org2,s in bib2.items():
                    org2 = org_mapping[org2]
                    f.write(f'\t --- {org2} ---\n')
                    f.write(f'\t\t {s}\n')

    with open(f'alloc_pairwise_{margin}','wb') as f:
        pickle.dump(alloc,f)
