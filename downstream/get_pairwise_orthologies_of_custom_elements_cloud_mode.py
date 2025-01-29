#! /usr/bin/env python3

from subprocess import run
import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
import pickle
from sys import argv
from pprint import pprint

if __name__ == "__main__":

    margin = int(argv[1])
    mode = 'cloud'

    org_mapping,chr_mapping = get_mapping()

    first = 1
    last = 0
    clouds = {}
    with open(f'out_only_relevant_{margin}_cloud/cloudAP.txt') as f:
        for line in f:
            if first == 1:
                first = 0
                continue
            cloud_id,ele1,ele2,x,y = line.strip().split()
            if last == 0:
                last = cloud_id
            elif cloud_id != last:
                customs = 0
                for ele11,ele22 in clouds[last]:
                    if 'CUSTOM' in ele11:
                        customs += 1
                if customs > len(clouds[last])/2:
                    del clouds[last]
                last = cloud_id
            if cloud_id not in clouds:
                clouds[cloud_id] = []
            clouds[cloud_id].append((ele1,ele2))

    customs = 0
    for ele1,ele2 in clouds[last]:
        if 'CUSTOM' in ele1:
            customs += 1
    if customs > len(clouds[last])/2:
        del clouds[last]

    #for cloud,eles in clouds.items():
    #    for ele1,ele2 in eles:
    #        if ele1=='318orgchr2CUSTOMele1517735to1518061' or ele2=='318orgchr2CUSTOMele1517735to1518061':
    #            print(cloud)
    #exit(0)

    alloc = {}

    for cloud,l in clouds.items():
        for ele1,ele2 in l:
            org1 = ele1.split('chr')[0]
            org2 = ele2.split('chr')[0]
            if org1 == org2 or 'CUSTOM' not in ele1 or 'CUSTOM' not in ele2:
                continue
            if org1 not in alloc:
                alloc[org1] = {}
            if org2 not in alloc:
                alloc[org2] = {}
            if ele1 not in alloc[org1]:
                alloc[org1][ele1] = {}
            if ele2 not in alloc[org2]:
                alloc[org2][ele2] = {}
            if org2 not in alloc[org1][ele1]:
                alloc[org1][ele1][org2] = set()
            if org1 not in alloc[org2][ele2]:
                alloc[org2][ele2][org1] = set()
            alloc[org1][ele1][org2].add(ele2)
            alloc[org2][ele2][org1].add(ele1)

    run(f'rm -rf pairwise_orthologies_{margin}_{mode} && mkdir -p pairwise_orthologies_{margin}_{mode}',shell=True) 
    for org,bib in alloc.items():
        org = org_mapping[org]
        with open(f'pairwise_orthologies_{margin}_{mode}/{org}','w') as f:
            for ele,bib2 in bib.items():
                f.write(f'------------------------ element: {ele} ------------------------\n')
                f.write(f'------------------------ aligns with ------------------------\n')
                for org2,s in bib2.items():
                    org2 = org_mapping[org2]
                    f.write(f'\t --- {org2} ---\n')
                    f.write(f'\t\t {s}\n')
    with open(f'alloc_pairwise_{margin}_{mode}','wb') as f:
        pickle.dump(alloc,f)
