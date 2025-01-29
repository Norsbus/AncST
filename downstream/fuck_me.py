#! /usr/bin/env python3

import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
org_mapping,chr_mapping = get_mapping()

org_families = {}

with open('458_orgs_with_family') as f:
    for line in f:
        if len(line.strip().split()) == 0:
            continue
        org,family = line.strip().split()
        org = org_mapping[org]
        org_families[org] = family
        if family not in org_families:
            org_families[family] = []
        org_families[family].append(org)

first = 1
last = 0
clouds = {}
with open(f'out_only_relevant_100000_cloud/cloudAP.txt') as f:
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

with open('clouds_with_halo.txt') as f:
    for line in f:
        cloud = line.strip()
        print('--------------------')
        print(cloud)
        for ele1,ele2 in clouds[cloud]:
            print(ele1,ele2)
            if ele1 == '318orgchr2CUSTOMele1517735to1518061':
                org2 = ele2.split('chr')[0]
                if org_families[org2] != 'Drosophilidae':
                    print(cloud)
                    print(org_families[org2])
            elif ele2 == '318orgchr2CUSTOMele1517735to1518061': 
                org1 = ele1.split('chr')[0]
                if org_families[org1] != 'Drosophilidae':
                    print(cloud)
                    print(org_families[org1])
