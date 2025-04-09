#! /usr/bin/env python3

import pickle,pathlib
from statistics import mean,stdev,median
from sys import argv

def aligned(org):
    
    orgs = []
    with open(f'{work_dir}/orgs','r') as f:
        for l in f:
            orgs.append(l.strip())
    
    size = {}
    for o in orgs:
        with open(f'{work_dir}/../utils/small_meta/{o}','rb') as f: 
            seqids,seqlen = pickle.load(f)
        size[o] = seqlen[-1]

    with open(f'{anchor_dir}/aligned/{org}','rb') as f:
        anchor_map = pickle.load(f)

    print('-----')
    print(org)
    aligned = {}
    scores = {}
    for o in orgs:
        if o == org:
            continue
        aligned[o] = 0
        scores[o] = []
    for i,d in anchor_map.items():
        for o,d2 in d['matches'].items():
            if o not in orgs:
                continue
            if d2['meta']['multiple matches out of tolerance range'] == 1 or d2['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                continue
            for j,bib_m in d2['matches'].items():
                scores[o].append(bib_m['match score'])
                aligned[o] += bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0]
    
    for o in orgs:
        if o == org:
            continue

        if len(scores[o]) > 0:
            print(f'{aligned[o]/size[org]*100:.2f}% are aligned to with {o}')
            print(f'avg score with {o}: {mean(scores[o]):.2f}')
            print(f'median score with {o}: {median(scores[o]):.2f}')


if __name__ == "__main__":
    org = argv[1]
    print(f'--- {org} ---')
    work_dir = argv[2]
    anchor_dir = argv[3]
    aligned(org)
