#! /usr/bin/env python3

from sys import argv
import pickle
from copy import deepcopy

anchor_dir = argv[1]

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())

for org in orgs:
    with open(f'{anchor_dir}/aligned/{org}','rb') as f:
        am = pickle.load(f)

    bib2 = {}

    for i,bib in am.items():
        if 'dup' in bib:
            continue
        bib2[i] = deepcopy(bib)
        for org2,mbib1 in bib['matches'].items():
            if 'dups_matches' in mbib1:
                del bib2[i]['matches'][org2]['dups_matches']

    with open(f'{anchor_dir}/aligned/{org}_without_dups','wb') as f:
        pickle.dump(bib2,f)
