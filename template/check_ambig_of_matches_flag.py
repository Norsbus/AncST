#! /usr/bin/env python3

import pickle
from copy import deepcopy

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())

am = {}
for org in orgs:
    with open(f'anchors/aligned/{org}','rb') as f:
        am[org] = pickle.load(f)

correct = 0
incorrect = 0

for org,bib in am.items():
    new = {}
    for i,bib2 in bib.items():
        new_entry = {'chromosome':bib2['chromosome'],'start':bib2['start'],'end':bib2['end'],'matches':{}} 
        new_matches = deepcopy(bib2['matches'])
        for org2,bib3 in bib2['matches'].items():
            if org2 not in orgs:
                continue
            if bib3['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                valid = 0
                for j,bib4 in bib3['matches'].items():
                    if am[org2][j]['matches'][org]['meta']['multiple matches out of tolerance range'] == 1:
                        valid = 1
                        break
                if valid == 0:
                    new_matches[org2]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 0
                    incorrect += 1
                else:
                    correct += 1
        new_entry['matches'] = new_matches
        new[i] = new_entry
    with open(f'aligned_corrected/{org}','wb') as f:
        pickle.dump(new,f)
print(correct,incorrect)
