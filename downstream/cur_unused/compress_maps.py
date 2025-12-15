#! /usr/bin/env python3

import pickle
from sys import argv

def check_orientation(match_bib):
    global cases
    global all_same
    oris = []
    for org,bib in match_bib['matches'].items():
        if bib['match is on other strand in other genome']:
            oris.append('reverse')
        else:
            oris.append('forward')
    if len(set(oris)) == 1:
        all_same += 1
    cases += 1
    f_count = oris.count('forward')
    r_count = oris.count('reverse')
    if f_count >= r_count:
        return('forward')
    else:
        return('reverse')
    return(False)

def check_forward_ambiguity(match_bib):
    if match_bib['meta']['multiple matches out of tolerance range'] == 1 or match_bib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
        return(True)
    return(False)

orgs = []
with open('orgs_all_matches') as f:
    for line in f:
        orgs.append(line.strip())

org = argv[1]
anchor_dir = argv[2]

new = {}
cases = 0
all_same = 0
with open(anchor_dir+f'/aligned/{org}','rb') as f:
    aligned = pickle.load(f)
for i,bib1 in aligned.items():
    valid = False
    new_entry = {'chromosome':bib1['chromosome'],'start':bib1['start'],'end':bib1['end'],'matches':{}}
    for org2,bib2 in bib1['matches'].items():
        # take out ones without matches 
        if len(bib2['matches'])==0:
            continue
        # take out ones with forward ambiguity
        if check_forward_ambiguity(bib2):
            continue
        valid = True
        # check orientation
        checked_orientation = check_orientation(bib2)
        new_entry['matches'][org2] = {}
        for j,bib3 in bib2['matches'].items():
            if bib3['match is on other strand in other genome']:
                ori = 'reverse'
            else:
                ori = 'forward'
            new_entry['matches'][org2][j] = (int(bib3['match score']),bib3[f'hit coordinates in (own) {org} candidate'],ori)
    if valid:
        new[i] = new_entry
    if i == 709101991:
        print(valid)
print(f'of {cases} alignments, {all_same} had the same orientation for all parts')
with open(f'compressed_maps/{org}','wb') as f:
    pickle.dump(new,f)
