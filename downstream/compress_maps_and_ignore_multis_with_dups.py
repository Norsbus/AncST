#! /usr/bin/env python3

import pickle
from sys import argv
from pprint import pprint
from os.path import isfile

def check_orientation(match_bib):
    global cases
    global all_same
    oris = []
    if 'matches' in match_bib:
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
    elif 'dups_matches' in match_bib:
        for org,bib in match_bib['dups_matches'].items():
            if org == 'syntenic':
                continue
            if bib['match is on other strand in other genome']:
                oris.append('reverse')
            else:
                oris.append('forward')
        if len(set(oris)) == 1:
            all_same += 1
        cases += 1
        f_count = oris.count('forward')
        r_count = oris.count('reverse')
    else:
        return(False)
    if f_count >= r_count:
        return('forward')
    else:
        return('reverse')
    return(False)

def check_forward_ambiguity(match_bib):
    if match_bib['meta']['multiple matches out of tolerance range'] == 1 or match_bib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
        return(True)
    return(False)

org = argv[1]
anchor_dir = argv[2]

new = {}
cases = 0
all_same = 0

if isfile(anchor_dir+f'/aligned/{org}_with_syn_eval'):
    ad = anchor_dir+f'/aligned/{org}_with_syn_eval'
else:
    ad = anchor_dir+f'/aligned/{org}'

with open(ad,'rb') as f:
    aligned = pickle.load(f)
for i,bib1 in aligned.items():
    valid = False
    new_entry = {'chromosome':bib1['chromosome'],'start':bib1['start'],'end':bib1['end'],'matches':{}}
    for org2,bib2 in bib1['matches'].items():
        # take out ones without matches 
        if 'matches' in bib2 and len(bib2['matches'])==0 and 'dups_matches' in bib2 and len(bib2['dups_matches'])==0:
            continue
        # take out ones with forward ambiguity
        if check_forward_ambiguity(bib2):
            continue
        valid = True
        # check orientation
        checked_orientation = check_orientation(bib2)
        new_entry['matches'][org2] = {}
        #if len(bib2['matches']) > 1:
        #    maxi = 0
        #    for j,bib3 in bib2['matches'].items():
        #        if int(bib3['match score']) > maxi:
        #            maxi = int(bib3['match score'])
        #            if bib3['match is on other strand in other genome']:
        #                ori = 'reverse'
        #                check_ori = True
        #            else:
        #                ori = 'forward'
        #                check_ori = False
        #    for j,bib3 in bib2['matches'].items():
        #        if int(bib3['match score']) == maxi and bib3['match is on other strand in other genome'] == check_ori:
        #            new_entry['matches'][org2][j] = (int(bib3['match score']),bib3[f'hit coordinates in (own) {org} candidate'],ori)
        #            break
        #else:
        if 'matches' in bib2:
            for j,bib3 in sorted(bib2['matches'].items()):
                if bib3['match is on other strand in other genome']:
                    ori = 'reverse'
                else:
                    ori = 'forward'
                new_entry['matches'][org2][j] = (int(bib3['match score']),bib3[f'hit coordinates in (own) {org} candidate'],ori)
                break
        if 'dups_matches' in bib2 and 'syntenic' in bib2['dups_matches'] and len(new_entry['matches'][org2]) == 0:
            syn = bib2['dups_matches']['syntenic']
            del bib2['dups_matches']['syntenic']
            for j,bib3 in sorted(bib2['dups_matches'].items()):
                if j not in syn:
                    continue
                pprint(bib3)
                if bib3['match is on other strand in other genome']:
                    ori = 'reverse'
                else:
                    ori = 'forward'
                new_entry['matches'][org2][j] = (int(bib3['match score']),bib3[f'hit coordinates in (own) {org} candidate'],ori)
                break
    if valid:
        new[i] = new_entry

print(f'of {cases} alignments, {all_same} had the same orientation for all parts')
with open(f'compressed_maps_multis_to_one/{org}','wb') as f:
    pickle.dump(new,f)
