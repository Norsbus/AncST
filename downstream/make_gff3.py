#! /usr/bin/env python3

from sys import argv
import pickle
from os.path import isfile

if len(argv) > 1:
    anchor_dir = argv[1]
else:
    anchor_dir = '../anchors/'

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())

all_anchors_gff = open('gff/all_anchors.gff3','w')
all_align_gff = open('gff/all_alignments.gff3','w')
for org in orgs:
    align_counter = 1
    if isfile(anchor_dir+f'/aligned/{org}_with_syn_eval'):
        ad = anchor_dir+f'/aligned/{org}_with_syn_eval'
    else:
        ad = anchor_dir+f'/aligned/{org}'
    with open(ad,'rb') as f:
        anchors = pickle.load(f)
    anchors_gff = open(f'gff/anchors_{org}.gff3','w')
    align_gff = open(f'gff/alignments_{org}.gff3','w')
    for i,bib in anchors.items():
        idd = f'{org}_{i}'
        chromo = bib['chromosome']
        start = bib['start']
        end = bib['end']
        any_match = 0
        if 'dup' in bib:
            dup = 'yes'
        else:
            dup = 'no'
        for org2,match_bib in bib['matches'].items():
            if ('matches' not in match_bib and ('dups_matches' not in match_bib or 'syntenic' not in match_bib['dups_matches'] or len(match_bib['dups_matches']['syntenic']) == 0)) or match_bib['meta']['multiple matches out of tolerance range'] == 1 or match_bib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                continue
            if 'matches' in match_bib:
                for j,m in match_bib['matches'].items():
                    if m['match is on other strand in other genome']:
                        strand = '-'
                    else:
                        strand = '+'
                    hitstart1,hitend1 = m[f'hit coordinates in (own) {org} candidate']
                    chr_hitstart1 = hitstart1 + start
                    chr_hitend1 = hitend1 + start
                    hitstart2,hitend2 = m[f'hit coordinates in {org2} candidate']
                    score = m['match score']
                    align_gff.write(f'{chromo} AncST alignment {chr_hitstart1} {chr_hitend1} {score} {strand} . ID=alignment{align_counter};Target={org2}_{j} {hitstart2} {hitend2}\n')
                    all_align_gff.write(f'{chromo} AncST alignment {chr_hitstart1} {chr_hitend1} {score} {strand} . ID=alignment{align_counter};Target={org2}_{j} {hitstart2} {hitend2}\n')
                    any_match = 1
                    align_counter += 1
            if 'dups_matches' in match_bib and 'syntenic' in match_bib['dups_matches']:
                js = match_bib['dups_matches']['syntenic']
                for j in js:
                    m = match_bib['dups_matches'][j]
                    if m['match is on other strand in other genome']:
                        strand = '-'
                    else:
                        strand = '+'
                    hitstart1,hitend1 = m[f'hit coordinates in (own) {org} candidate']
                    chr_hitstart1 = hitstart1 + start
                    chr_hitend1 = hitend1 + start
                    hitstart2,hitend2 = m[f'hit coordinates in {org2} candidate']
                    score = ['match score']
                    align_gff.write(f'{chromo} AncST alignment {chr_hitstart1} {chr_hitend1} {score} {strand} . ID=alignment{align_counter};Target={org2}_{j} {hitstart2} {hitend2}\n')
                    all_align_gff.write(f'{chromo} AncST alignment {chr_hitstart1} {chr_hitend1} {score} {strand} . ID=alignment{align_counter};Target={org2}_{j} {hitstart2} {hitend2}\n')
                    any_match = 1
                    align_counter += 1

        if any_match == 1:
            all_anchors_gff.write(f'{chromo}\tAncST\tDNA\t{start}\t{end}\t.\t+\t.\tID={idd};duplicate anchor={dup}\n') 
            anchors_gff.write(f'{chromo}\tAncST\tDNA\t{start}\t{end}\t.\t+\t.\tID={idd};duplicate anchor={dup}\n') 



    anchors_gff.close()
    align_gff.close()

all_anchors_gff.close()
all_align_gff.close()
