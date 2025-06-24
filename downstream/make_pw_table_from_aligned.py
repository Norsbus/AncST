#! /usr/bin/env python3


import pickle
import sys
from get_mapping_local import get_mapping

if __name__ == "__main__":

    anchor_dir = sys.argv[1]
    orgs = []
    with open('orgs','r') as f:
        for line in f:
            orgs.append(line.strip())

    org_mapping,chr_mapping = get_mapping()

    out = open('pairwise_alignments_table','w')

    done = {}
    c_id = 1
    score1against2 = 'NA'
    score2against1 = 'NA'

    candidates = {}
    for org in orgs:
        with open(f'{anchor_dir}/candidates//{org}','rb') as f:
            candidates[org] = pickle.load(f)

    for org in orgs:
        with open(f'{anchor_dir}/aligned/{org}','rb') as f:
            bib = pickle.load(f)
        done[org] = 1
        for i,bib2 in bib.items():
            chromo1 = bib2['chromosome'] 
            for org2,mbib in bib2['matches'].items():
                if mbib['meta']['multiple matches out of tolerance range'] == 1 or mbib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                    continue
                if org2 not in orgs or org2 in done:
                    continue
                for j,m in mbib['matches'].items():
                    chromo2 = candidates[org2][j]['chromosome']
                    start1,end1 = m[f'hit coordinates in (own) {org} candidate']
                    start1 += bib2['start']
                    end1 += bib2['start']
                    start2,end2 = m[f'hit coordinates in {org2} candidate']
                    start2 += candidates[org2][j]['start']
                    end2 += candidates[org2][j]['start']
                    score = m['match score']
                    ori = m['match is on other strand in other genome']
                    if ori == True:
                        ori = 'reverse'
                    else:
                        ori = 'forward'
                    out.write(f'{c_id}\t{org}\t{chromo1}\t{start1}\t{end1}\t{org2}\t{chromo2}\t{start2}\t{end2}\t{ori}\t{score}\t{score1against2}\t{score2against1}\n')
                    c_id += 1
                if 'dups_matches' in mbib and 'syntenic' in mbib['dups_matches']:
                    syn = mbib['dups_matches']['syntenic']
                    for j in syn:
                        m = mbib['dups_matches'][j]
                        chromo2 = candidates[org2][j]['chromosome']
                        start1,end1 = m[f'hit coordinates in (own) {org} candidate']
                        start1 += bib2['start']
                        end1 += bib2['start']
                        start2,end2 = m[f'hit coordinates in {org2} candidate']
                        start2 += candidates[org2][j]['start']
                        end2 += candidates[org2][j]['start']
                        score = m['match score']
                        ori = m['match is on other strand in other genome']
                        if ori == True:
                            ori = 'reverse'
                        else:
                            ori = 'forward'
                        out.write(f'{c_id}\t{org}\t{chromo1}\t{start1}\t{end1}\t{org2}\t{chromo2}\t{start2}\t{end2}\t{ori}\t{score}\t{score1against2}\t{score2against1}\n')
                        c_id += 1
    out.close()
