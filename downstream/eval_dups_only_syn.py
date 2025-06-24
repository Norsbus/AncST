#! /usr/bin/env python3

import pickle
from bisect import bisect_left
from pprint import pprint
from copy import deepcopy

if __name__ == "__main__":
    orgs = []
    with open('orgs') as f:
        for i in f:
            orgs.append(i.strip())
    for org in orgs:
        work_dir = '.'
        margin = tol = 100000
        with open(f'../anchors/aligned/{org}','rb') as f:
            am1 = pickle.load(f)
        with open(f'../anchors/aligned/{org}','rb') as f:
            am2 = pickle.load(f)
        iss = sorted(list(am2.keys()))
        new_bib = {}
        for i,bib in am1.items():
            new_bib[i] = deepcopy(bib)
            new_bib[i]['matches'] = {}
            for org2,bib2 in bib['matches'].items():
                new_bib[i]['matches'][org2] = {}
                new_bib[i]['matches'][org2]['meta'] = deepcopy(bib2['meta'])
                if 'matches' in bib2:
                    new_bib[i]['matches'][org2]['matches'] = deepcopy(bib2['matches'])
                if 'dups_matches' in bib2:
                    new_bib[i]['matches'][org2]['dups_matches'] = {}
                    local_iss = iss[bisect_left(iss,i-margin):min(len(iss), bisect_left(iss,i+(bib['end'] - bib['start'])+margin))]
                    syntenic = []
                    for j,m in bib2['dups_matches'].items():
                        syn_matches = [0,0]
                        ali_len = 0
                        for ii in local_iss:
                            if ii == i:
                                continue
                            if org2 in am2[ii]['matches'] and 'matches' in am2[ii]['matches'][org2] and am2[ii]['matches'][org2]['meta']['multiple matches out of tolerance range'] != 1:
                                for jj,mm in am2[ii]['matches'][org2]['matches'].items():
                                    if jj == j:
                                        continue
                                    if abs(jj-j) < tol:
                                        ali_len += mm['match score']
                                        if jj < j:
                                            syn_matches[0] += mm['match score']
                                        else:
                                            syn_matches[1] += mm['match score']
                        if syn_matches[0] >= 500 and syn_matches[1] >= 500:
                            syntenic.append(j)
                    if len(syntenic) > 0:
                        if 'syntenic' not in new_bib[i]['matches'][org2]['dups_matches']:
                            new_bib[i]['matches'][org2]['dups_matches']['syntenic'] = set()
                        for jjj in syntenic:
                            new_bib[i]['matches'][org2]['dups_matches'][jjj] = deepcopy(bib2['dups_matches'][jjj])
                            new_bib[i]['matches'][org2]['dups_matches']['syntenic'].add(jjj)
        with open(f'../anchors/aligned/{org}_with_syn_eval','wb') as f:
            pickle.dump(new_bib,f)
