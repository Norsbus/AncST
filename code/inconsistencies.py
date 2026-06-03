#!/usr/bin/env python3

from Bio import SeqIO
import pickle
from math import sqrt,floor
from bisect import bisect_left
from subprocess import run
from multiprocessing import Pool
import os,pathlib
from sys import argv
from subprocesses import blast,clasp


def inconsistencies(bib,files,org):

    for file in files:
        with open(f'{work_dir}/inconsistencies/{org}/{file}','rb') as f:
            incons = pickle.load(f)
        org2 = f'{file}'.strip()
        for i,s in incons.items():
            for j in s:
                # Explicit key checks instead of bare try/except
                if i not in bib:
                    print(f"WARNING: anchor {i} not in bib for {org} - possible data corruption")
                    continue
                # 'matches' at anchor level always exists
                if org2 not in bib[i]['matches']:
                    # This can happen if the match was removed - not necessarily corruption
                    continue
                if 'matches' not in bib[i]['matches'][org2]:
                    # Only has dups_matches, skip
                    continue
                # 'matches not considered...' always exists (initialized in parse_bcamm)
                not_considered_key = 'matches not considered upon applying stricter score criterion since there are consistency issues'

                if j in bib[i]['matches'][org2]['matches']:
                    bib[i]['matches'][org2][not_considered_key][j] = bib[i]['matches'][org2]['matches'][j]
                    del bib[i]['matches'][org2]['matches'][j]
                else:
                    if j not in bib[i]['matches'][org2][not_considered_key]:
                        print(f"WARNING: match {j} not in 'matches' or 'not considered' for anchor {i} to {org2} - should not happen")

    return(bib)

def mark_in_me(bib,files,org):

    for file in files:
        with open(f'{work_dir}/mark_in_others/{org}/{file}','rb') as f:
            incons = pickle.load(f)
        org2 = f'{file}'.strip()
        # important not to check if its zero (from collect_output.py some may) cos that could overwrite important
        # info from parsing with itself as focal species
        for i,match_i_bib in incons.items():
            for match_i,marks in match_i_bib.items():
                if marks[0] == 1 or marks[1] == 1:
                    if i not in bib:
                        print(f"WARNING: anchor {i} not in bib for {org} (marked by {org2} anchor {match_i}) - possible data corruption")
                        continue
                    # 'matches' at anchor level always exists
                    if org2 not in bib[i]['matches']:
                        # This anchor in org doesn't have matches to org2, but org2 had matches to it
                        # This can happen legitimately if the match wasn't reciprocal
                        print(f"INFO: anchor {i} in {org} has no matches to {org2}, but {org2} anchor {match_i} matched to it")
                        continue
                    # 'meta' always exists when org2 entry exists (created in parse_bcamm)
                    bib[i]['matches'][org2]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 1

    return(bib)

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir_candidates = root + '/anchors/candidates'
    anchor_dir_aligned = root + '/anchors/aligned'
    work_dir = argv[-1]
    
    org = argv[1]

    try:
        with open(anchor_dir_aligned + f'/{org}','rb') as f:
            bib = pickle.load(f)
    except:
        with open(anchor_dir_candidates + f'/{org}','rb') as f:
            bib = pickle.load(f)

    # Handle missing directories gracefully
    incons_dir = f'{work_dir}/inconsistencies/{org}'
    if os.path.exists(incons_dir):
        files = [name for name in os.listdir(incons_dir)]
    else:
        files = []

    ssew_bib = inconsistencies(bib,files,org)

    mark_dir = f'{work_dir}/mark_in_others/{org}'
    if os.path.exists(mark_dir):
        files = [name for name in os.listdir(mark_dir)]
    else:
        files = []

    new_bib = mark_in_me(bib,files,org)
    
    with open(anchor_dir_aligned + f'/{org}','wb') as f:
        pickle.dump(new_bib,f)
