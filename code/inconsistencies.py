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
                try:
                    if j in bib[i]['matches'][org2]['matches']:
                        bib[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues'][j] = bib[i]['matches'][org2]['matches'][j]
                        del bib[i]['matches'][org2]['matches'][j]
                    else:
                        if j not in bib[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues']:
                            print('should not happen')
                except:
                    pass

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
                    try:
                        bib[i]['matches'][org2]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 1
                    except:
                        print(f'could not set ambiguous flag for {org} anchor {i} matches in {org2}')

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
    
    files = [name for name in os.listdir(f'{work_dir}/inconsistencies/{org}')]

    ssew_bib = inconsistencies(bib,files,org)

    files = [name for name in os.listdir(f'{work_dir}/mark_in_others/{org}')]
    
    new_bib = mark_in_me(bib,files,org)
    
    with open(anchor_dir_aligned + f'/{org}','wb') as f:
        pickle.dump(new_bib,f)
