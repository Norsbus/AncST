#! /usr/bin/env python3

from sys import argv
from subprocess import run
import os
import multiprocessing as mp
from os.path import isfile

def exec_blast_clasp(org):
    
    hit_no = 1
    with open(f'temp_coords_{org}','w') as out:
        for orientation in ['forward','reverse']:
            with open(f'clasp_out_{orientation}/{org}') as f:
                for line in f:
                    if line[0] == '#':
                        continue
                    line = line.split()
                    orig_prot = line[1]
                    chromo = line[2]
                    start = line[5]
                    end = line[6]
                    score = float(line[-1])
                    out.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\thit-{hit_no}_{org}_{orig_prot}\t{score}\n')
                    hit_no += 1

    return(0)

if __name__ == "__main__":

    path = os.getcwd()

    orgs = []
    with open('orgs') as f:
        for line in f:
            org = line.strip()
            if isfile(path + f'/fastas/{org}/all.fasta'):
                orgs.append(org)

    with mp.Pool(processes=11) as pool:
        res = pool.map_async(exec_blast_clasp,orgs).get()
