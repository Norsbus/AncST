#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os

def thread_process(org):
    #run(f'rm -r ../genmap_indices/{org} ../genmap_out/{org}',shell=True)
    #return 0
    if not pathlib.Path(f'../genmap_indices/{org}').exists():
        print(f'making index for {org}')
        run(f'genmap index -F "../genomes/{org}.fasta" -I "../genmap_indices/{org}"',shell=True)
    else:
        print(f'{org} index exists')
    if not pathlib.Path(f'../genmap_out/{org}/17_2.freq16').exists():
        print(f'making map for {org}')
        run(f'mkdir -p ../genmap_out/{org}',shell=True)
        #run(f'genmap map -fl -K 13 -E 0 -I ../genmap_indices/{org} -O ../genmap_out/{org}/13_0 -r',shell=True)
        #run(f'genmap map -fl -K 15 -E 0 -I ../genmap_indices/{org} -O ../genmap_out/{org}/15_0 -r',shell=True)
        run(f'genmap map -fl -K 15 -E 1 -I ../genmap_indices/{org} -O ../genmap_out/{org}/15_1 -r',shell=True)
        run(f'genmap map -fl -K 17 -E 2 -I ../genmap_indices/{org} -O ../genmap_out/{org}/17_2 -r',shell=True)
        run(f'genmap map -fl -K 21 -E 4 -I ../genmap_indices/{org} -O ../genmap_out/{org}/21_4 -r',shell=True)
    else:
        print(f'{org} map exists')
    return 0

if __name__ == "__main__":
    os.environ['TMP_DIR'] = 'HOMEDIR/tmp/'
    orgs = set()
    run(f'mkdir -p logs',shell=True)
    with open('orgs_full') as f:
        for line in f:
            orgs.add(line.strip())
    for org in orgs:
        thread_process(org)
