#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
import multiprocessing as mp
import os

def thread_process(org):
    if not pathlib.Path(f'{work_dir}/../utils/macle_indices/{org}').exists():
        print(f'making macle index for {org}')
        cmd = f"""\
macle -s {work_dir}/../utils/genomes/{org}.fasta > {work_dir}/../utils/macle_indices/{org} && touch touch/{org}_macle_index_done\n\
"""
        run(cmd,shell=True)
    else:
        print(f'{org} index exists')
        #run(f'touch touch/{org}_macle_index_done',shell=True)
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = '/scr/k80san/karl//tmp'
    orgs = set()
    with open('compute_anchors_for') as f:
        for line in f:
            orgs.add(line.strip())
    #for org in orgs:
    #    thread_process(org)
    with mp.Pool(processes=9) as pool:
        res = pool.map_async(thread_process,orgs).get()
