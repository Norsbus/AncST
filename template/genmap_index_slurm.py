#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os

def thread_process(org):
    if not pathlib.Path(f'{work_dir}/../utils/genmap_indices/{org}').exists():
        print(f'making index for {org}')
        cmd = f"""\
genmap index -F {work_dir}/../utils/genomes/{org}.fasta -I {work_dir}/../utils/genmap_indices/{org} && touch touch/{org}_index_done\n\
"""
        run(cmd,shell=True)
    else:
        print(f'{org} index exists')
        #run(f'touch touch/{org}_index_done',shell=True)
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = '/scr/k80san/karl//tmp'
    orgs = set()
    with open('compute_anchors_for') as f:
        for line in f:
            orgs.add(line.strip())
    for org in orgs:
        thread_process(org)
