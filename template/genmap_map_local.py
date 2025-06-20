#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os
import pathlib

def thread_process(org,k_e):
    for k,e in k_e:
        if not pathlib.Path(f'{work_dir}/../utils/genmap_out/{org}/{k}_{e}.freq16').exists():
            print(f'making {k} {e} map for {org}')
            cmd = f"""\
genmap map -T 64 -fl -K {k} -E {e} -I {work_dir}/../utils/genmap_indices/{org} -O {work_dir}/../utils/genmap_out/{org}/{k}_{e} -r && touch touch/{org}_{k}_{e}_map_done\n\
"""
            run(cmd,shell=True)
        else:
            print(f'{org} map {k} {e} exists')
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = '/scr/k80san/karl/tmp'
    orgs = set()
    with open('compute_anchors_for') as f:
        for line in f:
            orgs.add(line.strip())
    k_e = {}
    with open(f'{work_dir}/genmap_params.txt') as f:
        for line in f:
            if '#' in line:
                continue
            line = line.strip().split()
            genome,k,e = line[:3]
            if genome not in k_e:
                k_e[genome] = set()
            k_e[genome].add((k,e))

    with open(f'{work_dir}/dups_params.txt') as f:
        for line in f:
            if '#' in line:
                continue
            line = line.strip().split()
            genome,k,e = line[:3]
            if genome not in k_e:
                k_e[genome] = set()
            k_e[genome].add((k,e))
            k,e = line[3:5]
            k_e[genome].add((k,e))
    for org in orgs:
        for k,e in k_e[org]:
            thread_process(org,[(k,e)])
