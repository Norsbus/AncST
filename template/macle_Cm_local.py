#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os
import pathlib

def thread_process(org,w,k):
    if not pathlib.Path(f'{work_dir}/../utils/macle_out/{org}/{w}_{k}.txt').exists():
        print(f'making {w} {k} map for {org}')
        cmd = f"""\
macle -i {work_dir}/../utils/macle_indices/{org} -w {w} -k {k} > {work_dir}/../utils/macle_out/{org}/{w}_{k}.txt && touch touch/{org}_{w}_{k}_macle_Cm_done\n\
"""
        run(cmd,shell=True)
    else:
        print(f'{org} macle {w} {k} exists')
        #run(f'touch touch/{org}_{w}_{k}_macle_Cm_done',shell=True)
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = '/scr/k80san/karl//tmp'
    orgs = set()
    with open('compute_anchors_for') as f:
        for line in f:
            orgs.add(line.strip())
    macle_paras = []
    with open('macle_params.txt') as f:
        for line in f:
            macle_paras.append(line.split())
    paras = {}
    for para in macle_paras:
        genome = para[0] 
        if genome not in paras:
            paras[genome] = set()
        len_lmers = int(para[1]) 
        interval = int(para[2]) 
        paras[genome].add((len_lmers,interval))
    for org in orgs:
        for len_lerms,interval in paras[org]:
            thread_process(org,len_lmers,interval)
