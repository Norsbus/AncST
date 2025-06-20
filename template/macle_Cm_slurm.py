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
#!/usr/bin/env bash\n\
#SBATCH --job-name=macle_Cm_{org}\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --mem=32000\n\
#SBATCH --time=3-00:00:00\n\
#SBATCH --error log/macle_Cm/stderr_map_{org}\n\
#SBATCH --output log/macle_Cm/stdout_map_{org}\n\
eval "$(conda shell.bash hook)"\n\
conda activate /homes/biertank/karl//miniconda3/envs/snakemake\n\
macle -i {work_dir}/../utils/macle_indices/{org} -w {w} -k {k} > {work_dir}/../utils/macle_out/{org}/{w}_{k}.txt && touch touch/{org}_{w}_{k}_macle_Cm_done\n\
"""
        with open('run.slurm','w') as f:
            f.write(cmd)
        run(f'sbatch run.slurm',shell=True)
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
