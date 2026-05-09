#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os
import pathlib

def thread_process(org,k_e):
    for k,e in k_e:
        if not pathlib.Path(f'../genmap_out/{org}/{k}_{e}.freq16').exists():
            print(f'making index for {org}')
            cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --job-name=genmap_{org}\n\
#SBATCH --cpus-per-task=30\n\
#SBATCH --mem=64000\n\
#SBATCH --time=3-00:00:00\n\
#SBATCH --error stderr_map_{org}\n\
#SBATCH --output stdout_map_{org}\n\
eval "$(conda shell.bash hook)"\n\
conda activate /homes/biertank/karl//miniconda3/envs/snakemake\n\
genmap map -fl -K {k} -E {e} -I {work_dir}/../utils/genmap_indices/{org} -O {work_dir}/../utils/genmap_out/{org}/{k}_{e} -r && touch touch/{org}_{k}_{e}_map_done\n\
"""
            with open('run.slurm','w') as f:
                f.write(cmd)
            run(f'sbatch run.slurm',shell=True)
        else:
            print(f'{org} map {k} {e} exists')
        return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = 'HOMEDIR/tmp'
    orgs = set()
    run(f'mkdir -p logs',shell=True)
    with open('orgs') as f:
        for line in f:
            orgs.add(line.strip())
    k_e = set()
    with open(f'{work_dir}/filter_params') as f:
        for line in f:
            if '#' in line:
                continue
            line = line.strip().split()
            k,e,l,i,p = line
            k_e.add((k,e))
    k_e = list(k_e)
    for org in orgs:
        thread_process(org,k_e)
