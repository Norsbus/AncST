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
#!/usr/bin/env bash\n\
#SBATCH --job-name=genmap_{org}\n\
#SBATCH --cpus-per-task=30\n\
#SBATCH --mem=64000\n\
#SBATCH --time=3-00:00:00\n\
#SBATCH --error log/gm_map/stderr_map_{org}\n\
#SBATCH --output log/gm_map/stdout_map_{org}\n\
eval "$(conda shell.bash hook)"\n\
conda activate /homes/biertank/karl//miniconda3/envs/snakemake\n\
genmap map -fl -K {k} -E {e} -I {work_dir}/../utils/genmap_indices/{org} -O {work_dir}/../utils/genmap_out/{org}/{k}_{e} -r && touch touch/{org}_{k}_{e}_map_done\n\
"""
            with open('run.slurm','w') as f:
                f.write(cmd)
            run(f'sbatch run.slurm',shell=True)
        else:
            print(f'{org} map {k} {e} exists')
            #run(f'touch touch/{org}_{k}_{e}_map_done',shell=True)
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = '/scr/k80san/karl/tmp'
    orgs = set()
    with open('compute_anchors_for') as f:
        for line in f:
            orgs.add(line.strip())
    k_e = {}
    if os.path.isfile(f'{work_dir}/genmap_params.txt'):
        with open(f'{work_dir}/genmap_params.txt') as f:
            for line in f:
                if '#' in line:
                    continue
                line = line.strip().split()
                genome,k,e = line[:3]
                if genome not in k_e:
                    k_e[genome] = set()
                k_e[genome].add((k,e))

    if os.path.isfile(f'{work_dir}/dups_params.txt'):
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
        if org in k_e:
            for k,e in k_e[org]:
                thread_process(org,[(k,e)])
