#! /usr/bin/env python3

import pathlib
from subprocess import run
from sys import argv
from multiprocessing import Pool,cpu_count
import os

def thread_process(org):
    if not pathlib.Path(f'../genmap_indices/{org}').exists():
        print(f'making index for {org}')
        cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --job-name=genmap_{org}\n\
#SBATCH --cpus-per-task=30\n\
#SBATCH --mem=64000\n\
#SBATCH --time=3-00:00:00\n\
#SBATCH --error stderr_index_{org}\n\
#SBATCH --output stdout_index_{org}\n\
eval "$(conda shell.bash hook)"\n\
conda activate CONDAHOMEDIR/miniconda3/envs/snakemake\n\
genmap index -F {work_dir}/../utils/genomes/{org}.fasta -I {work_dir}/../utils/genmap_indices/{org} && touch touch/{org}_index_done\n\
"""
        with open('run.slurm','w') as f:
            f.write(cmd)
        run(f'sbatch run.slurm',shell=True)
    else:
        print(f'{org} index exists')
    return 0

if __name__ == "__main__":
    work_dir = argv[1]
    os.environ['TMP_DIR'] = 'HOMEDIR/tmp'
    orgs = set()
    run(f'mkdir -p logs',shell=True)
    with open('orgs') as f:
        for line in f:
            orgs.add(line.strip())
    for org in orgs:
        thread_process(org)
