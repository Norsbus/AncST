#!/usr/bin/env python3

from subprocess import run,check_output
import os,pickle
from time import sleep
from sys import argv

orgs = []
with open(f"orgs","r") as f:
    for line in f:
        org = line.strip()
        orgs.append(org)

work_dir = argv[1]
code_dir = argv[2]
anchor_dir = argv[3]


for org in orgs:
    first = f'{code_dir}/checks.py {org} {work_dir} {anchor_dir}'
    cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --partition=main\n\
#SBATCH --job-name={org}_checks\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --mem=32000\n\
#SBATCH --time=1-00:00:00\n\
#SBATCH --error checks_log/{org}_err\n\
#SBATCH --output checks_log/{org}_out\n\
eval "$(conda shell.bash hook)"\n\
conda activate CONDAHOMEDIR/miniconda3/envs/snakemake\n\
{first}"""
    with open(f'slurm_scripts/myslurm_{org}','w') as f2:
        f2.write(cmd)
    run(f'sbatch slurm_scripts/myslurm_{org}',shell=True)
for org in orgs:
    first = f'{code_dir}/aligned_statistics.py {org} {work_dir} {anchor_dir}'
    cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --partition=main\n\
#SBATCH --job-name={org}_statistics\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --mem=32000\n\
#SBATCH --time=1-00:00:00\n\
#SBATCH --error statistics_log/{org}_err\n\
#SBATCH --output statistics_log/{org}_out\n\
eval "$(conda shell.bash hook)"\n\
conda activate CONDAHOMEDIR/miniconda3/envs/snakemake\n\
{first}"""
    with open(f'slurm_scripts/myslurm_{org}','w') as f2:
        f2.write(cmd)
    run(f'sbatch slurm_scripts/myslurm_{org}',shell=True)
