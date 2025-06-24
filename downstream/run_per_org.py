#!/usr/bin/env python3

from subprocess import run,check_output
import os
from time import sleep
import pathlib
from sys import argv

anchor_dir = argv[1]

with open('orgs_all_matches','r') as f:
    orgs = f.read().splitlines()

skip_some = False

for org in orgs:
    print(org)
    cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH -e slurm_log/{org}.err\n\
#SBATCH -o slurm_log/{org}.out\n\
#SBATCH --job-name={org}\n\
#SBATCH --cpus-per-task=5\n\
#SBATCH --mem=64000\n\
#SBATCH --time=7-00:00:00\n\
\n\
    ./compress_maps_and_ignore_multis_with_dups.py {org} {anchor_dir} && touch touch/compressed_{org}_done"""
    with open(f'slurm_job_scripts/myslurm_{org}','w') as f2:
        f2.write(cmd)
    if skip_some and counter_skip >= 0:
        run(f'sbatch slurm_job_scripts/myslurm_{org}',shell=True)
        counter_skip -= 1
        continue
    out = check_output(['sinfo','-o','%C'])
    idle = int(out.decode().split()[1].split('/')[1])
    while idle < 5:
        print(f'sleeping since only {idle} cores are idle and I want at least 10')
        out = check_output(['sinfo','-o','%C'])
        idle = int(out.decode().split()[1].split('/')[1])
        sleep(2)
    run(f'sbatch slurm_job_scripts/myslurm_{org}',shell=True)
    if idle > 10:
        skip_some = True
        counter_skip = 50
    sleep(1)

for org in orgs:
    while not os.path.exists(f'compressed_maps_multis_to_one/{org}'):
        sleep(10)
