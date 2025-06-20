#!/usr/bin/env python3

from subprocess import run,check_output
import os,pickle
from time import sleep
from sys import argv

#with open('orgs_tuple_to_bcamm','rb') as f:
#    otto = pickle.load(f)

work_dir = argv[1]
code_dir = argv[2]
anchor_dir = argv[3]

otto = []
with open('species_to_compare') as f:
    for line in f:
        org1,org2 = line.strip().split()
        orgs_tuple = org1 + '---' + org2
        if not os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'):
            orgs_tuple = org2 + '---' + org1
        otto.append(orgs_tuple)

for org_tuple in otto:
    first = f'{work_dir}/write_less_seqs.py {org_tuple} {work_dir} {code_dir} {anchor_dir} && touch touch/write_less_seqs_done_{org_tuple}'
    cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --job-name={org_tuple}\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --mem=64000\n\
#SBATCH --time=1-00:00:00\n\
#SBATCH --error logs/{org_tuple}_write_less.err\n\
#SBATCH --output logs/{org_tuple}_write_less.out\n\
eval "$(conda shell.bash hook)"\n\
conda activate /homes/biertank/karl//miniconda3/envs/snakemake\n\
{first}"""
    with open(f'slurm_scripts/myslurm_{org_tuple}_write_less','w') as f2:
        f2.write(cmd)
    sleep(1)
    out = check_output(['sinfo','-o','%C'])
    idle = int(out.decode().split()[1].split('/')[1])
    while idle < 5:
        print(f'sleeping since only {idle} cores are idle and I want at least 5')
        sleep(30)
        out = check_output(['sinfo','-o','%C'])
        idle = int(out.decode().split()[1].split('/')[1])
        print(f'idle (loop) : {idle}')
    run(f'sbatch slurm_scripts/myslurm_{org_tuple}_write_less',shell=True)
    sleep(1)
