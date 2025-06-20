#!/usr/bin/env python3

from subprocess import run,check_output
import os,pickle
from time import sleep
from sys import argv

#run('./make_directories.py',shell=True)

orgs = []
with open(f"orgs","r") as f:
    for line in f:
        org = line.strip()
        orgs.append(org)

work_dir = argv[1]
code_dir = argv[2]

orgs_tuples = []
files = []

for i,org1 in enumerate(orgs):
    for org2 in orgs[i+1:]:
        orgs_tuple = org1 + '---' + org2
        if not os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'):
            orgs_tuple = org2 + '---' + org1
        chunks = [name.split('forward.')[1] for name in os.listdir(f'{work_dir}/sequences_to_compare/{orgs_tuple}/{org2}/forward_split')]
        for f in chunks:
            orgs_tuples.append(orgs_tuple)
            files.append(f)

x = 50
for org_tuple,f in zip(orgs_tuples,files):
    first = f'{code_dir}/bcamm.py {org_tuple} {f} 11 {work_dir} && awk \'{{{{ if ($1 != "#" && $8>40) print }}}}\' {work_dir}/clasp_out_forward/{org_tuple}/{f} > {work_dir}/clasp_out_forward/{org_tuple}/{f}_awk && awk \'{{{{ if ($1 != "#" && $8>40) print }}}}\' {work_dir}/clasp_out_reverse/{org_tuple}/{f} > {work_dir}/clasp_out_reverse/{org_tuple}/{f}_awk'
    second = f'{code_dir}/parse_bcamm.py {org_tuple} {f} 11 {work_dir} && touch {work_dir}/touch/{org_tuple}_parse_bcamm_done_{f}'
    cmd = f"""\
#!/usr/bin/env bash\n\
#SBATCH --job-name={org_tuple}.{f}\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --mem=4000\n\
#SBATCH --time=1-00:00:00\n\
#SBATCH --error log/run_bcamm/{org_tuple}.{f}.global.err\n\
#SBATCH --output log/run_bcamm/{org_tuple}.{f}.global.out\n\
eval "$(conda shell.bash hook)"\n\
conda activate /homes/biertank/karl//miniconda3/envs/snakemake\n\
{first} && {second}"""
    with open(f'slurm_scripts/myslurm_{org}_{f}','w') as f2:
        f2.write(cmd)
    out = check_output(['sinfo','-o','%C'])
    idle = int(out.decode().split()[1].split('/')[1])
    while idle < 10:
        print(f'sleeping since only {idle} cores are idle and I want at least 10')
        sleep(1)
        out = check_output(['sinfo','-o','%C'])
        idle = int(out.decode().split()[1].split('/')[1])
        print(f'idle (loop) : {idle}')
    run(f'sbatch slurm_scripts/myslurm_{org}_{f}',shell=True)
    x -= 1
    if x < 0:
        x = 50
        sleep(10)
