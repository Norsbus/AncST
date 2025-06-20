#! /usr/bin/env python3

from subprocess import run,check_output
from sys import argv
import pathlib
from random import randint
from time import sleep
import os
import multiprocessing

no_cores = multiprocessing.cpu_count()

work_dir = pathlib.Path(__file__).parent.resolve()
root = str(pathlib.Path(__file__).parents[1])
code_dir = root + '/code'
anchor_dir = root + '/anchors'
#savespace = 'REMOTE/scr/k80san/karl/synteny/savespace'

#rando = randint(0,10000)
#while pathlib.Path(f'{savespace}/{rando}').exists():
#    rando = randint(0,10000)

orgs = []
with open(f"{work_dir}/orgs","r") as f:
    for line in f:
        org = line.strip()
        orgs.append(org)
k_e = set()
with open(f'{work_dir}/filter_params') as f:
    for line in f:
        if '#' in line:
            continue
        line = line.strip().split()
        k,e,l,i,p = line[:5]
        k_e.add((k,e,l,i,p))

# YOU NOW NEED FILTER_PARAMS with first 5 (k,e,l,i,p) defined
#input('are you happy with the part size in snakefile_get_pot and update_compared? i changed to 100000 in snakefile from 1000000. probably good to have a smaller number in update_compared cos there its not against whole genome but other candidate sets')
#input('have you changed to focus on specific chromosomes or general filtering?')
#input('have you checked if you are happy with the evalue for blast in subprocess?')
#input('have you changed clasp.x bin dir in ../code/subprocesses.py?')
# save anchors before new run

#run(f'cp -r {anchor_dir} {savespace}/save_anchor_dir/{rando}_starting_15_0',shell=True)

run(f'{work_dir}/make_directories.py',shell=True)
run(f'rm -f {work_dir}/.snakemake/log/*',shell=True)
caf = set()
with open('compute_anchors_for') as f:
    for line in f:
        caf.add(line.strip())

run(f'./genmap_index_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep index_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle != len(caf):
    sleep(30)
    out = check_output('ls touch | grep index_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])

run(f'./genmap_map_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep map_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle != len(caf)*len(k_e):
    sleep(30)
    out = check_output('ls touch | grep map_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])

run(f'./run_metadata.py {work_dir}',shell=True)
out = check_output('ls touch | grep metadata_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle != len(caf):
    sleep(30)
    out = check_output('ls touch | grep metadata_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])

no_new = {}
for org in orgs:
    no_new[org] = 0
for k,e,l,i,p in k_e:
    run(f'echo {k}_{e} >> save_params',shell=True)
    run(f'{work_dir}/make_directories.py',shell=True)
    run(f'echo {k} {e} {l} {i} {p} > params',shell=True)
    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_get_pot_windows -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
    for org in orgs:
        if pathlib.Path(f'{work_dir}/touch/no_new_candidates_{org}').exists():
            no_new[org] += 1
run(f'cp orig_params params',shell=True)
run(f'rm ../anchors/candidates/*',shell=True)
no = 0
for org in orgs:
    if no_new[org] == len(k_e):
        no += 1
if no == len(orgs):
    print('no new candidates...exiting')
    exit(1)
else:
     
    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete --keep-incomplete -s {code_dir}/snakefile_blast_clasp -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)

    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_anchor_candidates_and_prepare_comparison -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)

    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_compared -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)

    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_bcamm -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True) 

    run(f'rm -f {work_dir}/.snakemake/locks/*,shell=True',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_matches -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
    
    run(f'rm -f {work_dir}/.snakemake/locks/*,shell=True',shell=True)
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_checks -c {no_cores} --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
