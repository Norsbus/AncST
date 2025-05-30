#! /usr/bin/env python3

from subprocess import run,check_output
from sys import argv
import pathlib
from random import randint
from time import sleep,time
import os
import math

run('chmod +x *py',shell=True)
run('chmod +x ../code/*py',shell=True)

g_start = time()
execution_times = open('execution_times.txt','w')

work_dir = pathlib.Path(__file__).parent.resolve()
root = str(pathlib.Path(__file__).parents[1])
code_dir = root + '/code'
anchor_dir = root + '/anchors'
#savespace = 'REMOTE/scr/k80san/karl/synteny/savespace'
#
#rando = randint(0,10000)
#while pathlib.Path(f'{savespace}/{rando}').exists():
#    rando = randint(0,10000)
#
orgs = []
with open(f"{work_dir}/orgs","r") as f:
    for line in f:
        org = line.strip()
        orgs.append(org)
genome_sizes = {}
for org in orgs:
    genome_sizes[org] = os.path.getsize(f'../utils/genomes/{org}.fasta')
    with open(f'k_e_{org}','w') as f:
        k = math.ceil(math.log(genome_sizes[org],4))
        if k - math.log(genome_sizes[org],4) > 0.5:
            k -= 1
        f.write(f'k={k}\ne=0\n')
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

start = time()
run(f'{work_dir}/make_directories.py',shell=True)
run(f'rm -f {work_dir}/.snakemake/log/*',shell=True)
passed = time() - start

execution_times.write(f'make dirs: {passed}\n')


caf = set()
with open('compute_anchors_for') as f:
    for line in f:
        caf.add(line.strip())

start = time()
run(f'./genmap_index_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep index_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle != len(caf):
    sleep(30)
    out = check_output('ls touch | grep index_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'genmap indices: {passed}\n')

start = time()

run(f'./genmap_map_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep map_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle != len(caf)*len(k_e):
    sleep(30)
    out = check_output('ls touch | grep map_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'genmap maps: {passed}\n')

start = time()

run(f'./run_metadata.py {work_dir}',shell=True)
out = check_output('ls touch | grep metadata_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle != len(caf):
    sleep(30)
    out = check_output('ls touch | grep metadata_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'metadata: {passed}\n')


start_filter = time()
no_new = {}
for org in orgs:
    no_new[org] = 0
run(f'rm ../anchors/candidates/*',shell=True)
run(f'rm ../anchors/aligned/*',shell=True)
for k,e,l,i,p in k_e:
    run(f'echo {k}_{e} >> save_params',shell=True)
    #run(f'{work_dir}/make_directories.py',shell=True)
    run(f'echo {k} {e} {l} {i} {p} > params',shell=True)
    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    start = time()
    run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_get_pot_windows_macle -c 63 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
    passed = time() - start
    execution_times.write(f'filtering k={k} e={e}: {passed}\n')
    for org in orgs:
        if pathlib.Path(f'{work_dir}/touch/no_new_candidates_{org}').exists():
            no_new[org] += 1
passed = time() - start_filter

execution_times.write(f'total window filtering: {passed}\n')
run(f'cp orig_params params',shell=True)
#run(f'rm ../anchors/candidates/*',shell=True)
#run(f'rm ../anchors/aligned/*',shell=True)
no = 0
for org in orgs:
    if no_new[org] == len(k_e):
        no += 1
if no == len(orgs):
    print('no new candidates...exiting')
    exit(0)
#exit(0)
start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete --keep-incomplete -s {code_dir}/snakefile_blast_clasp -c 63 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'blasting candidates: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_anchor_candidates_and_prepare_comparison -c 63 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'reconciling candidates: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_compared -c 63 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'prepping pw comparisons: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_bcamm -c 63 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True) 
passed = time() - start
execution_times.write(f'pariwise comparison: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*,shell=True',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_matches -c 63 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'match making: {passed}\n')

#start = time() 
#run(f'./run_checks_multi.py {work_dir} {code_dir} {anchor_dir}',shell=True)
#passed = time() - start
#execution_times.write(f'running_checks: {passed}\n')

passed = time() - g_start

execution_times.write(f'total: {passed}\n')

execution_times.close()
