#! /usr/bin/env python3

from subprocess import run,check_output
from sys import argv
import pathlib
from random import randint
from time import sleep,time
import os

run('chmod +x *py',shell=True)
run('chmod +x ../code/*py',shell=True)

g_start = time()
execution_times = open('execution_times4.txt','w')

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

k_e = {}
with open(f'{work_dir}/genmap_params.txt') as f:
    for line in f:
        if '#' in line:
            continue
        line = line.strip().split()
        genome,k,e = line[:3]
        if genome not in k_e:
            k_e[genome] = set()
        k_e[genome].add((k,e))

# Dups parameters are optional
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

macle_paras = {}
macle_paras_no = {}

with open(f'{work_dir}/macle_params.txt') as f:
    for line in f:
        if '#' in line:
            continue
        line = line.strip().split()
        genome,w,p = line[:3]
        if genome not in macle_paras:
            macle_paras[genome] = set()
        macle_paras[genome].add((w,p))
        

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
done = 0
for o in caf:
    if os.path.isfile(f'../utils/macle_indices/{o}'):
        done += 1
run(f'./macle_index_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep macle_index_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
print(done,idle,len(caf))
while idle+done < len(caf):
    sleep(30)
    out = check_output('ls touch | grep macle_index_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'macle idx: {passed}\n')

start = time()
missing = 0
for o in caf:
    if o in macle_paras:
        for para in macle_paras[o]:
            if not os.path.isfile(f'../utils/macle_out/{o}/{para[0]}_{para[1]}.txt'):
                missing += 1
run(f'./macle_Cm_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep macle_Cm_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle < missing:
    sleep(30)
    out = check_output('ls touch | grep macle_Cm_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'macle Cm: {passed}\n')

start = time()
done = 0
for o in caf:
    if os.path.isdir(f'../utils/genmap_indices/{o}'):
        done += 1
run(f'./genmap_index_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep index_done  | wc -l',shell=True)
out2 = check_output('ls touch | grep macle_index_done  | wc -l',shell=True)
idle = int(out.decode().split()[0]) - int(out2.decode().split()[0])
while idle+done < len(caf):
    sleep(30)
    out = check_output('ls touch | grep index_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'genmap indices: {passed}\n')

start = time()

missing = 0
for o in caf:
    if o in k_e:
        for para in k_e[o]:
            if not os.path.isfile(f'../utils/genmap_out/{o}/{para[0]}_{para[1]}.freq16'):
                missing += 1
run(f'./genmap_map_slurm.py {work_dir}',shell=True)
out = check_output('ls touch | grep map_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle < missing:
    sleep(30)
    out = check_output('ls touch | grep map_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'genmap maps: {passed}\n')

start = time()

done = 0
for o in caf:
    if os.path.isdir(f'../utils/blastdbs/{o}.nsq') and os.path.isdir(f'../utils/small_meta/{o}') and os.path.isdir(f'../utils/metadata_genomes/{o}'):
        done += 1
run(f'./run_metadata.py {work_dir}',shell=True)
out = check_output('ls touch | grep metadata_done  | wc -l',shell=True)
idle = int(out.decode().split()[0])
while idle+done != len(orgs):
    sleep(30)
    out = check_output('ls touch | grep metadata_done  | wc -l',shell=True)
    idle = int(out.decode().split()[0])
passed = time() - start

execution_times.write(f'metadata: {passed}\n')
#
#
no_new = {}
for org in orgs:
    no_new[org] = 0

start = time()
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_get_pot_windows_macle_and_genmap_manual -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'filtering k={k} e={e}: {passed}\n')

for org in orgs:
    if pathlib.Path(f'{work_dir}/touch/no_new_candidates_{org}').exists():
        no_new[org] += 1

no = 0
for org in orgs:
    if no_new[org] == len(k_e):
        no += 1
if no == len(orgs):
    print('no new candidates...exiting')
    exit(0)

#for org in orgs:
#    if org not in caf:
#        run(f'cp ../save_anchor_candidates/{org} ../anchors/candidates/{org}',shell=True)

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete --keep-incomplete -s {code_dir}/snakefile_blast_clasp -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'blasting candidates: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_anchor_candidates_and_prepare_comparison -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'reconciling candidates: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete --keep-incomplete -s {code_dir}/snakefile_get_pot_dups -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'getting dups: {passed}\n')

run('./add_dups_to_anchors.py',shell=True)

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_compared -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'prepping pw comparisons: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_bcamm -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True) 
passed = time() - start
execution_times.write(f'pariwise comparison: {passed}\n')

start = time() 
run(f'rm -f {work_dir}/.snakemake/locks/*,shell=True',shell=True)
run(f'snakemake --scheduler greedy --conda-frontend conda --use-conda --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_matches -c 64 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
passed = time() - start
execution_times.write(f'match making: {passed}\n')

start = time() 
run(f'./run_checks_multi.py {work_dir} {code_dir} {anchor_dir}',shell=True)
passed = time() - start
execution_times.write(f'running_checks: {passed}\n')

run('./eval_dups_only_syn.py',shell=True)

passed = time() - g_start

execution_times.write(f'total: {passed}\n')

execution_times.close()
