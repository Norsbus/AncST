#! /usr/bin/env python3

from subprocess import run
from sys import argv
import pathlib
from random import randint

work_dir = pathlib.Path(__file__).parent.resolve()
root = str(pathlib.Path(__file__).parents[1])
code_dir = root + '/code'
anchor_dir = root + '/anchors'
savespace = 'HOMEDIR/savespace'

rando = randint(0,10000)
while pathlib.Path(f'{savespace}/{rando}').exists():
    rando = randint(0,10000)

orgs = []
with open(f"{work_dir}/orgs","r") as f:
    for line in f:
        org = line.strip()
        orgs.append(org)

# save anchors before new run

run(f'cp -r {anchor_dir} {savespace}/save_anchor_dir/{rando}',shell=True)

run(f'{work_dir}/make_directories.py',shell=True)
run(f'rm {work_dir}/.snakemake/log/*',shell=True)

for org in orgs:
    if not pathlib.Path(f'{root}/utils/genmap_indicesi/{org}').exists():
        run(f'genmap index -F "{root}/utils/genomes/{org}.fasta" -I "{root}/utils/genmap_indices/{org}" 1> {work_dir}/logs/genmap_index_out 2> {work_dir}/logs/genmap_index_err',shell=True)

run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
run(f'snakemake --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_get_pot_windows -c6 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
no = 0
for org in orgs:
    if pathlib.Path(f'{work_dir}/touch/no_new_candidates_{org}').exists():
        no += 1
if no == len(orgs):
    print('no new candidates...exiting')
    exit(1)
else:

    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --keep-incomplete --rerun-incomplete --keep-incomplete -s {code_dir}/snakefile_blast_clasp -c6 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)

    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_anchor_candidates_and_prepare_comparison -c6 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)
    
    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_compared_full -c6 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)

    run(f'rm -f {work_dir}/.snakemake/locks/*',shell=True)
    run(f'snakemake --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_bcamm -c6 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True) 

    run(f'rm -f {work_dir}/.snakemake/locks/*,shell=True',shell=True)
    run(f'snakemake --keep-incomplete --rerun-incomplete -s {code_dir}/snakefile_update_matches -c6 --config work_dir="{work_dir}" code_dir="{code_dir}" root_dir="{root}"',shell=True)

for org in orgs:
    run(f'{code_dir}/checks.py {org} 1> checks_log/{org}/out 2> checks_log/{org}/err {work_dir}',shell=True)
    run(f'{code_dir}/aligned_statistics.py {org} 1> statistics_log/{org}/out 2> statistics_log/{org}/err {work_dir}',shell=True)
