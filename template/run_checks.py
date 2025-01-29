#! /usr/bin/env python3

from subprocess import run
import pathlib

work_dir = pathlib.Path(__file__).parent.resolve()
root = str(pathlib.Path(__file__).parents[1])
code_dir = root + '/code'
anchor_dir = root + '/anchors/'

orgs = []

with open(f"{work_dir}/orgs","r") as f:
    for line in f:
        org = line.strip()
        orgs.append(org)

for org in orgs:
    run(f'{code_dir}/checks.py {org} 1> checks_log/{org}/out 2> checks_log/{org}/err {work_dir} {anchor_dir}',shell=True)
    run(f'{code_dir}/aligned_statistics.py {org} 1> statistics_log/{org}/out 2> statistics_log/{org}/err {work_dir} {anchor_dir}',shell=True)
