#! /usr/bin/env python3

import sys,os,ast,re
from get_mapping import get_mapping
from subprocess import run

org_mapping,chr_mapping = get_mapping()

orgs = []
with open('to_filter') as f:
    for line in f:
        if len(line) == 0 or '#' in line:
            continue
        orgs.append(line.strip())
orgs = [org_mapping[org] for org in orgs]

run('rm -rf ../filtered && mkdir -p ../filtered && cp -r ../* ../filtered',shell=True)

with open('which_files_to_produce_for_global') as f:
    files = f.read().splitlines()

for file in files:
    filtered = open('../filtered/'+file,'w')
    with open('../'+file) as f:
        for line in f:
            for org in orgs:
                pattern = '\D*' + org + '\D*'
                if re.match(pattern,line):
                    filtered.write(line)
                    break

    filtered.close()

with open('which_dirs_to_produce_for_global') as f:
    dirs = f.read().splitlines()

for dir in dirs:
    subdirs = os.listdir('../filtered/'+dir)
    for subdir in subdirs:
        valid = False
        for org in orgs:
            if org in subdir:
                valid = True
                break
        if valid == False:
            run(f'rm -r ../filtered/{dir}/{subdir}',shell=True)
