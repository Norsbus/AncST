#! /usr/bin/python3

from os import listdir
import sys
sys.path.append('./utils/')
from get_mapping import get_mapping
from pprint import pprint
    

#org_mapping,chr_mapping = get_mapping()


mcscanx_orgs = listdir('MCScanX/pairwise_orthologies')
for org in mcscanx_orgs:
    print(org)
    with open(f'MCScanX/pairwise_orthologies/{org}') as f:
        for line in f:
            if 'element:' in line:
                ele = line.split('element:')[1].split(' -')[0].strip()
                continue
            elif 'GCA' in line or 'GCF' in line:
                org2 = line.split('---')[1].split('---')[0].strip()
                continue
            elif '{' in line:
                line = line.strip()[2:-2]
                eles2 = [x.strip().replace('\'','').replace('\"','') for x in line.split(',')]
                continue

