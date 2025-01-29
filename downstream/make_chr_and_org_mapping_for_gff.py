#! /usr/bin/env python3

import pickle
from sys import argv
from statistics import mean,median,stdev
import re

anchor_dir = argv[1]

orgs = []
with open('orgs', 'r') as f:
    for line in f:
        orgs.append(line.strip())
mapping = {}
counter = 1

for org in orgs:

    mapping[org] = {'org_id': '','chromosomes': {}}
    org_alt = org.strip()
    #org_alt = re.sub(r'[^a-zA-Z0-9]', '',org_alt)
    #org_alt = re.sub(r'[^0-9]', '',org_alt)

    mapping[org]['org_id'] = str(counter)+'org'
    #mapping[org][str(counter)+'org'] = org
    counter += 1

#anchor_dir = '../../anchors_to_save/anchors_after_successful_only_new_ones_run_with_7_0_100_50_10_11_11_50_THEY_PASSED_CHECKS'
#anchor_dir = 'aligned_6_flies'
#anchor_dir = '../../anchors_plants_new'

for org in orgs:
    counter = 1
    #with open(anchor_dir+'/candidates/'+org, 'rb') as f:
    with open(anchor_dir+'/aligned/'+org, 'rb') as f:
        candidates = pickle.load(f)
    for i,bib in candidates.items():
        if bib['chromosome'] not in mapping[org]['chromosomes']:
            mapping[org]['chromosomes'][bib['chromosome']] = f'{mapping[org]["org_id"]}chr{counter}'
            #mapping[org]['chromosomes'][f'{mapping[org]["org_id"]}chr{counter}'] = bib['chromosome']
            counter += 1

with open('mapping','w') as f:
    for org,bib in mapping.items():
        f.write('#######################\n')
        f.write(f'###\tspecies\t###\n')
        f.write(org+'\t'+str(bib["org_id"])+'\n')
        f.write(f'###\tchromosomes\t###\n')
        for orig_chromo,other_chromo in bib['chromosomes'].items():
            f.write(orig_chromo+'\t'+other_chromo+'\n')
        f.write('#######################\n')

with open('mapping_pickle','wb') as f:
    mapping = pickle.dump(mapping,f)
