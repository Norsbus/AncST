#! /usr/bin/env python3

import pickle
from sys import argv

if __name__ == "__main__":

    anchor_dir = argv[1]
    orgs = []
    for line in open('orgs'):
        orgs.append(line.strip())
    
    all_orgs = set()

    for org in orgs:
        with open(anchor_dir+'/aligned/'+org,'rb') as f:
            am = pickle.load(f)
        for i,bib in am.items():
            for org2 in bib['matches']:
                all_orgs.add(org2)
    with open('orgs_all_matches','w') as f:
        for org in all_orgs:
            f.write(f'{org}\n')
