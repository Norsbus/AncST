#! /usr/bin/env python3

from subprocess import run
import pickle
from sys import argv
from os.path import isfile

def get_al(org):
    orgs = []
    res = {}
    with open(f'orgs','r') as f:
        for l in f:
            orgs.append(l.strip())

    size = {}
    for o in orgs:
        with open(f'../utils/small_meta/{o}','rb') as f: 
            seqids,seqlen = pickle.load(f)
        size[o] = seqlen[-1]

    if isfile(f'../utils/anchors/aligned/{org}_with_syn_eval'):
        ad = f'../utils/anchors/aligned/{org}_with_syn_eval'
    else:
        ad = f'../utils/anchors/aligned/{org}'
    with open(ad,'rb') as f:
        anchor_map = pickle.load(f)

    aligned = {}
    scores = {}
    for o in orgs:
        if o == org:
            continue
        aligned[o] = 0
        scores[o] = []
    for i,d in anchor_map.items():
        for o,d2 in d['matches'].items():
            if o not in orgs:
                continue
            if d2['meta']['multiple matches out of tolerance range'] == 1 or ['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                continue
            if 'matches' in d2:
                for j,bib_m in d2['matches'].items():
                    scores[o].append(bib_m['match score'])
                    aligned[o] += bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0]
            if 'dups_matches' in d2 and 'syntenic' in d2['dups_matches']:
                js = d2['dups_matches']['syntenic']
                for j in js:
                    bib_m = d2['dups_matches'][j]
                    scores[o].append(bib_m['match score'])
                    aligned[o] += bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0]

    for o in orgs:
        if o == org:
            continue

        if len(scores[o]) > 0:
            res[o] = aligned[o]/size[org]
    return res

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())
gres = {}
for org in orgs:
    gres[org] = get_al(org)

orders = []
last = orgs[0]
while len(orders) != len(orgs) - 1:
    maxi = [0,'']
    for o in orgs:
        if o not in gres[last] or gres[last] == 0:
            continue
        if gres[last][o] > maxi[0] and o not in [xxx[1] for xxx in orders] and o != orgs[0]:
            maxi = [gres[last][o],o] 
    orders.append((last,maxi[1]))
    last = maxi[1]
first = 1
with open('orders_orgs','w') as f:
    for ref,target in orders:
        f.write(f'{target}\t{ref}\n')
        run(f'./sc_mcdraw.py {ref} {target} {first} {argv[1]}',shell=True)
        first = 0
