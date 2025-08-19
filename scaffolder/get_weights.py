#! /usr/bin/env python3

from sys import argv
from subprocess import run
import pickle
from os.path import isfile

anchors_path = '../anchors/'

refs = []
with open('refs.txt','r') as f:
    for line in f:
        refs.append(line.strip())

target = argv[1]

weights = {}
max_w = 0
for ref in refs:
    aligned_len = 0
    with open(f'{anchors_path}/aligned/{ref}','rb') as f:
        anchor_map = pickle.load(f)
    for i,d in anchor_map.items():
        if target in d['matches']:
            if 'matches' in d['matches'][target]:
                if d['matches'][target]['meta']['multiple matches out of tolerance range'] == 1 or d['matches'][target]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                    continue
                for j,bib_m in d['matches'][target]['matches'].items():
                    aligned_len += bib_m[f'hit coordinates in {target} candidate'][1] - bib_m[f'hit coordinates in {target} candidate'][0]
    weights[ref] = aligned_len
    if aligned_len > max_w:
        max_w = aligned_len

while max_w > 1e5:
    print(max_w)
    for ref in weights:
        if int(weights[ref]) == int(max_w):
            max_w = int(max_w/10)
        weights[ref] = int(weights[ref]/10)
        if weights[ref] < 1:
            weights[ref] = 1


with open('ref_weights.txt','w') as f:
    for ref,aligned_len in weights.items():
        f.write(f'{ref}\t{int(aligned_len)}\n')
