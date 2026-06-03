#! /usr/bin/env python3

import pickle
import sys
import argparse
from os.path import isfile

def list_of_strings(arg):
    return arg.split(',')
parser = argparse.ArgumentParser()
parser.add_argument('--refs', type=list_of_strings)
parser.add_argument('--target', type=str)
args = parser.parse_args()
refs = args.refs
target_genome = target = args.target

ref_weights = {}
if isfile('ref_weights.txt'):
    with open(f'ref_weights.txt') as f:
        for line in f:
            ref,w = line.strip().split()
            w = int(float(w))
            ref_weights[ref] = w
nodes = {}
s_tot_v = set()

pairs = {}

for ref_genome in refs:
    ref_weight = ref_weights[ref_genome]
    first = 1
    with open(f'singles_out/blossom_v_infile_target_{target_genome}_ref_{ref_genome}.txt') as f:
        for line in f:
            if first == 1:
                first = 0
                continue
            n1,n2,s = line.strip().split()
            s = int(s)
            if int(n1) > int(n2):
                loc_n1 = n2
                loc_n2 = n1
            else:
                loc_n1 = n1
                loc_n2 = n2
            s_tot_v.add(loc_n1)
            s_tot_v.add(loc_n2)
            if loc_n1[:-1] not in pairs:
                pairs[loc_n1[:-1]] = set()
            pairs[loc_n1[:-1]].add(loc_n1)
            if loc_n2[:-1] not in pairs:
                pairs[loc_n2[:-1]] = set()
            pairs[loc_n2[:-1]].add(loc_n2)
            if loc_n1 not in nodes:
                nodes[loc_n1] = {}
            if loc_n2 not in nodes[loc_n1]:
                nodes[loc_n1][loc_n2] = s*ref_weight*-1
            else:
                nodes[loc_n1][loc_n2] -= s*ref_weight

count = 0
mapping = {}
mapping_for_reconstruction = {}

for p,m in pairs.items():
    if len(m) == 1:
        node = m.pop()
        if node[-1] == '1':
            n = node[:-1]+'2'
        else:
            n = node[:-1]+'1'
        s_tot_v.add(n)
        if n[-1] == '1':
            rec = n[:-1] + '_head'
        else:
            rec = n[:-1] + '_tail'
        mapping_for_reconstruction[count] = rec
        mapping[n] = count
        n = count
        count += 1


tot_v = len(s_tot_v)
with open(f'singles_out/blossom_v_infile_merged.txt','w') as f:
    f.write(f'{tot_v} fill_in_number_of_lines_-1\n')
    for n1,bib in nodes.items():
        if n1 not in mapping:
            if n1[-1] == '1':
                rec = n1[:-1] + '_head'
            else:
                rec = n1[:-1] + '_tail'
            mapping_for_reconstruction[count] = rec
            mapping[n1] = count
            n1 = count
            count += 1
        else:
            n1 = mapping[n1]
        for n2,s in bib.items():
            if n2 not in mapping:
                if n2[-1] == '1':
                    rec = n2[:-1] + '_head'
                else:
                    rec = n2[:-1] + '_tail'
                mapping_for_reconstruction[count] = rec
                mapping[n2] = count
                n2 = count
                count += 1
            else:
                n2 = mapping[n2]
            f.write(f'{n1} {n2} {s}\n')
    
    processed = {}
    s_tot_v = list(s_tot_v)
    for n1 in s_tot_v:
        for n2 in s_tot_v:
            if n1[:-1] == n2[:-1]:
                continue
            if int(n1) > int(n2):
                loc_n1 = n2
                loc_n2 = n1
            else:
                loc_n1 = n1
                loc_n2 = n2
            if (loc_n1,loc_n2) in processed:
                continue
            processed[(loc_n1,loc_n2)] = 1
            if loc_n1 in nodes and loc_n2 in nodes[loc_n1]:
                continue
            else:
                loc_n1 = mapping[loc_n1]
                loc_n2 = mapping[loc_n2]
                f.write(f'{loc_n1} {loc_n2} 0\n')

with open(f'singles_out/blossom_mapping.pickle','wb') as f:
    pickle.dump(mapping_for_reconstruction,f)
