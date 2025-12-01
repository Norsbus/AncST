#! /usr/bin/env python3

from sys import argv
from subprocess import run
import pickle
from os.path import isfile
import sys

def get_chr_l(org,chromo):
    idx = small_meta[org][0].index(chromo)
    if idx == 0:
        return(small_meta[org][1][0])
    else:
        return(small_meta[org][1][idx] - small_meta[org][1][idx-1])

ref_genome = argv[1]

small_meta_path = '../utils/small_meta/'
anchors_path = '../utils/anchors/aligned/'
anchors_path_cand = '../utils/anchors/candidates/'
small_meta = {}

target_genome = argv[2]
orgs = [ref_genome,target_genome]
for org in orgs:
    with open(f'{small_meta_path}/{org}','rb') as f:
        small_meta[org] = pickle.load(f)

targets = [target_genome]
chromo_mapping2 = {}

orders_others = {}
oris_others = {}
orders_others[target_genome] = {}
oris_others[target_genome] = {}

with open(anchors_path+f'/{ref_genome}','rb') as f:
    anchors = pickle.load(f)
with open(anchors_path_cand+f'/{target_genome}','rb') as f:
    t_anchors = pickle.load(f)

mapping_contig_chromo = {}

first = int(argv[3])
if first == 1:
    iss = sorted(list(anchors.keys()))
else:
    iss = []
    ref_order = []
    with open(f'orders/{ref_genome}') as f:
        for line in f:
            seq, ori = line.strip().split('\t')
            ref_order.append((seq,ori))
    for seq,ori in ref_order:
        tmp = []
        for i,bib in anchors.items():
            if bib['chromosome'] == seq:
                tmp.append(i)
        iss += sorted(tmp)

for i in iss:
    chromo = anchors[i]['chromosome'] 
    if get_chr_l(ref_genome,chromo) < int(argv[4]):
        continue
    if chromo not in orders_others[target_genome]:
        orders_others[target_genome][chromo] = []
        oris_others[target_genome][chromo] = {}
    if 'matches' in anchors[i] and target_genome in anchors[i]['matches'] and 'matches' in anchors[i]['matches'][target_genome] and anchors[i]['matches'][target_genome]['meta']['multiple matches out of tolerance range'] == 0:
        for j,match_bib in anchors[i]['matches'][target_genome]['matches'].items():
            chromo2 = t_anchors[j]['chromosome'] 
            if get_chr_l(target_genome,chromo2) < int(argv[4]):
                continue
            if chromo2 not in oris_others[target_genome][chromo]:
                oris_others[target_genome][chromo][chromo2] = []
            if chromo2 not in mapping_contig_chromo:
                mapping_contig_chromo[chromo2] = {}
            if chromo not in mapping_contig_chromo[chromo2]:
                mapping_contig_chromo[chromo2][chromo] = 0
            mapping_contig_chromo[chromo2][chromo] += 1
            oris = match_bib['match is on other strand in other genome']
            orders_others[target_genome][chromo].append(chromo2)
            if oris==0:
                oris_others[target_genome][chromo][chromo2].append('forward')
            else:
                oris_others[target_genome][chromo][chromo2].append('reverse')

f2 = open(f'singles_out/significantly_divergent_orientations_from_AncST_target_{target_genome}_ref_{ref_genome}.txt','w')
orientations = {}
with open(f'singles_out/orientations_and_hits_from_AncST_target_{target_genome}_ref_{ref_genome}.txt','w') as f:
    for org in oris_others:
        for chromo in oris_others[org]:
            for chromo2 in oris_others[org][chromo]:
                if oris_others[org][chromo][chromo2].count('forward') > oris_others[org][chromo][chromo2].count('reverse'):
                    f.write(f'{chromo2} - {chromo} - forward ({oris_others[org][chromo][chromo2].count("forward")} out of {len(oris_others[org][chromo][chromo2])})\n')
                    if chromo2 not in orientations:
                        orientations[chromo2] = {}
                    orientations[chromo2][chromo] = '+'
                else:
                    f.write(f'{chromo2} - {chromo} - reverse ({oris_others[org][chromo][chromo2].count("reverse")} out of {len(oris_others[org][chromo][chromo2])})\n')
                    if chromo2 not in orientations:
                        orientations[chromo2] = {}
                    orientations[chromo2][chromo] = '-'
                tot_counts = len(oris_others[org][chromo][chromo2])
                if oris_others[org][chromo][chromo2].count('forward') > tot_counts*0.1 and oris_others[org][chromo][chromo2].count('forward') < tot_counts*0.9:
                    f2.write(f'{chromo2} - {chromo} - for:  {oris_others[org][chromo][chromo2].count("forward")} , rev: {oris_others[org][chromo][chromo2].count("reverse")}\n')
f2.close()

if len(orders_others) == 0:
    print(f'no syn between {ref_genome} and {targets[0]}')
    exit(0)

simple_map = {}
for t_chromo,bib in mapping_contig_chromo.items():
    r_chromos_with_count = list(bib.items())
    simple_map[t_chromo] = r_chromos_with_count[[x[1] for x in r_chromos_with_count].index(max([x[1] for x in r_chromos_with_count]))]

g_taken = {}
written = {}
new_orders_others = {targets[0]:{}}
contig_appearance_counter = {}

for ref in small_meta[ref_genome][0]:
    if ref not in orders_others[targets[0]]:
        #print(f'{ref_genome} has chromo {ref} which doesnt have anchor alignments to any target chromsomes')
        continue
    new_orders_others[targets[0]][ref] = []
    c = orders_others[targets[0]][ref]
    anything_valid = 0
    new_scaffold_len = 0
    taken = []
    for x in c:
        if x not in contig_appearance_counter:
            contig_appearance_counter[x] = 0
        contig_appearance_counter[x] += 1
        if x not in written and simple_map[x][0] == ref and contig_appearance_counter[x] >= simple_map[x][1]/2:
            anything_valid = 1
            written[x] = 1
            new_orders_others[targets[0]][ref].append(x)
            taken.append(x)
    if anything_valid == 1:
        for x in taken:
            g_taken[x] = 1
               
orders_others = new_orders_others

with open(f'simple_maps/{target_genome}','wb') as f:
    pickle.dump(simple_map,f)

### make main output file with contig order
written = {}
with open(f'orders/{target_genome}','w') as f:
    ref_order = []
    if first == 1:
        ref_order = [(x,'') for x in small_meta[ref_genome][0]]
    else:
        with open(f'orders/{ref_genome}') as f2:
            for line in f2:
                seq, ori = line.strip().split('\t')
                ref_order.append((seq,ori))
    for ref,ori1 in ref_order:
        if ref not in new_orders_others[targets[0]]:
            continue
        c = new_orders_others[targets[0]][ref]
        if len(c) == 0:
            continue
        for x in c:
            if x not in written and simple_map[x][0] == ref:
                ori = orientations[x][ref]
                f.write(f'{x}\t{ori}\n')
                written[x] = 1
