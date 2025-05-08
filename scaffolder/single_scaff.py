#! /usr/bin/env python3

from sys import argv
from subprocess import run
import pickle
from os.path import isfile
import sys

ref_genome = argv[1]

small_meta_path = '../utils/small_meta/'
anchors_path = '../anchors/aligned/'
anchors_path_cand = '../anchors/candidates/'
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


mapping_contig_chromo = {}
with open('singles_out/contig_number_mapping.pickle','rb') as f:
    contig_number_mapping = pickle.load(f)

if target_genome == 'human_chopped':
    ref_genome_name = ref_genome+'_vs_human_chopped'
else:
    ref_genome_name = ref_genome

with open(anchors_path+f'/{ref_genome_name}','rb') as f:
    anchors = pickle.load(f)
with open(anchors_path_cand+f'/{target_genome}','rb') as f:
    t_anchors = pickle.load(f)

iss = sorted(list(anchors.keys()))
for i in iss:
    chromo = anchors[i]['chromosome'] 
    if chromo not in orders_others[target_genome]:
        orders_others[target_genome][chromo] = []
        oris_others[target_genome][chromo] = {}
    if target_genome in anchors[i]['matches'] and anchors[i]['matches'][target_genome]['meta']['multiple matches out of tolerance range'] == 0:
        for j,match_bib in anchors[i]['matches'][target_genome]['matches'].items():
            chromo2 = t_anchors[j]['chromosome'] 
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
        print(f'{ref_genome} has chromo {ref} which doesnt have anchor alignments to any target chromsomes')
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

limited = {}
with open(f'singles_out/divergently_aligned_contigs_from_AncST_ref_{ref_genome}.txt','w') as f:
    for contig,bib in mapping_contig_chromo.items():
        if len(bib) > 1:
            tot_count = sum(bib.values())
            good_matches = []
            for chromo,count in bib.items():
                if count >= tot_count*0.2:
                    good_matches.append((chromo,count))
            if len(good_matches) > 1:
                f.write(f'{contig}: [')
                for chromo,count in good_matches:
                    f.write(f'{chromo}: {count},')
                f.write(']\n')
            else:
                limited[contig] = 1
        else:
            limited[contig] = 1

### make main output file with contig order
count = 1
written = {}
with open(f'singles_out/contig_order_target_{target_genome}_ref_{ref_genome}.out','w') as f:
    for ref in small_meta[ref_genome][0]:
        if ref not in new_orders_others[targets[0]]:
            continue
        f.write(f'>Scaffold_{count}\n')
        c = new_orders_others[targets[0]][ref]
        for x in c:
            if x not in written and simple_map[x][0] == ref:
                ori = orientations[x][ref]
                f.write(f'{x}\t{ori}\n')
                written[x] = 1
        count += 1

# make Blossom V input
count = 1
written = {}
no_v = 0
no_e = 0
for ref in small_meta[ref_genome][0]:
    if ref not in new_orders_others[targets[0]]:
        continue
    c = new_orders_others[targets[0]][ref]
    loc_no_v = 0
    for enu,x in enumerate(c[:-1]):
        if simple_map[x][0] == ref and simple_map[c[enu+1]][0] == ref:
            loc_no_v += 2 # two for each contig
    # last contig is missing from count but thats alright since one would have to deduct 2 for the start and end anyway
    no_v += loc_no_v
    no_e += int(loc_no_v/2)
with open(f'singles_out/blossom_v_infile_target_{target_genome}_ref_{ref_genome}.txt','w') as f:
    f.write(f'{no_v} {no_e}\n')
    for ref in small_meta[ref_genome][0]:
        if ref not in new_orders_others[targets[0]]:
            continue
        c = new_orders_others[targets[0]][ref]
        for enu,x in enumerate(c[:-1]):
            if simple_map[x][0] == ref and simple_map[c[enu+1]][0] == ref and x in limited and c[enu+1] in limited:
                ori = orientations[x][ref]
                number1 = contig_number_mapping[x]
                head1 = number1 + '1'
                tail1 = number1 + '2'
                if ori == '-':
                    head1,tail1 = tail1,head1
                number2 = contig_number_mapping[c[enu+1]]
                head2 = number2 + '1'
                tail2 = number2 + '2'
                ori = orientations[c[enu+1]][ref]
                if ori == '-':
                    head2,tail2 = tail2,head2
                f.write(f'{tail1} {head2} 1\n')
        count += 1


