#! /usr/bin/env python

import networkx as nx
import pickle
import sys

target_genome = sys.argv[1]

with open(f'singles_out/blossom_mapping.pickle','rb') as f:
    mapping_for_reconstruction = pickle.load(f)

with open(f'singles_out/contig_number_mapping.pickle','rb') as f:
    contig_number_mapping = pickle.load(f)

original_edges = {}

with open(f'singles_out/blossom_v_infile_merged.txt') as f:
    first = 1
    for line in f:
        if first == 1:
            first = 0
            continue
        n1,n2,s = line.strip().split()
        if int(n1) > int(n2):
            loc_n1 = n2
            loc_n2 = n1
        else:
            loc_n1 = n1
            loc_n2 = n2
        original_edges[(loc_n1,loc_n2)] = int(s)

G = nx.Graph()
combine = set()
with open(f'blossom5-v2.05.src/out') as f:
    first = 1
    for line in f:
        if first == 1:
            first = 0
            continue
        n1,n2 = line.strip().split()
        if int(n1) > int(n2):
            loc_n1 = n2
            loc_n2 = n1
        else:
            loc_n1 = n1
            loc_n2 = n2
        weight = original_edges[(loc_n1,loc_n2)]
        loc_n1 = mapping_for_reconstruction[int(loc_n1)] # int because in merge....py count was usedas the variable which is int but then strings are obviously written to the file
        loc_n2 = mapping_for_reconstruction[int(loc_n2)]
        G.add_node(loc_n1)
        G.add_node(loc_n2)
        combine.add('_'.join(loc_n1.split('_')[:-1]))
        combine.add('_'.join(loc_n2.split('_')[:-1]))
        if weight != 0:
            G.add_edge(loc_n1,loc_n2,weight=weight)

for x in combine:
    G.add_edge(x+'_head',x+'_tail',weight='headtailedge')

sc = [x for x in nx.simple_cycles(G)]

while len(sc) > 0:
    for c in sc:
        c += [c[0]]
        weights = []
        for enu,n1 in enumerate(c[:-1]):
            n2 = c[enu+1]
            if G[n1][n2]['weight'] == 'headtailedge':
                continue
            else:
                weights.append(G[n1][n2]['weight'])
        min_weight = max(weights)
        for enu,n1 in enumerate(c[:-1]):
            n2 = c[enu+1]
            if G[n1][n2]['weight'] == min_weight:
                G.remove_edge(n1,n2)
    sc = [x for x in nx.simple_cycles(G)]


count = 1
degrees = {node:val for (node, val) in G.degree()}
g_taken = {}
with open(f'multi_out/scaffolds_numbers.out','w') as f:
    for cc in sorted(nx.connected_components(G), key = len, reverse=True):
        if len(cc) < 4:
            continue
        ends = []
        for node in cc:
            if degrees[node] == 1:
                ends.append(node)
        if len(ends) == 2:
            paths = [p for p in nx.all_simple_paths(G,ends[0],ends[1])]
            if len(paths) > 1:
                continue
            f.write(f'>scaffold_{count}\n')
            taken = []
            count += 1
            last_contig = contig = 'contig'+paths[0][0]
            if 'head' in contig:
                f.write(f'{"_".join(contig.split("_")[:-1])} +\n')
            else:
                f.write(f'{"_".join(contig.split("_")[:-1])} -\n')
            taken.append("_".join(contig.split("_")[:-1]))
            same = 1
            for contig in paths[0][1:]:
                contig = 'contig'+contig
                if same == 1:
                    same = 0
                    continue
                else:
                    same = 1
                if 'head' in contig:
                    f.write(f'{"_".join(contig.split("_")[:-1])} +\n')
                else:
                    f.write(f'{"_".join(contig.split("_")[:-1])} -\n')
                taken.append("_".join(contig.split("_")[:-1]))
                last_contig = contig
            for gt in taken:
                g_taken[gt] = 1
        else:
            continue

out2 = open(f'multi_out/scaffolds_names.out','w')
with open(f'multi_out/scaffolds_numbers.out','r') as f:
    for line in f:
        if '>scaffold' in line:
            out2.write(line)
        elif 'contig' in line:
            name = contig_number_mapping[line.split('contig')[1].strip().split()[0].strip()]
            if '+' in line[-5:]:
                out2.write(f'{name}+\n')
            else:
                out2.write(f'{name}-\n')
out2.close()
