#! /usr/bin/python3

from os import listdir
import sys
sys.path.append('./utils/')
from get_mapping import get_mapping
from pprint import pprint
    

org_mapping,chr_mapping = get_mapping()

groups_m = {}
elements_m = {}
group_id_m = 1

mcscanx_clusters = listdir('MCScanX/clusters')
for c in mcscanx_clusters:
    groups_m[group_id_m] = set()
    with open(f'MCScanX/clusters/{c}') as f:
        for line in f:
            org,name = line.strip().split()
            groups_m[group_id_m].add(name)
            if name not in elements_m:
                elements_m[name] = set()
            elements_m[name].add(group_id_m)
    group_id_m += 1

groups_c = {}
elements_c = {}
group_id_c = 1
custom_clusters = listdir('custom/clusters/50000')
for c in custom_clusters:
    groups_c[group_id_c] = set()
    with open(f'custom/clusters/50000/{c}') as f:
        for line in f:
            if 'species' in line:
                org = line.strip().split()[1]
            # element on chromosome seq453 start 2740 end 2792 strand 1
            if 'element on chromosome' in line:
                chromo = line.strip().split('element on chromosome')[1].split()[0]
                chromo = chr_mapping[org][chromo]
                strand = line.strip().split('element on chromosome')[1].split()[6]
                start = line.strip().split('element on chromosome')[1].split()[2]
                end = line.strip().split('element on chromosome')[1].split()[4]
                if strand == '1':
                    strand = '+'
                else:
                    strand = '-'
                name = f'{chromo}CUSTOMele{start}to{end}{strand}'
                if name not in elements_c:
                    elements_c[name] = set()
                groups_c[group_id_c].add(name)
                elements_c[name].add(group_id_c)
    group_id_c += 1

                    

groups_a = {}
elements_a = {}
group_id_a = 1
segments_genes = {}
segments = {}
segments_multiplicons = {}
genes_multiplicons = {}

first = 1
with open('i-ADHoRe/out_only_relevant/list_elements.txt') as f:
    for line in f:
        if first == 1:
            first = 0
            continue
        # id  segment gene    position    orientation
        id,segment,gene,position,orientation = line.strip().split()
        gene = gene+orientation
        if 'CUSTOM' in gene:
            if gene not in segments_genes:
                segments_genes[gene] = set()
            if segment not in segments:
                segments[segment] = set()
            segments_genes[gene].add((segment,position))
            segments[segment].add(gene)

first = 1
with open('i-ADHoRe/out_only_relevant/segments.txt') as f:
    for line in f:
        if first == 1:
            first = 0
            continue
        # id  multiplicon genome  list    first   last    order
        id,multiplicon,genome,list,first,last,order = line.strip().split()
        if id in segments:
            if multiplicon not in genes_multiplicons:
                genes_multiplicons[multiplicon] = set()
                segments_multiplicons[multiplicon] = set()
            segments_multiplicons[multiplicon].add(id)
            genes_multiplicons[multiplicon].update(segments[id])


orthologies_a = set()
for m,s in genes_multiplicons.items():
    pos = {}
    for gene in s:
        for segment,position in segments_genes[gene]:
            if segment in segments_multiplicons[m]:
                if gene in pos:
                    print(gene,pos[gene],position)
                    input()
                pos[gene] = position
    for gene,p in pos.items():
        for gene2,p2 in pos.items():
            if gene != gene2:# and p==p2:
                orthologies_a.add(tuple(sorted([gene,gene2])))
                elements_a[gene] = 1
                elements_a[gene2] = 1

# EVAL
pw_rel_m = set()
pw_left_mc = set()
pw_left_ma = set()
pw_rel_c = set()
pw_left_cm = set()
pw_left_ca = set()
pw_rel_a = orthologies_a
pw_left_am = set()
pw_left_ac = set()

for g,s in groups_m.items():
    for ele in s:
        org1 = ele.split('chr')[0]
        if ele not in elements_c:
            for ele2 in s:
                org2 = ele2.split('chr')[0]
                if ele != ele2 and org1 != org2:
                    pw_left_mc.add(tuple(sorted([ele,ele2])))
        if ele not in elements_a:
            for ele2 in s:
                org2 = ele2.split('chr')[0]
                if ele != ele2 and org1 != org2:
                    pw_left_ma.add(tuple(sorted([ele,ele2])))
        else:
            for ele2 in s:
                org2 = ele2.split('chr')[0]
                if ele != ele2 and org1 != org2:
                    pw_rel_m.add(tuple(sorted([ele,ele2])))
for g,s in groups_c.items():
    for ele in s:
        org1 = ele.split('chr')[0]
        if ele not in elements_m:
            for ele2 in s:
                org2 = ele2.split('chr')[0]
                if ele != ele2 and org1 != org2:
                    pw_left_cm.add(tuple(sorted([ele,ele2])))
        if ele not in elements_a:
            for ele2 in s:
                org2 = ele2.split('chr')[0]
                if ele != ele2 and org1 != org2:
                    pw_left_ca.add(tuple(sorted([ele,ele2])))
        else:
            for ele2 in s:
                org2 = ele2.split('chr')[0]
                if ele != ele2 and org1 != org2:
                    pw_rel_c.add(tuple(sorted([ele,ele2])))
for p in orthologies_a:
    org1 = p[0].split('chr')[0]
    org2 = p[1].split('chr')[0]
    if org1 == org2:
        continue
    if p[0] not in elements_m or p[1] not in elements_m:
        pw_left_am.add(tuple(sorted([ele,ele2])))
    if p[0] not in elements_c or p[1] not in elements_c:
        pw_left_ac.add(tuple(sorted([ele,ele2])))


print(len(pw_rel_m.intersection(pw_rel_c)))
print(len(pw_rel_m.intersection(pw_rel_a)))
print(len(pw_rel_c.intersection(pw_rel_a)))
print(len(pw_rel_m),len(pw_rel_c),len(orthologies_a))
print(len(pw_left_mc),len(pw_left_ma))
print(len(pw_left_cm),len(pw_left_ca))
print(len(pw_left_am),len(pw_left_ac))
print(len(elements_m),len(elements_c),len(elements_a))

#pprint(pw_rel_m - pw_rel_m.intersection(pw_rel_c) - pw_rel_m.intersection(pw_left_mc))
#pprint(pw_rel_c - pw_rel_m.intersection(pw_rel_c) - pw_rel_c.intersection(pw_left_cm))

#adhore_clusters = listdir('i-ADHoRe/
#for ele,ids in elements_m.items():
#    if len(ids) == 1:
#        for x in ids:
#            pprint(groups_m[x])
#
#for ele,ids in elements_c.items():
#    if len(ids) == 1:
#        for x in ids:
#            pprint(groups_c[x])
#for ele,ids in elements_a.items():
#    if len(ids) == 1:
#        for x in ids:
#            pprint(groups_a[x])
