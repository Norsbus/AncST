#! /usr/bin/env python3

import pickle
from sys import argv
from pprint import pprint

morgs = []
with open('orgs') as f:
    for line in f:
        morgs.append(line.strip())

with open('mapping_pickle','rb') as f:
    mapping = pickle.load(f)

chr_mapping = {}
org_mapping = {}
org_order = []
for org,bib in mapping.items(): 
    org_order.append(org)
    org_mapping[org] = bib['org_id']
    org_mapping[bib['org_id']] = org
    chr_mapping[org] = {}
    chr_mapping[bib['org_id']] = {}
    for chromo,other_chromo in bib['chromosomes'].items():
        chr_mapping[org][chromo] = other_chromo
        chr_mapping[org][other_chromo] = chromo
        chr_mapping[bib['org_id']][chromo] = other_chromo
        chr_mapping[bib['org_id']][other_chromo] = chromo

am = {}
for org in morgs:
    with open(f'compressed_maps_multis_to_one/{org}','rb') as f:
        am[org] = pickle.load(f)

coords = set()
alignments = {}
c = 0
for org,bib in am.items():
    for i,bib2 in bib.items():
        chromo = chr_mapping[org][bib2['chromosome']]
        start = bib2['start']
        end = bib2['end']
        ele_name = f'{chromo}ele{start}to{end}'
        for org2,bib3 in bib2['matches'].items():
            if org2 not in morgs:
                continue
            if len(bib3) > 1:
                print(org,i)
                pprint(bib3)
            for j,mbib in bib3.items():
                if j not in am[org2] or org not in am[org2][j]['matches'] or i not in am[org2][j]['matches'][org]:
                    c += 1
                    continue
                chromo2 = chr_mapping[org2][am[org2][j]['chromosome']]
                start2 = am[org2][j]['start']
                end2 = am[org2][j]['end']
                score = mbib[0]
                ele2_name = f'{chromo2}ele{start2}to{end2}'
                #if (ele2_name in ['246orgchr11ele5565230to5567530','246orgchr8ele14618700to14639000'] and ele_name == '151orgchr1ele10029750to10044000') or (ele2_name == '151orgchr1ele10029750to10044000' and ele_name in ['246orgchr11ele5565230to5567530','246orgchr8ele14618700to14639000']):
                #    print(org,i,org2,j)
                #    pprint(bib2)
                #    pprint(mbib)
                #    input()
                if org_order.index(org) < org_order.index(org2):
                    if ele_name not in alignments:
                        alignments[ele_name] = {}
                    if ele2_name not in alignments[ele_name]:
                        alignments[ele_name][ele2_name] = set()
                    alignments[ele_name][ele2_name].add(score)
                else:
                    if ele2_name not in alignments:
                        alignments[ele2_name] = {}
                    if ele_name not in alignments[ele2_name]:
                        alignments[ele2_name][ele_name] = set()
                    alignments[ele2_name][ele_name].add(score)
                coords.add(f'{chromo}\t{ele_name}\t{start}\t{end}\n')
                coords.add(f'{chromo2}\t{ele2_name}\t{start2}\t{end2}\n')

print(c)
with open('pairwise_multis_to_one.gff','w') as f:
    for ele in coords:
        f.write(ele)
with open('pairwise_multis_to_one.homology','w') as f:
    for ele_name,bib in alignments.items():
        for ele2_name,s in bib.items():
            f.write(f'{ele_name}\t{ele2_name}\t{max(s)}\n')

#alignments = sorted(list(alignments))
#to_del = set()
#to_add = []
#for enu,x in enumerate(alignments[:-1]):
#    e11,e21,s1 = x.strip().split('\t')
#    s1 = float(s1)
#    maxi = s1
#    double = 0
#    for y in alignments[enu+1:]:
#        e12,e22,s2 = x.strip().split('\t')
#        s2 = float(s2)
#        if e12 != e11 and e22 != e21:
#            break
#        else:
#            double = 1
#            if s1 > maxi:
#                maxi = s1
#            if s2 > maxi:
#                maxi = s2
#            to_del.add(x)
#            to_del.add(y)
#    if double == 1:
#        to_add.append(f'{e11}\t{e22}\t{maxi}\n')
#
#for x in to_del:
#    alignments.remove(x)

#with open('MCScanX.homology_pairwise','w') as f:
#    for l in alignments:
#        f.write(l)
    #for l in to_add:
    #    f.write(l)
