#! /usr/bin/env python3

import pickle
from statistics import mean,stdev,median
from numpy import quantile
from subprocess import run

def get_ref(org):
    global ref
    if org == ref:
        return 'True'
    else:
        return 'False'
    
if __name__ == "__main__":
    
    with open('ref_org') as f:
        for line in f:
            ref = line.strip()
    
    with open('../../462/mapping_pickle','rb') as f:
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

    orgs = []
    with open('orgs') as f:
        for line in f:
            orgs.append(line.strip())

    am = {}
    for org in orgs:
        with open(f'../../462/compressed_maps_multis_to_one/{org}','rb') as f:
            am[org] = pickle.load(f)
    print('LOADING DONE')
    coords = set()
    alignments = {}
    for org,bib in am.items():
        for i,bib2 in bib.items():
            chromo = chr_mapping[org][bib2['chromosome']]
            start = bib2['start']
            end = bib2['end']
            ele_name = f'{chromo}ele{start}to{end}'
            for org2,bib3 in bib2['matches'].items():
                if org2 not in orgs:
                    continue
                if len(bib3) > 1:
                    print(org,i)
                    pprint(bib3)
                for j,mbib in bib3.items():
                    if j not in am[org2] or org not in am[org2][j]['matches'] or i not in am[org2][j]['matches'][org]:
                        continue
                    chromo2 = chr_mapping[org2][am[org2][j]['chromosome']]
                    start2 = am[org2][j]['start']
                    end2 = am[org2][j]['end']
                    score = mbib[0]
                    length = mbib[1][1] - mbib[1][0]
                    ori = mbib[2]
                    ele2_name = f'{chromo2}ele{start2}to{end2}'
                    if org_order.index(org) < org_order.index(org2):
                        if ele_name not in alignments:
                            alignments[ele_name] = {}
                        if ele2_name not in alignments[ele_name]:
                            alignments[ele_name][ele2_name] = set()
                        alignments[ele_name][ele2_name].add((score,length,ori))
                    else:
                        if ele2_name not in alignments:
                            alignments[ele2_name] = {}
                        if ele_name not in alignments[ele2_name]:
                            alignments[ele2_name][ele_name] = set()
                        alignments[ele2_name][ele_name].add((score,length,ori))
                    coords.add(f'{chromo}\t{ele_name}\t{start}\t{end}\n')
                    coords.add(f'{chromo2}\t{ele2_name}\t{start2}\t{end2}\n')

    with open('MCScanX.gff','w') as f:
        for ele in coords:
            f.write(ele)
    with open('MCScanX.homology','w') as f:
        for ele_name,bib in alignments.items():
            for ele2_name,s in bib.items():
                maxi = max([x[0] for x in s])
                score = [x[0] for x in s if x[0] == maxi][0]
                f.write(f'{ele_name}\t{ele2_name}\t{score}\n')

    # for Smore

    run(f'mkdir -p smore_anchors',shell=True)
    for org in orgs:
        org = org_mapping[org]
        exec(f'smore_{org} = open("smore_anchors/{org}.bed","w+")')
    
    count = 0
    count_2 = 0
    anchor_lines = {}
    for org in orgs:
        org = org_mapping[org]
        anchor_lines[org] = set()
    for ele_name,bib in alignments.items():
        org = ele_name.split('ele')[0].split('chr')[0]
        chromo = ele_name.split('ele')[0]
        start = ele_name.split('ele')[1].split('to')[0]
        start = int(start)
        end = ele_name.split('ele')[1].split('to')[1]
        end = int(end)
        length = end - start
        length = str(length)
        for ele2_name,s in bib.items():
            org2 = ele2_name.split('ele')[0].split('chr')[0]
            chromo2 = ele2_name.split('ele')[0]
            start2 = ele2_name.split('ele')[1].split('to')[0]
            start2 = int(start2)
            end2 = ele2_name.split('ele')[1].split('to')[1]
            end2 = int(end2)
            length2 = end2 - start2
            length2 = str(length2)
            count += 1
            count_2 += 1
            maxi = max([x[0] for x in s])
            score = [x for x in s if x[0] == maxi][0][0]
            score = str(score)
            # length below would get actual max length of both alignments...but most likely smore originally worked on lengths o-f MAF blocks or sth like that, so thats why im putting actual candidate lengths
            #length = [x for x in s if x[0] == maxi][0][1]
            ori = [x for x in s if x[0] == maxi][0][2]
            if ori == 'forward':
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count) + '\t' + str(start) + '\t' + length + '\t' + '+' + '\t' + get_ref(org) + '\t' + score + '\n')
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count_2) + '\t' + str(start) + '\t' + length + '\t' + '-' + '\t' + get_ref(org) + '\t' + score + '\n')
                anchor_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count) + '\t' + str(start2) + '\t' + length2 + '\t' + '+' + '\t' + get_ref(org2) + '\t' + score + '\n')
                anchor_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count_2) + '\t' + str(start2) + '\t' + length2 + '\t' + '-' + '\t' + get_ref(org2) + '\t' + score + '\n')
            elif ori == 'reverse':
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count_2) + '\t' + str(start) + '\t' + length + '\t' + '+' + '\t' + get_ref(org) + '\t' + score + '\n')
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count) + '\t' + str(start) + '\t' + length + '\t' + '-' + '\t' + get_ref(org) + '\t' + score + '\n')
                anchor_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count_2) + '\t' + str(start2) + '\t' + length2 + '\t' + '+' + '\t' + get_ref(org2) + '\t' + score + '\n')
                anchor_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count) + '\t' + str(start2) + '\t' + length2 + '\t' + '-' + '\t' + get_ref(org2) + '\t' + score + '\n')


    for org,lines in anchor_lines.items():
        for line in lines:
            exec(f'smore_{org}.write(line)')

    # i-AdHoRe
    adhore_pairwise = open('adhore_pairwise.table','w')

    adhore_bib = {}
    for org in orgs:
        org = org_mapping[org]
        adhore_bib[org] = {}
    count = 0
    count_2 = 0
    anchor_lines = {}
    for org in orgs:
        anchor_lines[org] = set()
    for ele_name,bib in alignments.items():
        org = ele_name.split('ele')[0].split('chr')[0]
        chromo = ele_name.split('ele')[0]
        if chromo not in adhore_bib[org]:
            adhore_bib[org][chromo] = set()
        # randomly assigning forward orientation to the first elements which are randomly assigned order according to org_order
        strand = '+'
        adhore_bib[org][chromo].add(f'{ele_name}-')
        adhore_bib[org][chromo].add(f'{ele_name}+')
        start = ele_name.split('ele')[1].split('to')[0]
        start = int(start)
        end = ele_name.split('ele')[1].split('to')[1]
        end = int(end)
        length = end - start
        for ele2_name,s in bib.items():
            org2 = ele2_name.split('ele')[0].split('chr')[0]
            chromo2 = ele2_name.split('ele')[0]
            if chromo2 not in adhore_bib[org2]:
                adhore_bib[org2][chromo2] = set()
            start2 = ele2_name.split('ele')[1].split('to')[0]
            start2 = int(start2)
            end2 = ele2_name.split('ele')[1].split('to')[1]
            end2 = int(end2)
            length2 = end2 - start2
            count += 1
            count_2 += 1
            maxi = max([x[0] for x in s])
            score = [x for x in s if x[0] == maxi][0][0]
            # length below would get actual max length of both alignments...but most likely smore originally worked on lengths o-f MAF blocks or sth like that, so thats why im putting actual candidate lengths
            #length = [x for x in s if x[0] == maxi][0][1]
            ori = [x for x in s if x[0] == maxi][0][2]
            if ori == 'forward':
                strand = '+'
            else:
                strand = '-'
            adhore_bib[org2][chromo2].add(f'{ele2_name}+')
            adhore_bib[org2][chromo2].add(f'{ele2_name}-')
            adhore_pairwise.write(f'{ele_name}\t{ele2_name}\n')


    adhore_config = open('adhore_config','w')
    for org,bib in adhore_bib.items():
        run(f'mkdir -p adhore_gene_lists/{org}',shell=True)
        adhore_config.write(f'genome={org}\n')
        for chromo,l in bib.items():
            adhore_config.write(f'{chromo} adhore_gene_lists/{org}/{chromo}.lst\n')
            l = list(set(l))
            exec(f'file_{chromo} = open("adhore_gene_lists/{org}/{chromo}.lst","w")')
            for ele in l:
                ele = ele + '\n'
                exec(f'file_{chromo}.write(ele)')
            exec(f"file_{chromo}.close()")

    #for org in orgs:
    #    org = org_mapping[org]
    #    exec(f"smore_{org}.close()")

    template_adhore_config = """
blast_table=adhore_pairwise.table

output_path=out
cluster_type=hybrid

gap_size= 50
cluster_gap= 70
max_gaps_in_alignment=100
cloud_gap_size=30
cloud_cluster_gap=50

q_value=0.7
prob_cutoff=0.01
anchor_points=3
alignment_method=gg2
level_2_only=false
multiple_hypothesis_correction=FDR
visualizeGHM=true
visualizeAlignment=true
number_of_threads=1
    """
    adhore_config.write(template_adhore_config)
    adhore_config.close()
    adhore_pairwise.close()
