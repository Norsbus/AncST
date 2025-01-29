#! /usr/bin/env python3

import pickle
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

    orgs = []
    with open('orgs') as f:
        for line in f:
            orgs.append(line.strip())

    am = {}
    for org in orgs:
        with open(f'compressed_maps_multis_to_one/{org}','rb') as f:
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
                    hit_coords1 = mbib[1]
                    hit_coords2 = am[org2][j]['matches'][org][i][1]
                    ele2_name = f'{chromo2}ele{start2}to{end2}'
                    if org_order.index(org) < org_order.index(org2):
                        if ele_name not in alignments:
                            alignments[ele_name] = {}
                        if ele2_name not in alignments[ele_name]:
                            alignments[ele_name][ele2_name] = set()
                        alignments[ele_name][ele2_name].add((score,length,ori,(hit_coords1[0],hit_coords1[1]),(hit_coords2[0],hit_coords2[1])))
                    else:
                        if ele2_name not in alignments:
                            alignments[ele2_name] = {}
                        if ele_name not in alignments[ele2_name]:
                            alignments[ele2_name][ele_name] = set()
                        alignments[ele2_name][ele_name].add((score,length,ori,(hit_coords1[0],hit_coords1[1]),(hit_coords2[0],hit_coords2[1])))
                    coords.add(f'{chromo}\t{ele_name}\t{start}\t{end}\n')
                    coords.add(f'{chromo2}\t{ele2_name}\t{start2}\t{end2}\n')

    
    with open('MCScanX.gff','w+') as f:
        f.writelines(coords)

    mc_lines = set()

    di_lines = set()

    count = 0
    count_2 = 0
    smore_lines = {}
    for org in orgs:
        org = org_mapping[org]
        smore_lines[org] = set()
    
    adhore_bib = {}
    for org in orgs:
        org = org_mapping[org]
        adhore_bib[org] = {}
    adhore_pairwise = set()


    for ele_name,bib in alignments.items():
        org = ele_name.split('ele')[0].split('chr')[0]
        chromo = ele_name.split('ele')[0]
        if chromo not in adhore_bib[org]:
            adhore_bib[org][chromo] = set()
        start = ele_name.split('ele')[1].split('to')[0]
        start = int(start)
        end = ele_name.split('ele')[1].split('to')[1]
        end = int(end)
        adhore_bib[org][chromo].add((start,f'{ele_name}+'))
        adhore_bib[org][chromo].add((start,f'{ele_name}-'))
        # SMORE : length below would get actual max length of both alignments...but most likely smore originally worked on lengths o-f maf blocks or sth like that, so thats why im putting actual candidate lengths
        cand_length1 = end - start
        cand_length1 = str(cand_length1)
        for ele2_name,s in bib.items():
            count += 1
            count_2 += 1
            org2 = ele2_name.split('ele')[0].split('chr')[0]
            chromo2 = ele2_name.split('ele')[0]
            if chromo2 not in adhore_bib[org2]:
                adhore_bib[org2][chromo2] = set()
            start2 = ele2_name.split('ele')[1].split('to')[0]
            start2 = int(start2)
            end2 = ele2_name.split('ele')[1].split('to')[1]
            end2 = int(end2)
            adhore_bib[org2][chromo2].add((start2,f'{ele2_name}+'))
            adhore_bib[org2][chromo2].add((start2,f'{ele2_name}-'))
            maxi = max([x[0] for x in s])
            score = [x for x in s if x[0] == maxi][0][0]
            score = str(score)
            # SMORE : length below would get actual max length of both alignments...but most likely smore originally worked on lengths o-f maf blocks or sth like that, so thats why im putting actual candidate lengths
            cand_length2 = end2 - start2
            cand_length2 = str(cand_length2)
            ori = [x for x in s if x[0] == maxi][0][2]
            hit_start1,hit_end1 = [x for x in s if x[0] == maxi][0][3]
            hit_start1,hit_end1 = int(hit_start1),int(hit_end1)
            hit_length1 = hit_end1 - hit_start1
            hit_start2,hit_end2 = [x for x in s if x[0] == maxi][0][4]
            hit_start2,hit_end2 = int(hit_start2),int(hit_end2)
            hit_length2 = hit_end2 - hit_start2

            # DIALIGN : originally they work on ungapped alignment lengths, so closest i get here is longest actual gapped alignment length
            hit_length = max(hit_length1,hit_length2)
            
            mc_lines.add(f'{ele_name}\t{ele2_name}\t{score}\n')

            di_lines.add(f'{ele_name}\t{chromo}\t{start}\t{end}\t{hit_start1}\t{hit_end1}\t{ele2_name}\t{chromo2}\t{start2}\t{end2}\t{hit_start2}\t{hit_end2}\t{score}\t{hit_length}\t{ori}\n')
            
            if ori == 'forward':
                # for adhore
                strand = '+'
                smore_lines[org].add(chromo + '\t' + org + '_' + str(count) + '\t' + str(start) + '\t' + cand_length1 + '\t' + '+' + '\t' + get_ref(org) + '\t' + score + '\n')
                smore_lines[org].add(chromo + '\t' + org + '_' + str(count_2) + '\t' + str(start) + '\t' + cand_length1 + '\t' + '-' + '\t' + get_ref(org) + '\t' + score + '\n')
                smore_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count) + '\t' + str(start2) + '\t' + cand_length2 + '\t' + '+' + '\t' + get_ref(org2) + '\t' + score + '\n')
                smore_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count_2) + '\t' + str(start2) + '\t' + cand_length2 + '\t' + '-' + '\t' + get_ref(org2) + '\t' + score + '\n')
            elif ori == 'reverse':
                # for adhore
                strand = '-'
                smore_lines[org].add(chromo + '\t' + org + '_' + str(count_2) + '\t' + str(start) + '\t' + cand_length1 + '\t' + '+' + '\t' + get_ref(org) + '\t' + score + '\n')
                smore_lines[org].add(chromo + '\t' + org + '_' + str(count) + '\t' + str(start) + '\t' + cand_length1 + '\t' + '-' + '\t' + get_ref(org) + '\t' + score + '\n')
                smore_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count_2) + '\t' + str(start2) + '\t' + cand_length2 + '\t' + '+' + '\t' + get_ref(org2) + '\t' + score + '\n')
                smore_lines[org2].add(chromo2 + '\t' + org2 + '_' + str(count) + '\t' + str(start2) + '\t' + cand_length2 + '\t' + '-' + '\t' + get_ref(org2) + '\t' + score + '\n')

            adhore_pairwise.add(f'{ele_name}\t{ele2_name}\n')


    with open('MCScanX.homology','w+') as f:
        f.writelines(mc_lines)

    with open('dialign.homology','w+') as f:
        f.writelines(di_lines)

    run(f'rm -rf smore_anchors && mkdir -p smore_anchors',shell=True)
    for org in orgs:
        org = org_mapping[org]
        with open(f"smore_anchors/{org}.bed","w+") as f:
            f.writelines(smore_lines[org])
    
    with open('adhore_pairwise.table','w+') as f:
        f.writelines(adhore_pairwise)

    first = 1
    adhore_config = open('adhore_config','w+')
    for org,bib in adhore_bib.items():
        run(f'rm -rf adhore_gene_lists/{org} && mkdir -p adhore_gene_lists/{org}',shell=True)
        if first == 1:
            space = ''
            first = 0
        else:
            space = '\n'
        adhore_config.write(f'{space}genome={org}\n')
        for chromo,l in bib.items():
            adhore_config.write(f'{chromo} adhore_gene_lists/{org}/{chromo}.lst\n')
            l = [x[1]+'\n' for x in sorted(l)]
            with open(f"adhore_gene_lists/{org}/{chromo}.lst","w+") as f:
                f.writelines(l)

    template_adhore_config = """
blast_table=adhore_pairwise.table

output_path=out
cluster_type=collinear

gap_size=20
cluster_gap=20
max_gaps_in_alignment=40
tandem_gap=2

q_value=0.7
prob_cutoff=0.01
anchor_points=3
alignment_method=gg2
multiple_hypothesis_correction=FDR
visualizeGHM=false
visualizeAlignment=true
number_of_threads=1
    """
    adhore_config.write(template_adhore_config)
    adhore_config.close()
