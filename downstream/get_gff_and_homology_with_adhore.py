#! /usr/bin/env python3

import pickle
from statistics import mean,stdev,median
from numpy import quantile
from subprocess import run
from sys import argv

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
    for org,bib in mapping.items(): 
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
        with open('compressed_maps/'+org,'rb') as f:
            am[org] = pickle.load(f)
    
    gff_MCScanX = open('MCScanX.gff','w')
    MCScanX_homology = open('MCScanX.homology','w')
    dialign_homology = open('dialign.homology','w')
    adhore_pairwise = open('adhore_pairwise.table','w')
    adhore_families = open('adhore_families.table','w')

    adhore_bib = {}

    for org in orgs:
        org = org_mapping[org]
        exec(f'smore_{org} = open("smore_anchors/{org}.bed","w+")')

        adhore_bib[org] = {}



    with open(f'clusters','rb') as f:
        clusters = pickle.load(f)
    with open(f'i_bib','rb') as f:
        i_bib = pickle.load(f)
    count_2 = max(list(clusters.keys()))
    for count,cluster in clusters.items():
        count_2 += 1
        elements = []
        cluster_score = cluster['representative_score']
        cluster_length = cluster['representative_length']
        focal_orientation = cluster['focal_orientation']
        del cluster['representative_score']
        del cluster['representative_length']
        del cluster['focal_orientation']
        # for MCScanX
        for org,bib in cluster.items():
            s = bib['matches']
            chromo = am[org][list(s)[0]]['chromosome']
            if len(s) > 1:
                start = 1e15
                end = -1
                for ele in s:
                    if am[org][ele]['start'] < start:
                        start = am[org][ele]['start']
                    if am[org][ele]['end'] > end:
                        end = am[org][ele]['end']
            else:
                ele = list(s)[0]
                start = am[org][ele]['start']
                end = am[org][ele]['end']
            elements.append((org,chromo,start,end,bib['orientation']))
        for enu,ele in enumerate(elements[:-1]):
            org,chromo,start,end,ori = ele
            chromo = chr_mapping[org][chromo]
            org = org_mapping[org]
            ele_name = f'{chromo}ele{start}to{end}'
            gff_MCScanX.write(f'{chromo}\t{ele_name}\t{start}\t{end}\n')
            for ele2 in elements[enu+1:]:
                org2,chromo2,start2,end2,ori2 = ele2
                chromo2 = chr_mapping[org2][chromo2]
                org2 = org_mapping[org2]
                ele2_name = f'{chromo2}ele{start2}to{end2}'
                
                MCScanX_homology.write(f'{ele_name}\t{ele2_name}\t{cluster_score}\n')
        # add last
        org,chromo,start,end,ori = elements[-1]
        chromo = chr_mapping[org][chromo]
        org = org_mapping[org]
        ele_name = f'{chromo}ele{start}to{end}'
        gff_MCScanX.write(f'{chromo}\t{ele_name}\t{start}\t{end}\n')
        
        # for dialign
        
        #for org,bib in cluster.items():
        #    s = bib['matches']
        #    lines = set() 
        #    for j in s:
        #        chromo = am[org][j]['chromosome']
        #        chromo = chr_mapping[org][chromo]
        #        start = am[org][j]['start']
        #        end = am[org][j]['end']
        #        ele_name = f'{chromo}ele{start}to{end}'
        #        for org2,match_bib in am[org][j]['matches'].items():
        #            if org2 not in cluster:
        #                continue
        #            for k,score_hit_coords_ori in match_bib.items():
        #                if k not in cluster[org2]['matches'] or org not in am[org2][k]['matches'] or j not in am[org2][k]['matches'][org]:
        #                    continue
        #                chromo2 = am[org2][k]['chromosome']
        #                chromo2 = chr_mapping[org2][chromo2]
        #                start2 = am[org2][k]['start']
        #                end2 = am[org2][k]['end']
        #                ele2_name = f'{chromo2}ele{start2}to{end2}'
        #                score = int(score_hit_coords_ori[0])
        #                hit_start1,hit_end1 = score_hit_coords_ori[1]
        #                hit_start1,hit_end1 = int(hit_start1),int(hit_end1)
        #                length1 = hit_end1 - hit_start1
        #                hit_start2,hit_end2 = am[org2][k]['matches'][org][j][1]
        #                hit_start2,hit_end2 = int(hit_start2),int(hit_end2)
        #                length2 = hit_end2 - hit_start2
        #                ori = score_hit_coords_ori[2]
        #                length = max(length1,length2)
        #                lines.add((ele_name,chromo,start,end,hit_start1,hit_end1,ele2_name,chromo2,start2,end2,hit_start2,hit_end2,score,length,ori))
        #    for ele_name,chromo,start,end,hit_start1,hit_end1,ele2_name,chromo2,start2,end2,hit_start2,hit_end2,score,length,ori in lines:
        #        dialign_homology.write(f'{ele_name}\t{chromo}\t{start}\t{end}\t{hit_start1}\t{hit_end1}\t{ele2_name}\t{chromo2}\t{start2}\t{end2}\t{hit_start2}\t{hit_end2}\t{score}\t{length}\t{ori}\n')
    
        # for Smore
        anchor_lines = {}
        cluster_score = str(cluster_score)
        cluster_length = str(cluster_length)
        for org,bib in cluster.items():
            anchor_lines[org_mapping[org]] = set()
            s = bib['matches']
            chromo = am[org][list(s)[0]]['chromosome']
            if len(s) > 1:
                start = 1e15
                end = -1
                for ele in s:
                    if am[org][ele]['start'] < start:
                        start = am[org][ele]['start']
                    if am[org][ele]['end'] > end:
                        end = am[org][ele]['end']
            else:
                ele = list(s)[0]
                start = am[org][ele]['start']
                end = am[org][ele]['end']
            elements.append((org,chromo,start,end,bib['orientation']))
        for enu,ele in enumerate(elements[:-1]):
            org,chromo,start,end,ori = ele
            chromo = chr_mapping[org][chromo]
            org = org_mapping[org]
            if ori == 'forward':
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count) + '\t' + str(start) + '\t' + cluster_length + '\t' + '+' + '\t' + get_ref(org) + '\t' + cluster_score + '\n')
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count_2) + '\t' + str(start) + '\t' + cluster_length + '\t' + '-' + '\t' + get_ref(org) + '\t' + cluster_score + '\n')
            elif ori == 'reverse':
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count_2) + '\t' + str(start) + '\t' + cluster_length + '\t' + '+' + '\t' + get_ref(org) + '\t' + cluster_score + '\n')
                anchor_lines[org].add(chromo + '\t' + org + '_' + str(count) + '\t' + str(start) + '\t' + cluster_length + '\t' + '-' + '\t' + get_ref(org) + '\t' + cluster_score + '\n')

        for org,lines in anchor_lines.items():
            for line in lines:
                exec(f'smore_{org}.write(line)')

        # i-AdHoRe
    
        for org,bib in cluster.items():
            s = bib['matches']
            chromo = am[org][list(s)[0]]['chromosome']
            if chr_mapping[org][chromo] not in adhore_bib[org_mapping[org]]:
                adhore_bib[org_mapping[org]][chr_mapping[org][chromo]] = []
            if len(s) > 1:
                start = 1e15
                end = -1
                for ele in s:
                    if am[org][ele]['start'] < start:
                        start = am[org][ele]['start']
                    if am[org][ele]['end'] > end:
                        end = am[org][ele]['end']
            else:
                ele = list(s)[0]
                start = am[org][ele]['start']
                end = am[org][ele]['end']
            elements.append((org,chromo,start,end,bib['orientation']))
        for enu,ele in enumerate(elements[:-1]):
            org,chromo,start,end,ori = ele
            chromo = chr_mapping[org][chromo]
            org = org_mapping[org]
            ele_name = f'{chromo}ele{start}to{end}'
            adhore_families.write(f'{ele_name}\tfamily{count}\n')
            if ori == 'forward':
                strand = '+'
            else:
                strand = '-'
            adhore_bib[org][chromo].append(f'{ele_name}{strand}')
            for ele2 in elements[enu+1:]:
                org2,chromo2,start2,end2,ori2 = ele2
                chromo2 = chr_mapping[org2][chromo2]
                org2 = org_mapping[org2]
                ele2_name = f'{chromo2}ele{start2}to{end2}'
                adhore_pairwise.write(f'{ele_name}\t{ele2_name}\n')
        # add last
        org,chromo,start,end,ori = elements[-1]
        if ori == 'forward':
            strand = '+'
        else:
            strand = '-'
        chromo = chr_mapping[org][chromo]
        org = org_mapping[org]
        ele_name = f'{chromo}ele{start}to{end}'
        adhore_families.write(f'{ele_name}\tfamily{count}\n')
        adhore_bib[org][chromo].append(f'{ele_name}{strand}')


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
blast_table=adhore_families.table

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
table_type=family
multiple_hypothesis_correction=FDR
visualizeGHM=true
visualizeAlignment=true
number_of_threads=1
    """
    adhore_config.write(template_adhore_config)
    adhore_config.close()
                        
    gff_MCScanX.close()
    MCScanX_homology.close()
    dialign_homology.close()
    adhore_families.close()

    adhore_pairwise.close()
