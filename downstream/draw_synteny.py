#! /usr/bin/env python3


import sys
sys.path.append('../utils/')
from get_mapping import get_mapping

import pickle
from random import randint
from matplotlib.pyplot import cm
from pygenomeviz import GenomeViz
import numpy as np
from subprocess import run
from bisect import bisect_left,bisect_right

def get_gff():

    gff = {}
    with open('../MCScanX/MCScanX.gff') as f:
        for line in f:
            chromo,name,start,end = line.strip().split()
            gff[name] = (chromo,start,end)

    return gff

def get_homology(gff):

    clusters = {}
    alloc = {}
    cluster_id = 1

    homology = {}
    with open('../MCScanX/MCScanX.homology') as f:
        for line in f:
            name1,name2,score = line.strip().split()
            chromo1,start1,end1 = gff[name1]
            chromo2,start2,end2 = gff[name2]
            org1 = org_mapping[name1.split('ele')[0].split('chr')[0]]
            org2 = org_mapping[name2.split('ele')[0].split('chr')[0]]
            chromo1 = chr_mapping[org1][chromo1]
            chromo2 = chr_mapping[org2][chromo2]
            start1 = int(start1)
            start2 = int(start2)
            end1 = int(end1)
            end2 = int(end2)
            
            if org1 not in homology:
                homology[org1] = {}
            if chromo1 not in homology[org1]:
                homology[org1][chromo1] = {}
            if start1 not in homology[org1][chromo1]:
                homology[org1][chromo1][start1] = {'end':end1,'matches':{},'score':int(score)}
            homology[org1][chromo1][start1]['matches'][org2] = (chromo2,start2,int(score))

            if org2 not in homology:
                homology[org2] = {}
            if chromo2 not in homology[org2]:
                homology[org2][chromo2] = {}
            if start2 not in homology[org2][chromo2]:
                homology[org2][chromo2][start2] = {'end':end2,'matches':{},'score':int(score)}
            homology[org2][chromo2][start2]['matches'][org1] = (chromo1,start1,int(score))

            if org1 not in alloc:
                alloc[org1] = {}
            if chromo1 not in alloc[org1]:
                alloc[org1][chromo1] = {}
            if org2 not in alloc:
                alloc[org2] = {}
            if chromo2 not in alloc[org2]:
                alloc[org2][chromo2] = {}
            if start1 not in alloc[org1][chromo1] and start2 not in alloc[org2][chromo2]:
                clusters[cluster_id] = {}
                clusters[cluster_id][org1] = {}
                clusters[cluster_id][org2] = {}
                clusters[cluster_id][org1][chromo1] = set()
                clusters[cluster_id][org2][chromo2] = set()
                clusters[cluster_id][org1][chromo1].add(start1)
                clusters[cluster_id][org2][chromo2].add(start2)
                alloc[org1][chromo1][start1] = cluster_id
                alloc[org2][chromo2][start2] = cluster_id
                cluster_id += 1
            elif start1 not in alloc[org1][chromo1]:
                c_id = alloc[org2][chromo2][start2]
                if org1 not in clusters[c_id]:
                    clusters[c_id][org1] = {}
                if chromo1 not in clusters[c_id][org1]:
                    clusters[c_id][org1][chromo1] = set()
                clusters[c_id][org1][chromo1].add(start1)
                alloc[org1][chromo1][start1] = alloc[org2][chromo2][start2]
            elif start2 not in alloc[org2]:
                c_id = alloc[org1][chromo1][start1]
                if org2 not in clusters[c_id]:
                    clusters[c_id][org2] = {}
                if chromo2 not in clusters[c_id][org2]:
                    clusters[c_id][org2][chromo2] = set()
                clusters[c_id][org2][chromo2].add(start2)
                alloc[org2][chromo2][start2] = alloc[org1][chromo1][start1]
            elif alloc[org1][chromo1][start1] != alloc[org2][chromo2][start2]:
                to_del = alloc[org1][chromo1][start1]
                c_id = alloc[org2][chromo2][start2]
                for org11,bib in clusters[to_del].items():
                    if org11 not in clusters[c_id]:
                        clusters[c_id][org11] = {}
                    for chromo11,s in bib.items():
                        if chromo11 not in clusters[alloc[org2][start2]][org11]:
                            clusters[c_id][org11][chromo11] = set()
                        for x in s:
                            clusters[c_id][org11][chromo11].add(x)
                            alloc[org11][chromo11][x] = c_id
                del clusters[to_del]

    return homology,clusters,alloc

def main(pivot,pivot_coords):
    
    g_anchors = {}


    for p_start,p_end,p_ori,p_chromo in pivot_coords:

        p_start,p_end = int(p_start),int(p_end)

        print('+++++++++++++++++++++++++')
        print(f'element of {pivot}:') 
        print(f'p_chromo:{p_chromo} p_start:{p_start} p_end:{p_end} p_ori:{p_ori}')

        if p_chromo not in sorted_homology_starts[pivot]:
            print(f'no anchors on {p_chromo}...skippping')
            continue
        lower_bound = max(p_start-margin1,0)
        if lower_bound == 0:
            bigger_idx_start = 0
        else:
            bigger_idx_start = bisect_right(sorted_homology_starts[pivot][p_chromo],lower_bound)
        upper_bound = min(p_end+margin1,homology[pivot][p_chromo][sorted_homology_starts[pivot][p_chromo][-1]]['end'])
        if upper_bound == homology[pivot][p_chromo][sorted_homology_starts[pivot][p_chromo][-1]]['end']:
            smaller_idx_end = len(sorted_homology_starts[pivot][p_chromo])
        else:
            smaller_idx_end = bisect_left(sorted_homology_starts[pivot][p_chromo],upper_bound)
        while smaller_idx_end > 0 and homology[pivot][p_chromo][sorted_homology_starts[pivot][p_chromo][smaller_idx_end-1]]['end'] > p_end + margin1:
            smaller_idx_end -= 1
        anchors = sorted_homology_starts[pivot][p_chromo][bigger_idx_start:smaller_idx_end]
        if len(anchors) == 0:
            print('there are no suitable anchors around this elements')
            continue
        if p_chromo not in g_anchors:
            g_anchors[p_chromo] = {}
        if p_start not in g_anchors[p_chromo]:
            g_anchors[p_chromo][p_start] = set()
        g_anchors[p_chromo][p_start].update(anchors)

    return(g_anchors)

if __name__ == "__main__":
    
    margin1 = int(sys.argv[1])
    if len(sys.argv) > 2:
        margin2 = int(sys.argv[2])
    else:
        margin2 = margin1
  
    org_mapping,chr_mapping = get_mapping()

    sorted_starts = {}
    coords_bib = {}
    with open('../coords') as f:
        for line in f:
            line = line.strip().split()
            if len(line) < 5:
                continue
            org,chromo,start,end,ori = line
            start = int(start)
            end = int(end)
            if org not in coords_bib:
                coords_bib[org] = {}
                sorted_starts[org] = {}
            if chromo not in coords_bib[org]:
                coords_bib[org][chromo] = {}
                sorted_starts[org][chromo] = []
            coords_bib[org][chromo][start] = (start,end,ori,chromo)
            sorted_starts[org][chromo].append(start)
    for org,bib in sorted_starts.items():
        for chromo,starts in bib.items():
            sorted_starts[org][chromo] = sorted(starts)

    gff = get_gff()
    homology,clusters,alloc = get_homology(gff)

    sorted_homology_starts = {}
    
    anchors = {}

    for pivot in coords_bib:
        sorted_homology_starts[pivot] = {}
        for chromo,starts in homology[pivot].items():
            sorted_homology_starts[pivot][chromo] = sorted(list(starts.keys()))
        pivot_coords = []
        for chromo,starts in coords_bib[pivot].items():
            pivot_coords += list(starts.values())
        anchors[pivot] = main(pivot,pivot_coords)


    tol = 10000
    check = {}
    drawing_order = {}
    cluster_alloc = {}
    order_c = 1
    ele_draw_alloc = {}
    for org,bib in anchors.items():
        if org not in ele_draw_alloc:
            ele_draw_alloc[org] = {}
        for chromo,bib2 in bib.items():
            for ele_start,ass in bib2.items():
                for i in ass:
                    if alloc[org][chromo][i] not in check:
                        check[alloc[org][chromo][i]] = set()
                    if alloc[org][chromo][i] not in cluster_alloc:
                        cluster_alloc[alloc[org][chromo][i]] = set()
                    check[alloc[org][chromo][i]].add(i)
                    cluster_alloc[alloc[org][chromo][i]].add((org,ele_start))
    for c_id,members in cluster_alloc.items():
        order_c_alt = set()
        drawing_order[order_c] = set()
        for org,ele_start in members:
            if ele_start in ele_draw_alloc[org]:
                order_c_alt.add(ele_draw_alloc[org][ele_start])
        for o_c in order_c_alt:
            for org,ele_start in drawing_order[o_c]:
                ele_draw_alloc[org][ele_start] = order_c
                drawing_order[order_c].add((org,ele_start))
            del drawing_order[o_c]
        for org,ele_start in members:
            ele_draw_alloc[org][ele_start] = order_c
            drawing_order[order_c].add((org,ele_start))
        order_c += 1

    final_drawing_order = {}
    c_new = 1
    for c,members in drawing_order.items():
        if len(members) > 1:
            final_drawing_order[c_new] = members
            c_new += 1

    counter = {}
    no_colors = 0
    for c_id,members in check.items():
        alloc_times = len(members)
        if alloc_times not in counter:
            counter[alloc_times] = 0
        counter[alloc_times] += 1
        if alloc_times > 1:
            clusters[c_id]['color'] = 1
            no_colors += 1
    
    color = cm.rainbow(np.linspace(0, 1, no_colors))
    used = {}
    for c,bib in clusters.items():
        if 'color' in bib:
            ri = randint(0,len(color)-1)
            while ri in used:
                ri = randint(0,len(color)-1)
            used[ri] = True
            bib['color'] = color[ri]

    run(f'mkdir -p clusters && mkdir -p clusters/{margin1}',shell=True)
    run(f'mkdir -p images && mkdir -p images/{margin1}',shell=True)

    for c,members in final_drawing_order.items():
        gv = GenomeViz()
        out = open(f'clusters/{margin1}/cluster_{c}.txt','w')
        out.write(f'cluster {c}:')
        for org,ele_start in members:
            out.write('--------------------\n')
            out.write(f'species {org}\n')
            min_dis_to_gene = 1e15
            for chromo,bib2 in anchors[org].items():
                if ele_start not in bib2:
                    continue
                ass = bib2[ele_start]
                to_draw = []
                mini = 1e15
                maxi = 0

                start,end,ori,chro = coords_bib[org][chromo][ele_start]
                start_gene = start
                end_gene = end
                if ori == 'forward':
                    ori = 1
                else:
                    ori = -1
                if ori == -1:
                    turn = True
                else:
                    turn = False
                ori = 1
                to_draw.append([start, end,'red','arrow','gene',0,"center","center",ori])
                out.write(f'element on chromosome {chromo} start {start} end {end} strand {ori}\n')
                if start < mini:
                    mini = start
                if end > maxi:
                    maxi = end

                for i in ass:
                    start = i
                    if len(check[alloc[org][chromo][start]]) < 2:
                        continue
                    end = homology[org][chromo][start]['end']
                    score = homology[org][chromo][start]['score']
                    mid = start + (end-start)/2
                    start_score = mid - score/2
                    end_score = mid + score/2
                    if start_score >= start_gene and start_score <= end_gene:
                        continue
                    if end_score >= start_gene and end_score <= end_gene:
                        continue
                    if end_score < start_gene:
                        dis_to_gene = start_gene - end_score
                        start_score = (-1,dis_to_gene)
                    elif start_score > end_gene:
                        dis_to_gene = start_score - end_gene
                        start_score = (1,dis_to_gene)
                    else:
                        continue
                    if dis_to_gene < 1000:
                        continue

                    if dis_to_gene < min_dis_to_gene:
                        min_dis_to_gene = dis_to_gene

                    start_score = start_score 
                    to_draw.append([start_score[0], start_score[1], alloc[org][chromo][start],'bigbox',score])#,alloc[org][chromo][i]])
                    if start < mini:
                        mini = start
                    if end > maxi:
                        maxi = end
                name, genome_size = org + '_' + str(chromo) + '_' + str(ele_start), 2*margin1 + tol
                track = gv.add_feature_track(name,genome_size)
                out.write(f'has {len(to_draw)-1} anchors (one or more significant alignments) that are shared with something in this cluster\n')
                for ele in to_draw:
                    if len(ele) > 5:
                        ele[0] = genome_size/2 - end_gene + start_gene
                        ele[1] = genome_size/2 + end_gene - start_gene
                    else:
                        if ele[0] == -1:

                            new_end = genome_size/2 - end_gene + start_gene - ele[1]
                            new_start = new_end - ele[4]

                        else:

                            new_start = genome_size/2 - end_gene + start_gene + ele[1]
                            new_end = new_start + ele[4]
                        ele[0] = new_start
                        ele[1] = new_end

                        if turn:
                            temp = genome_size - ele[0]
                            ele[0] = genome_size - ele[1]
                            ele[1] = temp

                    if len(ele) > 5:
                        #track.add_feature(ele[0],ele[1],facecolor=ele[2],plotstyle=ele[3],label=ele[4],labelrotation=ele[5],labelvpos=ele[6],labelha=ele[7],strand=ele[8])
                        track.add_feature(ele[0],ele[1],facecolor=ele[2],plotstyle=ele[3],label=ele[4],labelrotation=ele[5],labelha=ele[7],strand=ele[8])
                    else:
                        track.add_feature(ele[0],ele[1],facecolor=clusters[ele[2]]['color'],plotstyle=ele[3],label='')
            out.write('--------------------\n')
        gv.savefig(f'images/{margin1}/cluster_{c}.svg')
        gv.savefig(f'images/{margin1}/cluster_{c}.png')
        out.close()
