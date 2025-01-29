#! /usr/bin/env python3

from pprint import pprint
import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
from Bio import SeqIO
import pickle
from random import randint
from matplotlib.pyplot import cm
from pygenomeviz import GenomeViz
import numpy as np
from subprocess import run
from bisect import bisect_left,bisect_right

def get_genome(genomes_dir,org):
    
    seqs = SeqIO.parse(f'{genomes_dir}/{org}.fasta', "fasta")
    s = {}
    for seq in seqs:
        s[seq.id] = seq

    return(s)

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
    genomes = {}
    maps = {}
    aligned = {}
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

    sorted_anchors = {}
    tol = 10000
    check = {}
    drawing_order = {}
    cluster_alloc = {}
    order_c = 1
    ele_draw_alloc = {}
    for org,bib in anchors.items():
        genomes[org] = get_genome('HOMEDIR/adrena/stable_synteny/utils/genomes',org)
        with open(f'../../../../anchors_to_save/andrena/candidates/{org}','rb') as f2:
            maps[org] = pickle.load(f2)
        with open(f'../../../../anchors_to_save/andrena/aligned/{org}','rb') as f2:
            aligned[org] = pickle.load(f2)
        sorted_anchors[org] = {}
        for i,bibus in maps[org].items():
            if bibus['chromosome'] not in sorted_anchors[org]:
                sorted_anchors[org][bibus['chromosome']] = []
            sorted_anchors[org][bibus['chromosome']].append((bibus['start'],bibus['end'],i))
        for chromo,l in sorted_anchors[org].items():
            sorted_anchors[org][chromo] = sorted(l)
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

    run(f'rm -r clusters && mkdir -p clusters && mkdir -p clusters/{margin1}',shell=True)
    run(f'rm -r images && mkdir -p images && mkdir -p images/{margin1}',shell=True)


    for c,members in final_drawing_order.items():
        gv = GenomeViz()
        out = open(f'clusters/{margin1}/cluster_{c}.txt','w')
        out.write(f'cluster {c}:')
        clusters_drawn = set()
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
                    clusters_drawn.add(alloc[org][chromo][start])
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
        member_orgs = [x[0] for x in members]
        run(f'mkdir -p clusters/{margin1}/alignments/',shell=True)
        #datentechnisch = open(f'clusters/{margin1}/succinct_alignments_cluster_{c}.txt','w')
        datentechnisch2 = open(f'clusters/{margin1}/alignments_cluster_{c}_local_and_with_scores.txt','w')
        #datentechnisch.write('# species1 chromo1 start1 end1 species2 chromo2 start2 end score\n')
        datentechnisch2.write('# species1 chromo1 start1 end1 species2 chromo2 start2 end score score_against_org1_genome score_against_org2_genome\n')
        with open(f'clusters/{margin1}/alignments_cluster_{c}.txt','w') as f:
            #f.write(f'-----------------------------------\n')
            f.write(f'cluster {c} alignments:\n')
            f.write('----------------------\n')
            run(f'mkdir -p clusters/{margin1}/alignments/cluster_{c}/',shell=True)
            for cc in clusters_drawn:
                run(f'mkdir -p clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}',shell=True)
                run(f'mkdir -p clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas',shell=True)
                run(f'mkdir -p clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/blast_out',shell=True)
                to_write = {}
                f.write(f'-----------------------------------\n')
                start_ends = set()
                to_blast_trans = set()
                no_s = {}
                for org,bib in clusters[cc].items():
                    if org == 'color':
                        continue
                    to_write[org] = []
                    no_s[org] = []
                    for chromo,bib2 in bib.items():
                        for start in bib2:
                            if org not in member_orgs:
                                pass
                                f.write(f'### element not in cluster ({org}\t{chromo}\t{start}\t{homology[org][chromo][start]["end"]}) ###\n')
                                idx_start = bisect_left(sorted_anchors[org][chromo],start,key=lambda x: x[0])
                                idx_end = bisect_right(sorted_anchors[org][chromo],homology[org][chromo][start]["end"],key=lambda x: x[1])
                                for start,end,orig_i in sorted_anchors[org][chromo][idx_start:idx_end]:
                                    to_blast_trans.add(org)
                                    seq = genomes[org][chromo][int(start):int(end)]
                                    seq.id = f'{org} {chromo} {start} {end} forward'
                                    seq.description = ''
                                    to_write[org].append((seq,orig_i))
                            else:
                                f.write(f'REGION: {org}\t{chromo}\t{start}\t{homology[org][chromo][start]["end"]}\n')
                                f.write(f'ORIGINAL ANCHORS (start - end):\n')
                                f.write('\t\t========\n')
                                idx_start = bisect_left(sorted_anchors[org][chromo],start,key=lambda x: x[0])
                                idx_end = bisect_right(sorted_anchors[org][chromo],homology[org][chromo][start]["end"],key=lambda x: x[1])
                                for start,end,orig_i in sorted_anchors[org][chromo][idx_start:idx_end]:
                                    f.write(f'\t\t{start}\t{end}\n')
                                    start_ends.add((org,chromo,start,end,orig_i))
                                    f.write('\t\t========\n')
                                    to_blast_trans.add(org)
                                    seq = genomes[org][chromo][int(start):int(end)]
                                    seq.id = f'{org} {chromo} {start} {end} forward'
                                    seq.description = ''
                                    to_write[org].append((seq,orig_i))
                start_ends = list(start_ends)
                done = {}
                for enu,se in enumerate(start_ends):
                    for se2 in start_ends[enu+1:]:
                        if se[0] == se2[0]:
                            continue
                        else:
                            if se2[0] not in aligned[se[0]][se[4]]['matches']:
                                print(f'org1={se[0]} anchor start={se[2]} and org2={se2[0]} with start={se2[2]} do not have a direct match')
                                print('lets see all the pw alignments in the cluster')
                                print('---------------------------------------------')
                                for orgus,bibus in clusters[cc].items():
                                    print(f'org = {orgus}')
                                    if orgus == 'color':
                                        continue
                                    for chromous,bib2us in bibus.items():
                                        print(f'chromosome = {chromous}')
                                        for startus in bib2us:
                                            endus = homology[orgus][chromous][startus]["end"]
                                            idx_start = bisect_left(sorted_anchors[orgus][chromous],startus,key=lambda x: x[0])
                                            idx_end = bisect_right(sorted_anchors[orgus][chromous],endus,key=lambda x: x[1])
                                            for start,end,orig_i in sorted_anchors[orgus][chromous][idx_start:idx_end]:
                                                print(f'anchor {orig_i} with matches:')
                                                for orgus2,matchbib in aligned[orgus][orig_i]['matches'].items():
                                                    print(f'matches with {orgus2}:')
                                                    for m,b in matchbib['matches'].items():
                                                        print(f'anchor 2 chromo = {maps[orgus2][m]["chromosome"]}')
                                                        print(f'anchor 2 start = {maps[orgus2][m]["start"]}')
                                                        print(f'anchor 2 end = {maps[orgus2][m]["end"]}')
                                                        print(f'anchor 2 score = {b["match score"]}')
                                                
                                print('---------------------------------------------')


                                pass
                            elif se2[4] not in aligned[se[0]][se[4]]['matches'][se2[0]]['matches']:
                                pass
                            else:
                                start2 = aligned[se[0]][se[4]]['matches'][se2[0]]['matches'][se2[4]][f'hit coordinates in {se2[0]} candidate'][0]
                                end2 = aligned[se[0]][se[4]]['matches'][se2[0]]['matches'][se2[4]][f'hit coordinates in {se2[0]} candidate'][1]
                                start1 = aligned[se[0]][se[4]]['matches'][se2[0]]['matches'][se2[4]][f'hit coordinates in (own) {se[0]} candidate'][0]
                                end1 = aligned[se[0]][se[4]]['matches'][se2[0]]['matches'][se2[4]][f'hit coordinates in (own) {se[0]} candidate'][1]
                                score = aligned[se[0]][se[4]]['matches'][se2[0]]['matches'][se2[4]]['match score']
                                start1 = int(start1)
                                end1 = int(end1)
                                start2 = int(start2)
                                end2 = int(end2)
                                start1 += se[2]
                                end1 += se[2]
                                start2 += se2[2]
                                end2 += se2[2]

                                #datentechnisch.write(f'{se[0]}\t{se[1]}\t{se[2]}\t{se[3]}\t{se2[0]}\t{se2[1]}\t{se2[2]}\t{se2[3]}\t{score}\n')
                                #HERE WAR ES EINFACH OHNE DAS ZEUG ZWISCHEN UND DEM NEUEN GLEICHEN
                                #datentechnisch2.write(f'{se[0]}\t{se[1]}\t{start1}\t{end1}\t{se2[0]}\t{se2[1]}\t{start2}\t{end2}\t{score}\n')

                                for seq,orig_i in to_write[se[0]]:
                                    if orig_i == se[4]:
                                        no_s[se[0]].append(orig_i)
                                        SeqIO.write(seq,f'clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{se[0]}_{orig_i}.fasta','fasta')
                                        break
                                run(f'blastn -query clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{se[0]}_{orig_i}.fasta -db HOMEDIR/adrena/stable_synteny/utils/blastdbs/{se2[0]} -outfmt 6 -evalue 1e-10 -out tmp_out',shell=True)
                                maxi = 0
                                with open('tmp_out') as tmp_out:
                                    for line in tmp_out:
                                        if float(line.strip().split()[-1]) > maxi:
                                            maxi = float(line.strip().split()[-1])

                                snd_1_2 = maxi

                                for seq,orig_i in to_write[se2[0]]:
                                    if orig_i == se2[4]:
                                        no_s[se2[0]].append(orig_i)
                                        SeqIO.write(seq,f'clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{se2[0]}_{orig_i}.fasta','fasta')
                                        break
                                run(f'blastn -query clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{se2[0]}_{orig_i}.fasta -db HOMEDIR/adrena/stable_synteny/utils/blastdbs/{se[0]} -outfmt 6 -evalue 1e-10 -out tmp_out',shell=True)
                                maxi = 0
                                with open('tmp_out') as tmp_out:
                                    for line in tmp_out:
                                        if float(line.strip().split()[-1]) > maxi:
                                            maxi = float(line.strip().split()[-1])

                                snd_2_1 = maxi

                                datentechnisch2.write(f'{se[0]}\t{se[1]}\t{start1}\t{end1}\t{se2[0]}\t{se2[1]}\t{start2}\t{end2}\t{score}\t{snd_1_2}\t{snd_2_1}\n')


                                done[se] = 1
                                done[se2] = 1
                written = 0
                for se in start_ends:
                    if se in done:
                        continue
                    else:
                        if written == 0:
                            written = 1
                            for oggo,seqs in to_write.items():
                                for seq,orig_i in seqs:
                                    if orig_i in no_s[oggo]:
                                        continue
                                    no_s[oggo].append(orig_i)
                                    SeqIO.write(seq,f'clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{oggo}_{orig_i}.fasta','fasta')

                        for org_to_blast_trans_against in to_blast_trans:
                            if org_to_blast_trans_against == se[0]:
                                continue
                            for no_seq in no_s[org_to_blast_trans_against]:
                                run(f'blastn -query clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{se[0]}_{se[4]}.fasta -subject clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{org_to_blast_trans_against}_{no_seq}.fasta -outfmt 6 -evalue 1e-10 -out clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/blast_out/{se[0]}_{se[4]}_{org_to_blast_trans_against}_{no_seq}',shell=True)
                        #datentechnisch2.write(f'# {se[0]}\t{se[1]}\t{se[2]}\t{se[3]}\ttransitive membership\tblast scores against:\n# ------\n')
                        lines = []
                        lines.append(f'# {se[0]}\t{se[1]}\t{se[2]}\t{se[3]}\ttransitive membership\tblast scores against (best when checked against all seqs from this species) + (best (up to) 3 lines of blast output like normally above when against whole other genome (to see how ambiguous/repetitive)):\n# ------\n')
                        maxi = {}
                        for org_to_blast_trans_against in to_blast_trans:
                            if org_to_blast_trans_against == se[0]:
                                continue
                            check_trans = 0
                            if org_to_blast_trans_against not in member_orgs:
                                #datentechnisch2.write(f'# (not in cluster of genes -> no genes near these matches or not clear enough evidence) \t')
                                lines.append(f'# (not in cluster of genes -> no genes near these matches or not clear enough evidence) \t')
                            #datentechnisch2.write(f'# {org_to_blast_trans_against}\t')
                            else:
                                check_trans = 1
                            if check_trans == 1:
                                for chromoos,startos in clusters[cc][org_to_blast_trans_against].items():
                                    for sos in startos:
                                        endos = homology[org_to_blast_trans_against][chromoos][sos]["end"] 
                                        lines.append(f'# {org_to_blast_trans_against}, syntenic region in this cluster:{chromoos}:{sos}-{endos}\t')
                            else:
                                lines.append(f'# {org_to_blast_trans_against}\t')
                            maxi[org_to_blast_trans_against] = 0
                            for no_seq in no_s[org_to_blast_trans_against]:
                                with open(f'clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/blast_out/{se[0]}_{se[4]}_{org_to_blast_trans_against}_{no_seq}') as blast_out:
                                    for line in blast_out:
                                        if float(line.split()[-1]) > maxi[org_to_blast_trans_against]:
                                            maxi[org_to_blast_trans_against] = float(line.split()[-1])
                            #datentechnisch2.write(f'# {maxi[org_to_blast_trans_against]} \n')
                            lines.append(f'# {maxi[org_to_blast_trans_against]} \n')
                            if check_trans == 1:
                                oggo = org_to_blast_trans_against
                                run(f'blastn -query clusters/{margin1}/alignments/cluster_{c}/anchors_{cc}/fastas/{se[0]}_{se[4]}.fasta -db HOMEDIR/adrena/stable_synteny/utils/blastdbs/{oggo} -outfmt 6 -evalue 1e-10 -out tmp_out',shell=True)
                                scores_lines = []
                                with open('tmp_out') as tmp_out:
                                    for line in tmp_out:
                                        scores_lines.append((float(line.strip().split()[-1]),line))
                                scores_lines = sorted(scores_lines, key=lambda tup: tup[0])
                                x = 3
                                for s,l in scores_lines:
                                    if x > 0:
                                        lines.append('# '+l)
                                        x -= 1
                                    else:
                                        break
                           
                        if sum(list(maxi.values())) == 0:
                            continue
                        else:
                            for line in lines:
                                datentechnisch2.write(line)

                        ### what does the orig anchor dict say

                        #print('#### original hit dict ####')

                        #pprint(aligned[se[0]][se[4]],stream=datentechnisch2)

                        #print('####') 

                        ###

                        datentechnisch2.write(f'# ------\n')

                f.write(f'-----------------------------------\n')
            #f.write(f'-----------------------------------\n')
        datentechnisch = open(f'succint_alignments_cluster_{c}.txt','w')
