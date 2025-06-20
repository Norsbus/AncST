#! /usr/bin/env python3

import pickle
from bisect import bisect_left,bisect_right
import multiprocessing as mp
from copy import deepcopy
import numpy as np
from dbscan1d.core import DBSCAN1D
import argparse
import networkx as nx
import markov_clustering as mc
import os

def which_chromo(seqlen,seqids,i):
    seqlen = sorted(seqlen)
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

def get_regions(org,iss,j):
    iss.sort()
    labels = dbs.fit_predict(iss)
    regions = {}
    for enu,ident in enumerate(labels):
        if ident == -1:
            continue
        else:
            if ident not in regions:
                regions[ident] = []
            regions[ident].append(iss[enu])
    res = {}
    for ident,region in regions.items():
        res[ident] = [min(region),max(region)]
    return org,res

def get_elements_and_syntenies_of_regions(global_regions):
    anchors_per_region = {}
    cont_eles = {}
    for org,bib in global_regions.items():
        anchors_per_region[org] = {}
        cont_eles[org] = {}
        for ident,t in sorted(bib.items(), key=lambda x: x[1]):
            anchors_per_region[org][ident] = {}
            start,end = t
            abs_s = start
            abs_e = end
            seqids,seqlen = small_meta[org]
            chromo = which_chromo(seqlen,seqids,start)
            chromo_end = which_chromo(seqlen,seqids,end)
            idx_l = seqids.index(chromo)
            if idx_l > 0:
                le = seqlen[idx_l-1]
            else:
                le = 0
            len_chr = seqlen[idx_l]
            start -= le
            start = max(0,start)
            end -= le
            end = min(len_chr,end)
            
            if org not in coords or chromo not in coords[org]:
                elements = ''
            else:
                elements = sorted(coords[org][chromo])
                starts,ends,names = zip(*elements)
                idx_first_ele = bisect_right(starts,start)
                idx_last_ele = bisect_left(ends,end)
                if idx_last_ele - idx_first_ele > 0:
                    elements = ';'.join(names[idx_first_ele:idx_last_ele])
                else:
                    elements = ''

            anchors = sorted(list(collect[org].keys()))
            idx_first_ele = bisect_right(anchors,abs_s)
            idx_last_ele = bisect_left(anchors,abs_e)
            if idx_last_ele - idx_first_ele > 0:
                for a in anchors[idx_first_ele:idx_last_ele]:
                    anchors_per_region[org][ident][a] = deepcopy(collect[org][a])
            cont_eles[org][ident] = elements

    pw_syn = {}
    G=nx.Graph()
    
    for org,bib in anchors_per_region.items():
        pw_syn[org] = {}
        for ident1,bib2 in bib.items():
            pw_syn[org][ident1] = {}
            for a,bib3 in bib2.items():
                for org2,lj in bib3.items():
                    if org2 in anchors_per_region:
                        for ident2,bib22 in anchors_per_region[org2].items():
                            for j in lj:
                                if j in bib22:
                                    if org2 not in pw_syn[org][ident1]:
                                        pw_syn[org][ident1][org2] = {}
                                    if ident2 not in pw_syn[org][ident1][org2]:
                                        pw_syn[org][ident1][org2][ident2] = 0
                                    pw_syn[org][ident1][org2][ident2] += 1

    G.add_nodes_from([f'{o}$$${x}' for o,x in pw_syn.items()])
    for org,bib in pw_syn.items():
        for id1,bib2 in bib.items():
            for org2,bib3 in bib2.items():
                for id2,no in bib3.items():
                    G.add_edge(f'{org}$$${id1}',f'{org2}$$${id2}',weight=no)

    G_trans = nx.transitive_closure(G,reflexive=None)
    trans_syn = {}

    for c in sorted(nx.connected_components(G_trans), key=len, reverse=True):
        c = list(c)
        if len(c) > 1:
            eles = []
            for i in c:
                org,ident = i.split('$$$')
                ident = int(ident)
                eles.append((org,ident))
            for ele in eles:
                for ele2 in eles:
                    if ele != ele2 and ele[0] != ele2[0]:
                        org1,ident1 = ele
                        org2,ident2 = ele2
                        if org1 in pw_syn and ident1 in pw_syn[org1] and org2 in pw_syn[org1][ident1] and ident2 in pw_syn[org1][ident1][org2]:
                            continue
                        if org1 not in trans_syn:
                            trans_syn[org1] = {}
                        if ident1 not in trans_syn[org1]:
                            trans_syn[org1][ident1] = {}
                        if org2 not in trans_syn[org1][ident1]:
                            trans_syn[org1][ident1][org2] = {}
                        trans_syn[org1][ident1][org2][ident2] = -42

    nodelist = []
    mapping = {}
    c = 0
    for node in G.nodes():
        nodelist.append(node)
        mapping[c] = node
        c += 1
    matrix = nx.to_scipy_sparse_array(G, nodelist=nodelist)
    result = mc.run_mcl(matrix)
    clusters_mcl = mc.get_clusters(result)

    pw_syn_mcl = {}

    for c in clusters_mcl:
        c = list(c)
        if len(c) > 1:
            eles = []
            for i in c:
                i = mapping[i]
                org,ident = i.split('$$$')
                ident = int(ident)
                eles.append((org,ident))
            for ele in eles:
                for ele2 in eles:
                    if ele != ele2 and ele[0] != ele2[0]:
                        org1,ident1 = ele
                        org2,ident2 = ele2
                        if org1 not in pw_syn_mcl:
                            pw_syn_mcl[org1] = {}
                        if ident1 not in pw_syn_mcl[org1]:
                            pw_syn_mcl[org1][ident1] = {}
                        if org2 not in pw_syn_mcl[org1][ident1]:
                            pw_syn_mcl[org1][ident1][org2] = {}
                        pw_syn_mcl[org1][ident1][org2][ident2] = -21

    return cont_eles,pw_syn,trans_syn,pw_syn_mcl

def get_new_regions_to_eval(org,global_reg,old_reg):
    if len(old_reg) == 0:
        return org,[(s,e) for s,e in global_reg.values()]
    new_regions_to_eval = []
    starts = sorted([x[0] for x in old_reg])
    ends = sorted([x[1] for x in old_reg])
    for ident,region in global_reg.items():
        s,e = region
        idx_s = bisect_left(starts,s)
        idx_e = bisect_right(ends,e)
        if idx_s != idx_e:
            continue
        else:
            if (idx_s == 0 and e < starts[0]) or idx_s == len(old_reg):
                new_regions_to_eval.append((s,e))
            else:
                if e < starts[idx_s] and starts[idx_s] - e > margin_to_be_new and s > ends[idx_s - 1] and s - ends[idx_s - 1] > margin_to_be_new:
                    new_regions_to_eval.append((s-1000,e+1000))

    return org,new_regions_to_eval

def get_coords(query_org,target_orgs,regions,first=0):
    if first == -1:
        regions = list(regions.values())
    with open(f'{aligned_path}/{query_org}','rb') as f:
        aligned = pickle.load(f)
    iss = sorted(aligned.keys())
    matching_anchors = {}
    for s,e in regions:
        local_iss = iss[bisect_left(iss,s):min(len(iss), bisect_left(iss,e))]
        for i in local_iss:
            if first == 0 and query_org in coords_to_exclude:
                valid = 1
                for sx,ex in coords_to_exclude[query_org]:
                    if i >= sx and i <= ex:
                        valid = 0
                        break
                    elif i < sx and sx - i > 100:
                        valid = 0
                        break
                    elif i > ex and i - ex > 100:
                        valid = 0
                        break
                if valid == 0:
                    continue
            bib = aligned[i]
            for org in target_orgs:
                if org in bib['matches'] and bib['matches'][org]['meta']['multiple matches out of tolerance range'] != 1:
                    matches = list(bib['matches'][org]['matches'].keys())
                    if i not in matching_anchors:
                        matching_anchors[i] = {}
                    if org not in matching_anchors[i]:
                        matching_anchors[i][org] = set()
                    matching_anchors[i][org].update(matches)
    return query_org,matching_anchors

def reconcile_collect(res,collect):
    for qo,ma in res:
        if qo not in collect:
            collect[qo] = {}
        for i,bib in ma.items():
            if i not in collect[qo]:
                collect[qo][i] = {}
            for to,sj in bib.items():
                if to not in collect[qo][i]:
                    collect[qo][i][to] = set()
                if to not in collect:
                    collect[to] = {}
                for j in sj:
                    if j not in collect[to]:
                        collect[to][j] = {}
                    if qo not in collect[to][j]:
                        collect[to][j][qo] = set()
                    collect[qo][i][to].add(j)
                    collect[to][j][qo].add(i)
    return collect 


if __name__ == "__main__":
    parser= argparse.ArgumentParser(description='Find Syntenic Regions From List of Genetic Elements of Interest',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--cores', type=int,nargs='?', default=1,help="number of cores to run with")
    parser.add_argument('--iter', type=int,nargs='?', default=3,help="number of iterations to find anchors and thus syntenic regions")
    parser.add_argument('--threshold_new_region', type=int,nargs='?', default=1e6,help="margin of syntenic regions, i.e.: number of nulceotides within newly added anchors are considered to be within a neighbouring syntenic region")
    parser.add_argument('--anchor_scope', type=int,nargs='?', default=5e5,help="number of nucleotides up-and downstream of input genetic elements to be considered for initial anchor search")
    parser.add_argument('--no_anchors_form_region', type=int,nargs='?', default=4,help="minimum number of anchors to form a'syntenic region'")
    parser.add_argument('--anchor_path', type=str,nargs='?', default='../utils/anchors/')
    parser.add_argument('--small_meta_path', type=str,nargs='?', default='../utils/small_meta/')
    parser.add_argument('--coords', type=str,nargs='?', default='coords',help="file with coordinates of elements of interest in format:species*tab*chr*tab*start*tab*end*tab*strand (+/-)*tab*[optional:name]")
    parser.add_argument('--coords_to_exclude', type=str,nargs='?',help="file with coordinates of additional elements to exclude from being considered as anchors (in same format as elements of interest: species*tab*chr*tab*start*tab*end*tab*strand(+/-)*tab*[optional:name]. if you happen to know your genomic environment it may be helpful to exclude known elements which should not be considered as synteny anchors since AncST produces <= 5 percent false positives. so if you e.g. want to be rather conservative, are in an environment with complex duplication histories, are working with rather badly assembled genomes, etc. you might want to exclude those elements. on the other hand, this will generally reduce the resolution of synteny anchors and in almost all cases it is fine to leave this option empty and just consider all anchors around since there are additional 'safety' checks and 'fishy' results will show as inconsistent and/or 'too broad' clustering of the regions considered")
    args = parser.parse_args()
    if args.cores:
        no_cores = int(args.cores)
    if args.iter:
        no_iter = int(args.iter)
    if args.threshold_new_region:
        margin_to_be_new = int(args.threshold_new_region)
    if args.anchor_scope:
        margin_anchors = int(args.anchor_scope)
    if args.no_anchors_form_region:
        no_anchors_form_region = int(args.no_anchors_form_region)
    if args.anchor_path:
        anchor_path = args.anchor_path
    if args.small_meta_path:
        small_meta_path = args.small_meta_path
    if args.coords:
        coords_file = args.coords
    
    aligned_path = f'{anchor_path}/aligned/'
    candidates_path = f'{anchor_path}/candidates/'

    orgs = []
    if os.path.isfile('orgs'):
        with open('orgs','r') as f:
            for line in f:
                orgs.append(line.strip())
    else:
        orgs = os.listdir(f'{anchor_path}/aligned/')

    target_orgs = orgs

    small_meta = {}
    for org in orgs:
        with open(f'{small_meta_path}/{org}','rb') as f:
            small_meta[org] = pickle.load(f)

    dbs = DBSCAN1D(eps=margin_anchors/3, min_samples=no_anchors_form_region)

    collect = {}
    coords = {}
    coords_to_exclude = {}
    
    name_count = 1
    with open(coords_file) as f:
        for line in f:
            line = line.strip().split()
            if len(line) == 6:
                org,chromo,start,end,ori,name = line
            elif len(line) == 5:
                org,chromo,start,end,ori = line
                name = f'element_on_line_{name_count}'
                name_count += 1
            else:
                print(line)
                print('please provide a coodinates file of form species(name of genome/anchors)\tchr\tstart\tend\tstrand (+/-)\t[optional:name]')
                exit(0)
            start = int(start)
            end = int(end)
            mid = int(start+(end-start)/2)
            seqids,seqlen = small_meta[org]
            idx_l = seqids.index(chromo)
            if idx_l > 0:
                l = seqlen[idx_l-1]
            else:
                l = 0
            len_chr = seqlen[idx_l]
            mid += l
            if org not in collect:
                collect[org] = {}
                coords[org] = {}
            collect[org][mid-10000] = {}
            collect[org][mid] = {}
            collect[org][mid+10000] = {}
            if chromo not in coords[org]:
                coords[org][chromo] = []
            coords[org][chromo].append((start,end,name))
            if org not in coords_to_exclude:
                coords_to_exclude[org] = []
            coords_to_exclude[org].append((start,end))

    if args.coords_to_exclude:
        coords_to_exclude_file = args.coords_to_exclude
        with open(coords_to_exclude_file) as f:
            for line in f:
                line = line.strip().split()
                org,chromo,start,end,ori_gene,name = line
                start = int(start)
                end = int(end)
                seqids,seqlen = small_meta[org]
                idx_l = seqids.index(chromo)
                if idx_l > 0:
                    l = seqlen[idx_l-1]
                else:
                    l = 0
                len_chr = seqlen[idx_l]
                start += l
                end += l
                if org not in coords_to_exclude:
                    coords_to_exclude[org] = []
                coords_to_exclude[org].append((start,end))

    initial_regions = {}
    for org,bib in collect.items():
        initial_regions[org] = []
        for mid in bib:
            initial_regions[org].append((max(0,mid-margin_anchors),mid+margin_anchors))

    print(f'iteration no 1')

    with mp.Pool(processes=no_cores) as p:
        res = p.starmap(get_coords,[(query_org,target_orgs,regions,1) for query_org,regions in initial_regions.items()])

    collect = reconcile_collect(res,collect)

    j = 1
    with mp.Pool(processes=no_cores) as p:
        res = p.starmap(get_regions,[(org,np.array(list(bib.keys())),j) for org,bib in collect.items()])

    global_regions = {}
    for org,regions in res:
        if len(regions) > 0:
            global_regions[org] = regions

    old_regions = deepcopy(initial_regions)
    for org in list(global_regions.keys()):
        if org not in old_regions:
            old_regions[org] = []

    with mp.Pool(processes=no_cores) as p:
        res = p.starmap(get_new_regions_to_eval,[(org,global_regions[org],old_regions[org]) for org in list(global_regions.keys())])
   
    new_regions = {}
    for org,regions in res:
        new_regions[org] = regions

    for j in range(2,no_iter+1):

        if len(new_regions) == 0:
            print(f'exiting before iteration no {j} since no new regions found')
            break

        print(f'iteration no {j}')

        with mp.Pool(processes=no_cores) as p:
            res = p.starmap(get_coords,[(query_org,target_orgs,regions) for query_org,regions in new_regions.items()])
        
        collect = reconcile_collect(res,collect)

        old_regions = {}
        for org,bib in global_regions.items():
            old_regions[org] = list(bib.values())

        with mp.Pool(processes=no_cores) as p:
            res = p.starmap(get_regions,[(org,np.array(list(bib.keys())),j) for org,bib in collect.items()])

        global_regions = {}
        for org,regions in res:
            if len(regions) > 0:
                global_regions[org] = regions

        for org in list(global_regions.keys()):
            if org not in old_regions:
                old_regions[org] = []

        with mp.Pool(processes=no_cores) as p:
            res = p.starmap(get_new_regions_to_eval,[(org,global_regions[org],old_regions[org]) for org in list(global_regions.keys())])
       
        new_regions = {}
        for org,regions in res:
            new_regions[org] = regions
    
    with mp.Pool(processes=no_cores) as p:
        res = p.starmap(get_coords,[(query_org,target_orgs,regions,-1) for query_org,regions in global_regions.items()])
    collect = reconcile_collect(res,collect)

    cont_eles,pw_syn,trans_syn,pw_syn_mcl = get_elements_and_syntenies_of_regions(global_regions)

    bp_no = 1

    with open(f'syntenic_regions','w') as f:
        for org,bib in global_regions.items():
            f.write(f'==============================================================================================================\n')
            f.write(f'{org}\n')
            for ident,t in sorted(bib.items(), key=lambda x: x[1]):
                start,end = t
                elements = cont_eles[org][ident]
                seqids,seqlen = small_meta[org]
                chromo = which_chromo(seqlen,seqids,start)
                chromo_end = which_chromo(seqlen,seqids,end)
                if chromo_end != chromo:
                    seqlen_chr_start = seqlen[seqids.index(chromo)]
                    save_end = end
                    end = seqlen_chr_start - 1
                idx_l = seqids.index(chromo)
                if idx_l > 0:
                    le = seqlen[idx_l-1]
                else:
                    le = 0
                len_chr = seqlen[idx_l]
                start -= le
                start = max(0,start)
                end -= le
                end = min(len_chr,end)
                if chromo_end != chromo:
                    elements1 = []
                    for x in coords[org][chromo]:
                        if x[0] >= start and x[1] <= end:
                            elements1.append(x[-1])
                    f.write('-----------------------------------------------------------------------------------------------------------------\n')
                    f.write(f'<<< {org}\tregion {ident}\t{chromo}\t{start}\t{end} (== end of chr) \tPOTENTIAL CHROMOSOMAL BREAKPOINT {bp_no} REGION 1 >>>\n')
                    f.write(f'\telements: {elements1}\n')
                    f.write(f'\tdirectly pairwise syntenic to [1 = after markov chain clustering (strength connections between regions == number of anchor matches); 0 = could not hold up MCL]:\n')
                    if org in pw_syn and ident in pw_syn[org]:
                        for org2,l in pw_syn[org][ident].items():
                            if len(l) == 0:
                                continue
                            l_mcl = []
                            if org in pw_syn_mcl and ident in pw_syn_mcl[org] and org2 in pw_syn_mcl[org][ident]:
                                l_mcl = pw_syn_mcl[org][ident][org2]
                            if len(l) > 1:
                                f.write(f'\t\t{org2} (potential breakpoints as multiple regions match - likely just artificial partition into multiple regions due to actual sequence insertions/deletions or program parameter/anchor availability/assembly artefacts):\n')
                            else:
                                f.write(f'\t\t{org2}:\n')
                            for reg2 in l:
                                if reg2 in l_mcl:
                                    f.write(f'\t\t\tregion {reg2} [1] \n')
                                else:
                                    f.write(f'\t\t\tregion {reg2} [0] \n')
                    f.write(f'\tadditionally transitively pairwise syntenic to:\n')
                    if org in trans_syn and ident in trans_syn[org]:
                        for org2,l in trans_syn[org][ident].items():
                            if len(l) == 0:
                                continue
                            if len(l) > 1:
                                f.write(f'\t\t{org2} (potential breakpoints as multiple regions match - likely just artificial partition into multiple regions due to actual sequence insertions/deletions or program parameter/anchor availability/assembly artefacts):\n')
                            else:
                                f.write(f'\t\t{org2}:\n')
                            for reg2 in l:
                                f.write(f'\t\t\tregion {reg2}\n')
                    f.write('-----------------------------------------------------------------------------------------------------------------\n')
                    start = 1
                    end = save_end
                    idx_l = seqids.index(chromo_end)
                    if idx_l > 0:
                        le = seqlen[idx_l-1]
                    else:
                        le = 0
                    len_chr = seqlen[idx_l]
                    end -= le
                    end = min(len_chr,end)
                    elements2 = []
                    for x in coords[org][chromo_end]:
                        if x[0] >= start and x[1] <= end:
                            elements2.append(x[-1])
                    f.write('-----------------------------------------------------------------------------------------------------------------\n')
                    f.write(f'<<< {org}\tregion {ident}\t{chromo_end}\t{start} (== start of chr) \t{end}\tPOTENTIAL CHROMOSOMAL BREAKPOINT {bp_no} REGION 2 >>>\n')
                    f.write(f'\telements: {elements2}\n')
                    f.write(f'\tdirectly pairwise syntenic to:\n')
                    f.write(f'\tdirectly pairwise syntenic to [1 = after markov chain clustering (strength connections between regions == number of anchor matches); 0 = could not hold up MCL]:\n')
                    if org in pw_syn and ident in pw_syn[org]:
                        for org2,l in pw_syn[org][ident].items():
                            if len(l) == 0:
                                continue
                            l_mcl = []
                            if org in pw_syn_mcl and ident in pw_syn_mcl[org] and org2 in pw_syn_mcl[org][ident]:
                                l_mcl = pw_syn_mcl[org][ident][org2]
                            if len(l) > 1:
                                f.write(f'\t\t{org2} (potential breakpoints as multiple regions match - likely just artificial partition into multiple regions due to actual sequence insertions/deletions or program parameter/anchor availability/assembly artefacts):\n')
                            else:
                                f.write(f'\t\t{org2}:\n')
                            for reg2 in l:
                                if reg2 in l_mcl:
                                    f.write(f'\t\t\tregion {reg2} [1] \n')
                                else:
                                    f.write(f'\t\t\tregion {reg2} [0] \n')
                    f.write(f'\tadditionally transitively pairwise syntenic to:\n')
                    if org in trans_syn and ident in trans_syn[org]:
                        for org2,l in trans_syn[org][ident].items():
                            if len(l) == 0:
                                continue
                            if len(l) > 1:
                                f.write(f'\t\t{org2} (potential breakpoints as multiple regions match - likely just artificial partition into multiple regions due to actual sequence insertions/deletions or program parameter/anchor availability/assembly artefacts):\n')
                            else:
                                f.write(f'\t\t{org2}:\n')
                            for reg2 in l:
                                f.write(f'\t\t\tregion {reg2}\n')
                    f.write('-----------------------------------------------------------------------------------------------------------------\n')
                    bp_no += 1
                else:
                    f.write('-----------------------------------------------------------------------------------------------------------------\n')
                    f.write(f'<<< {org}\tregion {ident}\t{chromo}\t{start}\t{end} >>>\n')
                    f.write(f'\telements: {elements}\n')
                    f.write(f'\tdirectly pairwise syntenic to [1 = after markov chain clustering (strength connections between regions == number of anchor matches); 0 = could not hold up MCL]:\n')
                    if org in pw_syn and ident in pw_syn[org]:
                        for org2,l in pw_syn[org][ident].items():
                            if len(l) == 0:
                                continue
                            l_mcl = []
                            if org in pw_syn_mcl and ident in pw_syn_mcl[org] and org2 in pw_syn_mcl[org][ident]:
                                l_mcl = pw_syn_mcl[org][ident][org2]
                            if len(l) > 1:
                                f.write(f'\t\t{org2} (potential breakpoints as multiple regions match - likely just artificial partition into multiple regions due to actual sequence insertions/deletions or program parameter/anchor availability/assembly artefacts):\n')
                            else:
                                f.write(f'\t\t{org2}:\n')
                            for reg2 in l:
                                if reg2 in l_mcl:
                                    f.write(f'\t\t\tregion {reg2} [1] \n')
                                else:
                                    f.write(f'\t\t\tregion {reg2} [0] \n')
                    f.write(f'\tadditionally transitively pairwise syntenic to:\n')
                    if org in trans_syn and ident in trans_syn[org]:
                        for org2,l in trans_syn[org][ident].items():
                            if len(l) == 0:
                                continue
                            if len(l) > 1:
                                f.write(f'\t\t{org2} (potential breakpoints as multiple regions match - likely just artificial partition into multiple regions due to actual sequence insertions/deletions or program parameter/anchor availability/assembly artefacts):\n')
                            else:
                                f.write(f'\t\t{org2}:\n')
                            for reg2 in l:
                                f.write(f'\t\t\tregion {reg2}\n')
                    f.write('-----------------------------------------------------------------------------------------------------------------\n')
            f.write(f'==============================================================================================================\n')
