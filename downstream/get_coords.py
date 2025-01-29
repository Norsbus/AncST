#! /usr/bin/env python3

import pickle
from bisect import bisect_left
import multiprocessing as mp

def which_chromo(seqlen,seqids,i):
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

def get_coords(pivot,orgs,old_coords):
    p_org,p_start,p_end = pivot
    with open(f'{aligned_path}/{p_org}','rb') as f:
        aligned = pickle.load(f)
    iss = sorted(aligned.keys())
    iss = iss[bisect_left(iss,p_start):bisect_left(iss,p_end)]
    for i in iss:
        bib = aligned[i]
        for org in orgs:
            if org in bib['matches']:
                matches = list(bib['matches'][org]['matches'].keys())
                if org in old_coords:
                    old_coords[org].update(matches)
                else:
                    old_coords[org] = set(matches)
    return old_coords

if __name__ == "__main__":
    aligned_path = 'HOMEDIR/andrena/stable_synteny/anchors/aligned/'
    candidates_path = 'HOMEDIR/andrena/stable_synteny/anchors/candidates/'
    small_meta_path = 'HOMEDIR/andrena/stable_synteny/utils/small_meta/'
    orgs = []
    small_meta = {}
    with open('orgs','r') as f:
        for line in f:
            orgs.append(line.strip())

    for org in orgs:
        with open(f'{small_meta_path}/{org}','rb') as f:
            small_meta[org] = pickle.load(f)

    no_iter = 3
    no_cores = 24

    margin_to_be_new = 50000
    margin_anchors = 50000

    elements = []
    with open('coords') as f:
        for line in f:
            org,chromo,start,end,ori = line.strip().split()
            start = int(start)
            end = int(end)
            seqids,seqlen = small_meta[org]
            idx_l = seqids.index(chromo)
            if idx_l > 0:
                l = seqlen[idx_l-1]
            else:
                l = 0
            l_end = seqlen[idx_l]
            start += l
            end += l
            elements.append((org,max(start-margin_anchors,l),min(end+margin_anchors,l_end)))
    
    old_coords = {}
        
    for i in range(no_iter):
        print(f'iter no {i} elements of length {len(elements)}')
        with mp.Pool(processes=no_cores) as p:
            res = p.starmap(get_coords,[(pivot,orgs,old_coords) for pivot in elements])
        elements = []
        for new_coords in res:
            for org,s in new_coords.items():
                if org not in old_coords:
                    old_coords[org] = set()
                for i in s:
                    seqdids,seqlen = small_meta[org]
                    chromo = which_chromo(seqlen,seqdids,i)
                    idx_l = seqids.index(chromo)
                    if idx_l > 0:
                        l = seqlen[idx_l-1]
                    else:
                        l = 0
                    l_end = seqlen[idx_l]
                    if i < min(old_coords[org]) and min(old_coords[org]) - i < margin_to_be_new:
                        elements.append((org,max(i-margin_anchors,l),min(i+margin_anchors,l_end)))
                    elif i > max(old_coords[org]) and i - max(old_coords[org]) < margin_to_be_new:
                        elements.append((max(org,i-margin_anchors,l),min(i+margin_anchors,l_end)))
                old_coords[org].update(s)
        elements = list(set(elements))
        pprint(old_coords)
        print('---')
    with open('syntenic_regions','w') as f:
        for org,s in old_coords.items():
            with open(f'{candidates_path}/{org}','rb') as f2:
                candidates = pickle.load(f2)
            to_write = {}
            for i in s:
                chromo = candidates[i]['chromosome']
                start = candidates[i]['start']
                end = candidates[i]['end']
                if chromo not in to_write:
                    to_write[chromo] = []
                to_write[chromo].append((start,end))
            for chromo,regions in to_write.items():
                f.write(f'{org}\t{chromo}\t{min([r[0] for r in regions])}\t{max([r[1] for r in regions])}\n')
