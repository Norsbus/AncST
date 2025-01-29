#! /usr/bin/env python3

from bisect import bisect_left
import sys
sys.path.append('../utils/')
from get_mapping import get_mapping

def overlap(start,end,coords_to_check):
    idx = bisect_left([x[0] for x in coords_to_check],start)
    if idx == len(coords_to_check):
        if coords_to_check[idx-1][1] > start:
            return True
        else:
            return False
    if end > coords_to_check[idx][0]:
        return True
    return False

def write_gff():

    elements = []
    take = {}
    out = open(f'MCScanX_plus_coords_only_relevant_{margin}.gff','w')
    with open('MCScanX.gff') as f:
        for line in f:
            chromo,name,start,end = line.strip().split()
            start = int(start)
            end = int(end)
            if chromo not in filtered_coords_per_chr or end < min_max_per_chr[chromo][0] or start > min_max_per_chr[chromo][1]:
                continue
            if chromo in filtered_coords_per_chr and overlap(start,end,filtered_coords_per_chr[chromo]):
                continue
            else:
                out.write(f'{chromo}\t{name}\t{start}\t{end}\n')
                if chromo not in take:
                    take[chromo] = []
                take[chromo].append(name)
    for ele in coords:
        chromo,id_seq,start,end = ele
        out.write(f'{chromo}\t{id_seq}\t{start}\t{end}\n')
        elements.append((id_seq,chromo))
    out.close()

    return elements,take

def write_homology(elements,take):

    homology = {}
    out = open(f'MCScanX_plus_coords_only_relevant_{margin}.homology','w')
    with open('MCScanX.homology') as f:
        for line in f:
            id1,id2,score = line.strip().split()
            chromo1 = id1.split('ele')[0]
            chromo2 = id2.split('ele')[0]
            if chromo1 not in take or id1 not in take[chromo1] or chromo2 not in take or id2 not in take[chromo2]:
                continue
            out.write(f'{id1}\t{id2}\t{score}\n')
    for enu,t in enumerate(elements):
        id1,org1 = t
        for t2 in elements[enu+1:]:
            id2,org2 = t2
            if org1 == org2:
                continue
            out.write(f'{id1}\t{id2}\t9999999999\n')
    out.close()

    return 0

if __name__ == "__main__":
    
    margin = int(sys.argv[1])
    org_mapping,chr_mapping = get_mapping()
   
    coords = []
    coords_per_chr = {}
    with open('../coords') as f:
        for line in f:
            line = line.strip().split()
            org,chromo,start,end,ori,hit_name = line
            org = org_mapping[org]
            try:
                chromo = chr_mapping[org][chromo]
            except:
                print(f'there are no anchors on {chromo}...')
                print(f'not putting element of {org_mapping[org]} on {chromo} starting at {start} ending at {end} in MCScanX files (will be neglected for analysis)')
                continue
            if ori == 'forward':
                strand = '+'
            elif ori == 'reverse':
                strand = '-'
            ele_name = f'{chromo}CUSTOMele{start}to{end}{strand}'
            coords.append((chromo,ele_name,int(start),int(end)))
            if chromo not in coords_per_chr:
                coords_per_chr[chromo] = []
            coords_per_chr[chromo].append((int(start),int(end),ele_name))
            
    for chromo in coords_per_chr:
        coords_per_chr[chromo].sort()

    filtered_coords_per_chr = {}
    to_del = set()
    for chromo,l in coords_per_chr.items():
        if len(l) == 1:
            start1,end1,ele_name1 = l[0]
            if chromo not in filtered_coords_per_chr:
                filtered_coords_per_chr[chromo] = [(start1,end1)]
        else:
            for enu,t in enumerate(l[:-1]):
                start1,end1,ele_name1 = t
                if (chromo,ele_name1,start1,end1) in to_del:
                    continue
                for start2,end2,ele_name2 in l[enu+1:]:
                    if end1 < start2:
                        if chromo in filtered_coords_per_chr:
                            filtered_coords_per_chr[chromo].append((start1,end1))
                        else:
                            filtered_coords_per_chr[chromo] = [(start1,end1)]
                        break
                    else:
                        to_del.add((chromo,ele_name2,start2,end2))
            if chromo not in filtered_coords_per_chr:
                lens = [x[1] - x[0] for x in l]
                idx = lens.index(max(lens))
                filtered_coords_per_chr[chromo] = l[idx:idx+1]
            if l[-2][1] < filtered_coords_per_chr[chromo][-1][0]:
                filtered_coords_per_chr[chromo].append((l[-1][0],l[-1][1]))
        filtered_coords_per_chr[chromo].sort()
    for ele in to_del:
        coords.remove(ele)
    
    min_max_per_chr = {}
    for chromo,l in filtered_coords_per_chr.items():
        min_max_per_chr[chromo] = (l[0][0]-margin,l[-1][1]+margin)

    elements,take = write_gff()
    write_homology(elements,take)
