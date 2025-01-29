#! /usr/bin/env python3


import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
from subprocess import run
from pprint import pprint

def merge(c1,c2):
    c = c1 + c2
    return(list(set(c)))

def get_coords():
    coords = {}
    with open('MCScanX_plus_coords_only_relevant.gff') as f:
        for line in f:
            chromo,name,start,end = line.strip().split()
            if 'CUSTOM' in line:
                coords[name] = org_mapping[chromo.split('chr')[0]]
    return(coords)


if __name__ == "__main__":

    org_mapping,chr_mapping = get_mapping()
    coords = get_coords()

    alloc = {}
    cluster_alloc = {}
    cluster_count = 1

    skip = {}
    lines = 0

    if len(sys.argv) > 1 and sys.argv[2] == 'exclude_tandem_regions':

        with open('MCScanX_plus_coords_only_relevant.collinearity','r') as f:
            for line in f:
                if 'Alignment' in line:
                    if lines != 0 and valid == 0:
                        for l in lines:
                            skip[l] = 1
                    lines = []
                    valid = 0
                else:
                    if line.count('CUSTOM') != 2:
                        valid = 1
                    else:
                        lines.append(line)
        if valid == 0:
            for l in lines:
                skip[l] = 1

    with open('MCScanX_plus_coords_only_relevant.collinearity','r') as f:
        for line in f:
            if 'CUSTOM' in line and line not in skip:
                line = line.split(':')[1].split()
                first = line[0]
                second = line[1]
                c_id_first = False
                if first in alloc:
                    c_id_first = alloc[first]
                c_id_second = False
                if second in alloc:
                    c_id_second = alloc[second]
                
                if c_id_first == False and c_id_second == False:
                    cluster_alloc[cluster_count] = [first,second]
                    alloc[first] = alloc[second] = cluster_count
                    cluster_count += 1
                    c_id_first = c_id_second = cluster_count
                elif c_id_second == False:
                    c_id_second = c_id_first
                    cluster_alloc[c_id_first].append(second)
                    alloc[second] = c_id_first
                elif c_id_first == False:
                    c_id_first = c_id_second
                    cluster_alloc[c_id_second].append(first)
                    alloc[first] = c_id_second
                elif c_id_first != c_id_second:
                    cluster_alloc[cluster_count] = merge(cluster_alloc[c_id_first],cluster_alloc[c_id_second])
                    for i in cluster_alloc[c_id_first]:
                        alloc[i] = cluster_count
                    for i in cluster_alloc[c_id_second]:
                        alloc[i] = cluster_count
                    del cluster_alloc[c_id_first]
                    del cluster_alloc[c_id_second]
                    cluster_count += 1

    run('rm -rf clusters && mkdir -p clusters',shell=True) 
    new_count = 1
    singletons = open('clusters/singletons','w')
    for c in cluster_alloc.values():
        if len(c) < 2:
            for i in c:
                singletons.write(coords[i]+'\t'+i+'\n')
        else:
            with open('clusters/cluster_'+str(new_count),'w') as f:
                for i in c:
                    f.write(coords[i]+'\t'+i+'\n')
            run(f'sort clusters/cluster_{new_count} -o clusters/cluster_{new_count}',shell=True)
            new_count += 1
    singletons.close()
