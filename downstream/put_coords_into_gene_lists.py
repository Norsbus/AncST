#! /usr/bin/env python3

import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
from bisect import bisect_left,bisect_right
from subprocess import run
from os.path import isfile

if __name__ == "__main__":

    margin = int(sys.argv[1])
    org_mapping,chr_mapping = get_mapping()

    relevant_anchors = {}
    limits_per_chr = {}
    eles = {}
    config_write = {}
    skip = {}
    with open('../coords') as f:
        for line in f:
            line = line.strip().split()
            org,chromo,start,end,ori,hit_name = line
            org_orig = org
            chromo_orig = chromo
            org = org_mapping[org]
            if org not in chr_mapping:
                print(f'skipping {org_orig} since it has no anchors')
                continue
            if chromo not in chr_mapping[org]:
                print(f'skipping {chromo} of {org_orig} since it has no anchors (ignoring all respective elements)')
                continue
            chromo_orig = chromo
            chromo = chr_mapping[org][chromo]
            if not isfile(f'adhore_gene_lists/{org}/{chromo}.lst'):
                print(f'skipping {chromo_orig} of {org_orig} since it has no anchors or significant alignment (ignoring all respective elements)')
                continue
            if chromo not in limits_per_chr:
                limits_per_chr[chromo] = []
                eles[chromo] = set()
            if ori == 'forward':
                strand = '+'
            elif ori == 'reverse':
                strand = '-'
            name = f'{chromo}CUSTOMele{start}to{end}{strand}'
            eles[chromo].add((int(start),int(end),name))
            limits_per_chr[chromo].append((max(int(start)-margin,0),int(end)+margin))
            
    for chromo in limits_per_chr:
        limits_per_chr[chromo].sort()
        org = chromo.split('chr')[0]
        if org not in config_write:
            config_write[org] = set()
        config_write[org].add(chromo)

    for chromo,l in limits_per_chr.items():
        org = chromo.split('chr')[0]
        new_l = []
        for enu,x1 in enumerate(l[:-1]):
            new_x = [x1[0],x1[1]]
            for x2 in l[enu+1:]:
                if x2[0] <= x1[1]:
                    new_x[1] = x2[1]
                else:
                    break
            new_l.append((new_x[0],new_x[1]))
        if len(new_l) == 0 or l[-1][1] != new_l[-1][1]:
            new_l.append(l[-1])

        limits_per_chr[chromo] = sorted(new_l)

        anchors = []
        with open(f'adhore_gene_lists/{org}/{chromo}.lst') as f:
            for line in f:
                name = line.strip()
                start = int(name.split('ele')[1].split('to')[0])
                end = int(name.split('ele')[1].split('to')[1][:-1])
                anchors.append((start,end,name))

        anchors = list(sorted(anchors))
        limited_anchors = set()
        # search for start and end of anchor list according to limits
        for limit_start,limit_end in new_l:
            idx_start = bisect_right(anchors,limit_start, key=lambda i: i[0])
            idx_end = bisect_left(anchors,limit_end, key=lambda i: i[0])
            limited_anchors.update(anchors[idx_start:idx_end])

        if len(limited_anchors) == 0:
            if org not in skip:
                skip[org] = {}
            skip[org][chromo] = 1
            print(f'there are no anchors on {chromo} of {org_mapping[org]} for the element starting at {start} ending at {end} in at a margin of {margin}. the element will be ignored for the analysis')
            continue
    
        for x in limited_anchors:
            relevant_anchors[x[2][:-1]] = 1
        limited_anchors = sorted(list(limited_anchors))
        # add coordinates and delete overlapping anchors
        for ele in eles[chromo]:
            start,end,name = ele
            idx_smaller = bisect_left(limited_anchors,start,key=lambda i: i[0]) - 1
            to_del = []
            if idx_smaller >= 0:
                a_start,a_end,a_name = limited_anchors[idx_smaller]
                while idx_smaller >= 0 and (start >= a_start and start <= a_end) or (end >= a_start and end <= a_end):
                    to_del.append(idx_smaller)
                    idx_smaller -= 1
                    if idx_smaller < 0:
                        break
                    a_start,a_end,a_name = limited_anchors[idx_smaller]
            
            idx_bigger = bisect_right(limited_anchors,start,key=lambda i: i[0])
            if idx_bigger < len(limited_anchors):
                a_start,a_end,a_name = limited_anchors[idx_bigger]
                while idx_bigger < len(limited_anchors) and (start >= a_start and start <= a_end) or (end >= a_start and end <= a_end):
                    to_del.append(idx_bigger)
                    idx_bigger += 1
                    if idx_bigger == len(limited_anchors):
                        break
                    a_start,a_end,a_name = limited_anchors[idx_bigger]

                for idx in sorted(to_del,reverse=True):
                    del limited_anchors[idx]

            idx_ele = bisect_left(limited_anchors,start,key=lambda i: i[0])
            limited_anchors.insert(idx_ele,(start,end,name))

        run(f'mkdir -p adhore_gene_lists_only_relevant_{margin} && mkdir -p adhore_gene_lists_only_relevant_{margin}/{org}',shell=True)
        with open(f'adhore_gene_lists_only_relevant_{margin}/{org}/{chromo}.lst','w+') as f:
                l = [x[2]+'\n' for x in sorted(limited_anchors)]
                f.writelines(l)

    out = open(f'adhore_pairwise.table_only_relevant_{margin}','w+')
    with open('adhore_pairwise.table') as f:
        for line in f:
            ele1,ele2 = line.strip().split()
            if ele1 in relevant_anchors and ele2 in relevant_anchors:
                out.write(line)

    for chromo,l in eles.items():
        org = chromo.split('chr')[0]
        for chromo2,l2 in eles.items():
            org2 = chromo2.split('chr')[0]
            if org == org2:
                continue
            for x in l:
                for x2 in l2:
                    out.write(f'{x[2][:-1]}\t{x2[2][:-1]}\n')
    out.close()

    adhore_config = open(f'adhore_config_only_relevant_{margin}','w+')
    first = 1
    for org,chromos in config_write.items():
        if first == 1:
            space = ''
            first = 0
        else:
            space = '\n'
        adhore_config.write(f'{space}genome={org}\n')
        for chromo in chromos:
            if org in skip and chromo in skip[org]:
                continue
            adhore_config.write(f'{chromo} adhore_gene_lists_only_relevant_{margin}/{org}/{chromo}.lst\n')
    template_adhore_config = f"""
blast_table=adhore_pairwise.table_only_relevant_{margin}

output_path=out_only_relevant_{margin}
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
