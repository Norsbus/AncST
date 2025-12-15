#! /usr/bin/env python3

from subprocess import run
import sys
sys.path.append('../../utils/')
from get_mapping import get_mapping
import pickle
from sys import argv

if __name__ == "__main__":

    margin = int(argv[1])
    if len(argv) > 2:
        mode = '_' + argv[2]
    else:
        mode = ''

    org_mapping,chr_mapping = get_mapping()
    segments_genes = {}
    segments_multiplicons = {}
    genes_multiplicons = {}

    ### check which segments dont only contain CUSTOM elements and are therefore valid ###
    first = 1
    segments = {}
    with open(f'out_only_relevant_{margin}{mode}/list_elements.txt') as f:
        for line in f:
            if first == 1:
                first = 0
                continue
            # id  segment gene    position    orientation
            id,segment,gene,position,orientation = line.strip().split()
            if segment not in segments:
                segments[segment] = set()
            segments[segment].add(gene)
    
    valid_segments = set()
    for segment,genes in segments.items():
        n_customs = 0
        for g in genes:
            if 'CUSTOM' in g:
                n_customs += 1
        if n_customs < len(genes) / 2:
            valid_segments.add(segment)

    first = 1
    segments = {}
    with open(f'out_only_relevant_{margin}{mode}/list_elements.txt') as f:
        for line in f:
            if first == 1:
                first = 0
                continue
            # id  segment gene    position    orientation
            id,segment,gene,position,orientation = line.strip().split()
            if segment not in valid_segments:
                continue
            gene = gene+orientation
            if 'CUSTOM' in gene:
                if gene not in segments_genes:
                    segments_genes[gene] = set()
                if segment not in segments:
                    segments[segment] = set()
                segments_genes[gene].add((segment,position))
                segments[segment].add(gene)

    first = 1
    with open(f'out_only_relevant_{margin}{mode}/segments.txt') as f:
        for line in f:
            if first == 1:
                first = 0
                continue
            # id  multiplicon genome  list    first   last    order
            id,multiplicon,genome,list,first,last,order = line.strip().split()
            if id in segments:
                if multiplicon not in genes_multiplicons:
                    genes_multiplicons[multiplicon] = set()
                    segments_multiplicons[multiplicon] = set()
                segments_multiplicons[multiplicon].add(id)
                genes_multiplicons[multiplicon].update(segments[id])

    alloc = {}
    for m,s in genes_multiplicons.items():
        pos = {}
        for gene in s:
            for segment,position in segments_genes[gene]:
                if segment in segments_multiplicons[m]:
                    if gene in pos:
                        print(gene,pos[gene],position)
                        input()
                    pos[gene] = position
        for gene,p in pos.items():
            org1 = gene.split('chr')[0]
            for gene2,p2 in pos.items():
                org2 = gene2.split('chr')[0]
                if gene != gene2 and p == p2 and org1 != org2:
                    if org1 not in alloc:
                        alloc[org1] = {}
                    if org2 not in alloc:
                        alloc[org2] = {}
                    if gene not in alloc[org1]:
                        alloc[org1][gene] = {}
                    if gene2 not in alloc[org2]:
                        alloc[org2][gene2] = {}
                    if org2 not in alloc[org1][gene]:
                        alloc[org1][gene][org2] = set()
                    if org1 not in alloc[org2][gene2]:
                        alloc[org2][gene2][org1] = set()
                    alloc[org1][gene][org2].add(gene2)
                    alloc[org2][gene2][org1].add(gene)
    
                    
    run(f'rm -rf pairwise_orthologies_{margin}{mode} && mkdir -p pairwise_orthologies_{margin}{mode}',shell=True) 
    for org,bib in alloc.items():
        org = org_mapping[org]
        with open(f'pairwise_orthologies_{margin}{mode}/{org}','w') as f:
            for ele,bib2 in bib.items():
                f.write(f'------------------------ element: {ele} ------------------------\n')
                f.write(f'------------------------ aligns with ------------------------\n')
                for org2,s in bib2.items():
                    org2 = org_mapping[org2]
                    f.write(f'\t --- {org2} ---\n')
                    f.write(f'\t\t {s}\n')
    with open(f'alloc_pairwise_{margin}{mode}','wb') as f:
        pickle.dump(alloc,f)
