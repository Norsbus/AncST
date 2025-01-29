#! /usr/bin/env python3


if __name__ == "__main__":

    segments_genes = {}
    segments = {}
    segments_multiplicons = {}
    genes_multiplicons = {}

    first = 1
    with open('out_only_relevant/list_elements.txt') as f:
        for line in f:
            if first == 1:
                first = 0
                continue
            # id  segment gene    position    orientation
            id,segment,gene,position,orientation = line.strip().split()
            gene = gene+orientation
            if 'CUSTOM' in gene:
                if gene not in segments_genes:
                    segments_genes[gene] = set()
                if segment not in segments:
                    segments[segment] = set()
                segments_genes[gene].add((segment,position))
                segments[segment].add(gene)

    first = 1
    with open('out_only_relevant/segments.txt') as f:
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

    orthologies = set()
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
            for gene2,p2 in pos.items():
                if gene != gene2 and p==p2:
                    orthologies.add((gene,gene2))

    alloc = {}

    run('rm -rf pairwise_orthologies && mkdir -p pairwise_orthologies',shell=True) 
    for org,bib in alloc.items():
        org = org_mapping[org]
        with open(f'pairwise_orthologies/{org}','w') as f:
            for ele,bib2 in bib.items():
                f.write(f'------------------------ element: {ele} ------------------------\n')
                f.write(f'------------------------ aligns with ------------------------\n')
                for org2,s in bib2.items():
                    org2 = org_mapping[org2]
                    f.write(f'\t --- {org2} ---\n')
                    f.write(f'\t\t {s}\n')
