#! /usr/bin/env python3

import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
import pickle
from pprint import pprint

if __name__ == "__main__":
    org_mapping,chr_mapping = get_mapping()
    origs = {}
    prot_mapping = {}
    with open('origs') as f:
        for line in f:
            if len(line.split()) != 6:
                continue
            org,chromo,start,end,ori,name = line.split()
            prot_name = name.split('_')[-1]
            org = org_mapping[org]
            chromo = chr_mapping[org][chromo]
            if ori == 'forward':
                ori = '+'
            else:
                ori = '-'
            ele = f'{chromo}CUSTOMele{start}to{end}{ori}'
            if org not in origs:
                origs[org] = {}
            if ele not in origs[org]:
                origs[org][ele] = 1
            prot_mapping[prot_name] = ele

    pprint(prot_mapping)

    alloc_to_halos = {}
    halos_where = {}
    for org,bib in origs.items():
        for ele in bib:
            halos_where[ele] = set()


    with open('alloc_pairwise','rb') as f:
        alloc_pairwise = pickle.load(f)


    for org,bib in alloc_pairwise.items():
        if org not in alloc_to_halos:
            alloc_to_halos[org] = {}
        for ele,bib2 in bib.items():
            for org2 in origs:
                if org2 not in bib2:
                    continue
                else:
                    for ele2 in bib2[org2]:
                        if ele2 not in origs[org2]:
                            print(f'{ele2} of {org2} not orig ele of d.mel?')
                            continue
                        else:
                            if ele not in alloc_to_halos[org]:
                                alloc_to_halos[org][ele] = set()
                            alloc_to_halos[org][ele].add(ele2)
                            halos_where[ele2].add(org)
            if ele in alloc_to_halos[org]:
                for orig_ele in alloc_to_halos[org][ele]:
                    for org2,s in bib2.items():
                        for ele2 in s:
                            if org2 not in alloc_to_halos:
                                alloc_to_halos[org2] = {}
                            if ele2 not in alloc_to_halos[org2]:
                                alloc_to_halos[org2][ele2] = set()
                            alloc_to_halos[org2][ele2].add(orig_ele)
                            halos_where[orig_ele].add(org2)

    with open('alloc_to_halos','wb') as f:
        pickle.dump(alloc_to_halos,f)

    org_families = {}
    with open('458_orgs_with_family') as f:
        for line in f:
            if len(line.strip().split()) == 0:
                continue
            org,family = line.strip().split()
            org = org_mapping[org]
            org_families[org] = family
            if family not in org_families:
                org_families[family] = []
            org_families[family].append(org)

    orig_eles_blast_hits = {}
    number_orgs_blast_hits = {}
    blast_hit_characterization = {}
    with open('../utils/coords_both_genomic_and_regional_nono') as f:
        for line in f:
            org,chromo,start,end,ori,original_hit = line.strip().split()
            blast_hit_characterization[original_hit] = (org,chromo,start,end,ori)
            org = org_mapping[org]
            try:
                chromo = chr_mapping[org][chromo]
            except:
                print(f'no anchors for {chromo} of {org}')
            family = org_families[org]
            prot_name = original_hit.split('_')[-1]
            orig_ele = prot_mapping[prot_name]
            if ori == 'forward':
                ori = '+'
            else:
                ori = '-'
            ele = f'{chromo}CUSTOMele{start}to{end}{ori}'
            if orig_ele not in orig_eles_blast_hits:
                orig_eles_blast_hits[orig_ele] = set()
                number_orgs_blast_hits[orig_ele] = {}
            if family not in number_orgs_blast_hits[orig_ele]:
                number_orgs_blast_hits[orig_ele][family] = set()
            number_orgs_blast_hits[orig_ele][family].add(org)
            if org in alloc_to_halos:
                if ele not in alloc_to_halos[org]:
                    orig_eles_blast_hits[orig_ele].add(family)
            else:
                orig_eles_blast_hits[orig_ele].add(family)


    families_origs = {}
    number_orgs = {}
    for orig_ele,l in halos_where.items():
        families_origs[orig_ele] = set()
        number_orgs[orig_ele] = {}
        for org in l:
            family = org_families[org]
            families_origs[orig_ele].add(family)
            if family not in number_orgs[orig_ele]:
                number_orgs[orig_ele][family] = 1
            else:
                number_orgs[orig_ele][family] += 1
    for orig_ele,s in families_origs.items():
        print(f'orig ele: {orig_ele} has orthologs in following families:\n')
        for f in s:
            if f in orig_eles_blast_hits[orig_ele]:
                print(f'\t{f}({number_orgs[orig_ele][f]}) ({len(number_orgs_blast_hits[orig_ele][f])})\n')
            else:
                print(f'\t{f}({number_orgs[orig_ele][f]})\n')
        for f in orig_eles_blast_hits[orig_ele]:
            if f not in s:
                print(f'BLAST HIT in {f} ({len(number_orgs_blast_hits[orig_ele][f])})\n')
