#! /usr/bin/env python3

import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
import pickle
from pprint import pprint
from os import listdir

if __name__ == "__main__":
    margin = int(sys.argv[1])
    if len(sys.argv) > 2:
        mode = '_' + sys.argv[2]
    else:
        mode = ''
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
                ori2 = '-'
            else:
                ori = '-'
                ori2 = '+'
            ele = f'{chromo}CUSTOMele{start}to{end}'
            if org not in origs:
                origs[org] = {}
            if ele not in origs[org]:
                origs[org][ele] = 1

            prot_mapping[prot_name] = ele
    
    small_meta = {}
    orgs = listdir('HOMEDIR/insects/stable_synteny/utils/small_meta/')
    for org in orgs:
        with open(f'HOMEDIR/insects/stable_synteny/utils/small_meta/{org}','rb') as f:
            #seqids,seqlen
            small_meta[org] = pickle.load(f)

    alloc_to_halos = {}
    halos_where = {}
    for org,bib in origs.items():
        for ele in bib:
            halos_where[ele] = set()

    with open(f'alloc_pairwise_{margin}{mode}','rb') as f:
        alloc_pairwise = pickle.load(f)
    
    #for ele,bib in alloc_pairwise['318org'].items():
    #    if ele[-1] in ['+','-']:
    #        ele = ele[:-1]
    #    if ele not in origs['318org']:
    #        print(f'{ele} not orig ele of d.mel?')
    #        continue
    #    for org2,set2 in bib.items():
    #        for ele2 in set2:
    #            halos_where[ele].add(org2)
    #            if org2 not in alloc_to_halos:
    #                alloc_to_halos[org2] = {}
    #            if ele2 not in alloc_to_halos[org2]:
    #                alloc_to_halos[org2][ele2] = set()
    #            alloc_to_halos[org2][ele2].add(ele)

    which_Syrphidae = set()
    #pprint(alloc_pairwise['200org']['200orgchr2CUSTOMele2888685to2888885'])
    #exit(0)
    for org,bib in alloc_pairwise.items():
        if org not in alloc_to_halos:
            alloc_to_halos[org] = {}
        for ele,bib2 in bib.items():
            for org2 in origs:
                if org2 not in bib2:
                    continue
                else:
                    for ele2 in bib2[org2]:
                        if ele2[-1] in ['+','-']:
                            ele2 = ele2[:-1]
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
                            if orig_ele == '318orgchr2CUSTOMele1517735to1518061':
                                if org2 == '427org':# and ele2 == '427orgchr13CUSTOMele7875180to7875311':
                                    found = 0
                                    new = set([(org,ele)])
                                    while found == 0:
                                        new_new = set()
                                        for org3,ele3 in new:
                                            if '318org' in alloc_pairwise[org3][ele3]:
                                                found = 1
                                                print(org,ele,org2,ele2,org3,ele3,orig_ele)
                                            else:
                                                for o,s in alloc_pairwise[org3][ele3].items():
                                                    for e in s:
                                                        new_new.add((o,e))
                                        new = new_new
    exit(0)
    #with open('alloc_to_halos','wb') as f:
    #    pickle.dump(alloc_to_halos,f)

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
    
    ### DEBUG ###
    #for org,bib in alloc_to_halos.items():
    #    #if org_families[org] == 'Drosophilidae':
    #    #    continue
    #    for ele,s in bib.items():
    #        print(f'{org} ele {ele} aligns with HALOS:')
    #        print(s)
    #        print('and aligns with other elements:')
    #        for org22,s22 in alloc_pairwise[org][ele].items():
    #            if org_families[org22] == 'Drosophilidae':
    #                continue
    #            print(org22,s22)
    #            for x in s22:
    #                if org22 in alloc_to_halos and x in alloc_to_halos[org22]:
    #                    print(alloc_to_halos[org22][x])
    #        input()
    ### DEBUG ###

    orig_eles_blast_hits = {}
    number_orgs_blast_hits = {}
    blast_hit_characterization = {}
    which_blast_hits = {}
    best_hits = {}
    with open('../utils/coords') as f:
        for line in f:
            org,chromo,start,end,ori,original_hit,hits,score = line.strip().split('\t')
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
            if '$' in hits:
                bh = [x for x in hits.split('$')]
                best_hits[original_hit] = bh
                #max_score = 0
                #for x in bh:
                #    h,score = x.split(':')
                #    score = float(score)
                #    if score > max_score:
                #        max_score = score
                #        best_hits[original_hit] = h
            else:
                best_hits[original_hit] = hits
            if orig_ele not in orig_eles_blast_hits:
                which_blast_hits[orig_ele] = {}
                orig_eles_blast_hits[orig_ele] = set()
            if orig_ele not in number_orgs_blast_hits:
                number_orgs_blast_hits[orig_ele] = {}
            if family not in number_orgs_blast_hits[orig_ele]:
                number_orgs_blast_hits[orig_ele][family] = set()
            number_orgs_blast_hits[orig_ele][family].add(org)
            if org in alloc_to_halos:
                if ele not in alloc_to_halos[org]:
                    if family not in which_blast_hits[orig_ele]:
                        which_blast_hits[orig_ele][family] = set()
                    which_blast_hits[orig_ele][family].add(original_hit)
                    orig_eles_blast_hits[orig_ele].add(family)
            else:
                if family not in which_blast_hits[orig_ele]:
                    which_blast_hits[orig_ele][family] = set()
                which_blast_hits[orig_ele][family].add(original_hit)
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
                number_orgs[orig_ele][family] = set()
            number_orgs[orig_ele][family].add(org)
            if family == 'Dolichopodidae' and orig_ele == '318orgchr2CUSTOMele1517735to1518061':
                which_Syrphidae.add(org)

    for org in which_Syrphidae:
        print('----')
        print(f'{org}:{org_mapping[org]}')
    exit(0)

    for orig_ele,s in families_origs.items():
        print(f'orig ele: {orig_ele} has orthologs in following families:\n')
        for f in s:
            if f in number_orgs_blast_hits[orig_ele]:
                print(f'\t{f}({len(number_orgs[orig_ele][f])}) ({len(number_orgs_blast_hits[orig_ele][f])})\n')
            else:
                print(f'\t{f}({len(number_orgs[orig_ele][f])})\n')
        if orig_ele in orig_eles_blast_hits:
            for f in orig_eles_blast_hits[orig_ele]:
                if f not in s:
                    print(f'BLAST HIT in {f} ({len(number_orgs_blast_hits[orig_ele][f])})\n')
                    for oh in which_blast_hits[orig_ele][f]:
                        if ':' in best_hits[oh]:
                            best_hits[oh] = best_hits[oh].split(':')[0]
                        org,chromo,start,end,ori = blast_hit_characterization[oh]
                        length = int(end) - int(start)
                        if length < 0:
                            print(f'{oh} {org} {chromo} {start} {end} {ori}')
                            input()
                        org = org_mapping[org]
                        if chromo not in chr_mapping[org]:
                            anchors = 'no'
                        else:
                            anchors = 'yes'
                        idx_sm = small_meta[org_mapping[org]][0].index(chromo)
                        if idx_sm == 0:
                            length_chr = small_meta[org_mapping[org]][1][idx_sm]
                        else:
                            length_chr = small_meta[org_mapping[org]][1][idx_sm] - small_meta[org_mapping[org]][1][idx_sm-1]
                            if ':' in best_hits[oh]:
                                best_hits[oh] = best_hits[oh].split(':')[0]
                        print(f'\tlen_hit: {length} org: {org} chr: {chromo} len_chr: {length_chr} anchors: {anchors} best_hit:{best_hits[oh]}') #{prot_mapping[best_hits[oh]]}')
