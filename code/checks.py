#! /usr/bin/env python3

from pprint import pprint
import pickle,pathlib,os
from sys import argv
from bisect import bisect_left
from statistics import mean,median,stdev

# just to check that no anchors are overlapping each other or chromosome ends

def check_correspondence_candidates_aligned_maps(org):

    with open(f'{root}/anchors/candidates/{org}','rb') as f_1:
        candidates = pickle.load(f_1)
    with open(f'{root}/anchors/aligned/{org}','rb') as f_1:
        aligned = pickle.load(f_1)
    for i in candidates:
        if i not in aligned:
            print(f'CASE 1: {i} in candidates but not aligned')
    for i in aligned:
        if i not in candidates:
            print(f'CASE 2: {i} in aligned but not in candidates')
    

    print('checks correspondence done')

    return(0)


def check_compared(org):

    other_orgs = []
    with open(f'{work_dir}/orgs','r') as f:
        for line in f:
            other_orgs.append(line.strip())
    other_orgs.remove(org)
    compared_maps = os.listdir(f'{root}/anchors/compared/{org}')
    for cmfile in compared_maps:
        with open(f'{root}/anchors/compared/{org}/{cmfile}','rb') as f_1:
            cm = pickle.load(f_1)
        org2 = cmfile.split('with_')[1].strip()
        if org2 not in other_orgs:
            continue    
        with open(f'{root}/anchors/candidates/{org2}','rb') as f_1:
            candidates = pickle.load(f_1)
        for i,end in cm.items():
            if i not in candidates:
                print(f'{i} of {org2} not in {org} {cmfile} compared bib')
            else:
                if end != candidates[i]['end']:
                    print(f'end {end} from update compared not the ones in candidates map of {org2}: candidates[i]["end"]')
    
    print('checks compard done')
            
    return(0)


def check_sanity(org):

    lens = []

    with open(f'{root}/anchors/aligned/{org}','rb') as f_1:
        am1 = pickle.load(f_1)

    with open(f'{root}/utils/small_meta/{org}','rb') as f: 
        seqids,seqlen = pickle.load(f)
    
    sorted_i = sorted(am1.keys())
    for ii,i in enumerate(sorted_i):
        start = i
        anchor = am1[i]
        length = anchor['end'] - anchor['start']
        lens.append(length)
        end = i + length
        if end < start:
            print(f'bigger start than end for {i}')
        idx = bisect_left(seqlen,start)
        if start == seqlen[idx]:
            idx += 1
        chromo_start = seqids[idx]
        idx = bisect_left(seqlen,end)
        if end == seqlen[idx]:
            idx += 1
        chromo_end = seqids[idx]
        if chromo_start != chromo_end:
            print(f'start chromosome is different than end chromosome for {i}')
        if ii < len(sorted_i) - 1:
            if end >= sorted_i[ii+1]:
                print(f'end of {i} ({end}) is bigger than next start for ({sorted_i[ii+1]}) -> anchors overlap')
    print(f'mean lengths anchors: {mean(lens)}')
    print(f'median lengths anchors: {median(lens)}')
    print(f'stdev lengths anchors: {stdev(lens)}')
    print('checks sanity done')
    return(0)
            
        

def check_aligned(org):
    with open(f'{root}/anchors/aligned/{org}','rb') as f:
        am1 = pickle.load(f)
    for i,abib in am1.items():
        for org2,bbib in abib['matches'].items():
            for j,cbib in bbib['matches'].items():
                if j in bbib['matches not considered upon applying stricter score criterion since there are consistency issues']:
                    if cbib == bbib['matches not considered upon applying stricter score criterion since there are consistency issues'][j]:
                        print(f'{j} in both match dicts with same match for i={i} org2={org2} j={j}')
                        print(f'match: {cbib}')
                        del am1[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues'][j]
    other_orgs = []
    with open(f'{work_dir}/orgs','r') as f:
        for line in f:
            other_orgs.append(line.strip())
    other_orgs.remove(org)

    for org2 in other_orgs:
        print(f'--- {org2} ---')
        with open(f'{root}/anchors/aligned/{org2}','rb') as f:
            am2 = pickle.load(f)
        for i,bib in am1.items():
            if org2 in bib['matches'].keys():
                for j in bib['matches'][org2]['matches']:
                    meta = bib['matches'][org2]['meta']
                    if j not in am2:
                        print(f'anchor {i} of {org} : match {j} not even in second {org2}\'s map:')
                        continue
                    if org not in list(am2[j]['matches'].keys()):
                        print(f'{org} not in second {org2}\'s maps matches:')
                        print(i)
                        pprint(am1[i]['matches'][org2])
                        #print(j)
                        #pprint(am2[j]['matches'])
                    elif i not in list(am2[j]['matches'][org]['matches'].keys()):
                        print(f'{i} not in second {org2}\'s maps matches for {org}:')
                        print(i)
                        pprint(am1[i]['matches'][org2])
                        print(j)
                        pprint(am2[j]['matches'][org])
                    else:
                        score1 = round(bib['matches'][org2]['matches'][j]['match score'])
                        score2 = round(am2[j]['matches'][org]['matches'][i]['match score'])
                        if score1 != score2:
                            print('match score wrong')
                            print(score1,score2)
                            print(i)
                            pprint(am1[i]['matches'][org2])
                            print(j)
                            pprint(am2[j]['matches'][org])
                        start11,end11 = bib['matches'][org2]['matches'][j][f'hit coordinates in (own) {org} candidate']
                        start12,end12 = am2[j]['matches'][org]['matches'][i][f'hit coordinates in {org} candidate']
                        if start11 != start12 or end11 != end12:
                            print('start/end coordinates wrong')
                            print(start11,end11)
                            print(start12,end12)
                            print(i)
                            pprint(am1[i])
                            print(j)
                            pprint(am2[j])
                        start21,end21 = bib['matches'][org2]['matches'][j][f'hit coordinates in {org2} candidate']
                        start22,end22 = am2[j]['matches'][org]['matches'][i][f'hit coordinates in (own) {org2} candidate']
                        if start21 != start22 or end21 != end22:
                            print('start/end coordinates wrong')
                            print(start21,end21)
                            print(start22,end22)
                            print(i)
                            pprint(am1[i])
                            print(j)
                            pprint(am2[j])
                        if (meta['multiple matches out of tolerance range'] == 1 and am2[j]['matches'][org]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 0) or (am2[j]['matches'][org]['meta']['multiple matches out of tolerance range'] == 1 and meta['matches have ambiguous matches (tolerance/chromosome out of range)'] == 0): 
                            print('multiple matches out of tolerance range wrong')
                            pprint(meta['multiple matches out of tolerance range'])
                            pprint(am2[j]['matches'][org]['meta']['multiple matches out of tolerance range'])
                            print(i)
                            pprint(am1[i]['matches'][org2])
                            print(j)
                            pprint(am2[j]['matches'][org])
                        if meta['multiple matches on different chromosomes'] != am2[j]['matches'][org]['meta']['multiple matches on different chromosomes']:
                            if (meta['multiple matches on different chromosomes'] == 1 and am2[j]['matches'][org]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] != 1) or\
                            (meta['matches have ambiguous matches (tolerance/chromosome out of range)'] != 1 and am2[j]['matches'][org]['meta']['multiple matches on different chromosomes'] == 1):
                                print('multiple matches on different chromosomes wrong')
                                pprint(meta['multiple matches on different chromosomes'])
                                pprint(am2[j]['matches'][org]['meta']['multiple matches on different chromosomes'])
                                print(i)
                                pprint(am1[i]['matches'][org2])
                                print(j)
                                pprint(am2[j]['matches'][org])
    print('checks matches done')
    return(0)

if __name__ == '__main__':

    org = argv[1]
    print(f'--- {org} ---')
    work_dir = argv[2]
    anchor_dir = argv[3]
    root = str(pathlib.Path(__file__).parents[1])
  
    check_correspondence_candidates_aligned_maps(org)
    check_sanity(org)
    #check_compared(org)
    check_aligned(org)
    print('checks done')
