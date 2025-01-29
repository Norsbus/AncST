#!/usr/bin/env python3

from Bio import SeqIO
import pickle
from math import sqrt,floor
from bisect import bisect_left
from subprocess import run
from multiprocessing import Pool
import os,pathlib
from sys import argv
from subprocesses import blast,clasp
from copy import deepcopy

def check_range_and_chromos(i_org,entry,bib,org2):
    tolerance = 10000

    chromos = set()

    distances = []
    
    for m in entry['matches']:
        try:
            chromos.add(bib[m]['chromosome'])
        except:
            print(f'returning entry {m} of {org2} in {i_org} (org1) since chromo was not found in anchor candidates map... (Error)')
            return(entry)
        start = bib[m]['start']
        end = bib[m]['end']
        distances.append((start,end))
    if len(chromos) <= 1:
        entry['meta']['multiple matches on different chromosomes'] = 0
        for m in entry['matches']:
            if org2 in mark_in_others:
                if m in mark_in_others[org2]:
                    if i_org in mark_in_others[org2][m]:
                        mark_in_others[org2][m][i_org][0] = 0
                    else:
                        mark_in_others[org2][m][i_org] = [0,0]
                else:
                    mark_in_others[org2][m] = {i_org:[0,0]}
            else:
                mark_in_others[org2] = {m:{i_org:[0,0]}}
    else:
        entry['meta']['multiple matches on different chromosomes'] = 1
        entry['meta']['multiple matches out of tolerance range'] = 1
        
        for m in entry['matches']:
            if org2 in mark_in_others:
                if m in mark_in_others[org2]:
                    if i_org in mark_in_others[org2][m]:
                        mark_in_others[org2][m][i_org] = [1,1]
                    else:
                        mark_in_others[org2][m][i_org] = [1,1]
                else:
                    mark_in_others[org2][m] = {i_org:[1,1]}
            else:
                mark_in_others[org2] = {m:{i_org:[1,1]}}
        
        # can return as different chromo means out of range too
        return(entry)

    distances = sorted(distances)
    dis = 0
    for i,d in enumerate(distances[:-1]):
        if distances[i+1][0] - d[1] > dis:
            dis = distances[i+1][0] - d[1]

    if dis <= tolerance:
        entry['meta']['multiple matches out of tolerance range'] = 0
        for m in entry['matches']:
            if org2 in mark_in_others:
                if m in mark_in_others[org2]:
                    if i_org in mark_in_others[org2][m]:
                        mark_in_others[org2][m][i_org] = [0,0]
                    else:
                        mark_in_others[org2][m][i_org] = [0,0]
                else:
                    mark_in_others[org2][m] = {i_org:[0,0]}
            else:
                mark_in_others[org2] = {m:{i_org:[0,0]}}
    else:
        entry['meta']['multiple matches out of tolerance range'] = 1

        for m in entry['matches']:
            if org2 in mark_in_others:
                if m in mark_in_others[org2]:
                    if i_org in mark_in_others[org2][m]:
                        mark_in_others[org2][m][i_org] = [0,1]
                    else:
                        mark_in_others[org2][m][i_org] = [0,1]
                else:
                    mark_in_others[org2][m] = {i_org:[0,1]}
            else:
                mark_in_others[org2] = {m:{i_org:[0,1]}}

    return(entry)

def check_score(score,org1,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2,strict = 0,max_score = 0):
    if len(bib[i_seq_org1]['regions']) == 0:
        threshold_1 = bib[i_seq_org1]['score']
    else:
        idx_in_scores1 = bisect_left(bib[i_seq_org1]['regions'],start1)
        if start1 == bib[i_seq_org1]['regions'][idx_in_scores1]:
            idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
        idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
        idx_in_scores2 = bisect_left(bib[i_seq_org1]['regions'],end1)
        # for end its okay if its the same...need to take the score for the one before anyway
        idx_in_scores2 -= 1
        maxi_score = 0
        
        for i in range(idx_in_scores1,idx_in_scores2+1):
            if bib[i_seq_org1]['scores_regions'][i] > maxi_score:
                maxi_score = bib[i_seq_org1]['scores_regions'][i]
        threshold_1 = maxi_score
        
    if len(bibs[org2][i_seq_org2]['regions']) == 0:
        threshold_2 = bibs[org2][i_seq_org2]['score']
    else:
        idx_in_scores1 = bisect_left(bibs[org2][i_seq_org2]['regions'],start2)
        if start2 == bibs[org2][i_seq_org2]['regions'][idx_in_scores1]:
            idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
        idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
        idx_in_scores2 = bisect_left(bibs[org2][i_seq_org2]['regions'],end2)
        # for end its okay if its the same...need to take the score for the one before anyway
        idx_in_scores2 -= 1
        maxi_score = 0

        for i in range(idx_in_scores1,idx_in_scores2+1):
            if bibs[org2][i_seq_org2]['scores_regions'][i] > maxi_score:
                maxi_score = bibs[org2][i_seq_org2]['scores_regions'][i]
        threshold_2 = maxi_score

    if strict == 0:
        if score > max(threshold_1,threshold_2) + sqrt(max(threshold_1,threshold_2)):
            return(True)
        else:
            return(False)
    else:
        if score > max(threshold_1,threshold_2)*2:
            # this is for cases where there is a significant hit somewhere which messes with an actually
            # much more significant hit somewhere else. so here criterion is that its not taken if the hits score 
            # is less than 1/10 of the max score (seen emprirically that it happens and especially for ones 
            # where there is no closest one found in own genome, hence scores of 40)
            if int(threshold_1) == int(40) and int(threshold_2) == int(40):
                if score < max_score/10:
                    return(False)
                else:
                    return(True)
            else:
                return(True)
        else:
            return(False)

def old_check_score(score,org1,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2,strict = 0,max_score = 0):
    
    idx_in_scores1 = bisect_left(bib[i_seq_org1]['regions'],start1)
    if start1 == bib[i_seq_org1]['regions'][idx_in_scores1]:
        idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
    idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
    idx_in_scores2 = bisect_left(bib[i_seq_org1]['regions'],end1)
    # for end its okay if its the same...need to take the score for the one before anyway
    idx_in_scores2 -= 1
    maxi_score = 0
    
    for i in range(idx_in_scores1,idx_in_scores2+1):
        if bib[i_seq_org1]['scores_regions'][i] > maxi_score:
            maxi_score = bib[i_seq_org1]['scores_regions'][i]
    threshold_1 = maxi_score
        
    idx_in_scores1 = bisect_left(bibs[org2][i_seq_org2]['regions'],start2)
    if start2 == bibs[org2][i_seq_org2]['regions'][idx_in_scores1]:
        idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
    idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
    idx_in_scores2 = bisect_left(bibs[org2][i_seq_org2]['regions'],end2)
    # for end its okay if its the same...need to take the score for the one before anyway
    idx_in_scores2 -= 1
    maxi_score = 0

    for i in range(idx_in_scores1,idx_in_scores2+1):
        if bibs[org2][i_seq_org2]['scores_regions'][i] > maxi_score:
            maxi_score = bibs[org2][i_seq_org2]['scores_regions'][i]
    threshold_2 = maxi_score

    if strict == 0:
        if score > max(threshold_1,threshold_2) + sqrt(max(threshold_1,threshold_2)):
            return(True)
        else:
            return(False)
    else:
        if score > max(threshold_1,threshold_2)*2:
            return(True)
        else:
            # this is for cases where there is a significant hit somewhere which messes with an actually
            # much more significant hit somewhere else. so here criterion is that its not taken if the hits score 
            # is less than 1/10 of the max score (seen emprirically that it happens and especially for ones 
            # where there is no closest one found in own genome, hence scores of 40)
            if threshold_1 == 40 and threshold_2 == 40 and score > max_score/10:
                return(True)
            else:
                return(False)
        
def reconcile(old_entry,new_entry,org2,i_seq_org1,inconsistencies):
    
    results_org1 = old_entry
     
    for i,i_bib in new_entry.items():

        if i in results_org1['matches']:
            if i_bib['match score'] > results_org1['matches'][i]['match score']:
                results_org1['matches'][i] = i_bib
        else:
            results_org1['matches'][i] = i_bib
    
    # probably a time-saving approximation would be to only check if not out of range already it has to be checked (like the if statement below)
    # if results_org1['meta']['multiple matches on different chromosomes'] == 0 or results_org1['meta']['multiple matches out of tolerance range'] == 0:
    # BUT: it could be that new candidate is very significant hit and should be considered due to relative score criteria defined in check_score
    results_org1 = check_range_and_chromos(i_seq_org1,results_org1,bibs[org2],org2)
    if results_org1['meta']['multiple matches on different chromosomes'] == 1 or results_org1['meta']['multiple matches out of tolerance range'] == 1:
        to_remove = []
        scores = []
        for i_seq_org2,m in results_org1['matches'].items():
            scores.append(m['match score'])
        max_score = max(scores)
        for i_seq_org2,m in results_org1['matches'].items():
            start2 = m[f'hit coordinates in {org2} candidate'][0]
            end2 = m[f'hit coordinates in {org2} candidate'][1]
            start1 = m[f'hit coordinates in (own) {org} candidate'][0]
            end1 = m[f'hit coordinates in (own) {org} candidate'][1]
            score = m['match score']
            if not check_score(score,org,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2,1,max_score):
                to_remove.append(i_seq_org2)
        if len(to_remove) > 0:
            results_org1['meta']['changed matches and metadata upon applying stricter score criterion'] = 1
            for m in to_remove:
                results_org1['matches not considered upon applying stricter score criterion since there are consistency issues'][m] = deepcopy(results_org1['matches'][m])
                if org2 in inconsistencies:
                    if m in inconsistencies[org2]:
                        inconsistencies[org2][m].add(i_seq_org1)
                    else:
                        inconsistencies[org2][m] = set([i_seq_org1])
                else:
                    inconsistencies[org2] = {m:set([i_seq_org1])}
                del results_org1['matches'][m]
                if len(results_org1['matches']) > 0:
                    results_org1 = check_range_and_chromos(i_seq_org1,results_org1,bibs[org2],org2)
                
    # check for order and invertedness in any case

    if len(results_org1['matches']) > 1:
        
        # inversion
        inv_first = 'first'
        for j,j_bib in results_org1['matches'].items():
            if inv_first == 'first':
                inv_first = j_bib['match is on other strand in other genome']
            else:
                if j_bib['match is on other strand in other genome'] != inv_first:
                    results_org1['meta']['multiple matches on different strands'] = 1

        # if inversion all the same, just check order. if inversion == True -> on other strand, it has to be exactly opposite order

        if results_org1['meta']['multiple matches on different strands'] == 0:

            idx_other = sorted(list(results_org1['matches'].keys()))
            idx_self = []
            for j,j_bib in results_org1['matches'].items():
                idx_self.append((j_bib[f'hit coordinates in (own) {org} candidate'][0],j))
            if inv_first:
                idx_self = [x[1] for x in sorted(idx_self,reverse=True)]
            else:
                idx_self = [x[1] for x in sorted(idx_self)]
            if idx_self != idx_other:
                results_org1['meta']['not all multiple matches in same order as hits'] = 1

        # if the inversions are different, order is hard to define

        else:
            results_org1['meta']['not all multiple matches in same order as hits'] = 1
    return results_org1,inconsistencies

def which_chromo(seqlen,seqids,i):
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

def collect(bib,files,org,inconsistencies):


    # first collect all the parsed bcamm files

    parsed = {}
    for file in files:
        with open(f'{work_dir}/parse_bcamm/{org}/{file}','rb') as f:
            temp = pickle.load(f)
        org2 = list(temp.keys())[0]
        temp = temp[org2]
        if org2 not in parsed:
            parsed[org2] = temp
        else:
            for i in temp:
                if i not in parsed[org2]:
                    parsed[org2][i] = temp[i]
                else:
                    for j,m in temp[i]['matches'].items():
                        if j in parsed[org2][i]['matches']:
                            if m['match score'] > parsed[org2][i]['matches'][j]['match score']:
                                parsed[org2][i]['matches'][j] = m
                        else:
                            parsed[org2][i]['matches'][j] = m

    ### PARALELIZABLE BY CONSIDERING PARIWISE ANCHOR MAPS AND MERGING THEM AFTER ###
    
    # matches_to_del will be shrunken (see below)

    try:

        with open(f'{work_dir}/to_del/{org}/matches_to_del','rb') as f:
            matches_to_del = pickle.load(f)

    except:
    
        matches_to_del = {}

    # now use parsed dict to collect all matches for all other orgs

    for org2,parsed_results in parsed.items():
        if org2 not in bibs:
            with open(anchor_dir_candidates + f'/{org2}','rb') as f:
                bibs[org2] = pickle.load(f)
        for i in parsed_results:
            if i in matches_to_del and org2 in matches_to_del[i]:
                del matches_to_del[i][org2]
            # below is a bit hacky
            # resuing reconcile function with old and new entries being the same (except that new is only considered as its matches)
            if org2 not in bib[i]['matches']:
                bib[i]['matches'][org2],inconsistencies = reconcile(parsed_results[i],parsed_results[i]['matches'],org2,i,inconsistencies)
            else:
                bib[i]['matches'][org2],inconsistencies = reconcile(bib[i]['matches'][org2],parsed_results[i]['matches'],org2,i,inconsistencies)

    # here checking again the matches for which one or more were deleted by filtering process of new windows
    # --> might have changed something about there ranges...
    # only need to do that for the ones not considered anyway because of new results (i expect that most will be in new results
    # as they will be covered by new candidates)
    for i,bib2 in matches_to_del.items():
        # did not delete empty ones above since the check would have to be performed only once here
        if len(matches_to_del[i]) == 0:
            continue
        for org2 in bib2:
            if org2 not in bibs:
                with open(anchor_dir_candidates + f'/{org2}','rb') as f:
                    bibs[org2] = pickle.load(f)
                # again hacky...reusing reconcile function
            bib[i]['matches'][org2],inconsistencies = reconcile(bib[i]['matches'][org2],bib[i]['matches'][org2]['matches'],org2,i,inconsistencies)

    # if you wanted to check all again

    #for i,abib in bib.items():
    #    for org2,bbib in abib['matches'].items():
    #        if org2 not in bibs:
    #            with open(anchor_dir_candidates + f'/{org2}','rb') as f:
    #                bibs[org2] = pickle.load(f)
    #        bib[i]['matches'][org2] = check_range_and_chromos(i,bbib,bibs[org2],org2)

    # checking if a match is in "considered" and "not considered" (has happened but probably along the way of changing code/maps as I cant see a reason why it would from this code)
    # in checks.py now but kept here if you wanna change things automatically (if its always a good match e.g., just delete the ones in "not considered")
    for i,abib in bib.items():
        for org2,bbib in abib['matches'].items():
            for j,cbib in bbib['matches'].items():
                if j in bbib['matches not considered upon applying stricter score criterion since there are consistency issues']:
                    if cbib == bbib['matches not considered upon applying stricter score criterion since there are consistency issues'][j]:
                        print(f'{j} in both match dicts with same match for i={i} org2={org2} j={j}')
                        print(f'match: {cbib}')
                        del bib[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues'][j]

    return(bib,inconsistencies)


if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir_candidates = root + '/anchors/candidates'
    anchor_dir_aligned = root + '/anchors/aligned'
    work_dir = argv[-1]
    
    org = argv[1]

    try:
        with open(anchor_dir_aligned + f'/{org}','rb') as f:
            bib = pickle.load(f)
    except:
        with open(anchor_dir_candidates + f'/{org}','rb') as f:
            bib = pickle.load(f)
    
    bibs = {}

    files = [name for name in os.listdir(f'{work_dir}/parse_bcamm/{org}')]

    # inconsistencies is for moving matches which were moved here to ones not considered because of stricter score criterion
    # mark in others is to mark the multi-matches out of range in other anchor matches -> there it may not correspond to the matches but still useful to have this meta info not having to check every time if its "a really safe anchor"
    inconsistencies,mark_in_others = {},{}
    new_bib,inconsistencies = collect(bib,files,org,inconsistencies)
    
    with open(anchor_dir_aligned + f'/{org}','wb') as f:
        pickle.dump(new_bib,f)
    
    for org2,bib in inconsistencies.items():
        if org2 == 'melanogaster':
            continue
        with open(f'{work_dir}/inconsistencies/{org2}/{org}','wb') as f:
            pickle.dump(bib,f)
    for org2,bib in mark_in_others.items():
        if org2 == 'melanogaster':
            continue
        with open(f'{work_dir}/mark_in_others/{org2}/{org}','wb') as f:
            pickle.dump(bib,f)
