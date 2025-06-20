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


def check_score(score,org1,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2,strict = 0,max_score = 0):
    if len(dict_org1[i_seq_org1]['regions']) == 0:
        threshold_1 = dict_org1[i_seq_org1]['score']
    else:
        idx_in_scores1 = bisect_left(dict_org1[i_seq_org1]['regions'],start1)
        if start1 == dict_org1[i_seq_org1]['regions'][idx_in_scores1]:
            idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
        idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
        idx_in_scores2 = bisect_left(dict_org1[i_seq_org1]['regions'],end1)
        # for end its okay if its the same...need to take the score for the one before anyway
        idx_in_scores2 -= 1
        maxi_score = 0
        
        for i in range(idx_in_scores1,idx_in_scores2+1):
            if dict_org1[i_seq_org1]['scores_regions'][i] > maxi_score:
                maxi_score = dict_org1[i_seq_org1]['scores_regions'][i]
        threshold_1 = maxi_score
        
    if len(dict_org2[i_seq_org2]['regions']) == 0:
        threshold_2 = dict_org2[i_seq_org2]['score']
    else:
        idx_in_scores1 = bisect_left(dict_org2[i_seq_org2]['regions'],start2)
        if start2 == dict_org2[i_seq_org2]['regions'][idx_in_scores1]:
            idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
        idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
        idx_in_scores2 = bisect_left(dict_org2[i_seq_org2]['regions'],end2)
        # for end its okay if its the same...need to take the score for the one before anyway
        idx_in_scores2 -= 1
        maxi_score = 0

        for i in range(idx_in_scores1,idx_in_scores2+1):
            if dict_org2[i_seq_org2]['scores_regions'][i] > maxi_score:
                maxi_score = dict_org2[i_seq_org2]['scores_regions'][i]
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
        
def parse_clasp_out(orgs_tuple,file,ws):
    
    org1 = orgs_tuple.split('---')[0]
    org2 = orgs_tuple.split('---')[1]
    
    print('parsing clasp_out of {} and {}'.format(org1,org2))
    
    results_org1 = {}
    results_org2 = {}
    check_again_org1 = set()
    check_again_org2 = set()

    for orientation in ['forward','reverse']:
    
        with open(f'{work_dir}/clasp_out_{orientation}/{orgs_tuple}/{file}_awk','r') as clasp_out:

            for line in clasp_out:
              
                line = line.split()
                if line[0] == '#':
                    continue
                
                i = line[1]
                i_seq_org2 = int(i.split('$')[1])
                end_seq_org2 = int(i.split('$')[-1])
                chromo2 = dict_org2[i_seq_org2]['chromosome']
                i = line[2]
                i_seq_org1 = int(i.split('$')[1])
                end_seq_org1 = int(i.split('$')[-1])
                chromo1 = dict_org1[i_seq_org1]['chromosome']

                score = float(line[7])
                
                start2 = int(line[3]) - 1 # blast starts at 1...i have indexed windows as well as regions scores starting at 0
                end2 = int(line[4]) - 1
                start1 = int(line[5]) - 1
                end1 = int(line[6]) - 1 
               
                if orientation == 'reverse':
                    inverted = True
                    # turn around hit coordinates if inverted in order for all of them to refer to forward strand
                    length_anchor2 = end_seq_org2 - i_seq_org2
                    len_hit = end2 - start2
                    start2 = length_anchor2 - end2
                    end2 = start2 + len_hit

                else:
                    inverted = False

                if 'dup' in dict_org1[i_seq_org1] and dict_org1[i_seq_org1]['dup'] == 1 or 'dup' in dict_org2[i_seq_org2] and dict_org2[i_seq_org2]['dup'] == 1:
                    dup = 1
                else:
                    dup = 0
                
                if dup == 1 or check_score(score,org1,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2):

                    len_seq_org1 = end_seq_org1 - i_seq_org1
                    len_seq_org2 = end_seq_org2 - i_seq_org2

                    # is it already in there
                    if i_seq_org1 not in results_org1:
                        results_org1[i_seq_org1] = {}
                        results_org1[i_seq_org1]['meta'] = {}

                        results_org1[i_seq_org1]['matches not considered upon applying stricter score criterion since there are consistency issues'] = {}
                        results_org1[i_seq_org1]['meta']['tolerance range for multiple matches'] = tolerance + dict_org1[i_seq_org1]['end'] - dict_org1[i_seq_org1]['start']
                        results_org1[i_seq_org1]['meta']['word size candidates matching'] = ws
                        results_org1[i_seq_org1]['meta']['multiple matches out of tolerance range'] = 0
                        results_org1[i_seq_org1]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 0
                        results_org1[i_seq_org1]['meta']['not all multiple matches in same order as hits'] = 0
                        results_org1[i_seq_org1]['meta']['multiple matches on different strands'] = 0
                        results_org1[i_seq_org1]['meta']['multiple matches on different chromosomes'] = 0
                        results_org1[i_seq_org1]['meta']['changed matches and metadata upon applying stricter score criterion'] = 0
                        if dup == 1:
                            results_org1[i_seq_org1]['dups_matches'] = {i_seq_org2:{}}
                            results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                            results_org1[i_seq_org1]['dups_matches'][i_seq_org2][f'hit coordinates in {org2} candidate'] = [start2,end2]
                            results_org1[i_seq_org1]['dups_matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate'] = [start1,end1]
                            results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match score'] = score
                        else:
                            results_org1[i_seq_org1]['matches'] = {i_seq_org2:{}}
                            results_org1[i_seq_org1]['matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                            results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in {org2} candidate'] = [start2,end2]
                            results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate'] = [start1,end1]
                            results_org1[i_seq_org1]['matches'][i_seq_org2]['match score'] = score
                    
                    else:
                        if dup == 1:
                            if 'dups_matches' not in results_org1[i_seq_org1]:
                                results_org1[i_seq_org1]['dups_matches'] = {}
                            if i_seq_org2 in results_org1[i_seq_org1]['dups_matches']:
                                score_there = results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match score']
                                if score_there < score:
                                    results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match score'] = score
                                    results_org1[i_seq_org1]['dups_matches'][i_seq_org2][f'hit coordinates in {org2} candidate']  = [start2,end2]
                                    results_org1[i_seq_org1]['dups_matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate']  = [start1,end1]
                                    results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                            else:

                                # append in any case (but collect starts and ends first)
                                
                                results_org1[i_seq_org1]['dups_matches'][i_seq_org2] = {}
                                results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                                results_org1[i_seq_org1]['dups_matches'][i_seq_org2][f'hit coordinates in {org2} candidate'] = [start2,end2]
                                results_org1[i_seq_org1]['dups_matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate'] = [start1,end1]
                                results_org1[i_seq_org1]['dups_matches'][i_seq_org2]['match score'] = score


                        else:
                            if 'matches' not in results_org1[i_seq_org1]:
                                results_org1[i_seq_org1]['matches'] = {}
                            if i_seq_org2 in results_org1[i_seq_org1]['matches']:
                                score_there = results_org1[i_seq_org1]['matches'][i_seq_org2]['match score']
                                if score_there < score:
                                    results_org1[i_seq_org1]['matches'][i_seq_org2]['match score'] = score
                                    results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in {org2} candidate']  = [start2,end2]
                                    results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate']  = [start1,end1]
                                    results_org1[i_seq_org1]['matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                            else:

                                # append in any case (but collect starts and ends first)
                                
                                results_org1[i_seq_org1]['matches'][i_seq_org2] = {}
                                results_org1[i_seq_org1]['matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                                results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in {org2} candidate'] = [start2,end2]
                                results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate'] = [start1,end1]
                                results_org1[i_seq_org1]['matches'][i_seq_org2]['match score'] = score


                    if i_seq_org2 not in results_org2:
                        results_org2[i_seq_org2] = {}
                        results_org2[i_seq_org2]['meta'] = {}
                        results_org2[i_seq_org2]['matches not considered upon applying stricter score criterion since there are consistency issues'] = {}
                        results_org2[i_seq_org2]['meta']['tolerance range for multiple matches'] = tolerance + dict_org2[i_seq_org2]['end'] - dict_org2[i_seq_org2]['start']
                        results_org2[i_seq_org2]['meta']['word size candidates matching'] = ws
                        results_org2[i_seq_org2]['meta']['multiple matches out of tolerance range'] = 0
                        results_org2[i_seq_org2]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 0
                        results_org2[i_seq_org2]['meta']['not all multiple matches in same order as hits'] = 0
                        results_org2[i_seq_org2]['meta']['multiple matches on different strands'] = 0
                        results_org2[i_seq_org2]['meta']['multiple matches on different chromosomes'] = 0
                        results_org2[i_seq_org2]['meta']['changed matches and metadata upon applying stricter score criterion'] = 0

                        if dup == 1:
                            results_org2[i_seq_org2]['dups_matches'] = {i_seq_org1:{}}
                            results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                            results_org2[i_seq_org2]['dups_matches'][i_seq_org1][f'hit coordinates in {org1} candidate'] = [start1,end1]
                            results_org2[i_seq_org2]['dups_matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate'] = [start2,end2]
                            results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match score'] = score

                        else:
                            results_org2[i_seq_org2]['matches'] = {i_seq_org1:{}}
                            results_org2[i_seq_org2]['matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                            results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in {org1} candidate'] = [start1,end1]
                            results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate'] = [start2,end2]
                            results_org2[i_seq_org2]['matches'][i_seq_org1]['match score'] = score
                    else:
                        if dup == 1:
                            if 'dups_matches' not in results_org2[i_seq_org2]:
                                results_org2[i_seq_org2]['dups_matches'] = {}
                            if i_seq_org1 in results_org2[i_seq_org2]['dups_matches']:
                                score_there = results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match score']
                                if score_there < score:
                                    results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match score'] = score
                                    results_org2[i_seq_org2]['dups_matches'][i_seq_org1][f'hit coordinates in {org1} candidate']  = [start1,end1]
                                    results_org2[i_seq_org2]['dups_matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate']  = [start2,end2]
                                    results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                            else:

                                results_org2[i_seq_org2]['dups_matches'][i_seq_org1] = {}
                                results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                                results_org2[i_seq_org2]['dups_matches'][i_seq_org1][f'hit coordinates in {org1} candidate'] = [start1,end1]
                                results_org2[i_seq_org2]['dups_matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate'] = [start2,end2]
                                results_org2[i_seq_org2]['dups_matches'][i_seq_org1]['match score'] = score
                        else:
                            if 'matches' not in results_org2[i_seq_org2]:
                                results_org2[i_seq_org2]['matches'] = {}
                            if i_seq_org1 in results_org2[i_seq_org2]['matches']:
                                score_there = results_org2[i_seq_org2]['matches'][i_seq_org1]['match score']
                                if score_there < score:
                                    results_org2[i_seq_org2]['matches'][i_seq_org1]['match score'] = score
                                    results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in {org1} candidate']  = [start1,end1]
                                    results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate']  = [start2,end2]
                                    results_org2[i_seq_org2]['matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                            else:

                                results_org2[i_seq_org2]['matches'][i_seq_org1] = {}
                                results_org2[i_seq_org2]['matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                                results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in {org1} candidate'] = [start1,end1]
                                results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate'] = [start2,end2]
                                results_org2[i_seq_org2]['matches'][i_seq_org1]['match score'] = score
                            

    return results_org1,results_org2 

def which_chromo(seqlen,seqids,i):
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir_candidates = root + '/anchors/candidates'
    work_dir = argv[-1]

    tolerance = 10000   
 
    orgs_tuple = argv[1]
    file = argv[2]
    ws = argv[3]

    org1 = orgs_tuple.split('---')[0]
    org2 = orgs_tuple.split('---')[1]
    
    with open(anchor_dir_candidates + f'/{org1}','rb') as f:
        dict_org1 = pickle.load(f)
    with open(anchor_dir_candidates + f'/{org2}','rb') as f:
        dict_org2 = pickle.load(f)

    results_org1,results_org2 = parse_clasp_out(orgs_tuple,file,ws)
    # just to know which org without parsing file names
    results_org1_new = {}
    results_org2_new = {}
    results_org1_new[org2] = results_org1
    results_org2_new[org1] = results_org2

    with open(f'{work_dir}/parse_bcamm/{org1}/{orgs_tuple}_{file}','wb') as f:
        pickle.dump(results_org1_new,f)
    with open(f'{work_dir}/parse_bcamm/{org2}/{orgs_tuple}_{file}','wb') as f:
        pickle.dump(results_org2_new,f)
