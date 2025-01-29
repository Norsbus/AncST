#!/usr/bin/env python3

from subprocess import run,check_output
import os,pickle
from time import sleep
from sys import argv
from subprocesses import blast,clasp
import pathlib
from Bio import SeqIO
from math import sqrt,floor
from bisect import bisect_left
import os,pathlib
import multiprocessing as mp

def check_score(score,org1,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2,strict = 0,max_score = 0):
    if len(dict_orgs[org1][i_seq_org1]['regions']) == 0:
        threshold_1 = dict_orgs[org1][i_seq_org1]['score']
    else:
        idx_in_scores1 = bisect_left(dict_orgs[org1][i_seq_org1]['regions'],start1)
        if start1 == dict_orgs[org1][i_seq_org1]['regions'][idx_in_scores1]:
            idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
        idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
        idx_in_scores2 = bisect_left(dict_orgs[org1][i_seq_org1]['regions'],end1)
        # for end its okay if its the same...need to take the score for the one before anyway
        idx_in_scores2 -= 1
        maxi_score = 0
        
        for i in range(idx_in_scores1,idx_in_scores2+1):
            if dict_orgs[org1][i_seq_org1]['scores_regions'][i] > maxi_score:
                maxi_score = dict_orgs[org1][i_seq_org1]['scores_regions'][i]
        threshold_1 = maxi_score
        
    if len(dict_orgs[org2][i_seq_org2]['regions']) == 0:
        threshold_2 = dict_orgs[org2][i_seq_org2]['score']
    else:
        idx_in_scores1 = bisect_left(dict_orgs[org2][i_seq_org2]['regions'],start2)
        if start2 == dict_orgs[org2][i_seq_org2]['regions'][idx_in_scores1]:
            idx_in_scores1 += 1 # +1 if the same because then the score should be the one of the "next" region
        idx_in_scores1 -= 1 # -1 since scores have one less members of list than regions
        idx_in_scores2 = bisect_left(dict_orgs[org2][i_seq_org2]['regions'],end2)
        # for end its okay if its the same...need to take the score for the one before anyway
        idx_in_scores2 -= 1
        maxi_score = 0

        for i in range(idx_in_scores1,idx_in_scores2+1):
            if dict_orgs[org2][i_seq_org2]['scores_regions'][i] > maxi_score:
                maxi_score = dict_orgs[org2][i_seq_org2]['scores_regions'][i]
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
            if threshold_1 == 20 and threshold_2 == 20 and score > max_score/10:
                return(True)
            else:
                return(False)
        
def parse_clasp_out(args):

    orgs_tuple,file,ws,work_dir = args.split('$$$')
    
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
                chromo2 = dict_orgs[org2][i_seq_org2]['chromosome']
                i = line[2]
                i_seq_org1 = int(i.split('$')[1])
                end_seq_org1 = int(i.split('$')[-1])
                chromo1 = dict_orgs[org1][i_seq_org1]['chromosome']

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
                
                if check_score(score,org1,i_seq_org1,start1,end1,org2,i_seq_org2,start2,end2):

                    to_remove_in_1 = []
                    to_remove_in_2 = []
                    
                    len_seq_org1 = end_seq_org1 - i_seq_org1
                    len_seq_org2 = end_seq_org2 - i_seq_org2
                   
                    # is it already in there
                    if i_seq_org1 not in results_org1:
                        results_org1[i_seq_org1] = {}
                        results_org1[i_seq_org1]['meta'] = {}
                        results_org1[i_seq_org1]['matches not considered upon applying stricter score criterion since there are consistency issues'] = {}
                        results_org1[i_seq_org1]['meta']['tolerance range for multiple matches'] = tolerance + dict_orgs[org1][i_seq_org1]['end'] - dict_orgs[org1][i_seq_org1]['start']
                        results_org1[i_seq_org1]['meta']['word size candidates matching'] = ws
                        results_org1[i_seq_org1]['meta']['multiple matches out of tolerance range'] = 0
                        results_org1[i_seq_org1]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 0
                        results_org1[i_seq_org1]['meta']['not all multiple matches in same order as hits'] = 0
                        results_org1[i_seq_org1]['meta']['multiple matches on different strands'] = 0
                        results_org1[i_seq_org1]['meta']['multiple matches on different chromosomes'] = 0
                        results_org1[i_seq_org1]['meta']['changed matches and metadata upon applying stricter score criterion'] = 0
    
                        
                        results_org1[i_seq_org1]['matches'] = {i_seq_org2:{}}
                        results_org1[i_seq_org1]['matches'][i_seq_org2]['match is on other strand in other genome'] = inverted
                        results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in {org2} candidate'] = [start2,end2]
                        results_org1[i_seq_org1]['matches'][i_seq_org2][f'hit coordinates in (own) {org1} candidate'] = [start1,end1]
                        results_org1[i_seq_org1]['matches'][i_seq_org2]['match score'] = score
                    
                    else:
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
                        results_org2[i_seq_org2]['meta']['tolerance range for multiple matches'] = tolerance + dict_orgs[org2][i_seq_org2]['end'] - dict_orgs[org2][i_seq_org2]['start']
                        results_org2[i_seq_org2]['meta']['word size candidates matching'] = ws
                        results_org2[i_seq_org2]['meta']['multiple matches out of tolerance range'] = 0
                        results_org2[i_seq_org2]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 0
                        results_org2[i_seq_org2]['meta']['not all multiple matches in same order as hits'] = 0
                        results_org2[i_seq_org2]['meta']['multiple matches on different strands'] = 0
                        results_org2[i_seq_org2]['meta']['multiple matches on different chromosomes'] = 0
                        results_org2[i_seq_org2]['meta']['changed matches and metadata upon applying stricter score criterion'] = 0

                        results_org2[i_seq_org2]['matches'] = {i_seq_org1:{}}
                        results_org2[i_seq_org2]['matches'][i_seq_org1]['match is on other strand in other genome'] = inverted
                        results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in {org1} candidate'] = [start1,end1]
                        results_org2[i_seq_org2]['matches'][i_seq_org1][f'hit coordinates in (own) {org2} candidate'] = [start2,end2]
                        results_org2[i_seq_org2]['matches'][i_seq_org1]['match score'] = score
                    else:
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

def blast_clasp_anchor_map_making(args):
    org_tuple,file,ws,work_dir = args.split('$$$')
    org1 = org_tuple.split('---')[0]
    org2 = org_tuple.split('---')[1]
    # blast forward
    blast(f'{work_dir}/blastdbs/anchor_candidates_'+org_tuple+'_forward',f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward_split/forward.{file}',f'{work_dir}/blast_out_forward/{org_tuple}/{file}',ws)
    # blast reverse 
    blast(f'{work_dir}/blastdbs/anchor_candidates_'+org_tuple+'_forward',f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse_split/reverse.{file}',f'{work_dir}/blast_out_reverse/{org_tuple}/{file}',ws)
    # clasp forward
    clasp(f'{work_dir}/blast_out_forward/{org_tuple}/{file}',f'{work_dir}/clasp_out_forward/{org_tuple}/{file}',2,0.1)
    # clasp reverse
    clasp(f'{work_dir}/blast_out_reverse/{org_tuple}/{file}',f'{work_dir}/clasp_out_reverse/{org_tuple}/{file}',2,0.1)
    print(f'FILE11 : {work_dir}/clasp_out_forward/{org_tuple}/{file}')
    print(f'FILE12 : {work_dir}/clasp_out_forward/{org_tuple}/{file}_awk')
    print(f'FILE21 : {work_dir}/clasp_out_reverse/{org_tuple}/{file}')
    print(f'FILE22 : {work_dir}/clasp_out_forward/{org_tuple}/{file}_awk')
    run(f'awk \'{{{{ if ($1 != "#" && $8>40) print }}}}\' {work_dir}/clasp_out_forward/{org_tuple}/{file} > {work_dir}/clasp_out_forward/{org_tuple}/{file}_awk && awk \'{{{{ if ($1 != "#" && $8>40) print }}}}\' {work_dir}/clasp_out_reverse/{org_tuple}/{file} > {work_dir}/clasp_out_reverse/{org_tuple}/{file}_awk',shell=True)

    return(0)


def which_chromo(seqlen,seqids,i):
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

def exe1(args):


    return 0

def exe(args):
    
    org_tuple,file,ws,work_dir = args.split('$$$')
    org1 = org_tuple.split('---')[0]
    org2 = org_tuple.split('---')[1]
    
    blast_clasp_anchor_map_making(args)

    results_org1,results_org2 = parse_clasp_out(args)
    # just to know which org without parsing file names
    results_org1_new = {}
    results_org2_new = {}
    results_org1_new[org2] = results_org1
    results_org2_new[org1] = results_org2

    with open(f'{work_dir}/parse_bcamm/{org1}/{org_tuple}_{file}','wb') as f:
        pickle.dump(results_org1_new,f)
    with open(f'{work_dir}/parse_bcamm/{org2}/{org_tuple}_{file}','wb') as f:
        pickle.dump(results_org2_new,f)

    return 0

if __name__ == "__main__":
    
    orgs = []
    with open(f"orgs","r") as f:
        for line in f:
            org = line.strip()
            orgs.append(org)

    anchor_dir_candidates = '../anchors/candidates'
    work_dir = '.'

    dict_orgs = {}
    for org in orgs:
        with open(anchor_dir_candidates+f'/{org}','rb') as f:
            dict_orgs[org] = pickle.load(f)

    tolerance = 10000   
 
    ws = 11


    orgs_tuples = []
    files = []

    for i,org1 in enumerate(orgs):
        for org2 in orgs[i+1:]:
            orgs_tuple = org1 + '---' + org2
            if not os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'):
                print('wroking')
                print(os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'))
                orgs_tuple = org2 + '---' + org1
                print(os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'))
                print(orgs_tuple)
                input()
            chunks = [name.split('forward.')[1] for name in os.listdir(f'{work_dir}/sequences_to_compare/{orgs_tuple}/{org2}/forward_split')]
            for f in chunks:
                if os.path.isfile(f"{work_dir}/touch/{orgs_tuple}_parse_bcamm_done_{f}"):
                    print(f'{org_tuple} {f} already done')
                    continue
                orgs_tuples.append(orgs_tuple)
                files.append(f)
    otto = [f'{org_tuple}$$${f}$$$11$$$.' for org_tuple,f in zip(orgs_tuples,files)]
    with mp.Pool(processes=12) as pool:
        res = pool.map_async(exe,otto).get()
    #for org_tuple,file in zip(orgs_tuples,files):
    #    run(f'awk \'{{{{ if ($1 != "#" && $8>40) print }}}}\' {work_dir}/clasp_out_forward/{org_tuple}/{file} > {work_dir}/clasp_out_forward/{org_tuple}/{file}_awk && awk \'{{{{ if ($1 != "#" && $8>40) print }}}}\' {work_dir}/clasp_out_reverse/{org_tuple}/{file} > {work_dir}/clasp_out_reverse/{org_tuple}/{file}_awk')
    #with mp.Pool(processes=12) as pool:
    #    res = pool.map_async(exe2,otto).get()
