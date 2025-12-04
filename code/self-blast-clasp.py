#! /usr/bin/env python3

import pathlib
from sys import argv
from subprocess import run
import os
from multiprocessing import Pool,current_process
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from bisect import bisect_left
import pickle
from subprocesses import blast,clasp
import re
from collections import Counter

_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
        return [int(text) if text.isdigit() else text.lower()
                            for text in re.split(_nsre, s)]

def make_new_regions(regions):

    new_regions = []
    starts = [(x[0][0],x[1],'start') for x in regions]
    ends = [(x[0][1],x[1],'end') for x in regions]
    se = sorted(starts+ends)
    score_counts = Counter()

    for j,i in enumerate(se[:-1]):
        if i[-1] == 'start':
            score_counts[i[1]] += 1
        else:
            score_counts[i[1]] -= 1
            if score_counts[i[1]] == 0:
                del score_counts[i[1]]
        if se[j+1][0] - i[0] == 0:
            if i[1] > se[j+1][1]:
                se[j+1][1] = i[1]
            continue

        if len(new_regions) > 0:
            if new_regions[-1][1] == max(score_counts):
                new_regions[-1][0][1] = se[j+1][0]
            else:
                new_regions.append([[i[0],se[j+1][0]],max(score_counts)])
        else:
            new_regions.append([[i[0],se[j+1][0]],max(score_counts)])

    return(new_regions)


def make_new_regions_old(regions):
    
    new_regions = []
    starts = [(x[0][0],x[1],'start') for x in regions]
    ends = [(x[0][1],x[1],'end') for x in regions]
    se = sorted(starts+ends)
    scores = []

    for j,i in enumerate(se[:-1]):
        if i[-1] == 'start':
            scores.append(i[1])
        else:
            scores.remove(i[1])
        if se[j+1][0] - i[0] == 0:
            if i[1] > se[j+1][1]:
                se[j+1][1] = i[1]
            continue

        if len(new_regions) > 0:
            if new_regions[-1][1] == max(scores):
                new_regions[-1][0][1] = se[j+1][0]
            else:
                new_regions.append([[i[0],se[j+1][0]],max(scores)])
        else:
            new_regions.append([[i[0],se[j+1][0]],max(scores)])
    
    return(new_regions)


def reconcile_regions(bib):
    if len(bib[0]) == 0:
        return(bib)
    elif len(bib[0]) == 1:
        return(bib)
    new_bib = [[],bib[1]]
    bib[0] = sorted(bib[0],key=lambda x: (x[0][0],x[0][1]))
    
    i = 0
    while i < len(bib[0]):
        
        regs = [bib[0][i]]
        cur_end = bib[0][i][0][1]
        i += 1

        while i < len(bib[0]) and cur_end >= bib[0][i][0][0]:
            if bib[0][i][0][1] > cur_end:
                cur_end = bib[0][i][0][1]
            regs.append(bib[0][i])
            i += 1
        if len(regs) > 1:
            regs = make_new_regions(regs)
        new_bib[0].extend(regs)
         
    return(new_bib)

overlap_tol = 0.5

def parse_clasp_out(clasp_out):
    with open(f'{root}/utils/small_meta/{org}','rb') as f: 
        seqids,seqlen = pickle.load(f)
    results = {}
    scores = []

    # parsing clasp.x output. need to get the second best chained hit

    for orientation in ['forward','reverse']:
        ff = f'{work_dir}/clasp_out_{orientation}/{org}/{clasp_out}'
        run('''awk '{{ if ($1 != "#" && $8>40) print }}' {} > {}_awk'''.format(ff,ff),shell=True)
        with open(ff+'_awk','r') as file:
            for line in file:
                line = line.split()
                name = line[1].split('$')
                i = int(name[1])
                i_end = int(name[-1])
                laenge = i_end - i
                score = float(line[7])
                start = int(line[3]) 
                end = int(line[4]) 
                laenge2 = end - start
                if orientation == 'reverse':
                    # the corrections are necessary because of blast starting pos...ive gone through this and its annoying so just trust (never)...
                    start = laenge - end + 1 
                    end = start + laenge2
                
                chromo_db = line[2]
                start_db = int(line[5])
                end_db = int(line[6])

                pos = seqids.index(chromo_db)
                l_chr = seqlen[pos-1]
                if pos == 0:
                    x = 0
                else:
                    x = 1
                abs_pos_start = start_db + x * l_chr - 1 # blast positions start/end with 1 and ae inclusive (== start should be -1, end not)
                abs_pos_end = end_db + x * l_chr
                
                if i not in results:
                    results[i] = [[],i_end]

                if orientation == 'forward':

                    if abs_pos_start <= i and abs_pos_end >= i_end:
                        continue
                    elif abs_pos_start >= i and abs_pos_start <= i_end:
                        if abs_pos_end <= i_end:
                            continue
                        elif (i_end - abs_pos_start)/laenge2 > overlap_tol:
                            continue
                    elif abs_pos_end >= i and abs_pos_end <= i_end:
                        if (abs_pos_end - i)/laenge2 > overlap_tol:
                            continue

                results[i][0].append([[start,end],score])
        
    idx = []
    
    for i,bib in results.items():
        bib = reconcile_regions(bib)
        scores_i = [x[1] for x in bib[0]]
        if len(scores_i) > 0 and max(scores_i) > bib[1]-i:
            print('highly similar to sth else:',i,bib)
            #continue
        idx.append([(i,bib[1]),bib[0]]) 

        # criterion for not taking it is that compared to its length it has at least as many points
    
    # i have seen a case where blast doesnt find its own sequence..need tto check and put them with high threshold cos they are conspicuous
    # they are probably repeats...to make it consistent and have all candidates in final bib, i add them anyway
   
    orig_idx = []
    with open(f'{work_dir}/split_out_{orientation}/{org}/{clasp_out}','r') as f:
        for line in f:
            if line[0] == '>':
                line = line.split('$')
                orig_idx.append([int(line[1]),int(line[-1].split()[0])])

    for i,end in orig_idx:
        if i not in results:
            # THESE CASES ARE DUE TO MASKING (can be checked by running blastn with "-dust np" option
            print(f'sequence start: {i} - end: {end} (abs) not found in own genome...adding it with arbitrary (but relatively high) threshold of its length * 2')
            idx.append([(i,end),[[[0,end-i],(end-i)*2]]])

    idx = sorted(idx)
    
    return(idx)



def exec_blast_clasp(blast_word_sizes,split_out_file):
        
    for ws in blast_word_sizes:

        for orientation in ['forward','reverse']:
            
            blast(f'{root}/utils/blastdbs/{org}',f'{work_dir}/split_out_{orientation}/{org}/{split_out_file}',f'{work_dir}/blast_out_{orientation}/{org}/{split_out_file}',ws)
            
            clasp(f'{work_dir}/blast_out_{orientation}/{org}/{split_out_file}',f'{work_dir}/clasp_out_{orientation}/{org}/{split_out_file}')

        idx = parse_clasp_out(split_out_file)
        with open(f'{work_dir}/indices/{org}/{split_out_file}_indices','wb') as f:
            pickle.dump(idx,f)

    return(0)


if __name__ == "__main__":
    
    infile = argv[1]
    org = '.'.join(infile.split('.')[:-2])
    
    work_dir = argv[-1]
    root = str(pathlib.Path(__file__).parents[1])

    blast_word_sizes = [int(argv[2])]

    exec_blast_clasp(blast_word_sizes,infile)
