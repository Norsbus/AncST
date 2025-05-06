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
from pprint import pprint
from copy import deepcopy
from time import time

_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
        return [int(text) if text.isdigit() else text.lower()
                            for text in re.split(_nsre, s)]

overlap_tol = 0.5

def filter_blast_file(f,orientation):
    skip = {}
    # kmer$452068600$NC_079137.1$11383464$452069600	NC_079137.1	100.00	1000	0	0	1	1000	11383465	11384464	0.0	1847
    newlines = {}
    ignore_candidates = {}
    with open(f) as ff:
        print(f'{f} opened')
        for line in ff:
            #print(line)
            candidate,chromo_db,x,y,z,a,cstart,cend,start_db,end_db,e,s = line.strip().split()
            candidate = candidate.strip()
            chromo_db = chromo_db.strip()
            start_db = int(start_db)
            end_db = int(end_db)
            if candidate not in newlines:
                newlines[candidate] = []
            newlines[candidate].append(line)
            e = float(e)
            s = float(s)
            name = candidate.split('$')
            cchromo = name[2]
            i = int(name[1])
            i_end = int(name[-1])
            laenge = i_end - i
            cstart = int(cstart)
            cend = int(cend)
            laenge2 = int(cend - cstart)
            if 'reverse' in f:
                # the corrections are necessary because of blast starting pos...ive gone through this and its annoying so just trust (never)...
                cstart = laenge - cend + 1 
                cend = cstart + laenge2

            pos = seqids.index(chromo_db)
            l_chr = seqlen[pos-1]
            if pos == 0:
                x = 0
            else:
                x = 1
            abs_pos_start = start_db + x * l_chr - 1 # blast positions start/end with 1 and ae inclusive (== start should be -1, end not)
            abs_pos_end = end_db + x * l_chr
            
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

            if end_db - start_db > laenge/4 and s > (laenge/4)*1.5:
                ignore_candidates[candidate] = 1
                print(f'debug filter blast: ignoring {candidate} because of line:')
                print(line)
    print('hihi')
    pprint(len(newlines))
    print(ignore_candidates.keys())
    return ignore_candidates,newlines

def make_new_regions(regions):
    
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


def reconcile_regions_old(bib):
    print('into rec regions')
    pprint(bib)
    if len(bib[0]) == 0:
        return(bib)
    elif len(bib[0]) == 1:
        return(bib)
    new_bib = [[],bib[1]]
    bib[0] = sorted(bib[0],key=lambda x: (x[0][0],x[0][1]))
    print('processing rec region')
    pprint(bib[0]) 
    i = 0
    while i < len(bib[0]):
        print(i)
        
        regs = [bib[0][i]]
        cur_end = bib[0][i][0][1]
        i += 1

        while i < len(bib[0]) and cur_end >= bib[0][i][0][0]:
            if bib[0][i][0][1] > cur_end:
                cur_end = bib[0][i][0][1]
            regs.append(bib[0][i])
            i += 1
        print(regs)
        if len(regs) > 1:
            print('making new regions')
            regs = make_new_regions(regs)
        new_bib[0].extend(regs)
         
    return(new_bib)

def reconcile_regions_incomplete(abs_start,bib):
    print('into rec regions')
    pprint(bib)
    if len(bib[0]) == 0:
        return(bib)
    elif len(bib[0]) == 1:
        return(bib)
    new_bib = [[],bib[1]]
    bib[0] = sorted(bib[0],key=lambda x: (x[1]),reverse=True)
    print('sorted bib acc to scores')
    pprint(bib[0])
    to_process = deepcopy(bib[0])
    while len(to_process) > 0:
        for start,end,score in bib[0]:
            cur_start = start
            cur_end = end
            add = [start,end]
            contained = 0
            for start2,end2,score2 in new_bib[0]:
                if (start,end,score) == (start2,end2,score2):
                    continue
                if cur_start > start2 and cur_start < end2 and cur_end > end2:
                    cur_start = end2 + 1
                    break
                elif cur_end > start2 and cur_and < end2 and cur_start < start2:
                    cur_end = start2 - 1
                if start >= start2 and start <= end2 and end >= start2 and end <= end2:
                    contained = 1
                    break
            if contained == 1:
                continue
            else:
                if cur_end - cur_start > 0:
                    new_bib[0].append(([cur_start,cur_end],score))

def reconcile_regions(abs_start,bib):
    print('into rec regions')
    pprint(bib)
    if len(bib[0]) == 0:
        return(bib)
    elif len(bib[0]) == 1:
        return(bib)
    new_bib = [[],bib[1]]
    tuple_bib = [(t[0],t[1],s) for t,s in bib[0]]
    without_duplicates = list(set(tuple_bib))
    bib[0] = [([st,e],s) for st,e,s in without_duplicates]
    bib[0] = sorted(bib[0],key=lambda x: (x[1]),reverse=True)
    print('sorted bib acc to scores')
    pprint(bib[0])
    iss = []
    for i in range(1,bib[1]-abs_start+1):
        limit = 40
        for t,score in bib[0]:
            start,end= t
            if i >= start and i <= end:
                limit = score
                break
        iss.append(limit)

    #print('iss')
    #print(iss)
    cur_start = 1
    cur_limit = iss[0]
    for i,limit in enumerate(iss):
        if limit != cur_limit:
            new_bib[0].append([[cur_start,i],cur_limit])
            cur_limit = limit
            cur_start = i+1

    if len(new_bib[0]) == 0 or new_bib[0][-1][0][1] != len(iss)+1:
        new_bib[0].append([[cur_start,len(iss)+1],limit])
    print('final bib')
    pprint(new_bib)
    
    return(new_bib)


def parse_clasp_out(clasp_out,ignore_candidates,suf):
    results = {}
    scores = []

    iss_ignore_candidates = {}
    for i in ignore_candidates['forward']:
        iss_ignore_candidates[int(i.split('$')[1])] = 1
    for i in ignore_candidates['reverse']:
        iss_ignore_candidates[int(i.split('$')[1])] = 1
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
        bib = reconcile_regions(i,bib)
        scores_i = [x[1] for x in bib[0]]
        if len(scores_i) > 0 and max(scores_i) > bib[1]-i:
            print('highly similar to sth else:',i,bib)
            continue
        idx.append([(i,bib[1]),bib[0]]) 

        # criterion for not taking it is that compared to its length it has at least as many points
    
    # i have seen a case where blast doesnt find its own sequence..need tto check and put them with high threshold cos they are conspicuous
    # they are probably repeats...to make it consistent and have all candidates in final bib, i add them anyway
   
    orig_idx = []
    with open(f'{work_dir}/split_out_{orientation}{suf}/{org}/{clasp_out}','r') as f:
        for line in f:
            if line[0] == '>':
                line = line.split('$')
                orig_idx.append([int(line[1]),int(line[-1].split()[0])])
    
    print(iss_ignore_candidates.keys())
    for i,end in orig_idx:
        if i not in results and i not in iss_ignore_candidates:
            # THESE CASES ARE DUE TO MASKING (can be checked by running blastn with "-dust np" option
            print(f'sequence start: {i} - end: {end} (abs) not found in own genome...adding it with arbitrary (but relatively high) threshold of its length * 2')
            idx.append([(i,end),[[[0,end-i],(end-i)*2]]])

    idx = sorted(idx)
    
    return(idx)



def exec_blast_clasp(blast_word_sizes,split_out_file,suf=''):
        
    for ws in blast_word_sizes:

        ignore_candidates = {}
        newlines_filtered = {}
        for orientation in ['forward','reverse']:

            starttime = time()
            ret = blast(f'{root}/utils/blastdbs/{org}',f'{work_dir}/split_out_{orientation}{suf}/{org}/{split_out_file}',f'{work_dir}/blast_out_{orientation}/{org}/{split_out_file}',ws,suf)
            print(f'blast took {time() - starttime}')
            if ret == 1 or ret == -1:
                print('hihi 11')
                return(11)

            starttime = time()
            ignore_candidates[orientation],newlines_filtered[orientation] = filter_blast_file(f'{work_dir}/blast_out_{orientation}/{org}/{split_out_file}',orientation)
            print(f'filtering took {time() - starttime}')

        for orientation in ['forward','reverse']:
            with open(f'{work_dir}/blast_out_{orientation}/{org}/{split_out_file}_filtered','w') as ff:
                for candidate,lines in newlines_filtered[orientation].items():
                    if candidate in ignore_candidates['forward'] or candidate in ignore_candidates['reverse']:
                        continue
                    else:
                        ff.writelines(lines)

        for orientation in ['forward','reverse']:
            starttime = time()
            ret = clasp(f'{work_dir}/blast_out_{orientation}/{org}/{split_out_file}_filtered',f'{work_dir}/clasp_out_{orientation}/{org}/{split_out_file}',suf,2,0.1)
            print(f'clasp took {time() - starttime}')
            if ret == 1 or ret == -1:
                print('hihi 12')
                return(12)

        starttime = time()
        idx = parse_clasp_out(split_out_file,ignore_candidates,suf)
        print(f'parsing took {time() - starttime}')
        with open(f'{work_dir}/indices/{org}/{split_out_file}_indices','wb') as f:
            pickle.dump(idx,f)

    return(0)


if __name__ == "__main__":
    g_starttime = time()
    
    infile = argv[1]
    org = '.'.join(infile.split('.')[:-2])
    
    work_dir = argv[-1]
    root = str(pathlib.Path(__file__).parents[1])
    with open(f'{root}/utils/small_meta/{org}','rb') as f: 
        seqids,seqlen = pickle.load(f)

    blast_word_sizes = [int(argv[2])]

    ret = exec_blast_clasp(blast_word_sizes,infile)
    if ret == 11 or ret == 12:
        new_files = []
        counter = 1
        record_ids = {}
        for record in SeqIO.parse(f'{work_dir}/split_out_forward/{org}/{infile}', "fasta"):
            SeqIO.write(record, f"{work_dir}/split_out_forward_add/{org}/add-{infile}-{counter}.fasta", "fasta")
            new_files.append(f'add-{infile}-{counter}.fasta')
            record_ids[record.id] = counter
            counter += 1
        for record in SeqIO.parse(f'{work_dir}/split_out_reverse/{org}/{infile}', "fasta"):
            name = record_ids[record.id]
            SeqIO.write(record, f"{work_dir}/split_out_reverse_add/{org}/add-{infile}-{name}.fasta", "fasta")
        for f in new_files:
            ret = exec_blast_clasp(blast_word_sizes,f,suf='_add')
            if ret == 11 or ret == 12:
                print(f'FINAL FAIL because of blast/clasp time out ({ret}): {work_dir}/split_out_reverse/{org}/{f}')
    print(f'everything took {time() - g_starttime}')
