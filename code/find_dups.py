#!/usr/bin/env python3

import numpy as np
from pprint import pprint
import pickle
from bisect import bisect_left
from sys import argv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import Popen,run,PIPE,TimeoutExpired
import resource
from os.path import isfile

def overlaps(x,y):
    if x[0] >= y[0] and x[0] <= y[1]:
        return True
    if x[1] >= y[0] and x[1] <= y[1]:
        return True
    if y[0] >= x[0] and y[0] <= x[1]:
        return True
    if y[1] >= x[0] and y[1] <= x[1]:
        return True
    return False

def find_overlapping(idx1,idx2,window_size):

    pos1 = []
    pos2 = []

    for i in idx1:
        pos1.append((i,i+window_size))
    for i in idx2:
        pos2.append((i,i+window_size))

    pos1 = sorted(pos1)
    mstarts = [xxx[0] for xxx in pos1]
    mends = [xxx[1] for xxx in pos1]
    pos2 = sorted(pos2)
    gstarts = [xxx[0] for xxx in pos2]
    gends = [xxx[1] for xxx in pos2]
    idx = set()
    for ms,me in pos1:
        bigger_end_than_start_idx = bisect_left(gends,ms)
        smaller_start_than_end_idx = max(bisect_left(gstarts,me) - 1,0)
        if smaller_start_than_end_idx - bigger_end_than_start_idx < 1:
            continue
        for gs,ge in pos2[max(0,bigger_end_than_start_idx-10):min(smaller_start_than_end_idx+10,len(pos2))]:
            t1 = (ms,me)
            t2 = (gs,ge)
            if overlaps(t1,t2):
                idx.add(ms)
                idx.add(gs)
    return(idx)

def which_chromo(seqlen,seqids,i):
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

def find_next_valid(kmer_counts,window_size):
    i = 0
    while i < len(kmer_counts) - window_size+1:
        if 0 in kmer_counts[i:i+window_size]:
            # take last found 0 and add one
            i = i + np.where(kmer_counts[i:i+window_size] == 0)[0][-1] + 1
        else:
            return(i)
    
    return('no valid')

def get_dups(kmer_counts,kmer_counts2,window_size,pitch,old_dups,thres11,thres12,thres21,thres22):

    dups = []

    size_array_for_zero_check = 10000
    j = 0

    while j < len(kmer_counts)-window_size+1:

        for t in old_dups:
            if overlaps(t,(j,j+window_size)):
                j+=pitch
                continue
        
        # get next valid window_size-mer
        
        if 0 in kmer_counts[j:j+window_size]:
            res = find_next_valid(kmer_counts[j:j+size_array_for_zero_check],window_size) 
            while res == 'no valid' and j < len(kmer_counts)-window_size+1:
                j += size_array_for_zero_check
                res = find_next_valid(kmer_counts[j:j+size_array_for_zero_check],window_size) 
            # if end of array and no valid window_size-mer
            if res =='no valid':
                return(dups)
            else:
                j+=res
        counts = kmer_counts[j:j+window_size]
        tot = 0
        for i in range(2,thres11):
            tot += np.count_nonzero(counts == i)
        tot2 = 0
        counts2 = kmer_counts2[j:j+window_size]
        tot2 = (counts2 > thres21).sum()
        if tot >= len(counts)*thres12 and tot2 < len(counts2)*thres22:
            dups.append(j)

        j += pitch

    return(set(dups))

def make_k_mers(org,idx):
    
    idx = list(set(idx))

    records = []
    records_rev = []

    for i,end in idx:
        
        pos = bisect_left(seqlen,i)
        # since seqlen has indices of chromosmoe ends + 1, if something starts exactly at first nucleotide it will be wrong
        if i == seqlen[pos]:
            pos += 1
        ident = seqids[pos]
        # check again for chromosome breaks although this should largely be achieved automatically since genmap has count 0 for ends which are not valid
        # could still be (and happens) when old and new ones are combined are lie between chromosome
        # should be dealt with in combine_and_filter tho (just a check because new ones should never cover different chromosomes)
        ident2 = seqids[bisect_left(seqlen,end)]
        if end == seqlen[pos]:
            pos += 1
        if ident != ident2:
            continue
            
        # get the i with respect to chromosomes
        if pos == 0:
            x = 0
        else:
            x = 1

        i_alt = i - x * seqlen[pos-1]

        records.append(SeqRecord(s[i:end],id="kmer${}${}${}${}".format(i,ident,i_alt,end)))
        records_rev.append(SeqRecord(s[i:end].reverse_complement(),id="kmer${}${}${}${}".format(i,ident,i_alt,end)))
    
    SeqIO.write(records, f"dups_windows/{org}.fasta", "fasta")
    SeqIO.write(records_rev, f"dups_windows_rev/{org}.fasta", "fasta")

    return 0 

if __name__ == "__main__":

    org = argv[1]
    work_dir = argv[-1]
    
    old_dups = [] 
    #if isfile(f'{work_dir}/../anchors/candidates' + f'/{org}'):
    #    with open(f'{work_dir}/../anchors/candidates'+ f'/{org}','rb') as f:
    #        old = pickle.load(f)
    #    for i,bib in old.items():
    #        if 'dup' in bib:
    #            old_dups.append((i,i+(bib['end']-bib['start'])))

    dups = {}

    dups_paras = []
    with open(f"{work_dir}/dups_params.txt","r") as f:
        line = f.read()
        if org in line:
            dups_paras.append(line.strip().split()[1:])

    # If organism not in dups_params.txt, create empty dups file and exit
    if not dups_paras:
        print(f'No dups parameters found for {org}, creating empty dups file')
        bib = {}
        with open(f'{work_dir}/dups/{org}','wb') as f:
            pickle.dump(bib,f)
        exit(0)

    #k1: k of k-mer counts for initial identification of candidate windows<br>
    #e1: errors allowed in k1-mer counting<br>
    #k2: k of k-mer counts for exclusion of windows according to x2 and y2<br>
    #e2: errors allowed in k2-mer counting<br>
    #w: size of initially overlapping windows<br>
    #p: pitch of w-windows<br>
    #x1: how many counts of initial windows are allowed (will be a range [2,x1])<br>
    #y1: percentage of counts in candidate windows adhering to count values defined with x1<br>
    #x2: how many counts of k2,e2-mers are regarded as indicative of bad window (and all bigger than x2)<br>
    #y2: percentage of counts in candidate windows from which candidate is discarded if x2 count bigger than or equal to y2<br>
    for paras in dups_paras:
        k1 = int(paras[0])
        e1 = int(paras[1])
        k2 = int(paras[2])
        e2 = int(paras[3])
        window_size = int(paras[4])
        pitch = int(paras[5])
        thres11 = int(paras[6])
        thres12 = float(paras[7])/100
        thres21 = int(paras[8])
        thres22 = float(paras[9])/100

        kmer_counts = np.fromfile(f'../utils/genmap_out/{org}/{k1}_{e1}.freq16', dtype=np.uint16)
        kmer_counts2 = np.fromfile(f'../utils/genmap_out/{org}/{k2}_{e2}.freq16', dtype=np.uint16)

        # if 'N's are plenty use this instead of get_dups function to take windows with 'N's (not done in get_dups function). old code - not adapted fo multiple parameter sets
        #dups = []
        #for i in range(250,len(kmer_counts) - 250,100):
        #    counts = kmer_counts[i-250:i+250]
        #    tot = 0
        #    for j in range(2,6):
        #        tot += np.count_nonzero(counts == j)
        #    if tot > len(counts)*thres:
        #        dups.append(i-250)
        
        dups[tuple(paras)] = get_dups(kmer_counts,kmer_counts2,window_size,pitch,old_dups,thres11,thres12,thres21,thres22)

    idx = dups[tuple(dups_paras[0])]
    for paras in dups_paras[1:]:
        
        idx = find_overlapping(idx,dups[f'{tuple(paras)}'],window_size)

    with open(f'{work_dir}/../utils/metadata_genomes/{org}','rb') as f:
        seqids,seqlen,s = pickle.load(f)

    idx = sorted(list(idx))
    new_idx = []
    i = 0
    while i < len(idx)-1:
        idx_start = idx[i]
        idx_end = idx[i] + window_size
        pos = bisect_left(seqlen,idx_start)
        # if short chromosomes should be excluded
        #if pos == 0:
        #    if seqlen[0] < 1000000:
        #        i += 1
        #        continue
        #else:
        #    if seqlen[pos] - seqlen[pos-1] < 1000000:
        #        i += 1
        #        continue
        if idx_start == seqlen[pos]:
            pos += 1
        ident1 = seqids[pos]
        pos = bisect_left(seqlen,idx[i+1]+window_size)
        if idx[i+1]+window_size == seqlen[pos]:
            pos += 1
        ident2 = seqids[pos]
        while i < len(idx)-1 and (idx_end >= idx[i+1] or idx[i+1] - idx_end <= pitch) and ident1 == ident2:
            idx_end = idx[i+1] + window_size
            i += 1
            # have to check again...pre-calculating ident2 for checks in while loops which uses i+1 but when used here it could be out of range
            if i < len(idx)-1:
                pos = bisect_left(seqlen,idx[i+1]+window_size)
                if idx[i+1]+window_size == seqlen[pos]:
                    pos += 1
                ident2 = seqids[pos]
        new_idx.append((idx_start,idx_end))
        i += 1

    if len(new_idx) > 0 and len(idx) > 0:
        if new_idx[-1][1] != idx[-1] + window_size:
            new_idx.append((idx[-1],idx[-1]+window_size))

    bib = {}
    for i,end in new_idx:
        chromo = which_chromo(seqlen,seqids,i)
        idx_l = seqids.index(chromo)
        if idx_l > 0:
            l = seqlen[idx_l-1]
        else:
            l = 0
        chromo_start = i - l
        chromo_end = end - l
        bib[i] = {\
                    'dup':1,\
                    'chromosome':chromo,\
                    'start':chromo_start,\
                    'end':chromo_end,\
                    'matches':{}\
                    }

    with open(f'{work_dir}/dups/{org}','wb') as f:
        pickle.dump(bib,f)
