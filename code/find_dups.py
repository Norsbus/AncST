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
import pathlib

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

def find_overlapping(regions1, regions2):
    """
    Find overlapping regions between two sets.
    Each region is a (start, end, pitch) tuple.
    Returns set of (start, end, pitch) tuples that overlap.
    """
    # Extract (start, end) pairs for overlap checking, preserve full tuples
    pos1 = [(start, end, pitch) for start, end, pitch in regions1]
    pos2 = [(start, end, pitch) for start, end, pitch in regions2]

    pos1 = sorted(pos1)
    mstarts = [xxx[0] for xxx in pos1]
    mends = [xxx[1] for xxx in pos1]
    pos2 = sorted(pos2)
    gstarts = [xxx[0] for xxx in pos2]
    gends = [xxx[1] for xxx in pos2]

    overlapping = set()
    for ms, me, pitch1 in pos1:
        bigger_end_than_start_idx = bisect_left(gends, ms)
        smaller_start_than_end_idx = max(bisect_left(gstarts, me) - 1, 0)
        if smaller_start_than_end_idx - bigger_end_than_start_idx < 1:
            continue
        for gs, ge, pitch2 in pos2[max(0, bigger_end_than_start_idx-10):min(smaller_start_than_end_idx+10, len(pos2))]:
            t1 = (ms, me)
            t2 = (gs, ge)
            if overlaps(t1, t2):
                overlapping.add((ms, me, pitch1))
                overlapping.add((gs, ge, pitch2))
    return overlapping

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
                return(set([(d, d + window_size, pitch) for d in dups]))
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

    return(set([(d, d + window_size, pitch) for d in dups]))

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
    root = str(pathlib.Path(__file__).parents[1])

    old_dups = [] 
    #if isfile(f'{work_dir}/../anchors/candidates' + f'/{org}'):
    #    with open(f'{work_dir}/../anchors/candidates'+ f'/{org}','rb') as f:
    #        old = pickle.load(f)
    #    for i,bib in old.items():
    #        if 'dup' in bib:
    #            old_dups.append((i,i+(bib['end']-bib['start'])))

    dups = {}
    best_dups = []  # flag=0 sets
    check_dups = []  # flag=1 sets

    dups_paras = []
    with open(f"{work_dir}/dups_params.txt","r") as f:
        for line in f:
            if org in line:
                dups_paras.append(line.strip().split()[1:])

    # If organism not in dups_params.txt, create empty dups file and exit
    if not dups_paras:
        print(f'No dups parameters found for {org}, creating empty dups file')
        bib = {}
        with open(f'{work_dir}/dups/{org}','wb') as f:
            pickle.dump(bib,f)
        exit(0)

    #k1: k of k-mer counts for initial identification of candidate windows
    #e1: errors allowed in k1-mer counting
    #k2: k of k-mer counts for exclusion of windows according to x2 and y2
    #e2: errors allowed in k2-mer counting
    #w: size of initially overlapping windows
    #p: pitch of w-windows
    #x1: how many counts of initial windows are allowed (will be a range [2,x1])
    #y1: percentage of counts in candidate windows adhering to count values defined with x1
    #x2: how many counts of k2,e2-mers are regarded as indicative of bad window (and all bigger than x2)
    #y2: percentage of counts in candidate windows from which candidate is discarded if x2 count bigger than or equal to y2
    #flag: 0=include all, 1=only if overlaps with any other set
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
        flag = int(paras[10]) if len(paras) > 10 else 0  # Default to 0 for backward compatibility

        kmer_counts = np.fromfile(f'{root}/utils/genmap_out/{org}/{k1}_{e1}.freq16', dtype=np.uint16)
        kmer_counts2 = np.fromfile(f'{root}/utils/genmap_out/{org}/{k2}_{e2}.freq16', dtype=np.uint16)

        # Get dups for this parameter set - returns set of (start, end, pitch) tuples
        regions = get_dups(kmer_counts, kmer_counts2, window_size, pitch, old_dups,
                          thres11, thres12, thres21, thres22)

        # Categorize by flag
        if flag == 0:
            best_dups.append(regions)
        elif flag == 1:
            check_dups.append(regions)
        else:
            print(f'ERROR: Invalid flag {flag} for {org} in dups_params.txt')
            exit(1)

    # Combine: all flag=0 + flag=1 that overlap with ANY set
    final_regions = set()

    # Add all flag=0 regions unconditionally
    for regions in best_dups:
        final_regions.update(regions)

    # For flag=1: add regions that overlap with ANY other set (including other flag=1)
    all_sets = best_dups + check_dups
    for i, regions1 in enumerate(check_dups):
        for j, regions2 in enumerate(all_sets):
            # Skip comparing set with itself
            if i == j and regions2 in check_dups:
                continue
            overlapping = find_overlapping(regions1, regions2)
            final_regions.update(overlapping)

    regions = final_regions

    with open(f'{root}/utils/metadata_genomes/{org}','rb') as f:
        seqids,seqlen,s = pickle.load(f)

    # Sort regions by start position, keeping (start, end, pitch) tuples
    regions = sorted(list(regions), key=lambda x: x[0])
    new_idx = []
    i = 0
    while i < len(regions)-1:
        idx_start, idx_end, pitch_i = regions[i]
        pos = bisect_left(seqlen, idx_start)
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

        # Get next region
        next_start, next_end, next_pitch = regions[i+1]
        pos = bisect_left(seqlen, next_end)
        if next_end == seqlen[pos]:
            pos += 1
        ident2 = seqids[pos]

        # Merge regions that overlap or are close (using longest pitch)
        merge_pitch = max(pitch_i, next_pitch)
        while i < len(regions)-1 and (idx_end >= next_start or next_start - idx_end <= merge_pitch) and ident1 == ident2:
            idx_end = max(idx_end, next_end)  # Extend to furthest end
            i += 1
            # Update for next iteration
            if i < len(regions)-1:
                next_start, next_end, next_pitch = regions[i+1]
                merge_pitch = max(pitch_i, next_pitch)
                pos = bisect_left(seqlen, next_end)
                if next_end == seqlen[pos]:
                    pos += 1
                ident2 = seqids[pos]
        new_idx.append((idx_start, idx_end))
        i += 1

    # Handle last region
    if len(new_idx) > 0 and len(regions) > 0:
        last_start, last_end, last_pitch = regions[-1]
        if last_start > new_idx[-1][1]:
            new_idx.append((last_start, last_end))

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
