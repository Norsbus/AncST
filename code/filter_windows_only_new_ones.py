#!/usr/bin/env python3

import os
import numpy as np
import pickle,json
import multiprocessing as mp
from subprocess import run
from Bio import SeqIO
from sys import argv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from bisect import bisect_left
import pathlib

def get_idents(window):
    start,end = window
    pos = bisect_left(seqlen,start)
    if start == seqlen[pos]:
        pos += 1
    ident1 = seqids[pos]
    pos = bisect_left(seqlen,end)
    if end == seqlen[pos]:
        pos += 1
    ident2 = seqids[pos]
    return(ident1,ident2)


def make_k_mers(org,idx,k,e,len_window,interval,percentile):
    
    idx = list(set(idx))

    len_lmers = int(len_window)
    interval = int(interval)
    
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
            # prints just to have look
            print(f'for window with start {i} and end {end} the start chromosome {ident} is not the same as end chromsome {ident2}')
            print('kmer_counts:')
            print(kmer_counts[i-13:end])
            print(kmer_counts[i])
            print(f'{seqlen}')
            print(f'len kmers: {len(kmer_counts)}')
            print(f'len seq: {len(s)}')
            
            continue

            
        # get the i with respect to chromosomes
        if pos == 0:
            x = 0
        else:
            x = 1

        i_alt = i - x * seqlen[pos-1]

        records.append(SeqRecord(s[i:end],id="kmer${}${}${}${}".format(i,ident,i_alt,end)))
        records_rev.append(SeqRecord(s[i:end].reverse_complement(),id="kmer${}${}${}${}".format(i,ident,i_alt,end)))
    
    SeqIO.write(records, f"{work_dir}//windows_low_kmers_fastas/{org}.fasta", "fasta")
    SeqIO.write(records_rev, f"{work_dir}/windows_low_kmers_fastas_rev/{org}.fasta", "fasta")

    return(0)


# function to find next index where l-mer starts without Ns in fasta
# Ns are represented as 0 in kmer_counts

def find_next_valid(kmer_counts_local,l):
    i = 0
    while i < len(kmer_counts_local) - l+1:
        if 0 in kmer_counts_local[i:i+l]:
            # take last found 0 and add one
            i = i + np.where(kmer_counts_local[i:i+l] == 0)[0][-1] + 1
        else:
            return(i)
    
    return('no valid')

# function to calculate sums of consecutive kmer_counts

def get_sums(args):
    a,z,l,interval = args
    global kmer_counts
    kmer_counts_local = kmer_counts[a:z]
    sums = []
    indices = []
    size_array_for_zero_check = 10000
    j = 0

    while j < len(kmer_counts_local)-l+1:
        
        # get next valid l-mer

        if 0 in kmer_counts_local[j:j+l]:
            res = find_next_valid(kmer_counts_local[j:j+size_array_for_zero_check],l) 
            while res == 'no valid' and j < len(kmer_counts_local)-l+1:
                j += size_array_for_zero_check
                res = find_next_valid(kmer_counts_local[j:j+size_array_for_zero_check],l) 
            # if end of array and no valid l-mer
            if res =='no valid':
                return(sums,indices)
            else:
                j+=res
        sums.append(np.sum(kmer_counts_local[j:j+l]))
        indices.append(j+a)
        j += interval
    return(sums,indices)

def get_low_count_windows(org,genmap_out_file,out_file_name,k,e,len_window,interval,percentile):
    l = int(len_window)

    percentile = int(percentile)

    interval  = int(interval)

    global kmer_counts
    kmer_counts = np.fromfile(genmap_out_file, dtype=np.uint16)
    
    max_concur_proc = int(argv[4])

    #sums_global = []
    #indices_global = []
    #pool = mp.Pool(max_concur_proc)
    #size_chunks = int(len(kmer_counts)/max_concur_proc)
    #args = []
    #indices_chunks = []
    #no_chunks = int(np.ceil(len(kmer_counts)/size_chunks))
    #for i in range(no_chunks):
    #    if i*size_chunks+size_chunks >= len(kmer_counts):
    #        args.append((i*size_chunks,len(kmer_counts),l,interval))
    #    else:
    #        args.append((i*size_chunks,i*size_chunks+size_chunks,l,interval))
    #        # add chunk for in between, otherwise those kmer counts not considered
    #        args.append((i*size_chunks+size_chunks-l+1,i*size_chunks+size_chunks+l-1,l,interval))
    #for sums,indices in pool.map(get_sums,args):
    #    sums_global += sums
    #    indices_global += indices


    ### CHANGED TO 1 THREAD. otherwise overlapping parts would have to be checked so that no matter which threads the same windows come out...
    # only became a "problem" when implementing a version where adding new windows/anchor candidates is supposed to be consistent with the recent runs,
    # so that different thread numbers dont influence the windows...would only be slightly different and it doesnt really matter (definitely not for single runs)
    # but to be consistent it does now...not so much of a speed gain until now with those smaller genomes anyway and will probably never be the bottleneck
    # otherwise just implement checks to make sure windows are always the same when pipeline is run with same parameters
    sums_global,indices_global = get_sums((0,len(kmer_counts),l,interval))
    perc = np.percentile(sums_global,percentile)
    idx = np.where(sums_global < perc)[0]
    final_idx = [indices_global[i] for i in idx]
    final_idx = sorted(final_idx)

    return(final_idx)

def combine_and_filter_new_ones(idx,bib=False):

    idx = sorted(idx)

    new_idx = []

    i = 0

    if percentile < 31:

        while i < len(idx)-1:
            idx_start = idx[i]
            idx_end = idx[i] + len_lmers
            pos = bisect_left(seqlen,idx_start)
            if idx_start == seqlen[pos]:
                pos += 1
            ident1 = seqids[pos]
            pos = bisect_left(seqlen,idx[i+1])
            if idx[i+1] == seqlen[pos]:
                pos += 1
            ident2 = seqids[pos]
            while i < len(idx)-1 and (idx_end >= idx[i+1] or idx[i+1] - idx_end <= interval) and ident1 == ident2:
                idx_end = idx[i+1] + len_lmers
                i += 1
                # have to check again...pre-calculating ident2 for checks in while loops which uses i+1 but when used here it could be out of range
                if i < len(idx)-1:
                    pos = bisect_left(seqlen,idx[i+1])
                    if idx[i+1] == seqlen[pos]:
                        pos += 1
                    ident2 = seqids[pos]
            new_idx.append((idx_start,idx_end))
            i += 1

        if len(new_idx) > 0 and len(idx) > 0:
            if new_idx[-1][1] != idx[-1] + len_lmers:
                new_idx.append((idx[-1],idx[-1]+len_lmers))

    else:
        for i in idx:
            new_idx.append((i,i+len_lmers))

    if not bib:
        # for checks/debugging
        for start,end in new_idx:

            print('===============')
            print('orig',start,end)
        return(new_idx)

    new_ones = []

    old_idx = sorted(list(bib.keys()))
    old_idx_end = [i+bib[i]['end']-bib[i]['start'] for i in old_idx]

    for start,end in new_idx:
        bigger_idx = bisect_left(old_idx,start)
        if bigger_idx == len(old_idx):
            if old_idx_end[-1] < start:
                new_ones.append((start,end))
            continue
        if bigger_idx == 0:
            if end < old_idx[0]:
                new_ones.append((start,end))
            continue
        if start == old_idx[bigger_idx]:
            continue
        if old_idx_end[bigger_idx-1] >= start:
            continue
        if end >= old_idx[bigger_idx] or end >= old_idx_end[bigger_idx]:
            continue
        if old_idx_end[bigger_idx - 1] < start:
            new_ones.append((start,end))
        new_ones.append((start,end))

    return(new_ones)


if __name__ == "__main__":

    org = argv[1]
    infile = argv[2]
    outfile = argv[3]

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]

    max_concur_proc = int(argv[4])

    k,e,len_lmers,interval,percentile = int(argv[5]),int(argv[6]),int(argv[7]),int(argv[8]),int(argv[9])

    ### for logging ###

    print(f'filtering windows with parameters: k={k}, e={e}, len_lmers={len_lmers}, interval={interval}, percentile={percentile}')
    
    idx = get_low_count_windows(org,infile,outfile,k,e,len_lmers,interval,percentile)
    
    with open(root + '/utils/metadata_genomes/{}'.format(org),'rb') as f:
        seqids,seqlen,s = pickle.load(f)

    bib = False

    try:
        with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
            bib = pickle.load(f)
    except:
        pass

    if not bib:
        idx = combine_and_filter_new_ones(idx)
    else:
        idx = combine_and_filter_new_ones(idx,bib)
        if len(idx) == 0:
            print('no new candidates')
            run(f'touch {work_dir}/touch/no_new_candidates_{org}',shell=True)
        to_del = [] 
        with open(f'{work_dir}/to_del/{org}/to_del','wb') as f:
            pickle.dump(to_del, f)
    # outfile is indices directory/{org}
    with open(outfile,'wb') as f:
        pickle.dump(idx, f)

    make_k_mers(org,idx,k,e,len_lmers,interval,percentile)
