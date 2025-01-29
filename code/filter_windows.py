#!/usr/bin/env python3

import os
import pathlib
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

def get_low_count_windows_multiple_k_e(org,genmap_out_file,out_file_name,k_e,len_window,interval,percentile):
    l = int(len_window)

    percentile = int(percentile)

    interval  = int(interval)

    global kmer_counts
    kmer_counts = []

    for k,e in k_e:
        
        genmap_out_file = f'{root}/utils/genmap_out/{org}/{k}_{e}.freq16'

        kmer_counts_x = np.fromfile(genmap_out_file, dtype=np.uint16)
        print(kmer_counts_x[:100])
        kmer_counts_x = kmer_counts_x/np.sum(kmer_counts_x)
        print(kmer_counts_x[:100])
        if len(kmer_counts) == 0:
            kmer_counts = kmer_counts_x
        else:
            kmer_counts += kmer_counts_x
        print(kmer_counts[:100])
    
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

    cases_overlaps_single = 0
    cases_overlaps_double = 0
    cases_distances_single = 0
    cases_distances_double = 0
    cases_in_there = 0
    cases_in_there_but_longer = 0
    contains_old = 0
    contained_old = 0
    contains_new = 0
    two_sided_dual = 0
    lens_overlaps = []
    lens_distances = []

    new_ones = []

    # corresponds to 0-1 coverage of new one...has to be smaller than x to be incorporated and blasted again
    overlap_tolerance = 0.3
    distance_tolerance = 0
    print(f'=== overlap_tolerance = {overlap_tolerance} ===')
    ### CAUTION WITH DISTANCE TOLERANCE: can lead to the possibly unexpected behaviour that running again with same parameters 
    # yields more candidates because of closing the distance gaps iteratively... ###
    print(f'=== distance_tolerance = {distance_tolerance} ===')

    old_idx = sorted(list(bib.keys()))
    old_idx_end = [i+bib[i]['end']-bib[i]['start'] for i in old_idx]

    to_del = set()
    to_del_from_new_idx = set()

    for start,end in new_idx:

        left_check = True
        right_check= True

        # bunch of checks for chromosome breaks but really here i want to see if it starts or ends at chromsome ends the checks shall be turned off
        # it should still be checked below since when using some overlap or distance criterion it will still take it although it shouldnt
        pos = bisect_left(seqlen,start)
        # check if start is at chromsome start
        if start == seqlen[pos]:
            pos += 1
            left_check = False
        ident1 = seqids[pos]
        pos = bisect_left(seqlen,end)
        if end == seqlen[pos]:
            print('this should also never happen')
            pos += 1
        # check if end is at chromsome end
        elif end == seqlen[pos] - 1:
            print('this should also never happen, too (since kmer counts stop k - 1 bases before end)')
            right_check = False
        ident2 = seqids[pos]

        new = False
        
        print('===============')
        print('orig',start,end)
        orig_start = start
        orig_end = end
        # just a check...what im really checking is above that for winodws starting/ending at chromsome ends the checks are turned off
        if ident1 != ident2:
            print('SHOULD NOT HAPPEN FOR NEW WINDOWS')
            print(ident1,ident2)
            continue

        # this is just for counter
        overlaps = 0
        distances = 0

        # check all cases for which candidate is already covered
        # if the ones which are close are considered (distance cases here),
        # it does not make sense to consider overlap cases separately but
        # i keep the cases for future reference/changes

        if start in old_idx:
            cases_in_there += 1
            if not right_check:
                continue
            idx_start_in_old_idx = old_idx.index(start)
            if end > old_idx_end[idx_start_in_old_idx]:
                print(old_idx[idx_start_in_old_idx],old_idx_end[idx_start_in_old_idx])
                try:
                    print(old_idx[idx_start_in_old_idx-1],old_idx_end[idx_start_in_old_idx-1])
                except:
                    pass
                try:
                    print(old_idx[idx_start_in_old_idx+1],old_idx_end[idx_start_in_old_idx+1])
                except:
                    pass
                print('ets see')
                print(old_idx[idx_start_in_old_idx-3:idx_start_in_old_idx+3])
                print(old_idx_end[idx_start_in_old_idx-3:idx_start_in_old_idx+3])
                new = True
                cases_in_there_but_longer += 1
                to_del2 = set()
                print(idx_start_in_old_idx)
                for i in range(idx_start_in_old_idx+1,len(old_idx)):
                    pos = bisect_left(seqlen,old_idx[i])
                    if old_idx[i] == seqlen[pos]:
                        pos += 1
                    ident2 = seqids[pos]
                    if (end < old_idx[i] and old_idx[i] - end >= distance_tolerance) or ident1 != ident2:
                        break
                    else:
                        to_del2.add(i)
                # update end as well
                if end > old_idx_end[i-1]:
                    old_idx_end[idx_start_in_old_idx] = end
                else:
                    end = old_idx_end[i-1]
                    old_idx_end[idx_start_in_old_idx] = old_idx_end[i-1]
                to_del2 = sorted(list(to_del2),reverse=True)
                if len(to_del2) > 0:
                    print(to_del2)
                    for i in to_del2:
                        print(old_idx[i],old_idx_end[i])
                        del old_idx[i]
                        del old_idx_end[i]
                
                print('ets see2')
                print(old_idx[idx_start_in_old_idx-3:idx_start_in_old_idx+3])
                print(old_idx_end[idx_start_in_old_idx-3:idx_start_in_old_idx+3])

                print('case where idx was in there already but new one is bigger (may have contained one or more old ones)')
                print(start,end)
            
            if new:
                new_ones.append((start,end))
            continue
        
        else:
            bigger_idx = bisect_left(old_idx,start)
            
            # new one starts before all others
            # only one side has to be checked, can be continued
            if bigger_idx == 0:
                if not right_check:
                    continue
                pos = bisect_left(seqlen,old_idx[0])
                if old_idx[0] == seqlen[pos]:
                    pos += 1
                ident2 = seqids[pos]
                # new end overlaps with first old one
                if end >= old_idx[0]:
                    # contains old one
                    if end >= old_idx_end[0]:
                        new = True
                        contains_old += 1
                        contained_old += 1
                        old_idx[bigger_idx] = start
                        to_del2 = set()
                        for i in range(bigger_idx+1,len(old_idx)):
                            pos = bisect_left(seqlen,old_idx[i])
                            if old_idx[i] == seqlen[pos]:
                                pos += 1
                            ident2 = seqids[pos]
                            if (end < old_idx[i] and old_idx[i] - end >= distance_tolerance) or ident1 != ident2:
                                break
                            else:
                                to_del2.add(i)
                        # update end as well
                        if end >= old_idx_end[i-1]:
                            old_idx_end[bigger_idx] = end
                        else:
                            end = old_idx_end[i-1]
                            old_idx_end[bigger_idx] = old_idx_end[i-1]

                        to_del2 = sorted(list(to_del2),reverse=True)
                        if len(to_del2) > 0:
                            for i in to_del2:
                                del old_idx[i]
                                del old_idx_end[i]

                        print('case where new is smallest and contains old ones')
                        print(start,end)
                        print(old_idx[bigger_idx],old_idx_end[bigger_idx])
                    
                    elif (end - old_idx[0])/(end-start) <= overlap_tolerance and ident1 == ident2:
                        new = True
                        old_idx[0] = start
                        lens_overlaps.append(end-old_idx[0])
                        overlaps += 1
                        end = old_idx_end[0]
                    else:
                        if ident1 != ident2:
                            print('ident1 != ident2 no 1',ident1,ident2)
                        else:
                            to_del_from_new_idx.add((orig_start,orig_end))
                    
                else:
                    distance = old_idx[0] - end
                    # new one is close to first old one
                    if distance <= distance_tolerance and ident1 == ident2:
                        new = True
                        old_idx[0] = start
                        lens_distances.append(distance)
                        distances += 1
                        end = old_idx_end[0]
                print('case where new is smallest and anything could have happened(same or adjusted start/ends)')
                print(start,end)
                print(old_idx[bigger_idx],old_idx_end[bigger_idx])
                if new:
                    new_ones.append((start,end))
                continue
                    

            # new one starts after all others
            # only one side has to be checked, can be continued
            elif bigger_idx == len(old_idx):
                if not left_check:
                    continue
                pos = bisect_left(seqlen,old_idx[-1])
                if old_idx[-1] == seqlen[pos]:
                    pos += 1
                ident2 = seqids[pos]
                # old end overlaps with new start
                if old_idx_end[-1] >= start and ident1 == ident2:
                    # contains new one
                    if old_idx_end[-1] >= end:
                        contains_new += 1
                        continue
                    elif (old_idx_end[-1] - start)/(end-start) <= overlap_tolerance:
                        new = True
                        old_idx_end[-1] = end
                        lens_overlaps.append(old_idx_end[-1]-start)
                        overlaps += 1
                        start = old_idx[-1]
                    else:
                        if ident1 != ident2:
                            print('ident1 != ident2 no 2',ident1,ident2)
                        else:
                            to_del_from_new_idx.add((orig_start,orig_end))
                
                else:
                    distance = start - old_idx_end[-1]
                    # old end is close to new start
                    if distance <= distance_tolerance and ident1 == ident2:
                        new = True
                        old_idx_end[-1] = end
                        lens_distances.append(distance)
                        distances += 1
                        start = old_idx[-1]
                print('case where new is biggest and anything could have happened(same or adjusted start/ends')
                print(start,end)
                print(old_idx[bigger_idx-1],old_idx_end[bigger_idx-1])
                if new:
                    new_ones.append((start,end))
                continue
                   
            else:
                smaller_idx = bigger_idx - 1
                print('orig smaller/bigger old start-ends')
                print(old_idx[smaller_idx],old_idx_end[smaller_idx])
                print(old_idx[bigger_idx],old_idx_end[bigger_idx])
                # overlap upstreams
                pos = bisect_left(seqlen,old_idx[smaller_idx])
                if old_idx[smaller_idx] == seqlen[pos]:
                    pos += 1
                ident2 = seqids[pos]
                if old_idx_end[smaller_idx] >= start and left_check and ident1 == ident2:
                    # contains new one
                    if old_idx_end[smaller_idx] >= end:
                        contains_new += 1
                        continue
                    overlap = old_idx_end[smaller_idx] - start
                    if overlap/(end-start) <= overlap_tolerance:
                        new = True
                        old_idx_end[smaller_idx] = end
                        lens_overlaps.append(overlap)
                        overlaps += 1
                        start = old_idx[smaller_idx]
                    else:
                        if ident1 != ident2:
                            print('ident1 != ident2 no 3',ident1,ident2)
                        else:
                            to_del_from_new_idx.add((orig_start,orig_end))
                elif left_check and ident1 == ident2:
                    distance = start - old_idx_end[smaller_idx]
                    # small distance upstream
                    if distance <= distance_tolerance:
                        new = True
                        old_idx_end[smaller_idx] = end
                        lens_distances.append(distance)
                        distances += 1
                        start = old_idx[smaller_idx]
                
                # overlap downstream 
                pos = bisect_left(seqlen,old_idx[bigger_idx])
                if old_idx[bigger_idx] == seqlen[pos]:
                    pos += 1
                ident2 = seqids[pos]
                if end >= old_idx[bigger_idx] and right_check and ident1 == ident2:
                    overlap = end - old_idx[bigger_idx]
                    # contains old one
                    if end >= old_idx_end[bigger_idx]:
                        print('ets see')
                        print(old_idx[bigger_idx-3:bigger_idx+3])
                        print(old_idx_end[bigger_idx-3:bigger_idx+3])
                        new = True
                        contains_old += 1
                        contained_old += 1
                        old_idx[bigger_idx] = start
                        to_del2 = set()
                        for i in range(bigger_idx+1,len(old_idx)):
                            pos = bisect_left(seqlen,old_idx[i])
                            if old_idx[i] == seqlen[pos]:
                                pos += 1
                            ident2 = seqids[pos]
                            if (end < old_idx[i] and old_idx[i] - end >= distance_tolerance) or ident1 != ident2:
                                break
                            else:
                                to_del2.add(i)
                        if end > old_idx_end[i-1]:
                            print('case1')
                            old_idx_end[bigger_idx] = end
                        else:
                            print('case2')
                            end = old_idx_end[i-1]
                            old_idx_end[bigger_idx] = old_idx_end[i-1]
                        to_del2 = sorted(list(to_del2),reverse=True)
                        if len(to_del2) > 0:
                            print(len(to_del2))
                            for i in to_del2:
                                del old_idx[i]
                                del old_idx_end[i]

                        print('case where new is overlapping downstream and contains old ones')
                        print(start,end)
                        print(old_idx[bigger_idx-3:bigger_idx+3])
                        print(old_idx_end[bigger_idx-3:bigger_idx+3])
                    
                    elif overlap/(end-start) <= overlap_tolerance and ident1 == ident2:
                        new = True
                        old_idx[bigger_idx] = start
                        lens_overlaps.append(overlap)
                        overlaps += 1
                        end = old_idx_end[bigger_idx]
                    else:
                        if ident1 != ident2:
                            print('ident1 != ident2 no 4',ident1,ident2)
                        else:
                            to_del_from_new_idx.add((orig_start,orig_end))
                
                elif right_check and ident1 == ident2:
                    distance = old_idx[bigger_idx] - end
                    # small distance downstream
                    if distance <= distance_tolerance:
                        new = True
                        old_idx[bigger_idx] = start
                        lens_distances.append(distance)
                        distances += 1
                        end = old_idx_end[bigger_idx]

                if new:
                    new_ones.append((start,end))
                    print('new smaller/bigger')
                    print(old_idx[smaller_idx],old_idx_end[smaller_idx])
                    print(old_idx[bigger_idx],old_idx_end[bigger_idx])
                    print('new one',start,end)
                print('general case where anything could have happened to start/end')
                print(start,end)
                print(old_idx[bigger_idx-3:bigger_idx+3])
                print(old_idx_end[bigger_idx-3:bigger_idx+3])

                # final check to make sure old_idx is good
                # for cases in which both sides get changed
                ### THIS IS CRUCIAL ###
                #if (distances == 1 and overlaps == 1) or distances == 2 or overlaps == 2:
                #    old_idx[smaller_idx] = start
                #    old_idx_end[smaller_idx] = end
                #    del old_idx[bigger_idx]
                #    del old_idx_end[bigger_idx]

                # more robust way to check than with the counters above (independent of counters but harder to grasp immediately)
                if old_idx_end[smaller_idx] >= old_idx[bigger_idx]:
                    # old way, but actually i think its wrong because why not just take the start of smaller and end of bigger...
                    #old_idx[smaller_idx] = start
                    #old_idx_end[smaller_idx] = end
                    old_idx_end[smaller_idx] = old_idx_end[bigger_idx]
                    del old_idx[bigger_idx]
                    del old_idx_end[bigger_idx]
                print('after del')
                print(old_idx[bigger_idx-3:bigger_idx+3])
                print(old_idx_end[bigger_idx-3:bigger_idx+3])

                if overlaps == 1:
                    cases_overlaps_single += 1
                elif overlaps == 2:
                    cases_overlaps_double += 1
                if distances == 1:
                    cases_distances_single += 1
                elif distances == 2:
                    cases_distances_double += 1
                if distances == 1 and overlaps == 1:
                    two_sided_dual += 1

    # get elements to delete:
    # the indices in old_idx (which is now actually new,it was kept updated)
    # should correspond to bib. -> delete all in bib which are not in old_idx

    # also have to check ends (for which starts did not change)

    for i in bib:
        if i not in old_idx:
            to_del.add(i)
        else:
            end = i+bib[i]['end']-bib[i]['start']
            if end != old_idx_end[old_idx.index(i)]:
                to_del.add(i)

    # old_idx is actually new ones (except for really new ones which dont overlap), so might as well use it to collect them
    # otherwise if adding them each iteration there will need to be checks if start
    # and end have already occured (later iteration may change new ones from before)

    new_ones = []

    for i,j in enumerate(old_idx):
        if j not in bib:
            new_ones.append((j,old_idx_end[i]))
        else:
            if old_idx_end[i] != j + bib[j]['end'] - bib[j]['start']:
                new_ones.append((j,old_idx_end[i]))

    # delete ones that are kinda new but overlapped considerably

    for i in to_del_from_new_idx:
        new_idx.remove(i)
    # now add ones which are not covered
    for start,end in new_idx:
        if start not in old_idx:
            bigger_idx = bisect_left(old_idx,start)
            if bigger_idx == 0:
                new_ones.append((start,end))
                continue
            smaller_idx = bigger_idx - 1
            # if start is covered (or end) it is in there (or should be) -> if not, add it
            if start > old_idx_end[smaller_idx]:
                new_ones.append((start,end))
            # if it is supposedly in there, do checks if its fully covered
            else:
                if bigger_idx == 0:
                    # if it is the smallest element, it should be in old_idx directly (as done a few lines above)
                    print(start,end)
                    print(old_idx[bigger_idx],old_idx_end[bigger_idx])
                    print('WRONG 1')
                elif bigger_idx == len(old_idx):
                    if old_idx_end[-1] < end:
                        print(start,end)
                        print(old_idx[smaller_idx],old_idx_end[smaller_idx])
                        print('WRONG 2')
                else:
                    if not (start > old_idx[smaller_idx] and end <= old_idx_end[smaller_idx]):
                        print(start,end)
                        print(old_idx[smaller_idx],old_idx_end[smaller_idx])
                        print(old_idx[bigger_idx],old_idx_end[bigger_idx])
                        print('WRONG 3')

    test = sorted(old_idx.copy())
    test_end = sorted(old_idx_end.copy())

    print('huhu')
    print(test==old_idx)
    print(test_end==old_idx_end)

    for j,i in enumerate(old_idx[:-1]):
        if i > old_idx[j+1]:
            print(i,old_idx[j+1])
            print('damn')
        if i == old_idx[j+1]:
            print(i,i,old_idx[j+1])
            print('damn2')
    for j,i in enumerate(old_idx_end[:-1]):
        if i > old_idx_end[j+1]:
            print(old_idx[j],i,old_idx[j+1],old_idx_end[j+1])
            print('2damn')
        if i == old_idx_end[j+1]:
            print('2damn2')

    print(f'number of old candidates that were deleted: {len(to_del)} (of {len(old_idx)})')
    print(f'{contains_old} new ones contain {contained_old} old ones')
    print(f'{contains_new} new ones were contained in old ones')
    print(f'number of cases were the same index was already covered: {cases_in_there} and how many of the new ones were longer: {cases_in_there_but_longer}')
    print(f'number of cases where a potential new one overlapped on one side: {cases_overlaps_single}')
    print(f'number of cases where a potential new one overlapped on both sides: {cases_overlaps_double}')
    print(f'number of cases where both sides were still changed but one was overlap/one distance: {two_sided_dual}')
    if len(lens_overlaps) == 0:
        print('no overlaps')
    else:
        print(f'mean lengths of overlaps: {np.mean(lens_overlaps)}')
    print(f'number of cases were a potential new one close on one side: {cases_distances_single}')
    print(f'number of cases were a potential new one close on both sides: {cases_distances_double}')
    if len(lens_distances) == 0:
        print('no close ones')
    else:
        print(f'mean lengths of distances: {np.mean(lens_distances)}')

    if len(new_ones) == 0:
        print('no new ones')
    else:
        for i,j in new_ones:
            if j-i > 20000:
                print('bigger than 20000:',i,j)

    return(new_ones,to_del)


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

    k_e = []
    
    if pathlib.Path(f'{work_dir}/multiple_k_e').exists():
        with open(f'{work_dir}/multiple_k_e') as f:
            for line in f:
                if '#' in line:
                    continue
                line = line.strip().split()
                if len(line) != 2:
                    continue
                k,e = line
                k_e.append((k,e))

    k_e = list(set(k_e))


    if len(k_e) == 0:
        print(f'filtering windows with parameters: k={k}, e={e}, len_lmers={len_lmers}, interval={interval}, percentile={percentile}')
        idx = get_low_count_windows(org,infile,outfile,k,e,len_lmers,interval,percentile)
    else:
        print(f'filtering windows with parameters: len_lmers={len_lmers}, interval={interval}, percentile={percentile}')
        print(f' and multiple k and e: {k_e}')
        idx = get_low_count_windows_multiple_k_e(org,infile,outfile,k_e,len_lmers,interval,percentile)
        
    
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
        idx,to_del = combine_and_filter_new_ones(idx,bib)
        if len(idx) == 0:
            print('no new candidates')
            run(f'touch {work_dir}/touch/no_new_candidates_{org}',shell=True)
        
        with open(f'{work_dir}/to_del/{org}/to_del','wb') as f:
            pickle.dump(to_del, f)
    # outfile is indices directory/{org}
    with open(outfile,'wb') as f:
        pickle.dump(idx, f)

    make_k_mers(org,idx,k,e,len_lmers,interval,percentile)
