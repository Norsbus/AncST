#!/usr/bin/env python3

from pprint import pprint
import os
import math
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

def get_Cm_values(w,i):
    global_dict = {}
    global_dict[org] = []
    Cms = {}
    with open(f'{work_dir}/../utils/macle_out/{org}/{w}_{i}.txt') as f:
        for line in f:
            Cm = float(line.split()[-1].strip())
            if Cm not in Cms:
                Cms[Cm] = 0
            Cms[Cm] += 1
    s_unique = sorted(list(set(list(Cms.keys()))),reverse=True)
    tot = 0
    for Cm in s_unique:
        tot += Cms[Cm]
        global_dict[org].append((Cm,tot))
    return(global_dict)

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
        # could still be (and happens) when old and new ones are combined and lie between chromosome
        # should be dealt with in combine_and_filter tho (just a check because new ones should never cover different chromosomes)
        ident2 = seqids[bisect_left(seqlen,end)]
        if end == seqlen[pos]:
            pos += 1
        if ident != ident2:
            # prints just to have look
            print(f'for window with start {i} and end {end} the start chromosome {ident} is not the same as end chromsome {ident2}')
            print('lens chromosomes (abs ends):')
            print(f'{seqlen}')
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

def get_low_count_windows(org,genmap_out_file,out_file_name,k,e,len_window,interval,percentile,wanted_genome_size):

    l = int(len_window)

    percentile = int(percentile)

    interval  = int(interval)

    global kmer_counts
    kmer_counts = np.fromfile(genmap_out_file, dtype=np.uint16)
    
    max_concur_proc = int(argv[5])

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
    
    lll = sorted(list(zip(sums_global,indices_global)))

    how_many = int(wanted_genome_size/len_lmers)
    final_idx = set([xxx[1] for xxx in lll[:how_many]])
    last = lll[how_many-1][0]
    for xxx in lll[how_many:]:
        if xxx[0] <= last:
            final_idx.add(xxx[1])

    return(final_idx)

def get_low_count_windows_multiple_k_e(org,genmap_out_file,out_file_name,k_e,len_window,interval,percentile):
    l = int(len_window)

    percentile = int(percentile)

    interval  = int(interval)

    global kmer_counts
    kmer_counts = []

    for k,e in k_e:
        
        genmap_out_file = f'{work_dir}/../utils/genmap_out/{org}/{k}_{e}.freq16'

        kmer_counts_x = np.fromfile(genmap_out_file, dtype=np.uint16)
        kmer_counts_x = kmer_counts_x/np.sum(kmer_counts_x)
        if len(kmer_counts) == 0:
            kmer_counts = kmer_counts_x
        else:
            kmer_counts += kmer_counts_x
    
    max_concur_proc = int(argv[5])

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

def get_idx_from_macle(macle_file,half_the_window_size,pitch,percentile,global_Cm,wanted_genome_size):
    starts = []
    Cms = []
    lll = []
    with open(f'{work_dir}/../utils/macle_out/{org}/{half_the_window_size*2}_{pitch}.txt') as f:
        for line in f:
            chromo,midpoint,Cm = line.strip().split()
            Cm = float(Cm)
            midpoint = int(midpoint)
            idx_chr = seqids.index(chromo)
            if idx_chr == 0:
                le = 0
            else:
                le = seqlen[idx_chr-1]
            # macle output cuts half_the_window_sizes at end of chromosomes and uses "rest" as new half_the_window_size start
            if midpoint < half_the_window_size:
                abs_start = le
            else:
                abs_start = midpoint - half_the_window_size + le
            starts.append(abs_start)
            Cms.append(Cm)
            lll.append((Cm,abs_start))

    lll = sorted(lll,reverse=True)
    s_unique = sorted(list(set([xxx[0] for xxx in lll])),reverse=True)
    unique_d = {}
    for xxx in lll:
        if xxx[0] not in unique_d:
            unique_d[xxx[0]] = []
        unique_d[xxx[0]].append(xxx[1])
    final_idx = set([xxx[1] for xxx in lll if xxx[0] >= global_Cm])

    return(final_idx)


def combine_and_filter_new_ones(idx,bib=False):

    new_idx = []

    i = 0

    if percentile < 31:

        while i < len(idx)-1:
            print(f'processing {idx[i]}')
            idx_start = idx[i]
            if idx_start in genmap_info:
                len_lmers,interval = genmap_info[idx_start]
            elif idx_start in macle_info:
                len_lmers,interval = macle_info[idx_start]
            else:
                print(f'err {idx_start} not from macle nor genmap...?...')
            idx_end = idx[i] + len_lmers
            pos = bisect_left(seqlen,idx_start)
            if idx_start == seqlen[pos]:
                pos += 1
            if idx_end > seqlen[pos] - 1:
                # means that its last window on chromosome and midpoint + interval is bigger than chromosome
                print(f'last, not appending because too short:{idx_start} - {idx_end}')
                i += 1
                pos_tmp = bisect_left(seqlen,idx[i])
                if idx[i] == seqlen[pos_tmp]:
                    pos_tmp += 1
                ident_tmp = seqids[pos_tmp]
                while ident_tmp == ident1 and i < len(idx)-1:
                    i += 1
                    pos_tmp = bisect_left(seqlen,idx[i])
                    if idx[i] == seqlen[pos_tmp]:
                        pos_tmp += 1
                    ident_tmp = seqids[pos_tmp]
                continue
            ident1 = seqids[pos]
            pos = bisect_left(seqlen,idx[i+1])
            if idx[i+1] == seqlen[pos]:
                pos += 1
            ident2 = seqids[pos]
            if ident1 != ident2:
                # means that first one is exactly end of chromosme and next one already on next chromosome
                print(f'first one exactly end of chromosome, appending {idx_start} - {idx_end}')
                new_idx.append((idx_start,idx_end))
                i += 1
                pos_tmp = bisect_left(seqlen,idx[i])
                if idx[i] == seqlen[pos_tmp]:
                    pos_tmp += 1
                ident_tmp = seqids[pos_tmp]
                while ident_tmp == ident1 and i < len(idx)-1:
                    i += 1
                    pos_tmp = bisect_left(seqlen,idx[i])
                    if idx[i] == seqlen[pos_tmp]:
                        pos_tmp += 1
                    ident_tmp = seqids[pos_tmp]
                continue
            end2 = idx[i+1]+len_lmers
            if end2 > seqlen[pos] - 1:
                # means that next one is at end of chromosome. i += 2 cos 2 processed since i and i + 1 combined
                end2 = seqlen[pos] - 1
                print(f'second one is last one, appending {idx_start} - {end2}')
                new_idx.append((idx_start,end2))
                i += 2
                pos_tmp = bisect_left(seqlen,idx[i])
                if idx[i] == seqlen[pos_tmp]:
                    pos_tmp += 1
                ident_tmp = seqids[pos_tmp]
                while ident_tmp == ident1 and i < len(idx)-1:
                    i += 1
                    pos_tmp = bisect_left(seqlen,idx[i])
                    if idx[i] == seqlen[pos_tmp]:
                        pos_tmp += 1
                    ident_tmp = seqids[pos_tmp]
                continue
            while i < len(idx)-1 and (idx_end >= idx[i+1] or idx[i+1] - idx_end <= 100) and ident1 == ident2:
                next_i = idx[i+1]
                if next_i in genmap_info:
                    len_lmers,interval = genmap_info[next_i]
                elif next_i in macle_info:
                    len_lmers,interval = macle_info[next_i]
                else:
                    print(f'err {next_i} not from macle nor genmap...?...')
                idx_end = idx[i+1] + len_lmers
                i += 1
                if idx_end > seqlen[seqids.index(ident1)] - 1:
                    idx_end = seqlen[seqids.index(ident1)] - 1
                    pos_tmp = bisect_left(seqlen,idx[i])
                    if idx[i] == seqlen[pos_tmp]:
                        pos_tmp += 1
                    ident_tmp = seqids[pos_tmp]
                    while ident_tmp == ident1 and i < len(idx)-1:
                        i += 1
                        pos_tmp = bisect_left(seqlen,idx[i])
                        if idx[i] == seqlen[pos_tmp]:
                            pos_tmp += 1
                        ident_tmp = seqids[pos_tmp]
                    break
                # have to check again...pre-calculating ident2 for checks in while loops which uses i+1 but when used here it could be out of range
                if i < len(idx)-1:
                    pos = bisect_left(seqlen,idx[i+1])
                    if idx[i+1] == seqlen[pos]:
                        pos += 1
                    ident2 = seqids[pos]
                else:
                    break
            pos = bisect_left(seqlen,idx_end)
            if idx_end == seqlen[pos]:
                pos += 1
            ident2 = seqids[pos]
            if ident1 != ident2:
                print('why 1 diff chrs',idx_start,idx_end,ident1,ident2,seqlen[seqids.index(ident1)],seqlen[seqids.index(ident2)])
                continue
            print(f'normal one, appending {idx_start} - {idx_end}')
            new_idx.append((idx_start,idx_end))
            i += 1

        if len(new_idx) > 0 and len(idx) > 0:
            if new_idx[-1][1] != idx[-1] + len_lmers:
                idx_start = idx[-1]
                ###
                pos = bisect_left(seqlen,idx[-1])
                if idx[-1] == seqlen[pos]:
                    pos += 1
                ident1 = seqids[pos]
                if idx[-1] + len_lmers > seqlen[pos] - 1:
                    idx_end = seqlen[pos] - 1

                else:
                    idx_end = idx[-1] + len_lmers
                # if it was last, might actually be in there, so check and process only otherwise
                if idx_end != new_idx[-1][1]:
                    pos = bisect_left(seqlen,idx_end)
                    if idx_end == seqlen[pos]:
                        pos += 1
                    ident2 = seqids[pos]
                    if ident1 != ident2:
                        print('why 2 diff chrs',idx[-1],idx[-1]+len_lmers,ident1,ident2,seqlen[seqids.index(ident1)],seqlen[seqids.index(ident2)])
                    ###
                    print(f'last normal one, appending {idx_start} - {idx_end}')
                    new_idx.append((idx_start,idx_end))

    else:
        for i in idx:
            if idx_start in genmap_info:
                len_lmers,interval = genmap_info[idx_start]
            elif idx_start in macle_info:
                len_lmers,interval = macle_info[idx_start]
            else:
                print(f'err {idx_start} not from macle nor genmap...?...')
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
                if i < 0:
                    i = 0
                elif i >= len(old_idx_end):
                    i = len(old_idx_end) - 1
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
    infile_genmap = argv[2]
    infile_macle = argv[3]
    outfile = argv[4]

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]

    genome_size = os.path.getsize(f'{work_dir}/../utils/genomes/{org}.fasta')

    max_concur_proc = int(argv[5])

    k,e,len_lmers,interval,percentile = int(argv[6]),int(argv[7]),int(argv[8]),int(argv[9]),int(argv[10])

    with open(f'{work_dir}/../utils/metadata_genomes/{org}','rb') as f:
        seqids,seqlen,s = pickle.load(f)

    ### macle ###

    macle_paras = []
    if os.path.isfile('macle_params.txt'):
        with open('macle_params.txt') as f:
            for line in f:
                if org in line:
                    macle_paras.append(line.strip().split()[1:])
    
    idx1 = []
    best1 = []
    macle_info = {}
    for paras in macle_paras:
        len_lmers = int(paras[0]) 
        interval = int(paras[1]) 
        wanted_genome_size = int(int(paras[2])/100*genome_size)

        Cm_values = get_Cm_values(len_lmers,interval)
        Cm_value = 1e6
        tot = 0
        iii = 0
        while tot < wanted_genome_size:
            Cm_value = Cm_values[org][iii][0]
            tot += Cm_values[org][iii][1]*len_lmers
            iii+=1

        x = get_idx_from_macle(infile_macle,int(len_lmers/2),int(interval),percentile,Cm_value,wanted_genome_size)
        for iii in x:
            macle_info[iii] = [len_lmers,interval]
        if int(paras[3]) == 0:
            best1.append(x)
        elif int(paras[3]) == 1:
            idx1.append(x)
        else:
            print('wrong macle paras',macle_paras)

    ### genmap ###

    genmap_paras = []
    if os.path.isfile('genmap_params.txt'):
        with open('genmap_params.txt') as f:
            for line in f:
                if org in line:
                    genmap_paras.append(line.strip().split()[1:])

    idx2 = []
    best2 = []
    genmap_info = {}
    for paras in genmap_paras:

        k = int(paras[0])
        e = int(paras[1])
        len_lmers = int(paras[2])
        interval = int(paras[3])
        wanted_genome_size = int(int(paras[4])/100*genome_size)
        if wanted_genome_size == 0:
            continue
        infile = f'{work_dir}/../utils/genmap_out/{org}/{k}_{e}.freq16'
        x = get_low_count_windows(org,infile,outfile,k,e,len_lmers,interval,percentile,wanted_genome_size)
        for iii in x:
            genmap_info[iii] = [len_lmers,interval]
        if int(paras[5]) == 0:
            best2.append(x)
        elif int(paras[5]) == 1:
            idx2.append(x)
        else:
            print('wrong genmap paras',genmap_paras)

    idx = set()

    for s1 in idx1:
        for s2 in idx2:
            pos_macle = []
            pos_genmap = []
            for i in s1:
                pos_macle.append((i,i+macle_info[i][0]))
            for i in s2:
                pos_genmap.append((i,i+genmap_info[i][0]))

            pos_macle = sorted(pos_macle)
            mstarts = [xxx[0] for xxx in pos_macle]
            mends = [xxx[1] for xxx in pos_macle]
            pos_genmap = sorted(pos_genmap)
            gstarts = [xxx[0] for xxx in pos_genmap]
            gends = [xxx[1] for xxx in pos_genmap]
            for ms,me in pos_macle:
                bigger_end_than_start_idx = bisect_left(gends,ms)
                bigger_start_than_end_idx = max(bisect_left(gstarts,me),0)
                if bigger_start_than_end_idx - bigger_end_than_start_idx < 1:
                    continue
                for gs,ge in pos_genmap[max(0,bigger_end_than_start_idx-10):min(bigger_start_than_end_idx+10,len(pos_genmap))]:
                    t1 = (ms,me)
                    t2 = (gs,ge)
                    if overlaps(t1,t2):
                        idx.add(ms)
                        idx.add(gs)

    for s1 in best1:
        for i in s1:
            idx.add(i)
    for s2 in best2:
        for i in s2:
            idx.add(i)

    idx = sorted(list(idx))
    
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
    # chromo check
    to_del = set()
    for t in idx:
        pos = bisect_left(seqlen,t[0])
        if t[0] == seqlen[pos]:
            pos += 1
        ident1 = seqids[pos]
        pos = bisect_left(seqlen,t[1])
        if t[1] == seqlen[pos]:
            pos += 1
        ident2 = seqids[pos]
        if ident1 != ident2:
            print(f'ident diff : {t} {ident1} {ident2}')
            to_del.add(t)
    for i in sorted(list(to_del),reverse=True):
        idx.remove(i)


    with open(outfile,'wb') as f:
        pickle.dump(idx, f)
    
    for enu,t in enumerate(idx[:-1]):
        if idx[enu+1][0] <= t[1]:
            print(f'overlappping: {t} {idx[enu+1]}')
    tot_len = sum([x[1] - x[0] for x in idx])
    print(f'{len(idx)} of tot len {tot_len} non-overlapping candidates are processed')

    make_k_mers(org,idx,k,e,len_lmers,interval,percentile)
