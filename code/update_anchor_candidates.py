#! /usr/bin/env python3

import pickle
from sys import argv
import pathlib,os
from bisect import bisect_left
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
from datasketch import MinHash,MinHashLSHForest
import more_itertools

def update_anchors(bib,bib_aligned):

    new_ones = {}
    
    with open(f'{root}/utils/metadata_genomes/{org}','rb') as f:
        seqids,seqlen,s = pickle.load(f)

    with open(f'{work_dir}/indices/{org}/indices_global',"rb") as f:
        idx = pickle.load(f)

    idx = sorted(idx) 

    for i in idx:

        
        if i[0][0] in bib:
            print('should not happen as they shouldve been deleted (still continuing since the old one should then be overwritten)')
        
        chromo = which_chromo(seqlen,seqids,i[0][0])
        idx_l = seqids.index(chromo)
        if idx_l > 0:
            l = seqlen[idx_l-1]
        else:
            l = 0
        chromo_start = i[0][0] - l
        chromo_end = i[0][1] - l

        regions = [1]
        scores_regions = []
       
        i[1] = sorted(i[1],key=lambda x: x[0][0])

        for s_e,score in i[1]:
            if s_e[0] < regions[-1]:
                pass
                #print('wtf')
            elif s_e[0] == regions[-1]:
                pass
                #print('region same start as last end')
            else:
                scores_regions.append(40)
                regions.append(s_e[0])

            scores_regions.append(score)
            regions.append(s_e[1])

        if regions[-1] != i[0][1]-i[0][0]:
            scores_regions.append(40)
            regions.append(i[0][1]-i[0][0])
        

        bib[i[0][0]] = {\
                        #'kmer_profile':kmer_counts[i[0][0]:i[0][1]],\
                        'regions':regions,\
                        'scores_regions':scores_regions,\
                        'chromosome':chromo,\
                        'start':chromo_start,\
                        'end':chromo_end,\
			'blast_word_size': ws,\
                        'matches':{}\
                        }
        
        bib_aligned[i[0][0]] = {\
                        #'kmer_profile':kmer_counts[i[0][0]:i[0][1]],\
                        'regions':regions,\
                        'scores_regions':scores_regions,\
                        'chromosome':chromo,\
                        'start':chromo_start,\
                        'end':chromo_end,\
                        'blast_word_size': ws,\
                        'matches':{}\
                        }

        new_ones[i[0][0]] = bib[i[0][0]]

    return bib,bib_aligned,new_ones

def init_dict(org):

    bib = {}
    with open(f'{root}/utils/metadata_genomes/{org}','rb') as f:
        seqids,seqlen,s = pickle.load(f)

    with open(f'{work_dir}/indices/{org}/indices_global',"rb") as f:
        idx = pickle.load(f)

    idx = sorted(idx) 
    
    for i in idx:
        chromo = which_chromo(seqlen,seqids,i[0][0])
        idx_l = seqids.index(chromo)
        if idx_l > 0:
            l = seqlen[idx_l-1]
        else:
            l = 0
        chromo_start = i[0][0] - l
        chromo_end = i[0][1] - l

        regions = [1]
        scores_regions = []
       
        i[1] = sorted(i[1],key=lambda x: x[0][0])
        
        for s_e,score in i[1]:
            if s_e[0] < regions[-1]:
                pass
                #print('wtf')
            elif s_e[0] == regions[-1]:
                pass
                #print('region same start as last end')
            else:
                scores_regions.append(40)
                regions.append(s_e[0])

            scores_regions.append(score)
            regions.append(s_e[1])

        if regions[-1] != i[0][1]-i[0][0]:
            scores_regions.append(40)
            regions.append(i[0][1]-i[0][0])
        

        bib[i[0][0]] = {\
                        #'kmer_profile':kmer_counts[i[0][0]:i[0][1]],\
                        'regions':regions,\
                        'scores_regions':scores_regions,\
                        'chromosome':chromo,\
                        'start':chromo_start,\
                        'end':chromo_end,\
                        'blast_word_size': ws,\
                        'matches':{}\
                        }
    
    return bib

def which_chromo(seqlen,seqids,i):
    pos = bisect_left(seqlen,i)
    if i == seqlen[pos]:
        pos += 1
    return(seqids[pos])

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]

    org = argv[1]
    
    with open(f"{work_dir}/params","r") as f:
        params = f.readline().strip()
        params_list = params.split()
        ws = params_list[5]
        kmer_k = params_list[0]
        kmer_e = params_list[1]

    #enmap_out_file = f'{root}/utils/genmap_out/{org}/{kmer_k}_{kmer_e}.freq16'
    #kmer_counts = np.fromfile(genmap_out_file, dtype=np.uint16)
    
    global_idx = []
    indices = [f for f in os.listdir(f'{work_dir}/indices/{org}') if 'indices' in f and org in f and 'part' in f]
    try:
        indices += [f for f in os.listdir(f'{work_dir}/old_indices/{org}') if 'indices' in f and org in f and 'part' in f]
    except:
        print(f'no old indices for {org}')
    for index in indices:
        with open(f'{work_dir}/indices/{org}/' + index,'rb') as f:
            idx = pickle.load(f)
        global_idx += idx
    with open(f'{work_dir}/indices/{org}/indices_global','wb') as f:
        pickle.dump(global_idx,f)
    if len(global_idx) == 0:
        print(f'returning since most likely {org}\'s anchors were not updated and len(global_idx) == 0')
        exit(0) 
    # if first run, init bib, otherwise mostly there would be an laigned map as well (cos it would have run for a set of orgs
    # but in testing i wanted to run for only one without craeting an aligned map so theres another exception to update the candidates without updateing aligned)
    try:
        with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
            bib = pickle.load(f)
        
        with open(anchor_dir + '/aligned' + f'/{org}','rb') as f:
            bib_aligned = pickle.load(f)

        bib,bib_aligned,new_ones = update_anchors(bib,bib_aligned)
        with open(anchor_dir + '/candidates' + f'/{org}','wb') as f:
            pickle.dump(bib,f)
        with open(anchor_dir + '/aligned' + f'/{org}','wb') as f:
            pickle.dump(bib_aligned,f)
        
    except:
        try:
            with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
                bib = pickle.load(f)

            bib_aligned = {}
            bib,bib_aligned,new_ones = update_anchors(bib,bib_aligned)

        except:
        
            bib = init_dict(org)
    
        with open(anchor_dir + '/candidates' + f'/{org}','wb') as f:
            pickle.dump(bib,f)
