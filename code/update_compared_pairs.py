#! /usr/bin/env python3

import pickle
import pathlib
from sys import argv,getsizeof
from bisect import bisect_left
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from subprocess import run

def write_fastas(org,outpath,idx,bib):
    with open(f'{root}/utils/metadata_genomes/{org}','rb') as f: 
        seqids,seqlen,s = pickle.load(f)
    
    records = []
    records_rev = []

    for i in idx:
        
        i_start = i
        i_end = i + bib[i]['end'] - bib[i]['start']
        
        pos = bisect_left(seqlen,i_start)
        if i_start == seqlen[pos]:
            pos += 1
        # here i deleted a debug thing...if new error -> it was there before
        ident = seqids[pos]
        # get the i with respect to chromosomes
        if pos == 0:
            x = 0
        else:
            x = 1

        i_alt = i_start - x * seqlen[pos-1]

        records.append(SeqRecord(s[i_start:i_end],id="kmer${}${}${}${}".format(i_start,ident,i_alt,i_end)))
        records_rev.append(SeqRecord(s[i_start:i_end].reverse_complement(),id="kmer${}${}${}${}".format(i_start,ident,i_alt,i_end)))
        
    SeqIO.write(records, outpath+'/forward.fasta', "fasta")
    SeqIO.write(records_rev, outpath+'/reverse.fasta', "fasta")

    return(0)

def exe():

    #to_compare1 = list()
    #to_compare2 = list()
    #
    with open(anchor_dir + '/candidates' + f'/{org1}','rb') as f:
        bib1 = pickle.load(f)
    iss1 = list(bib1.keys())
    #try: 
    #    with open(anchor_dir + '/compared' + f'/{org1}/with_{org2}','rb') as f:
    #        compared1 = pickle.load(f)
    #except:
    #    compared1 = {}

    with open(anchor_dir + '/candidates' + f'/{org2}','rb') as f:
        bib2 = pickle.load(f)
    iss2 = list(bib2.keys())
    #try: 
    #    with open(anchor_dir + '/compared' + f'/{org2}/with_{org1}','rb') as f:
    #        compared2 = pickle.load(f)
    #except:
    #    compared2 = {}

    #for j in iss2:
    #    if j not in compared1:
    #        compared1[j] = bib2[j]['end']
    #        to_compare2.append(j)
    #    else:
    #        if compared1[j] != bib2[j]['end']:
    #            compared1[j] = bib2[j]['end']
    #            to_compare2.append(j)
    #for i in iss1:
    #    if i not in compared2:
    #        compared2[i] = bib1[i]['end']
    #        to_compare1.append(i)
    #    else:
    #        if compared2[i] != bib1[i]['end']:
    #            compared2[i] = bib1[i]['end']
    #            to_compare1.append(i)

    to_compare1 = iss1
    to_compare2 = iss2

    print(f'number of {org1} candidates to be compared to {org2}: {len(to_compare1)}')
    print(f'number of {org2} candidates to be compared to {org1}: {len(to_compare2)}')

    #with open(anchor_dir + '/compared' + f'/{org2}/with_{org1}','wb') as f:
    #    pickle.dump(compared2,f)
    # 
    #with open(anchor_dir + '/compared' + f'/{org1}/with_{org2}','wb') as f:
    #    pickle.dump(compared1,f)
    
    write_fastas(org1,f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org1}',to_compare1,bib1)
    run(f'makeblastdb -in "{work_dir}/sequences_to_compare/{org1}---{org2}/{org1}/forward.fasta" -out {work_dir}/blastdbs/anchor_candidates_{org1}---{org2}_forward -dbtype nucl',shell=True)
    write_fastas(org2,f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org2}',to_compare2,bib2)
    run(f'fasta-splitter --part-size 10000000 --out-dir {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward_split {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward.fasta && fasta-splitter --part-size 10000000 --out-dir {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse_split {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse.fasta',shell=True)

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]

    org1 = argv[1].split('---')[0]
    org2 = argv[1].split('---')[1]

    exe()
