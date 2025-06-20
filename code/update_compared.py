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

    with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
        bib1 = pickle.load(f)
    iss1 = list(bib1.keys())

    to_compare1 = iss1

    write_fastas(org,f'{work_dir}/sequences_to_compare/{org}',to_compare1,bib1)
    run(f'makeblastdb -in "{work_dir}/sequences_to_compare/{org}/forward.fasta" -out {work_dir}/blastdbs/anchor_candidates_{org}_forward -dbtype nucl',shell=True)
    run(f'fasta-splitter --part-size 10000000 --out-dir {work_dir}/sequences_to_compare/{org}/forward_split {work_dir}/sequences_to_compare/{org}/forward.fasta && fasta-splitter --part-size 10000000 --out-dir {work_dir}/sequences_to_compare/{org}/reverse_split {work_dir}/sequences_to_compare/{org}/reverse.fasta',shell=True)

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]

    org = argv[1].strip()

    exe()
