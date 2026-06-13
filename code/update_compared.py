#! /usr/bin/env python3

import pickle
import pathlib
from sys import argv,getsizeof
from bisect import bisect_left
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from subprocess import run
from subprocesses import get_split_size,set_config_dir

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

        records.append(SeqRecord(Seq(s[i_start:i_end]),id="kmer${}${}${}${}".format(i_start,ident,i_alt,i_end)))
        records_rev.append(SeqRecord(Seq(s[i_start:i_end]).reverse_complement(),id="kmer${}${}${}${}".format(i_start,ident,i_alt,i_end)))
        
    SeqIO.write(records, outpath+'/forward.fasta', "fasta")
    SeqIO.write(records_rev, outpath+'/reverse.fasta', "fasta")

    return(0)

def exe():

    with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
        bib1 = pickle.load(f)
    iss1 = list(bib1.keys())

    to_compare1 = iss1

    write_fastas(org,f'{work_dir}/sequences_to_compare/{org}',to_compare1,bib1)
    # Build makeblastdb command as list (safe from shell injection)
    makeblastdb_cmd = [
        'makeblastdb',
        '-in', f'{work_dir}/sequences_to_compare/{org}/forward.fasta',
        '-out', f'{work_dir}/blastdbs/anchor_candidates_{org}_forward',
        '-dbtype', 'nucl'
    ]
    run(makeblastdb_cmd)
    pairwise_part_size = get_split_size('pairwise')
    # Split fasta-splitter && chain into two separate calls (mimics && behavior)
    fasta_split_forward = [
        'fasta-splitter',
        '--part-size', str(pairwise_part_size),
        '--out-dir', f'{work_dir}/sequences_to_compare/{org}/forward_split',
        f'{work_dir}/sequences_to_compare/{org}/forward.fasta'
    ]
    result = run(fasta_split_forward)
    if result.returncode == 0:
        fasta_split_reverse = [
            'fasta-splitter',
            '--part-size', str(pairwise_part_size),
            '--out-dir', f'{work_dir}/sequences_to_compare/{org}/reverse_split',
            f'{work_dir}/sequences_to_compare/{org}/reverse.fasta'
        ]
        run(fasta_split_reverse)

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]
    set_config_dir(work_dir)

    org = argv[1].strip()

    exe()
