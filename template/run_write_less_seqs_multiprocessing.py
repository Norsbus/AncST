#!/usr/bin/env python3

from subprocess import run,check_output
import os,pickle
from time import sleep
import pathlib
from sys import argv,getsizeof
from bisect import bisect_left
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import multiprocessing as mp

def write_fastas(org,outpath,idx,bib):
    with open(f'{work_dir}/../utils/metadata_genomes/{org}','rb') as f: 
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

def exe(org_tuple):
    org1 = org_tuple.split('---')[0]
    org2 = org_tuple.split('---')[1]
    #print('gi')

    #orgs = []
    #with open(f"orgs","r") as f:
    #    for line in f:
    #        org = line.strip()
    #        orgs.append(org)

    #orgs_tuples = []
    #files = []
    #orgs_tuples_done = set()
    #for i,org1 in enumerate(orgs):
    #    for org2 in orgs[i+1:]:
    #        orgs_tuple = org1 + '---' + org2
    #        if not os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'):
    #            orgs_tuple = org2 + '---' + org1
    #        chunks = [name.split('forward.')[1] for name in os.listdir(f'{work_dir}/sequences_to_compare/{orgs_tuple}/{org2}/forward_split')]
    #        for f in chunks:
    #            if os.path.isfile(f'{work_dir}/touch/{orgs_tuple}_parse_bcamm_done_{f}'):
    #                orgs_tuples_done.add(orgs_tuple)
    #                continue
    #            orgs_tuples.append(orgs_tuple)
    #            files.append(f)
    ##print(orgs_tuples_done)
    #orgs_tuples = set(orgs_tuples)
    ##print(orgs_tuples.intersection(orgs_tuples_done))
    ##for org_tuple in orgs_tuples_done:
    ##    if org_tuple in orgs_tuples:
    ##        print(org_tuple)
    ##print(len(files))
    #with open('orgs_tuple_to_bcamm','wb') as f:
    #    pickle.dump(orgs_tuples,f)
    #exit(0)

    to_compare1 = list()
    to_compare2 = list()
    
    with open(anchor_dir + '/candidates' + f'/{org1}','rb') as f:
        bib1 = pickle.load(f)
    iss1 = list(bib1.keys())

    with open(anchor_dir + '/candidates' + f'/{org2}','rb') as f:
        bib2 = pickle.load(f)
    iss2 = list(bib2.keys())

    for j in iss2:
        l_40 = 0
        for jj,s in enumerate(bib2[j]['scores_regions']):
            if s <= 40:
                l_40 += bib2[j]['regions'][jj+1] - bib2[j]['regions'][jj]
        if l_40 >= (bib2[j]['end']-bib2[j]['start'])*0.8:
            continue
        to_compare2.append(j)
    for i in iss1:
        #l_40 = 0
        #for ii,s in enumerate(bib1[i]['scores_regions']):
        #    if s <= 40:
        #        l_40 += bib1[i]['regions'][ii+1] - bib1[i]['regions'][ii]
        #if l_40 >= (bib1[i]['end']-bib1[i]['start'])*0.8:
        #    continue
        to_compare1.append(i)

    print(f'number of {org1} candidates to be compared to {org2}: {len(to_compare1)} of {len(iss1)}')
    print(f'number of {org2} candidates to be compared to {org1}: {len(to_compare2)} of {len(iss2)}')
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward.fasta',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse.fasta',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward_split/*',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse_split/*',shell=True)
    write_fastas(org1,f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org1}',to_compare1,bib1)
    run(f'makeblastdb -in "{work_dir}/sequences_to_compare/{org1}---{org2}/{org1}/forward.fasta" -out {work_dir}/blastdbs/anchor_candidates_{org1}---{org2}_forward -dbtype nucl',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward.fasta',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse.fasta',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward_split/*',shell=True)
    run(f'rm {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse_split/*',shell=True)
    write_fastas(org2,f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org2}',to_compare2,bib2)
    run(f'fasta-splitter --part-size 10000000 --out-dir {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward_split {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward.fasta && fasta-splitter --part-size 10000000 --out-dir {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse_split {work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse.fasta',shell=True)
    return 0

#with open('orgs_tuple_to_bcamm','rb') as f:
#    otto = pickle.load(f)

work_dir = '/scr/k80san/karl//insects/stable_synteny/template'
code_dir = '/scr/k80san/karl//insects/stable_synteny/code'
anchor_dir = '/scr/k80san/karl//insects/stable_synteny/anchors'
orgs12 = []
flies = []
morgs = []
with open('methodenpaper_orgs') as f:
    for line in f:
        morgs.append(line.strip())
with open('orgs12') as f:
    for line in f:
        orgs12.append(line.strip())
with open('58orgs') as f:
    for line in f:
        flies.append(line.strip())
otto = []
for org1 in orgs12:
    for org2 in flies:
        orgs_tuple = org1 + '---' + org2
        if not os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'):
            orgs_tuple = org2 + '---' + org1
        #if not os.path.isfile(f'/scr/k80san/karl//insects/stable_synteny/template/blastdbs/anchor_candidates_{orgs_tuple}_forward.nsq'):
        otto.append(orgs_tuple)
print(len(otto))
for i,org1 in enumerate(orgs12[:-1]):
    for org2 in orgs12[i+1:]:
        if org1 in morgs or org2 in morgs:
            continue
        orgs_tuple = org1 + '---' + org2
        if not os.path.isdir(f'{work_dir}/clasp_out_forward/{orgs_tuple}'):
            orgs_tuple = org2 + '---' + org1
        #if not os.path.isfile(f'/scr/k80san/karl//insects/stable_synteny/template/blastdbs/anchor_candidates_{orgs_tuple}_forward.nsq'):
        otto.append(orgs_tuple)
print(len(otto))

with mp.Pool(processes=58) as pool:
    res = pool.map_async(exe,otto).get()
