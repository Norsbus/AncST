#! /usr/bin/env python3

from subprocess import run
from os.path import isfile
from os import getcwd
import os
from Bio import SeqIO
import pickle
from sys import argv
import multiprocessing as mp
from os.path import isfile

path = getcwd()

def makeblastdb(org):
    path = getcwd()
    if not (isfile(path + "/../blastdbs/{}.nin".format(org)) and isfile(path + "/../blastdbs/{}.nhr".format(org)) and isfile(path + "/../blastdbs/{}.nsq".format(org))):
        run("makeblastdb -in ../genomes/{}.fasta -out ../blastdbs/{} -dbtype nucl".format(org,org).split())

def make_metadata(org):
    
    path = getcwd()
    if not isfile(path + "/../metadata_genomes/{}".format(org)):
        
        seqs = SeqIO.parse('../genomes/{}.fasta'.format(org), "fasta")
        s = ''
        seqids = []
        seqlen = []
        tot = 0

        for seq in seqs:
            s += seq.seq
            seqids.append(seq.id)
            tot += len(seq.seq)
            seqlen.append(tot)

        with open('../metadata_genomes/{}'.format(org), 'wb') as f:
            pickle.dump((seqids, seqlen, s), f)

    return(0)

def make_small_metadata(org):
    
    path = getcwd()
    if not isfile(path + "/../small_meta/{}".format(org)):
        
        seqs = SeqIO.parse('../genomes/{}.fasta'.format(org), "fasta")
        seqids = []
        seqlen = []
        tot = 0
        for seq in seqs:
            seqids.append(seq.id)
            tot += len(seq.seq)
            seqlen.append(tot)

        with open('../small_meta/{}'.format(org), 'wb') as f:
            pickle.dump((seqids, seqlen), f)

    return(0)

def exec_funcs(org):
    makeblastdb(org)
    make_metadata(org)
    make_small_metadata(org)
    return 0

if __name__ == "__main__":
    if len(argv) > 1:
        orgs = [argv[1]]
    else:
        orgs = []
        with open('orgs','r') as f:
            for line in f:
                orgs.append(line.strip())

    with mp.Pool(processes=11) as pool:
        res = pool.map_async(exec_funcs,orgs).get()
