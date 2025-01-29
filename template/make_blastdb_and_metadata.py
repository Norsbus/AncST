#! /usr/bin/env python3

from subprocess import run
from os.path import isfile
from os import getcwd
import os
from Bio import SeqIO
import pickle
from sys import argv

def makeblastdb(org):
    path = getcwd()
    if not (isfile(path + "/../utils/blastdbs/{}.nin".format(org)) and isfile(path + "/../utils/blastdbs/{}.nhr".format(org)) and isfile(path + "/../utils/blastdbs/{}.nsq".format(org))):
        run("makeblastdb -in ../utils/genomes/{}.fasta -out ../utils/blastdbs/{} -dbtype nucl".format(org,org).split())

def make_metadata(org):
    
    path = getcwd()
    if not isfile(path + "/../utils/metadata_genomes/{}".format(org)):
        
        seqs = SeqIO.parse('../utils/genomes/{}.fasta'.format(org), "fasta")
        s = ''
        seqids = []
        seqlen = []
        tot = 0

        for seq in seqs:
            s += seq.seq
            seqids.append(seq.id)
            tot += len(seq.seq)
            seqlen.append(tot)

        with open('../utils/metadata_genomes/{}'.format(org), 'wb') as f:
            pickle.dump((seqids, seqlen, s), f)

    return(0)

def make_small_metadata(org):
    
    path = getcwd()
    if not isfile(path + "/../utils/small_meta/{}".format(org)):
        
        seqs = SeqIO.parse('../utils/genomes/{}.fasta'.format(org), "fasta")
        seqids = []
        seqlen = []
        tot = 0
        for seq in seqs:
            seqids.append(seq.id)
            tot += len(seq.seq)
            seqlen.append(tot)

        with open('../utils/small_meta/{}'.format(org), 'wb') as f:
            pickle.dump((seqids, seqlen), f)

    return(0)

if __name__ == "__main__":
    if len(argv) > 1:
        orgs = [argv[1]]
    else:
        orgs = []
        with open('orgs','r') as f:
            for line in f:
                orgs.append(line.strip())
    for org in orgs:
        makeblastdb(org)
        make_metadata(org)
        make_small_metadata(org)
