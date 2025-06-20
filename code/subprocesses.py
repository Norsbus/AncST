#! /usr/bin/env python3

from subprocess import run
from os.path import isfile
from os import getcwd
import os
from Bio import SeqIO
import pickle

def makeblastdb(org):
    path = getcwd()
    if not (isfile(path + "/blastdbs/{}.nin".format(org)) and isfile(path + "/blastdbs/{}.nhr".format(org)) and isfile(path + "/blastdbs/{}.nsq".format(org))):
        run("makeblastdb -in OTHER/scr/k80san/karl/genomes/{}.fasta -out blastdbs/{} -dbtype nucl".format(org,org).split())

def blast(db,query,outfile,word_size):
    cmd = 'blastn -db {} -query {} -strand plus -outfmt 6 -evalue 1e-3 -word_size {} -out {}'.format(db,query,word_size,outfile).split()
    run(cmd)
    return(0)

def clasp(infile,outfile,l=0.5,e=0):
    bin_dir = '../bin'
    cmd = '{}/clasp.x -m -i {} -c 7 8 9 10 12 -C 1 2 -l {} -e {} -o {}'.format(bin_dir,infile,l,e,outfile).split()
    run(cmd)
    return(0)

def make_metadata(org):
    
    path = getcwd()
    if not isfile(path + "/../utils/metadata_genomes/{}".format(org)):
        
        seqs = SeqIO.parse('OTHER/scr/k80san/karl/genomes/{}.fasta'.format(org), "fasta")
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
