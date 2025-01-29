#! /usr/bin/env python3

from subprocess import run
from sys import argv
from os.path import isfile

orgs = []
with open(argv[1]) as f:
    for line in f:
        orgs.append(line.strip())

for org in orgs:
    
    if not isfile("genomes/{}.fasta".format(org)):
        print(org,'genome')
    if not isfile("blastdbs/{}.nsq".format(org)):
        print(org,'blastdb')
    if not isfile("genmap_indices/{}/index.rev.lf.drv".format(org)):
        print(org,'index')
    if not isfile("small_meta/{}".format(org)):
        print(org,'small_meta')
    if not isfile("metadata_genomes/{}".format(org)):
        print(org,'big_meta')
