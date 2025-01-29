#! /usr/bin/env python3

from subprocess import run
from os.path import isfile
from os import getcwd

path = getcwd()

with open('orgs','r') as f:
    for line in f:
        if len(line.split()) > 0:
            accession = line.strip()
            if not (isfile(path + "/../genomes/{}.fasta".format(accession))):
                run('datasets download genome accession {} --filename ../genomes/{}.zip --include genome'.format(accession,accession), shell=True)
