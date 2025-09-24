#! /usr/bin/env python3

from subprocess import run
from os.path import isfile
from os import getcwd
import multiprocessing as mp

path = getcwd()


def download(org):
    if not (isfile(path + "/../genomes/{}.fasta".format(org))):
        run('datasets download genome accession {} --filename ../genomes/{}.zip --include genome'.format(org,org), shell=True)

orgs = []
with open('orgs','r') as f:
    for line in f:
        if len(line.split()) > 0:
            accession = line.strip()
            orgs.append(accession)

with mp.Pool(processes=1) as pool:
    res = pool.map_async(download,orgs).get()
