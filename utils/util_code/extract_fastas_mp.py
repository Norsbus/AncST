#! /usr/bin/env python3

from subprocess import run
import multiprocessing as mp
from os.path import isfile
from os import getcwd

path = getcwd()

def extract(org):
    run(f'unzip ../genomes/{org}.zip -d ../genomes/{org}_unzipped && rm ../genomes/{org}.zip', shell=True)
    run(f'mv ../genomes/{org}_unzipped/ncbi_dataset/data/{org}/*.fna ../genomes/{org}.fasta && rm -r ../genomes/{org}_unzipped', shell=True)

orgs = []
with open('orgs','r') as f:
    for line in f:
        if len(line.split()) > 0:
            accession = line.strip()
            if not (isfile(path + "/../genomes/{}.fasta".format(accession))):
                orgs.append(accession)

with mp.Pool(processes=30) as pool:
    res = pool.map_async(extract,orgs).get()
