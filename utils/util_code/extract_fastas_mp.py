#! /usr/bin/env python3

from subprocess import run
import multiprocessing as mp
from os.path import isfile
from os import getcwd
import sys

# Accept command-line arguments: accessions_file genomes_dir
if len(sys.argv) >= 3:
    accessions_file = sys.argv[1]
    genomes_dir = sys.argv[2]
else:
    # Fallback to old behavior
    accessions_file = 'orgs'
    genomes_dir = '../genomes'

path = getcwd()

def extract(org):
    run(f'unzip {genomes_dir}/{org}.zip -d {genomes_dir}/{org}_unzipped && rm {genomes_dir}/{org}.zip', shell=True)
    run(f'mv {genomes_dir}/{org}_unzipped/ncbi_dataset/data/{org}/*.fna {genomes_dir}/{org}.fasta && rm -r {genomes_dir}/{org}_unzipped', shell=True)

if __name__ == '__main__':
    orgs = []
    with open(accessions_file, 'r') as f:
        for line in f:
            if len(line.split()) > 0:
                accession = line.strip()
                if not isfile(f"{genomes_dir}/{accession}.fasta"):
                    orgs.append(accession)

    with mp.Pool(processes=1) as pool:
        res = pool.map_async(extract,orgs).get()
