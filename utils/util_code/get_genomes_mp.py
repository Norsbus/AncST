#! /usr/bin/env python3

from subprocess import run
from os.path import isfile
from os import getcwd
import multiprocessing as mp
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

# The conda environment should be activated by the wrapper script
# so 'datasets' should be available in PATH
DATASETS_CMD = 'datasets'


def download(org):
    output_path = f"{genomes_dir}/{org}.fasta"
    zip_path = f"{genomes_dir}/{org}.zip"
    if not isfile(output_path) and not isfile(zip_path):
        # Use shell=False with explicit args for better error handling and PATH resolution
        run([DATASETS_CMD, 'download', 'genome', 'accession', org, '--filename', zip_path, '--include', 'genome'])

if __name__ == '__main__':
    orgs = []
    with open(accessions_file, 'r') as f:
        for line in f:
            if len(line.split()) > 0:
                accession = line.strip()
                orgs.append(accession)

    with mp.Pool(processes=1) as pool:
        res = pool.map_async(download,orgs).get()
