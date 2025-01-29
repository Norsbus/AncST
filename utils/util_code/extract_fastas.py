#! /usr/bin/env python3

from subprocess import run

with open('accesssions','r') as f:
    for line in f:
        if len(line.split()) > 0:
            accession = line.strip()
        run(f'unzip ../genomes/{accession}.zip -d ../genomes/{accession}_unzipped && rm ../genomes/{accession}.zip', shell=True)
        run(f'mv ../genomes/{accession}_unzipped/ncbi_dataset/data/{accession}/*.fna ../genomes/{accession}.fasta && rm -r ../genomes/{accession}_unzipped', shell=True)
