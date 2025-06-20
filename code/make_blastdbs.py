#! /usr/bin/env python3

from subprocess import run

def makeblastdb(org):
    run("makeblastdb -in OTHER/scr/k80san/karl/genomes/{org}.fasta -out blastdbs/{org} -dbtype nucl")

