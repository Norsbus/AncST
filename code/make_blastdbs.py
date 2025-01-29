#! /usr/bin/env python3

from subprocess import run

def makeblastdb(org):
    run("makeblastdb -in OTHERHOMEDIRgenomes/{org}.fasta -out blastdbs/{org} -dbtype nucl")

