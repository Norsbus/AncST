#! /usr/bin/env python3

import os
from subprocess import run
import math

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())

genome_sizes = {}
automatic_k_e = {}
for org in orgs:
    genome_sizes[org] = os.path.getsize(f'/scr/k80san/karl/scaffolder/AncST/utils/genomes/{org}.fasta')
    k = math.ceil(math.log(genome_sizes[org],4))
    if k - math.log(genome_sizes[org],4) > 0.5:
        k -= 1
    automatic_k_e[org] = (k,0)

run('rm -f /scr/k80san/karl/scaffolder/AncST/template/genmap_params.txt && touch /scr/k80san/karl/scaffolder/AncST/template/genmap_params.txt',shell=True)
run('rm -f /scr/k80san/karl/scaffolder/AncST/template/macle_params.txt && touch /scr/k80san/karl/scaffolder/AncST/template/macle_params.txt',shell=True)
run('rm -f /scr/k80san/karl/scaffolder/AncST/template/dups_params.txt && touch /scr/k80san/karl/scaffolder/AncST/template/dups_params.txt',shell=True)
for org in orgs:
    run(f'echo "{org} {automatic_k_e[org][0]} {automatic_k_e[org][1]} 300 50 100 1" >> /scr/k80san/karl/scaffolder/AncST/template/genmap_params.txt',shell=True)
    run(f'echo "{org} {automatic_k_e[org][0]} {automatic_k_e[org][1]} 300 50 66 0" >> /scr/k80san/karl/scaffolder/AncST/template/genmap_params.txt',shell=True)
    run(f'echo "{org} 1000 50 100 1" >> /scr/k80san/karl/scaffolder/AncST/template/macle_params.txt',shell=True)
    run(f'echo "{org} 1000 50 33 0" >> /scr/k80san/karl/scaffolder/AncST/template/macle_params.txt',shell=True)
    run(f'echo "{org} 50 4 21 2 300 50 16 90 10 10" >> /scr/k80san/karl/scaffolder/AncST/template/dups_params.txt',shell=True)

