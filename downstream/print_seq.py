#! /usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from sys import argv
import pickle

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())
genomes = {}
aligned = {}
def get_genome(genomes_dir,org):
    
    seqs_without_margin = SeqIO.parse(f'{genomes_dir}/{org}.fasta', "fasta")
    s = {}
    for seq in seqs_without_margin:
        s[seq.id] = seq

    return(s)
for org in orgs:
    genomes[org] = get_genome('HOMEDIR/methodenpaper/test_pipeline/stable_synteny/utils/genomes/',org)
    with open(f'anchors_from_methodenpaper_test_pipeline/aligned/{org}','rb') as f2:
        aligned[org] = pickle.load(f2)

with open('clusters','rb') as f:
    clusters = pickle.load(f)
for c_id,clu in clusters.items():
    out=open('coords_test_dialign_ungapped','w')
    print(c_id)
    seqs = []
    for org,bib in clu.items():
        if org not in orgs:
            continue
        se = [min(bib['matches']),max(bib['matches'])]
        chromo = aligned[org][se[0]]['chromosome']
        start = aligned[org][se[0]]['start']
        end = aligned[org][se[1]]['end']
        ori = bib['orientation']
        out.write(f'{org}\t{chromo}\t{start-10000}\t{end+10000}\t{ori}\n')
        seqs.append(SeqRecord(genomes[org][chromo][int(start)-10000:int(end)+10000].seq,id='unkown'))
    #SeqIO.write(seqs,f'{c_id}', "fasta")
    out.close()
    exit(0)
    input()
