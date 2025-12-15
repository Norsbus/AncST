#! /usr/bin/env python3

from Bio import SeqIO
import pickle,os
from subprocess import run
import multiprocessing as mp

def get_genome(genomes_dir,org):
    
    seqs_without_margin = SeqIO.parse(f'{genomes_dir}/{org}.fasta', "fasta")
    s = {}
    for seq in seqs_without_margin:
        s[seq.id] = seq

    return(s)

def check_cluster(c_id):
    print(c_id)
    c = clusters[c_id]
    new = {}
    ranges = {}
    for org,bib in c.items():
        if org not in orgs:
            continue
        ranges[org] = [min(bib['matches']),max(bib['matches'])]
    seqs_without_margin = []
    for org,r in ranges.items():
        if r[0] == r[1]:
            chromo = aligned[org][r[0]]['chromosome']
            start = aligned[org][r[0]]['start']
            end = aligned[org][r[0]]['end']
        else:
            chromo = aligned[org][r[0]]['chromosome']
            start = aligned[org][r[0]]['start']
            end = aligned[org][r[1]]['end']
        seq = genomes[org][chromo][int(start):int(end)]
        seq.id = f'{org} {chromo} {start} {end} forward'
        seq.description = ''
        seqs_without_margin.append(seq)
    SeqIO.write(seqs_without_margin,f'seqs/{c_id}.fasta','fasta')
    return c_id,new

if __name__ == "__main__":
   

    genomes = {}
    aligned = {}
    with open('clusters','rb') as f:
        clusters = pickle.load(f)

    orgs = []
    with open('orgs','r') as f:
        for line in f:
            orgs.append(line.strip())
    for org in orgs:
        genomes[org] = get_genome('/scr/k80san/karl//methodenpaper/test_pipeline/stable_synteny/utils/genomes/',org)
        with open(f'anchors_from_methodenpaper_test_pipeline/aligned/{org}','rb') as f2:
            aligned[org] = pickle.load(f2)
    ids = []
    path = os.getcwd()
    for c_id in clusters:
        if os.path.isfile(path+f'/seqs/{c_id}.fasta'):
            continue
        else:
            ids.append(c_id)
    with mp.Pool(processes=24) as pool:
        res = pool.map_async(check_cluster,ids).get()
