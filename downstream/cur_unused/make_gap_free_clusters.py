#! /usr/bin/env python3

from Bio import SeqIO
import pickle,os
from subprocess import run
import multiprocessing as mp
import sys
sys.path.append('./out/utils/')
from get_mapping import get_mapping

def ungap(c_id):
    print(c_id)
    iss = {}
    lines = []
    for org,bib in clusters[c_id].items():
        if org not in orgs:
            continue
        iss[org] = min(bib['matches'])
    if not os.path.isfile(path+f'/blast_out_ungapped/{c_id}_both.txt'):
        run(f'blastn -ungapped -query seqs/{c_id}.fasta -subject seqs/{c_id}.fasta -outfmt 6 -word_size 4 -out blast_out_ungapped/{c_id}_both.txt',shell=True)
    else:
        print(f'cluster {c_id} already blasted')
    with open(f'blast_out_ungapped/{c_id}_both.txt','r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            if line[0] == line[1]:
                continue
            org1 = line[0]
            org2 = line[1]
            ione = iss[org1]
            itwo = iss[org2]
            start1 = int(line[6])
            end1 = int(line[7])
            start2 = int(line[8])
            end2 = int(line[9])
            if start2 > end2:
                start2,end2=end2,start2
                rev = 1
            else:
                rev = 0
            
            chromo1 = am[org1][ione]['chromosome']
            chromo1 = chr_mapping[org1][chromo1]
            absstart1 = am[org1][ione]['start']
            absend1 = am[org1][ione]['end']
            ele_name = f'{chromo1}ele{absstart1}to{absend1}'
            chromo2 = am[org2][itwo]['chromosome']
            chromo2 = chr_mapping[org2][chromo2]
            absstart2 = am[org2][itwo]['start']
            absend2 = am[org2][itwo]['end']
            ele2_name = f'{chromo2}ele{absstart2}to{absend2}'
            score = float(line[11])
            hit_start1,hit_end1 = int(start1),int(end1)
            length1 = hit_end1 - hit_start1
            hit_start2,hit_end2 = int(start2),int(end2)
            length2 = hit_end2 - hit_start2
            if rev == 0:
                ori = 'forward'
            else:
                ori = 'reverse'
            length = max(length1,length2)
            
            lines.append(f'{ele_name}\t{chromo1}\t{absstart1}\t{absend1}\t{hit_start1}\t{hit_end1}\t{ele2_name}\t{chromo2}\t{absstart2}\t{absend2}\t{hit_start2}\t{hit_end2}\t{score}\t{length}\t{ori}\n')
    return(lines)


if __name__ == "__main__":
    
    with open('clusters','rb') as f:
        clusters = pickle.load(f)

    org_mapping,chr_mapping = get_mapping()
    path = os.getcwd()
    lines = []
    orgs = []
    with open('orgs','r') as f:
        for line in f:
            orgs.append(line.strip())
    am = {}
    for org in orgs:
        with open(f'anchors_from_methodenpaper_test_pipeline/aligned/{org}','rb') as f2:
            am[org] = pickle.load(f2)
    with mp.Pool(processes=24) as pool:
        res = pool.map_async(ungap,list(clusters.keys())).get()
    with open('dialign_homology','w') as f:
        for lines in res:
            for line in lines:
                f.write(line)
