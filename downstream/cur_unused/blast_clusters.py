#! /usr/bin/env python3

from Bio import SeqIO
import pickle
from subprocess import run
import multiprocessing as mp

def blast_clasp_short(c_id):
    run(f'blastn -query seqs/{c_id}.fasta -subject seqs/{c_id}.fasta -outfmt 6 -word_size 4 -out blast_out/{c_id}_both.txt',shell=True)
    forward = []
    reverse = []
    with open(f'blast_out/{c_id}_both.txt','r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            if line[0] == line[1]:
                continue
            start2 = int(line[8])
            end2 = int(line[9])
            if start2 > end2:
                line[8],line[9] = line[9],line[8]
                reverse.append('\t'.join(line))
            else:
                forward.append('\t'.join(line))
    with open(f'blast_out/{c_id}_forward.txt','w') as f:
        f.write('\n'.join(forward))
    with open(f'blast_out/{c_id}_reverse.txt','w') as f:
        f.write('\n'.join(reverse))
    
    run(f'clasp.x -m -i blast_out/{c_id}_forward.txt -c 7 8 9 10 12 -C 1 2 -l 0.5 -e 0 -o clasp_out/{c_id}_forward',shell=True)
    run(f'clasp.x -m -i blast_out/{c_id}_reverse.txt -c 7 8 9 10 12 -C 1 2 -l 0.5 -e 0 -o clasp_out/{c_id}_reverse',shell=True)

    best_lines = {}
    best_scores = {}
    best_lens = {}
    with open(f'clasp_out/{c_id}_forward','r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            seq1 = line[1]
            seq2 = line[2]
            if seq1 not in best_lines:
                best_lines[seq1] = {}
                best_scores[seq1] = {}
                best_lens[seq1] = {}
            if seq2 not in best_scores[seq1]:
                best_scores[seq1][seq2] = 0
            if seq2 not in best_lines:
                best_lines[seq2] = {}
                best_scores[seq2] = {}
                best_lens[seq2] = {}
            if seq1 not in best_scores[seq2]:
                best_scores[seq2][seq1] = 0
            score = float(line[7])
            start1 = int(line[3])
            end1 = int(line[4])
            start2 = int(line[5])
            end2 = int(line[6])
            if score > best_scores[seq1][seq2]:
                best_scores[seq1][seq2] = score
                best_lines[seq1][seq2] = line
                best_lens[seq1][seq2] = end1-start1

            if score > best_scores[seq2][seq1]:
                best_scores[seq2][seq1] = score
                best_lines[seq2][seq1] = line
                best_lens[seq2][seq1] = end2-start2

    with open(f'clasp_out/{c_id}_reverse','r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split()
            seq1 = line[1]
            seq2 = line[2]
            if seq1 not in best_lines:
                best_lines[seq1] = {}
                best_scores[seq1] = {}
                best_lens[seq1] = {}
            if seq2 not in best_scores[seq1]:
                best_scores[seq1][seq2] = 0
            if seq2 not in best_lines:
                best_lines[seq2] = {}
                best_scores[seq2] = {}
                best_lens[seq2] = {}
            if seq1 not in best_scores[seq2]:
                best_scores[seq2][seq1] = 0
            score = float(line[7])
            start1 = int(line[3])
            end1 = int(line[4])
            start2 = int(line[5])
            end2 = int(line[6])
            if score > best_scores[seq1][seq2]:
                best_scores[seq1][seq2] = score
                best_lines[seq1][seq2] = line
                best_lens[seq1][seq2] = end1-start1

            if score > best_scores[seq2][seq1]:
                best_scores[seq2][seq1] = score
                best_lines[seq2][seq1] = line
                best_lens[seq2][seq1] = end2-start2
    return best_scores,best_lens,best_lines

def check_cluster(c_id):
    print(c_id)
    best_scores,best_lens,best_lines = blast_clasp_short(c_id)
    res = {'best_scores':best_scores,'best_lens':best_lens,'best_lines':best_lines}
    return(c_id,res)

if __name__ == "__main__":
    
    with open('clusters','rb') as f:
        clusters = pickle.load(f)

    orgs = []
    with open('orgs','r') as f:
        for line in f:
            orgs.append(line.strip())
    #for c_id in clusters:
    #    new = check_cluster(c_id)
    with mp.Pool(processes=24) as pool:
        res = pool.map_async(check_cluster,list(clusters.keys())).get()
    collect = {} 
    for c_id,new in res:
        collect[c_id] = new
    for c_id,new in collect.items():
        clusters[c_id]['region_blast_clasp_results'] = new
    with open('clusters_with_region_blast_clasp_results','wb') as f:
        pickle.dump(clusters,f)
