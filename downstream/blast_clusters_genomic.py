#! /usr/bin/env python3

from Bio import SeqIO
import pickle,os
from subprocess import run
import multiprocessing as mp
from pprint import pprint

def check_cluster(c_id):
    print(c_id)
    orgs = [org for org in list(clusters[c_id].keys()) if org in orgs_glob]
    #orgs.remove('representative_score')
    #orgs.remove('representative_length')
    #orgs.remove('focal_orientation')
    res = {}
    for org in orgs:
        mini = 1e15
        maxi = 0
        for i in clusters[c_id][org]['matches']:
            if i < mini:
                mini = i
            if i > maxi:
                maxi = i
        anchor_region_chromo = compressed[org][mini]['chromosome']
        anchor_region_start = compressed[org][mini]['start']
        anchor_region_end = compressed[org][maxi]['end']
        #if not os.path.isfile(path+f'/blast_out/{c_id}_{org}_both.txt') or not os.path.isfile(path+f'/blast_out/{c_id}_{org}_forward.txt') or not os.path.isfile(path+f'/blast_out/{c_id}_{org}_reverse.txt'):
        #    run(f'blastn -query seqs/{c_id}.fasta -db {db_path}/{org} -outfmt 6 -word_size 11 -evalue 1e-3 -out blast_out/{c_id}_{org}_both.txt',shell=True)
        #    forward = []
        #    reverse = []
        #    with open(f'blast_out/{c_id}_{org}_both.txt','r') as f:
        #        for line in f:
        #            if line.startswith('#'):
        #                continue
        #            line = line.strip().split()
        #            if line[0] == org:
        #                continue
        #            start2 = int(line[8])
        #            end2 = int(line[9])
        #            if start2 > end2:
        #                line[8],line[9] = line[9],line[8]
        #                reverse.append('\t'.join(line))
        #            else:
        #                forward.append('\t'.join(line))
        #    with open(f'blast_out/{c_id}_{org}_forward.txt','w') as f:
        #        f.write('\n'.join(forward))
        #    with open(f'blast_out/{c_id}_{org}_reverse.txt','w') as f:
        #        f.write('\n'.join(reverse))
        #else:
        #    print(f'{c_id}: blast already done')
        #
        #if not os.path.isfile(path+f'/clasp_out/{c_id}_{org}_forward'):
        #    run(f'clasp.x -m -i blast_out/{c_id}_{org}_forward.txt -c 7 8 9 10 12 -C 1 2 -l 0.5 -e 0 -o clasp_out/{c_id}_{org}_forward',shell=True)
        #else:
        #    print(f'{c_id}: clasp_for already done')
        #if not os.path.isfile(path+f'/clasp_out/{c_id}_{org}_reverse'):
        #    run(f'clasp.x -m -i blast_out/{c_id}_{org}_reverse.txt -c 7 8 9 10 12 -C 1 2 -l 0.5 -e 0 -o clasp_out/{c_id}_{org}_reverse',shell=True)
        #else:
        #    print(f'{c_id}: clasp_rev already done')

        best_inside = {}
        best_outside = {}
        for org2 in orgs:
            if org2 == org:
                continue
            best_outside[org2] = [0,'']
            best_inside[org2] = [0,'']

        with open(f'clasp_out/{c_id}_{org}_forward','r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split()
                try:
                    seq1 = line[1]
                    seq2 = line[2]
                    score = float(line[7])
                    start1 = int(line[3])
                    end1 = int(line[4])
                    start2 = int(line[5])
                    end2 = int(line[6])
                except:
                    out.write(f'{c_id}_{org}_forward\n')
                    out.write(f'{line}\n')
                    out.write('--------\n')
                    continue
                if seq2 == anchor_region_chromo and ((start2 >= anchor_region_start-margin and start2 <= anchor_region_end+margin) or (end2 >= anchor_region_start-margin and end2 <= anchor_region_end+margin) or (start2 < anchor_region_start and end2 > anchor_region_end)):
                    if score > best_inside[seq1][0]:
                        best_inside[seq1][0] = score
                        best_inside[seq1][1] = line
                else:
                    if score > best_outside[seq1][0]:
                        best_outside[seq1][0] = score
                        best_outside[seq1][1] = line


        with open(f'clasp_out/{c_id}_{org}_reverse','r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split()
                try:
                    seq1 = line[1]
                    seq2 = line[2]
                    score = float(line[7])
                    start1 = int(line[3])
                    end1 = int(line[4])
                    start2 = int(line[5])
                    end2 = int(line[6])
                except:
                    out.write(f'{c_id}_{org}_reverse\n')
                    out.write(f'{line}\n')
                    out.write('--------\n')
                    continue
                if seq2 == anchor_region_chromo and ((start2 >= anchor_region_start-margin and start2 <= anchor_region_end+margin) or (end2 >= anchor_region_start-margin and end2 <= anchor_region_end+margin) or (start2 < anchor_region_start and end2 > anchor_region_end)):
                    if score > best_inside[seq1][0]:
                        best_inside[seq1][0] = score
                        best_inside[seq1][1] = line
                else:
                    if score > best_outside[seq1][0]:
                        best_outside[seq1][0] = score
                        best_outside[seq1][1] = line
        res[org] = best_outside,best_inside

    return(c_id,res)

if __name__ == "__main__":

    db_path = 'HOMEDIR/methodenpaper/test_pipeline/stable_synteny/utils/blastdbs/'
    
    with open('clusters_with_region_blast_clasp_results','rb') as f:
        clusters = pickle.load(f)

    orgs_glob = []
    with open('orgs','r') as f:
        for line in f:
            orgs_glob.append(line.strip())
    
    compressed = {}
    for org in orgs_glob:
        with open('compressed_maps/'+org,'rb') as f:
            compressed[org] = pickle.load(f)
    path = os.getcwd()
    margin = 20000
    out = open('debug_bcg.txt','w')

    #for c_id in clusters:
    #    new = check_cluster(c_id)
    with mp.Pool(processes=16) as pool:
        res = pool.map_async(check_cluster,list(clusters.keys())).get()
        #res = pool.map_async(check_cluster,[1]).get()
    collect = {} 
    for c_id,new in res:
        collect[c_id] = new
    for c_id,new in collect.items():
        clusters[c_id]['genomic_blast'] = new
    with open('clusters_with_region_blast_clasp_results_and_genomic_blast_lenient_region_coords_20000','wb') as f:
        pickle.dump(clusters,f)
    out.close()
