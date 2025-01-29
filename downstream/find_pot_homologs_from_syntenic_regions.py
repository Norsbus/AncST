#! /usr/bin/env python3

from sys import argv
from subprocess import run
import os
import multiprocessing as mp
from os.path import isfile

def tblastn(subject,query,outfile):
    cmd = 'tblastn -subject {} -query {} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -evalue 1e-3 -word_size 2 -out {}'.format(subject,query,outfile)
    run(cmd,shell=True)
    return(0)

def blastn(subject,query,outfile):
    cmd = 'blastn -subject {} -query {} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -evalue 1e-3 -word_size 9 -out {}'.format(subject,query,outfile)
    run(cmd,shell=True)
    return(0)

def clasp(infile,outfile,l=0.5,e=0):
    cmd = './clasp.x -m -i {} -c 7 8 9 10 12 -C 1 2 -l {} -e {} -o {}'.format(infile,l,e,outfile).split()
    run(cmd)
    return(0)

def exec_blast_clasp(org):
    
    if argv[2] == 'prot':
        tblastn(f'fastas/{org}/all.fasta',f'{protein}',f'blast_out_both/{org}_regions')
    else:
        blastn(f'fastas/{org}/all.fasta',f'{protein}',f'blast_out_both/{org}_regions')

    
    with open(f'blast_out_both/{org}_regions','r') as f:
        newlines_forward = []
        newlines_reverse = []
        for line in f:
            line = line.split()
            if line[0] == '#' or len(line) == 0:
                continue
            if int(line[8]) > int(line[9]):
                line[8],line[9] = line[9],line[8]
                newlines_reverse.append('\t'.join(line))
            else:
                newlines_forward.append('\t'.join(line))

    with open(f'blast_out_forward/{org}_regions','w') as f:
        f.write('\n'.join(newlines_forward))
    with open(f'blast_out_reverse/{org}_regions','w') as f:
        f.write('\n'.join(newlines_reverse))

    for orientation in ['forward','reverse']:
    
        clasp(f'blast_out_{orientation}/{org}_regions',f'clasp_out_{orientation}/{org}_regions')

    hit_no = 1
    with open(f'temp_coords_{org}_from_regions','w') as out:
        for orientation in ['forward','reverse']:
            with open(f'clasp_out_{orientation}/{org}_regions') as f:
                for line in f:
                    if line[0] == '#':
                        continue
                    line = line.split()
                    orig_prot = line[1]
                    id_region = line[2]
                    split_id = id_region.split('_')
                    chromo = split_id[2]
                    try:
                        shift = int(split_id[3])
                    except:
                        shift = int(split_id[4])
                    start = int(line[5]) + shift
                    end = int(line[6]) + shift
                    out.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\thit-{hit_no}_{org}_{orig_prot}\n')
                    hit_no += 1

    return(0)

if __name__ == "__main__":

    protein = argv[1]
    path = os.getcwd()

    run('mkdir -p blast_out_both',shell=True)   
    run('mkdir -p blast_out_forward',shell=True)   
    run('mkdir -p blast_out_reverse',shell=True)   
    run('mkdir -p clasp_out_forward',shell=True)   
    run('mkdir -p clasp_out_reverse',shell=True)   
    orgs = []
    with open('orgs') as f:
        for line in f:
            org = line.strip()
            if isfile(path + f'/fastas/{org}/all.fasta'):
                orgs.append(org)

    with mp.Pool(processes=11) as pool:
        res = pool.map_async(exec_blast_clasp,orgs).get()
