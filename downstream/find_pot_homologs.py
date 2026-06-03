#! /usr/bin/env python3

from sys import argv
from subprocess import run
import os
import multiprocessing as mp

def tblastn(db,query,outfile):
    cmd = 'tblastn -db {} -query {} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -evalue 1e-3 -word_size 7 -out {}'.format(db,query,outfile)
    run(cmd,shell=True)
    return(0)

def blastn(db,query,outfile):
    cmd = 'blastn -db {} -query {} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -evalue 1e-3 -word_size 9 -out {}'.format(db,query,outfile)
    run(cmd,shell=True)
    return(0)

def clasp(infile,outfile,l=0.5,e=0):
    cmd = './clasp.x -m -i {} -c 7 8 9 10 12 -C 1 2 -l {} -e {} -o {}'.format(infile,l,e,outfile).split()
    run(cmd)
    return(0)

def exec_blast_clasp(org):
    
    if argv[3] == 'prot':
        tblastn(f'blastdbs/{org}',f'{protein}',f'blast_out_both/{org}')
    else:
        blastn(f'blastdbs/{org}',f'{protein}',f'blast_out_both/{org}')

    
    with open(f'blast_out_both/{org}','r') as f:
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

    with open(f'blast_out_forward/{org}','w') as f:
        f.write('\n'.join(newlines_forward))
    with open(f'blast_out_reverse/{org}','w') as f:
        f.write('\n'.join(newlines_reverse))

    for orientation in ['forward','reverse']:
    
        clasp(f'blast_out_{orientation}/{org}',f'clasp_out_{orientation}/{org}')

    hit_no = 1
    with open(f'temp_coords_{org}','w') as out:
        for orientation in ['forward','reverse']:
            with open(f'clasp_out_{orientation}/{org}') as f:
                for line in f:
                    if line[0] == '#':
                        continue
                    line = line.split()
                    orig_prot = line[1]
                    chromo = line[2]
                    start = line[5]
                    end = line[6]
                    out.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\thit-{hit_no}_{org}_{orig_prot}\n')
                    hit_no += 1

    return(0)

def mkbdb(org):
    if not os.path.isfile(path + f"/blastdbs/{org}.nsq"):
        run(f'makeblastdb -dbtype nucl -in {genomes}/{org}.fasta -out blastdbs/{org}',shell=True)

if __name__ == "__main__":
    
    protein = argv[1]
    genomes = argv[2]
    path = os.getcwd()

    run('mkdir -p blastdbs',shell=True)   
    run('mkdir -p blast_out_both',shell=True)   
    run('mkdir -p blast_out_forward',shell=True)   
    run('mkdir -p blast_out_reverse',shell=True)   
    run('mkdir -p clasp_out_forward',shell=True)   
    run('mkdir -p clasp_out_reverse',shell=True)   
    orgs = []
    with open('orgs') as f:
        for line in f:
            org = line.strip()
            orgs.append(org)

    with mp.Pool(processes=11) as pool:
        res = pool.map_async(mkbdb,orgs).get()
            
    with mp.Pool(processes=11) as pool:
        res = pool.map_async(exec_blast_clasp,orgs).get()
