#! /usr/bin/env python3

from subprocess import run,Popen,TimeoutExpired,PIPE
from os.path import isfile
from os import getcwd
import os
from Bio import SeqIO
import pickle
import resource

def setlimits():

    #if suf == '' or suf == '_add':
    #    resource.setrlimit(rsrc, (16000000000, hard))
    #else:
    #    resource.setrlimit(rsrc, (32000000000, hard))
    MAX_VIRTUAL_MEMORY = 16 * 1024 * 1024 * 1024 # 16 MB
    resource.setrlimit(resource.RLIMIT_AS, (MAX_VIRTUAL_MEMORY, resource.RLIM_INFINITY))

    return 0

def blast(db,query,outfile,word_size,suf):
    if suf == '':
        TO = 1200
    elif suf == '_add':
        TO = 420
    elif suf == 'bcamm':
        TO = 604800
    cmd = 'blastn -db {} -query {} -strand plus -outfmt 6 -evalue 1e-3 -word_size {} -out {}'.format(db,query,word_size,outfile).split()
    #run(cmd)
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, preexec_fn=setlimits)
    try:
        outs, errs = proc.communicate(timeout=TO)
    except TimeoutExpired:
        print(f'blast {query} TimeoutExpired')
        proc.kill()
        outs, errs = proc.communicate()
        return(1)
    if proc.returncode != 0:
        print(f'blast {query} probably process killed my MEM limit')
        return(-1)
    return(0)

def clasp(infile,outfile,suf,l=0.5,e=0):
    if suf == '':
        TO = 900
    elif suf == '_add':
        TO = 300
    elif suf == 'bcamm':
        TO = 604800
    bin_dir = '../bin'
    cmd = '{}/clasp.x -m -i {} -c 7 8 9 10 12 -C 1 2 -l {} -e {} -o {}'.format(bin_dir,infile,l,e,outfile).split()
    print(cmd)
    #run(cmd)
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, preexec_fn=setlimits)
    try:
        outs, errs = proc.communicate(timeout=TO)
    except TimeoutExpired:
        print(f'clasp {infile} TimeoutExpired')
        proc.kill()
        outs, errs = proc.communicate()
        return(1)
    if proc.returncode != 0:
        print(f'clasp {infile} probably process killed my MEM limit')
        return(-1)
    return(0)
