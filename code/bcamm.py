#! /usr/bin/env python3

from sys import argv
from subprocesses import blast,clasp,get_mem_limit_bytes,run_blast_clasp_with_retry
import pathlib
import os

def blast_clasp_anchor_map_making(orgs_tuple,file,ws):
    org1 = orgs_tuple.split('---')[0]
    org2 = orgs_tuple.split('---')[1]
    mem_limit = get_mem_limit_bytes()
    db = f'{work_dir}/blastdbs/anchor_candidates_{org1}_forward'
    for orient in ['forward','reverse']:
        other = 'reverse' if orient == 'forward' else 'forward'
        run_blast_clasp_with_retry(
            db,
            f'{work_dir}/sequences_to_compare/{org2}/{orient}_split/{orient}.{file}',
            f'{work_dir}/blast_out_{orient}/{orgs_tuple}/{file}',
            f'{work_dir}/clasp_out_{orient}/{orgs_tuple}/{file}',
            ws, mem_limit, clasp_mode='pairwise',
            rev_fasta_path=f'{work_dir}/sequences_to_compare/{org2}/{other}_split/{other}.{file}',
            orientation=orient)
    # check for one-sided failures (valid for pairwise, but log for visibility)
    for orient in ['forward','reverse']:
        cf = f'{work_dir}/clasp_out_{orient}/{orgs_tuple}/{file}'
        if os.path.exists(cf):
            with open(cf) as f:
                if not any(not l.startswith('#') for l in f):
                    print(f'[WARN] {orient} clasp for {orgs_tuple}/{file} has no data lines (all-fail)')
        else:
            print(f'[ERROR] {orient} clasp output missing: {cf}')
    return(0)

if __name__ == "__main__":
    
    work_dir = argv[-1]
    root = str(pathlib.Path(__file__).parents[1])
    
    blast_clasp_anchor_map_making(argv[1],argv[2],argv[3])
