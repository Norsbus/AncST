#! /usr/bin/env python3

from sys import argv
from subprocesses import blast,clasp
import pathlib

def blast_clasp_anchor_map_making(orgs_tuple,file,ws):
    org1 = orgs_tuple.split('---')[0]
    org2 = orgs_tuple.split('---')[1]
    # blast forward
    blast(f'{work_dir}/blastdbs/anchor_candidates_'+orgs_tuple+'_forward',f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/forward_split/forward.{file}',f'{work_dir}/blast_out_forward/{orgs_tuple}/{file}',ws)
    # blast reverse 
    blast(f'{work_dir}/blastdbs/anchor_candidates_'+orgs_tuple+'_forward',f'{work_dir}/sequences_to_compare/{org1}---{org2}/{org2}/reverse_split/reverse.{file}',f'{work_dir}/blast_out_reverse/{orgs_tuple}/{file}',ws)
    # clasp forward
    clasp(f'{work_dir}/blast_out_forward/{orgs_tuple}/{file}',f'{work_dir}/clasp_out_forward/{orgs_tuple}/{file}',2,0.1)
    # clasp reverse
    clasp(f'{work_dir}/blast_out_reverse/{orgs_tuple}/{file}',f'{work_dir}/clasp_out_reverse/{orgs_tuple}/{file}',2,0.1)

    return(0)

if __name__ == "__main__":
    
    work_dir = argv[-1]
    root = str(pathlib.Path(__file__).parents[1])
    
    blast_clasp_anchor_map_making(argv[1],argv[2],argv[3])
