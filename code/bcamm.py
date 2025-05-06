#! /usr/bin/env python3

from sys import argv
from subprocesses import blast,clasp
import pathlib
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def blast_clasp_anchor_map_making(orgs_tuple,file,ws,suf=''):
    org1 = orgs_tuple.split('---')[0]
    org2 = orgs_tuple.split('---')[1]
    # blast forward
    ret1 = blast(f'{work_dir}/blastdbs/anchor_candidates_'+org1+'_forward',f'{work_dir}/sequences_to_compare/{org2}/forward_split{suf}/forward.{file}',f'{work_dir}/blast_out_forward/{orgs_tuple}/{file}',ws,'bcamm')
    # blast reverse 
    ret2 = blast(f'{work_dir}/blastdbs/anchor_candidates_'+org1+'_forward',f'{work_dir}/sequences_to_compare/{org2}/reverse_split{suf}/reverse.{file}',f'{work_dir}/blast_out_reverse/{orgs_tuple}/{file}',ws,'bcamm')
    if ret1 != 0 or ret2 != 0:
        return(1)
    # clasp forward
    ret1 = clasp(f'{work_dir}/blast_out_forward/{orgs_tuple}/{file}',f'{work_dir}/clasp_out_forward/{orgs_tuple}/{file}','bcamm',2,0.1)
    # clasp reverse
    ret2 = clasp(f'{work_dir}/blast_out_reverse/{orgs_tuple}/{file}',f'{work_dir}/clasp_out_reverse/{orgs_tuple}/{file}','bcamm',2,0.1)
    if ret1 != 0 or ret2 != 0:
        return(1)
    return(0)

if __name__ == "__main__":
    
    work_dir = argv[-1]
    root = str(pathlib.Path(__file__).parents[1])
    
    ret = blast_clasp_anchor_map_making(argv[1],argv[2],argv[3])
    orgs_tuple = argv[1]
    file = argv[2]
    org1 = orgs_tuple.split('---')[0]
    org2 = orgs_tuple.split('---')[1]
    if ret != 0:
        new_files = []
        counter = 1
        record_ids = {}
        for record in SeqIO.parse(f'{work_dir}/sequences_to_compare/{org2}/forward_split/forward.{file}', "fasta"):
            SeqIO.write(record, f"{work_dir}/sequences_to_compare/{org2}/forward_split_add/forward.add-{counter}-{file}", "fasta")
            new_files.append(f'add-{counter}-{file}')
            record_ids[record.id] = counter
            counter += 1
        for record in SeqIO.parse(f'{work_dir}/sequences_to_compare/{org2}/reverse_split/reverse.{file}', "fasta"):
            name = record_ids[record.id]
            SeqIO.write(record, f"{work_dir}/sequences_to_compare/{org2}/reverse_split_add/reverse.add-{name}-{file}", "fasta")
        for f in new_files:
            ret = blast_clasp_anchor_map_making(argv[1],f,argv[3])
            if ret != 0:
                print(f'FINAL FAIL because of blast/clasp time out or MEM limit ({ret}): orig file {argv[2]} single file {f}')
