#! /usr/bin/env python3

import os,json,sys,shutil
from subprocess import run
from sys import argv
import pathlib

work_dir = pathlib.Path(__file__).parent.resolve()


def make_directories(directories):

    for directory in directories:
        directory = f'{work_dir}/'+directory
        if not os.path.exists(directory):
            os.mkdir(directory)

    return 0


def remove_directories(directories):
    
    for directory in directories:
        directory = f'{work_dir}/'+directory
        try:
            shutil.rmtree(directory)
        except FileNotFoundError as e:
            pass
            #logger.info(e)

    return 0

if __name__ == "__main__":
    
    try:
        remove_directories(['logs'])
    except:
        print('cannot delete old logs since they id not exist...making now')

    make_directories(['logs'])
    
    orgs = []
    with open(f"{work_dir}/orgs","r") as f:
        for line in f:
            orgs.append(line.strip())


    with open(f'{work_dir}/config.json','r') as f:

        config = json.loads(f.read())
        
        if "directories" not in config:
            pass
            #logger.error('directories not specified in config.json')
            #sys.exit()
        else:
            d=config["directories"]
            remove_directories(d)
            make_directories(d)
            per_org = ['split_out_forward','split_out_reverse','blast_out_forward','blast_out_reverse','clasp_out_forward','clasp_out_reverse','indices','to_del','inconsistencies','parse_bcamm','mark_in_others','checks_log','statistics_log']
            to_make = []
            for org in orgs:
                for d in per_org:
                    to_make.append(f'{d}/{org}')
            remove_directories(to_make)
            make_directories(to_make)

        orgs_combo = []
        
        for i,org1 in enumerate(orgs):
            orgs_combo += [(org1,org2) for org2 in orgs[i+1:]]
        
        to_make = []
        for d in ['blast_out_forward','blast_out_reverse','clasp_out_forward','clasp_out_reverse','sequences_to_compare']:
            for org1,org2 in orgs_combo:
                to_make.append(f'{d}/{org1}---{org2}')
        remove_directories(to_make)
        make_directories(to_make)

        to_make = []
        for d in ['sequences_to_compare']:
            for org in orgs:
                to_make.append(f'{d}/{org}')
            #for org1,org2 in orgs_combo:
            #    to_make.append(f'{d}/{org1}---{org2}/{org1}')
            #    to_make.append(f'{d}/{org1}---{org2}/{org2}')
        remove_directories(to_make)
        make_directories(to_make)
        to_make = []
        for d in ['sequences_to_compare']:
            for org in orgs:
                to_make.append(f'{d}/{org}/forward_split')
                to_make.append(f'{d}/{org}/reverse_split')
            #for org1,org2 in orgs_combo:
            #    to_make.append(f'{d}/{org1}---{org2}/{org1}/forward_split')
            #    to_make.append(f'{d}/{org1}---{org2}/{org1}/reverse_split')
            #    to_make.append(f'{d}/{org1}---{org2}/{org2}/forward_split')
            #    to_make.append(f'{d}/{org1}---{org2}/{org2}/reverse_split')
        remove_directories(to_make)
        make_directories(to_make)
