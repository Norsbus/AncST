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
    
    orgs = []
    with open(f"{work_dir}/orgs","r") as f:
        for line in f:
            orgs.append(line.strip())


    with open(f'{work_dir}/config.json','r') as f:

        config = json.loads(f.read())
        
        if "directories" not in config:
            print('no valid config file')
            pass
        else:
            d=config["directories"]
            remove_directories(d)
            make_directories(d)
            make_directories(["out/utils","out/custom","out/MCScanX","out/dialign","out/Smore"])
