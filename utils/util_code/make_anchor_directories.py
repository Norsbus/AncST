#! /usr/bin/env python3

import os,json,sys,shutil
from subprocess import run
from sys import argv
import pathlib

work_dir = pathlib.Path(__file__).parent.resolve()


def make_directories(directories):

    for directory in directories:
        directory = f'../../anchors/'+directory
        if not os.path.exists(directory):
            os.mkdir(directory)
            print(f'making {directory} dir')
        else:
            pass

    return 0


if __name__ == "__main__":
    
    orgs = []
    with open(f"orgs","r") as f:
        for line in f:
            orgs.append(line.strip())

    with open(f'{work_dir}/config.json','r') as f:

        config = json.loads(f.read())
        
        if "directories" not in config:
            pass
        else:
            ds=config["directories"]
            to_make = []
            for d in ds:
                to_make.append(f'{d}')
            make_directories(to_make)
            per_org = ['compared']
            to_make = []
            for d in per_org:
                for org in orgs:
                    to_make.append(f'{d}/{org}')
            make_directories(to_make)
