#!/usr/bin/env python3

from subprocess import run,check_output
import os
from time import sleep
import pathlib
from sys import argv

anchor_dir = argv[1]

with open('orgs_all_matches','r') as f:
    orgs = f.read().splitlines()

skip_some = False

for org in orgs:
    print(org)
    cmd = f"./compress_maps_and_ignore_multis.py {org} {anchor_dir} && touch touch/compressed_{org}_done"
    run(cmd.split())
