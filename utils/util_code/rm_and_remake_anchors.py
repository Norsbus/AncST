#! /usr/bin/env python3

from subprocess import run

run('rm -r ../../anchors', shell=True)
run('mkdir ../../anchors', shell=True)
run('./make_anchor_directories.py', shell=True)
