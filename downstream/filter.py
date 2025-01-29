#! /usr/bin/env python3
from sys import argv
out = open('edit.homology','w')
with open(argv[1]) as f:
    for line in f:
        ele1,ele2,score = line.split()
        org1 = ele1.split('chr')[0]
        org2 = ele2.split('chr')[0]
        if org1 == org2:
            continue
        else:
            out.write(line)
out.close()
