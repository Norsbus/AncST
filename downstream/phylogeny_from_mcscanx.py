#! /usr/bin/env python3

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceMatrix
import numpy as np
import pickle
import sys
sys.path.append('./utils/')
from get_mapping import get_mapping

def get_scores():
    scores = {}
    with open('MCScanX.homology') as f:
        for line in f:
            ele1,ele2,score = line.strip().split('\t')
            score = int(float(score))
            scores[(ele1,ele2)] = score
            scores[(ele2,ele1)] = score
    return scores

def get_matrix_from_mcscanx_match_counts():

    org_mapping,chr_mapping = get_mapping()

    orgs = []
    with open('orgs') as f:
        for line in f:
            orgs.append(line.strip())

    counter = {}

    homology = get_scores()
    with open(f'MCScanX.collinearity','r') as f:
        for line in f:
            if 'Alignment' in line or '#' in line:
                continue
            else:
                line = line.strip().split(':')[1].split()
                first = line[0].strip()
                second = line[1].strip()
                org1 = first.split('chr')[0]
                org2 = second.split('chr')[0]
                if (org1,org2) not in counter:
                    counter[(org1,org2)] = 0
                if (org2,org1) not in counter:
                    counter[(org2,org1)] = 0
                if org1 not in counter:
                    counter[org1] = set()
                if org2 not in counter:
                    counter[org2] = set()
                counter[(org1,org2)] += homology[(first,second)]
                counter[(org2,org1)] += homology[(first,second)]
                counter[org1].add((first,second,homology[(first,second)]))
                counter[org2].add((second,first,homology[(first,second)]))

    score_counter = {}
    for org,s in counter.items():
        if org[-3:] != 'org':
            continue
        else:
            score_counter[org] = sum([x[2] for x in s])

    matrix = [[]]
    for enu,org in enumerate(orgs):
        if enu == 0:
            continue
        loc_dists = [0 for i in range(enu)]
        for i in range(0,enu):
            loc_dists[i] = 1-(counter[(org_mapping[orgs[enu]],org_mapping[orgs[i]])]/((score_counter[org_mapping[orgs[enu]]]+score_counter[org_mapping[orgs[i]]])/2))
        matrix.append(loc_dists)

    return matrix,orgs

if __name__ == "__main__":

    matrix,orgs = get_matrix_from_mcscanx_match_counts()

    for i in range(len(matrix)):
        matrix[i].append(0)

    distMatrix = DistanceMatrix(orgs, matrix)

    constructor = DistanceTreeConstructor()

    UPGMATree = constructor.upgma(distMatrix)

    #Phylo.draw(UPGMATree)

    #Phylo.draw_ascii(UPGMATree)

    NJTree = constructor.nj(distMatrix)

    #Phylo.draw(NJTree)

    #Phylo.draw_ascii(NJTree)

    from Bio import Phylo

    Phylo.write(UPGMATree, "UPGMATree.nwk", "newick")
    Phylo.write(NJTree, "NJTree.nwk", "newick")
