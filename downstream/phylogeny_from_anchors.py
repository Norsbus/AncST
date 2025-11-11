# Importing necessary libraries
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceMatrix
import numpy as np


def get_matrix_from_anchors():

    orgs = []
    with open('orgs') as f:
        for line in f:
            orgs.append(line.strip())

    ams = {}
    for org in orgs:
        path = 'compressed_maps_multis_to_one/{}'.format(org)
        with open(path, 'rb') as f:
            ams[org] = pickle.load(f)

    lens = {}
    for org in orgs:
        path = '../utils/small_meta/{}'.format(org)
        with open(path, 'rb') as f:
            x = pickle.load(f)
            lens[org] = x[1][-1]

    matrix = [[]]
    for enu, org in enumerate(orgs):
        if enu == 0:
            continue
        loc_dists = [0 for i in range(enu)]
        for i in range(0, enu):
            org2 = orgs[i]
            if lens[org] < lens[org2]:
                dist_org = org
                other = org2
            else:
                dist_org = org2
                other = org
            for j, bib in ams[dist_org].items():
                if other not in bib['matches']:
                    continue
                else:
                    for k, t in bib['matches'][other].items():
                        loc_dists[i] += (t[1][1] - t[1][0])
        dists = []
        for d in loc_dists:
            dists.append(1 - (d / lens[dist_org]))
        matrix.append(dists)
    return matrix, orgs

if __name__ == "__main__":

    matrix,orgs = get_matrix_from_anchors()

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
