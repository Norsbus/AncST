#! /usr/bin/env python3

import networkx as nx
import pickle
G=nx.Graph()

orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())
for org in orgs:
    with open(f'REMOTEHOMEDIRsynteny/400_sth/anchors/candidates/{org}','rb') as f:
        am = pickle.load(f)
    iss = [f'{org}_{i}' for i in list(am.keys())]
    G.add_nodes_from(iss)
print('added nodes')
done = {}
for org in orgs:
    print(org)
    done[org] = 1
    with open(f'REMOTEHOMEDIRsynteny/400_sth/anchors/aligned/{org}','rb') as f:
        am = pickle.load(f)
    for i,bib in am.items():
        for org2,mbib in bib['matches'].items():
            if org2 in done:
                continue
            if mbib['meta']['multiple matches out of tolerance range'] == 1 or mbib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                continue
            for j,mmbib in mbib['matches'].items():
                G.add_edge(f'{org}_{i}', f'{org2}_{j}', weight=mmbib['match score'])
with open('graph.pickle','wb') as f:
    pickle.dump(G,f)
