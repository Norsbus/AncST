#!/usr/bin/env python3

from sys import argv
import pathlib,pickle

if __name__ == "__main__":

    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]
    org = argv[1]
    orgs = []
    to_del = {}
    matches_to_del = {}

    # exit if first run
    try:
        with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
            candidates = pickle.load(f)
    except:
        exit(0)
   
    try:
        with open(anchor_dir + '/aligned' + f'/{org}','rb') as f:
            aligned = pickle.load(f)
    except:
        aligned = False
    try:
        with open(anchor_dir + '/hashes' + f'/{org}','rb') as f:
            hashes = pickle.load(f)
    except:
        hashes = False
    
    with open(f"{work_dir}/orgs","r") as f:
        for line in f:
            orgs.append(line.strip())
    for o in orgs:
        try:
            with open(f'{work_dir}/to_del/{o}/to_del','rb') as f:
                to_del[o] = pickle.load(f)
        except:
            to_del[o] = []

    for i in to_del[org]:
        if aligned:
            if i not in aligned:
                print(f'{i} not in {org}\'s aligned map when looping to get matches')
                continue
            for org2,bib in aligned[i]['matches'].items():
                if org2 not in to_del:
                    try:
                        with open(f'{work_dir}/to_del/{org2}/to_del','rb') as f:
                            to_del[org2] = pickle.load(f)
                    except:
                        to_del[org2] = []
                s = set()
                for j in bib['matches not considered upon applying stricter score criterion since there are consistency issues']:
                    if j not in to_del[org2]:
                        s.add(j)
                for j in bib['matches']:
                    if j not in to_del[org2]:
                        s.add(j)
                if len(s) > 0:
                    for j in s:
                        # should not be initialised for orgs list in direcotory as matches can be from before for different set of orgs
                        if org2 not in matches_to_del:
                            matches_to_del[org2] = {}
                        if j not in matches_to_del[org2]:
                            matches_to_del[org2][j] = {org:set([i])}
                        else:
                            if org not in matches_to_del[org2][j]:
                                matches_to_del[org2][j][org] = set([i])
                            else:
                                matches_to_del[org2][j][org].add(i)
    for i in to_del[org]:
        if aligned:
            try:
                del aligned[i]
            except:
                print(f'{i} not in {org}\'s aligned map when trying to delete {i}')
        if hashes:
            try:
                del hashes[i]
            except:
                print(f'{i} not in {org}\'s hashes when trying to delete {i}')
        try:
            del candidates[i] 
        except:
            print(f'{i} not in {org}\'s candidates map when trying to delete {i}')

    for o in orgs:
        try:
            with open(anchor_dir + '/compared' + f'/{o}'+f'/with_{org}','rb') as f:
                compared = pickle.load(f)
            for i in to_del[org]:
                del compared[i]
            with open(anchor_dir + '/compared' + f'/{o}'+f'/with_{org}','wb') as f:
                pickle.dump(compared,f)
        except:
            continue

    with open(anchor_dir + '/candidates' + f'/{org}','wb') as f:
        pickle.dump(candidates,f)
    if aligned:
        with open(anchor_dir + '/aligned' + f'/{org}','wb') as f:
            pickle.dump(aligned,f)
    if hashes:
        with open(anchor_dir + '/hashes' + f'/{org}','wb') as f:
            pickle.dump(hashes,f)

    if aligned:
        for org2,bib in matches_to_del.items():
            try:
                with open(f'{work_dir}/to_del/{org2}/from_{org}','wb') as f:
                    pickle.dump(bib,f) 
            except:
                print(f'could not open {work_dir}/to_del/{org2}/from_{org}')
