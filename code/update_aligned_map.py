#!/usr/bin/env python3

from sys import argv
import pathlib,pickle,os

if __name__ == "__main__":

    org = argv[1]
    
    code_dir = pathlib.Path(__file__).parent.resolve()
    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'
    work_dir = argv[-1]
    # exit if first run
    try:
        with open(anchor_dir + '/candidates' + f'/{org}','rb') as f:
            bib = pickle.load(f)
        with open(anchor_dir + '/aligned' + f'/{org}','rb') as f:
            aligned = pickle.load(f)
    except:
        exit(0)
    del bib

    # collect it for all orgs to be used in collect_output.py:
    # there it needs to be checked again if the deleted ones have impact on ranges/metadata
    # better to do at the end cos most will be checked anyway since new ones will cover those and be checked

    matches_to_del_global = {}
    files = [name for name in os.listdir(f'{work_dir}/to_del/{org}') if 'from' in name]
    for file in files:
        with open(f'{work_dir}/to_del/{org}/{file}','rb') as f:
            matches_to_del = pickle.load(f)
        for i,bib in matches_to_del.items():
            if i not in matches_to_del_global:
                matches_to_del_global[i] = bib
            else:
                for org2,s in bib.items():
                    if org2 not in matches_to_del_global[i]:
                        matches_to_del_global[i][org2] = s
                    else:
                        matches_to_del_global[i][org2].union(s)
            # could include it in loop above but i dont for clarity
            for org2,s in bib.items():
                for j in s:
                    if j in aligned[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues']:
                        del aligned[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues'][j]
                        if j in aligned[i]['matches'][org2]['matches']:
                            print(f'WARNING: match {j} of {org2} was in "matches" and "matches not considered" of i={i}')
                    if j in aligned[i]['matches'][org2]['matches']:
                        del aligned[i]['matches'][org2]['matches'][j]
                    else:
                        print('should not happen')

    with open(f'{work_dir}/to_del/{org}/matches_to_del','wb') as f:
        pickle.dump(matches_to_del_global,f)
    with open(anchor_dir + '/aligned' + f'/{org}','wb') as f:
        pickle.dump(aligned,f)
