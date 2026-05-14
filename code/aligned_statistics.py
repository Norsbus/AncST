#! /usr/bin/env python3

import pickle,pathlib
from statistics import mean,stdev,median
from sys import argv
from finalize_aligned import check_eval_syn_reciprocity

def aligned(org):

    orgs = []
    with open(f'{work_dir}/orgs','r') as f:
        for l in f:
            orgs.append(l.strip())

    size = {}
    for o in orgs:
        with open(f'{work_dir}/../utils/small_meta/{o}','rb') as f:
            seqids,seqlen = pickle.load(f)
        size[o] = seqlen[-1]

    with open(f'{anchor_dir}/aligned/{org}','rb') as f:
        anchor_map = pickle.load(f)

    # match the in-memory reciprocity flag from finalize_succinct
    check_eval_syn_reciprocity(org, anchor_map, anchor_dir, work_dir)

    print('-----')
    print(org)
    aligned = {}
    scores = {}
    for o in orgs:
        if o == org:
            continue
        aligned[o] = 0
        scores[o] = []
    for i,d in anchor_map.items():
        for o,d2 in d['matches'].items():
            if o not in orgs:
                continue
            # 'meta' always exists when org2 entry exists (created in parse_bcamm)
            if d2['meta']['multiple matches out of tolerance range'] == 1 or d2['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                continue

            # At least one of 'matches' or 'dups_matches' must exist (parse_bcamm always creates one)
            has_matches = 'matches' in d2
            has_dups = 'dups_matches' in d2

            if not has_matches and not has_dups:
                raise KeyError(f"Data corruption: Neither 'matches' nor 'dups_matches' key found for anchor={i}, org={o}. Keys present: {list(d2.keys())}")

            # Process regular matches
            if has_matches:
                for j,bib_m in d2['matches'].items():
                    scores[o].append(int(bib_m['match score']))
                    aligned[o] += int(bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0])

            # Process dups_matches - only syntenic entries
            if has_dups:
                dups_dict = d2['dups_matches']
                # Check for 'syntenic' key (expected in _with_syn_eval files)
                if 'syntenic' in dups_dict:
                    syntenic_set = dups_dict['syntenic']
                else:
                    # No syntenic filtering - include all (older format)
                    syntenic_set = None

                for j,bib_m in dups_dict.items():
                    if j == 'syntenic':  # Skip the 'syntenic' key itself
                        continue
                    # If syntenic set exists, only include matches in that set
                    if syntenic_set is not None and j not in syntenic_set:
                        continue
                    scores[o].append(int(bib_m['match score']))
                    aligned[o] += int(bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0])
    
    for o in orgs:
        if o == org:
            continue

        if len(scores[o]) > 0:
            print(f'{aligned[o]/size[org]*100:.2f}% are aligned to with {o}')
            print(f'avg score with {o}: {mean(scores[o]):.2f}')
            print(f'median score with {o}: {median(scores[o]):.2f}')


if __name__ == "__main__":
    org = argv[1]
    print(f'--- {org} ---')
    work_dir = argv[2]
    anchor_dir = argv[3]
    aligned(org)
