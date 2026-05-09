#! /usr/bin/env python3

"""
Generate plotting order for draw_verbose.py based on alignment coverage.

This script calculates alignment coverage between organisms and chromosomes,
then generates a plotting_order file with (org, chromo) pairs ordered by
alignment relationships. Uses a greedy approach: start with first organism,
add its chromosomes, then move to the organism with highest alignment coverage.

Usage:
    ./generate_plotting_order.py

Reads from:
    - orgs: list of organism IDs
    - anchors/aligned/{org}: aligned anchor data for each organism
    - small_meta/{org}: chromosome metadata for each organism

Writes to:
    - plotting_order: chromosome-level plotting order (org chromo pairs)
"""

import pickle
import os
from os.path import isfile

def get_alignment_coverage(org, orgs, maps_dir='anchors/aligned', meta_dir='small_meta'):
    """
    Calculate alignment coverage from one organism to all others.

    Returns dict mapping org -> total aligned length (as fraction of org's genome)
    """
    res = {}

    # Load genome sizes
    size = {}
    for o in orgs:
        meta_found = False
        meta_paths_tried = []
        for meta_path in [meta_dir, f'../../utils/{meta_dir}', f'../../../utils/{meta_dir}']:
            meta_file = f'{meta_path}/{o}'
            meta_paths_tried.append(meta_file)
            try:
                with open(meta_file, 'rb') as f:
                    seqids, seqlen = pickle.load(f)
                size[o] = seqlen[-1]
                meta_found = True
                break
            except FileNotFoundError:
                continue
        if not meta_found:
            print(f'Warning: metadata for {o} not found in: {", ".join(meta_paths_tried)}')
            continue

    # Load alignment map for this organism
    map_found = False
    map_paths_tried = []
    for maps_path in [maps_dir, f'../../utils/{maps_dir}', f'../../../utils/{maps_dir}', f'../../utils/anchors/aligned', f'../../../utils/anchors/aligned']:
        map_file = f'{maps_path}/{org}'
        map_paths_tried.append(map_file)
        try:
            with open(map_file, 'rb') as f:
                anchor_map = pickle.load(f)
            map_found = True
            break
        except FileNotFoundError:
            continue
    if not map_found:
        print(f'Warning: alignment data for {org} not found in: {", ".join(map_paths_tried)}')
        return res

    # Calculate aligned lengths
    aligned = {}
    scores = {}
    for o in orgs:
        if o == org:
            continue
        aligned[o] = 0
        scores[o] = []

    for i, d in anchor_map.items():
        for o, d2 in d['matches'].items():
            if o not in orgs:
                continue
            if d2['meta']['multiple matches out of tolerance range'] == 1 or d2['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
                continue
            if 'matches' in d2:
                for j, bib_m in d2['matches'].items():
                    scores[o].append(bib_m['match score'])
                    aligned[o] += bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0]
            if 'dups_matches' in d2 and 'syntenic' in d2['dups_matches']:
                js = d2['dups_matches']['syntenic']
                for j in js:
                    bib_m = d2['dups_matches'][j]
                    scores[o].append(bib_m['match score'])
                    aligned[o] += bib_m[f'hit coordinates in (own) {org} candidate'][1] - bib_m[f'hit coordinates in (own) {org} candidate'][0]

    for o in orgs:
        if o == org:
            continue
        if len(scores[o]) > 0 and org in size:
            res[o] = aligned[o] / size[org]

    return res


def get_chromosomes_for_org(org, maps_dir='anchors/aligned'):
    """
    Get list of chromosomes for an organism that have alignment data.

    Returns list of chromosome IDs.
    """
    anchor_map = None
    map_paths_tried = []
    for maps_path in [maps_dir, f'../../utils/{maps_dir}', f'../../../utils/{maps_dir}', f'../../utils/anchors/aligned', f'../../../utils/anchors/aligned']:
        map_file = f'{maps_path}/{org}'
        map_paths_tried.append(map_file)
        try:
            with open(map_file, 'rb') as f:
                anchor_map = pickle.load(f)
            break
        except FileNotFoundError:
            continue

    if anchor_map is None:
        print(f'Warning: alignment data for {org} not found in: {", ".join(map_paths_tried)}')
        return []

    chromosomes = set()
    for i, d in anchor_map.items():
        if 'chromosome' in d:
            chromosomes.add(d['chromosome'])

    return sorted(list(chromosomes))


if __name__ == "__main__":
    # Read organism list - try multiple locations
    orgs = []
    orgs_paths = ['orgs', '../../utils/orgs', '../../../utils/orgs']
    orgs_file_found = False

    for orgs_path in orgs_paths:
        if os.path.isfile(orgs_path):
            with open(orgs_path) as f:
                for line in f:
                    orgs.append(line.strip())
            orgs_file_found = True
            print(f'Loaded {len(orgs)} organisms from {orgs_path}')
            break

    # Fallback: auto-detect from pairwise_alignments_table
    if not orgs_file_found or len(orgs) == 0:
        print('Warning: orgs file not found, auto-detecting organisms from pairwise_alignments_table...')
        orgs_set = set()
        if os.path.isfile('pairwise_alignments_table'):
            with open('pairwise_alignments_table') as f:
                for line in f:
                    if line.strip():
                        cols = line.strip().split('\t')
                        if len(cols) >= 2:
                            orgs_set.add(cols[0])
                            orgs_set.add(cols[1])
        orgs = sorted(list(orgs_set))
        print(f'Auto-detected {len(orgs)} organisms from pairwise_alignments_table')

    if len(orgs) == 0:
        print('Error: no organisms found')
        exit(1)

    # Calculate pairwise alignment coverage between all organisms
    print('Calculating alignment coverage between organisms...')
    gres = {}
    for org in orgs:
        gres[org] = get_alignment_coverage(org, orgs)

    # Greedily order organisms by alignment coverage
    # Start with first organism, then pick organism with highest coverage
    orders = []
    last = orgs[0]
    while len(orders) != len(orgs) - 1:
        maxi = [0, '']
        for o in orgs:
            if o not in gres[last] or gres[last][o] == 0:
                continue
            if gres[last][o] > maxi[0] and o not in [xxx[1] for xxx in orders] and o != orgs[0]:
                maxi = [gres[last][o], o]
        if maxi[1] == '':
            # No more organisms with alignments, add remaining orgs
            for o in orgs:
                if o not in [xxx[0] for xxx in orders] and o not in [xxx[1] for xxx in orders] and o != orgs[0]:
                    maxi = [0, o]
                    break
        if maxi[1] != '':
            orders.append((last, maxi[1]))
            last = maxi[1]
        else:
            break

    # Create chromosome-level plotting order
    # For each organism in order, add all its chromosomes
    plotting_order = []
    orgs_in_order = [orgs[0]]
    for ref, target in orders:
        if target not in orgs_in_order:
            orgs_in_order.append(target)

    for org in orgs_in_order:
        chromosomes = get_chromosomes_for_org(org)
        for chromo in chromosomes:
            plotting_order.append((org, chromo))

    # Write plotting order to file
    with open('plotting_order', 'w') as f:
        for org, chromo in plotting_order:
            f.write(f'{org} {chromo}\n')

    print(f'Generated plotting_order with {len(plotting_order)} chromosome tracks')
    print(f'Organism order: {" -> ".join(orgs_in_order)}')
