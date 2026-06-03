#!/usr/bin/env python3
"""
Unified aligned map update script - Phase 3a
Combines update_aligned_1 + update_aligned_2 logic

Can be called in two stages for Snakemake barrier:
  Stage 1: Track matches and write notifications
  Stage 2: Read notifications and clean aligned map

Usage:
  update_aligned.py <org> <work_dir> stage1  # Write notifications
  update_aligned.py <org> <work_dir> stage2  # Read notifications and clean

Full flow:
1. Stage 1 (per genome):
   - Read to_del/{org}/to_del (complete record from Phase 2)
   - Load/create aligned/{org}
   - Track matches in deleted candidates
   - Write to_del/{org2}/from_{org} notifications
   - Delete candidates from aligned

2. BARRIER: Wait for all orgs to complete stage1

3. Stage 2 (per genome):
   - Read all to_del/{org}/from_* notifications
   - Remove references to other orgs' deleted candidates
   - Write final clean aligned/{org}
"""

import pickle
from sys import argv
import pathlib
import os

def which_chromo(seqlen, seqids, i):
    """Get chromosome name for genomic position."""
    from bisect import bisect_left
    pos = bisect_left(seqlen, i)
    if i == seqlen[pos]:
        pos += 1
    return seqids[pos]

def stage1_write_notifications(org, work_dir, root):
    """
    Stage 1: Track matches and write cross-genome notifications.

    Does NOT require aligned map to exist (first run handling).
    If aligned exists: tracks matches and deletes from aligned.
    If aligned doesn't exist: creates empty structure from candidates.
    """

    anchor_dir = root + '/anchors'

    # Load organism list
    orgs = []
    with open(f"{work_dir}/orgs", "r") as f:
        for line in f:
            orgs.append(line.strip())

    # Load to_del (complete record from Phase 2: filter_windows + dups)
    try:
        with open(f'{work_dir}/to_del/{org}/to_del', 'rb') as f:
            to_del_org = pickle.load(f)
    except:
        to_del_org = []
        print(f'No to_del file for {org} (first run or no deletions)')

    print(f'Processing {len(to_del_org)} candidates to delete for {org}')

    # Try to load existing aligned map
    aligned_exists = False
    try:
        with open(anchor_dir + f'/aligned/{org}', 'rb') as f:
            aligned = pickle.load(f)
        aligned_exists = True
        print(f'Loaded existing aligned map with {len(aligned)} candidates')
    except:
        # First run or aligned doesn't exist yet - create from candidates
        print(f'No existing aligned map for {org} (first run)')
        try:
            with open(anchor_dir + f'/candidates/{org}', 'rb') as f:
                candidates = pickle.load(f)
            # Create aligned with same structure as candidates
            aligned = {}
            for idx, data in candidates.items():
                aligned[idx] = data.copy()
            print(f'Created aligned map from {len(aligned)} candidates')
            aligned_exists = True  # We created it
        except:
            print(f'ERROR: No candidates file for {org}')
            return

    # Track matches to notify other orgs
    matches_to_del = {}

    # For each candidate to delete, track its matches
    for i in to_del_org:
        if i not in aligned:
            print(f'{i} not in {org}\'s aligned map (already deleted or never existed)')
            continue

        # Check matches with other organisms
        if 'matches' in aligned[i]:
            for org2, bib in aligned[i]['matches'].items():
                # Collect matches that still exist in org2
                # (skip ones that org2 also deleted)
                s = set()

                # 'matches not considered...' always exists (initialized in parse_bcamm)
                for j in bib['matches not considered upon applying stricter score criterion since there are consistency issues']:
                    s.add(j)

                # Check regular matches (optional - may only have dups_matches)
                if 'matches' in bib:
                    for j in bib['matches']:
                        s.add(j)

                if len(s) > 0:
                    for j in s:
                        # Build notification: "org deleted i, which matched org2's j"
                        if org2 not in matches_to_del:
                            matches_to_del[org2] = {}
                        if j not in matches_to_del[org2]:
                            matches_to_del[org2][j] = {org: set([i])}
                        else:
                            if org not in matches_to_del[org2][j]:
                                matches_to_del[org2][j][org] = set([i])
                            else:
                                matches_to_del[org2][j][org].add(i)

    # Delete candidates from aligned
    deleted_count = 0
    for i in to_del_org:
        if i in aligned:
            del aligned[i]
            deleted_count += 1

    print(f'Deleted {deleted_count} candidates from aligned map')

    # Write updated aligned map
    # Directories created by setup_directories.py -> make_anchor_directories.py
    aligned_path = anchor_dir + f'/aligned/{org}'
    print(f'Writing aligned map to: {aligned_path}')
    print(f'  Directory exists: {os.path.isdir(anchor_dir + "/aligned")}')
    print(f'  Writing {len(aligned)} candidates')

    with open(aligned_path, 'wb') as f:
        pickle.dump(aligned, f)
        f.flush()  # Flush Python buffer
        os.fsync(f.fileno())  # Sync to disk

    # Verify write succeeded
    if os.path.exists(aligned_path):
        size = os.path.getsize(aligned_path)
        print(f'  SUCCESS: Wrote {size} bytes to {aligned_path}')
    else:
        print(f'  ERROR: File not found after write!')

    # Write cross-genome notifications
    for org2, bib in matches_to_del.items():
        try:
            with open(f'{work_dir}/to_del/{org2}/from_{org}', 'wb') as f:
                pickle.dump(bib, f)
        except:
            print(f'WARNING: Could not write notification to {org2}/from_{org}')

    print(f'Stage 1 complete: wrote notifications to {len(matches_to_del)} organisms')

def stage2_clean_references(org, work_dir, root):
    """
    Stage 2: Read notifications from other orgs and clean aligned map.

    Removes references to candidates deleted in other genomes.
    """

    anchor_dir = root + '/anchors'

    # Load aligned map
    try:
        with open(anchor_dir + f'/aligned/{org}', 'rb') as f:
            aligned = pickle.load(f)
    except:
        print(f'ERROR: No aligned map for {org} in stage2')
        return

    # Collect all notifications from other organisms
    matches_to_del_global = {}

    try:
        files = [name for name in os.listdir(f'{work_dir}/to_del/{org}') if 'from_' in name]
    except:
        files = []
        print(f'No notifications directory for {org}')

    print(f'Processing {len(files)} notifications from other organisms')

    for file in files:
        try:
            with open(f'{work_dir}/to_del/{org}/{file}', 'rb') as f:
                matches_to_del = pickle.load(f)

            # Merge notifications
            for i, bib in matches_to_del.items():
                if i not in matches_to_del_global:
                    matches_to_del_global[i] = bib
                else:
                    for org2, s in bib.items():
                        if org2 not in matches_to_del_global[i]:
                            matches_to_del_global[i][org2] = s
                        else:
                            matches_to_del_global[i][org2] = matches_to_del_global[i][org2].union(s)
        except Exception as e:
            print(f'ERROR reading notification {file}: {e}')
            continue

    # Clean aligned map
    cleaned_count = 0
    for i, bib in matches_to_del_global.items():
        if i not in aligned:
            print(f'WARNING: {i} from notifications not in aligned map')
            continue

        for org2, s in bib.items():
            if org2 not in aligned[i]['matches']:
                continue

            for j in s:
                # Remove from matches not considered - key always exists (initialized in parse_bcamm)
                if j in aligned[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues']:
                    del aligned[i]['matches'][org2]['matches not considered upon applying stricter score criterion since there are consistency issues'][j]
                    cleaned_count += 1

                # Remove from regular matches (optional - may only have dups_matches)
                if 'matches' in aligned[i]['matches'][org2]:
                    if j in aligned[i]['matches'][org2]['matches']:
                        del aligned[i]['matches'][org2]['matches'][j]
                        cleaned_count += 1

                # Bug #2 fix: Remove from dups_matches (optional key)
                if 'dups_matches' in aligned[i]['matches'][org2]:
                    if j in aligned[i]['matches'][org2]['dups_matches']:
                        del aligned[i]['matches'][org2]['dups_matches'][j]
                        cleaned_count += 1
                    # Also remove from syntenic set if present
                    if 'syntenic' in aligned[i]['matches'][org2]['dups_matches']:
                        aligned[i]['matches'][org2]['dups_matches']['syntenic'].discard(j)

    print(f'Cleaned {cleaned_count} match references')

    # Write final cleaned aligned map
    with open(anchor_dir + f'/aligned/{org}', 'wb') as f:
        pickle.dump(aligned, f)

    # Write matches_to_del for potential use by collect_output
    with open(f'{work_dir}/to_del/{org}/matches_to_del', 'wb') as f:
        pickle.dump(matches_to_del_global, f)

    print(f'Stage 2 complete: final aligned map has {len(aligned)} candidates')

if __name__ == "__main__":

    if len(argv) < 4:
        print('Usage: update_aligned.py <org> <work_dir> <stage1|stage2>')
        exit(1)

    org = argv[1]
    work_dir = argv[2]
    stage = argv[3]

    root = str(pathlib.Path(__file__).parents[1])

    print(f'======= {org} =======')

    if stage == 'stage1':
        print('Phase 3a Stage 1: Writing cross-genome notifications')
        stage1_write_notifications(org, work_dir, root)
    elif stage == 'stage2':
        print('Phase 3a Stage 2: Cleaning aligned map')
        stage2_clean_references(org, work_dir, root)
    else:
        print(f'ERROR: Unknown stage "{stage}". Use stage1 or stage2')
        exit(1)
