#!/usr/bin/env python3
"""
Standalone converter: aligned anchor pickles -> aligned_succinct pickles.

Reproduces the exact conversion logic from code/finalize_aligned.py:
  1. Reciprocity check on syntenic dups (in-memory flag setting only)
  2. Create succinct format with filtered/compressed match tuples

SAFETY:
  - All aligned pickle files are opened READ-ONLY ('rb'), never written to.
  - In-memory modifications (reciprocity flags) never touch disk.
  - Output goes exclusively to --succinct directory.
  - --succinct must differ from --aligned (enforced).

PREREQUISITE:
  The aligned pickles must have already been through eval_syn (syntenic
  marking). This is the normal state of production aligned/ files after
  the pipeline's finalize_aligned eval_syn step.

Usage:
  python convert_aligned_to_succinct.py --aligned /path/to/anchors/aligned \\
                                        --succinct /path/to/output \\
                                        --cores 4
"""

import sys
import os
import argparse
import pickle
from pathlib import Path
from multiprocessing import Pool


# Exact copies of functions from code/finalize_aligned.py

def check_forward_ambiguity(match_bib):
    """Check if match has forward ambiguity (should be filtered out)."""
    if match_bib['meta']['multiple matches out of tolerance range'] == 1 or \
       match_bib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
        return True
    return False


def check_eval_syn_reciprocity(org, aligned_map, aligned_dir, all_orgs, suffix=''):
    """
    Check that syntenic dup markings are reciprocal across organisms.

    For each anchor with syntenic dups to org2, verify that the partner's
    aligned map also marks the reverse as syntenic. If ANY syntenic dup in
    a (anchor, org2) relationship is NOT reciprocal, set the ambiguous flag
    on the org2_data['meta'] IN MEMORY ONLY.

    All partner file access is READ-ONLY.

    Returns number of (anchor, org2) relationships flagged as ambiguous.
    """
    other_orgs = [o for o in all_orgs if o != org]
    flagged_count = 0

    for org2 in other_orgs:
        # Quick scan: does own map have ANY syntenic dups for org2?
        has_syntenic_for_org2 = False
        for anchor_idx, anchor_data in aligned_map.items():
            if org2 not in anchor_data['matches']:
                continue
            org2_data = anchor_data['matches'][org2]
            if ('dups_matches' in org2_data and
                    'syntenic' in org2_data['dups_matches'] and
                    len(org2_data['dups_matches']['syntenic']) > 0):
                has_syntenic_for_org2 = True
                break

        if not has_syntenic_for_org2:
            continue

        # Load partner's aligned map (READ-ONLY)
        partner_path = Path(aligned_dir) / (org2 + suffix)
        if not partner_path.exists():
            print(f"  WARNING: Partner aligned map not found: {partner_path}", file=sys.stderr)
            continue

        with open(partner_path, 'rb') as f:
            partner_map = pickle.load(f)

        # Check each anchor with syntenic dups to org2
        for anchor_idx, anchor_data in aligned_map.items():
            if org2 not in anchor_data['matches']:
                continue

            org2_data = anchor_data['matches'][org2]
            if ('dups_matches' not in org2_data or
                    'syntenic' not in org2_data['dups_matches'] or
                    len(org2_data['dups_matches']['syntenic']) == 0):
                continue

            syntenic_set = org2_data['dups_matches']['syntenic']

            # Check each position j in the syntenic set
            relationship_nonreciprocal = False
            for j in syntenic_set:
                if j not in partner_map:
                    relationship_nonreciprocal = True
                    break

                if org not in partner_map[j]['matches']:
                    relationship_nonreciprocal = True
                    break

                partner_org_data = partner_map[j]['matches'][org]

                if 'dups_matches' not in partner_org_data:
                    relationship_nonreciprocal = True
                    break

                partner_dups = partner_org_data['dups_matches']

                if 'syntenic' not in partner_dups:
                    relationship_nonreciprocal = True
                    break

                if anchor_idx not in partner_dups['syntenic']:
                    relationship_nonreciprocal = True
                    break

            if relationship_nonreciprocal:
                # IN MEMORY ONLY — never written back to disk
                org2_data['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 1
                flagged_count += 1

        # Free partner map
        del partner_map

    return flagged_count


def create_succinct(org, aligned_map):
    """
    Create succinct format from aligned map.
    Exact replica of the logic in code/finalize_aligned.py.

    - Only includes anchors with at least one valid match
    - Only includes first match per org2 (sorted by position)
    - For dups_matches, only includes syntenic ones and only if no regular matches
    - Each match becomes tuple: (score_int, (start, end), orientation_str)
    """
    succinct = {}

    for i, bib1 in aligned_map.items():
        valid = False
        new_entry = {
            'chromosome': bib1['chromosome'],
            'start': bib1['start'],
            'end': bib1['end'],
            'matches': {}
        }

        for org2, bib2 in bib1['matches'].items():
            has_matches = 'matches' in bib2 and len(bib2['matches']) > 0
            has_syntenic_dups = ('dups_matches' in bib2 and
                                 'syntenic' in bib2['dups_matches'] and
                                 len(bib2['dups_matches']['syntenic']) > 0)

            if not has_matches and not has_syntenic_dups:
                continue

            if check_forward_ambiguity(bib2):
                continue

            valid = True
            new_entry['matches'][org2] = {}

            # Try regular matches first — take only the first one (sorted by position)
            if has_matches:
                for j, bib3 in sorted(bib2['matches'].items()):
                    if bib3['match is on other strand in other genome']:
                        ori = 'reverse'
                    else:
                        ori = 'forward'
                    new_entry['matches'][org2][j] = (
                        int(bib3['match score']),
                        tuple(bib3[f'hit coordinates in (own) {org} candidate']),
                        ori
                    )
                    break  # Only take first match

            # If no regular matches, try syntenic dups
            if len(new_entry['matches'][org2]) == 0 and has_syntenic_dups:
                syn = bib2['dups_matches']['syntenic']
                dups_items = [(k, v) for k, v in bib2['dups_matches'].items() if k != 'syntenic']
                for j, bib3 in sorted(dups_items):
                    if j not in syn:
                        continue
                    if bib3['match is on other strand in other genome']:
                        ori = 'reverse'
                    else:
                        ori = 'forward'
                    new_entry['matches'][org2][j] = (
                        int(bib3['match score']),
                        tuple(bib3[f'hit coordinates in (own) {org} candidate']),
                        ori
                    )
                    break  # Only take first syntenic dup

        if valid:
            succinct[i] = new_entry

    return succinct


# Worker + main

def convert_one(args):
    """Convert a single organism's aligned pickle to succinct format."""
    org, aligned_dir, succinct_dir, all_orgs, suffix = args

    aligned_path = Path(aligned_dir) / (org + suffix)
    succinct_path = Path(succinct_dir) / org

    # Load aligned map (READ-ONLY from disk)
    with open(aligned_path, 'rb') as f:
        aligned_map = pickle.load(f)

    # Reciprocity check (modifies IN-MEMORY copy only, reads partner files read-only)
    flagged = check_eval_syn_reciprocity(org, aligned_map, aligned_dir, all_orgs, suffix)

    # Create succinct from in-memory modified map
    succinct = create_succinct(org, aligned_map)

    # Write succinct to output directory (never touches aligned dir)
    with open(succinct_path, 'wb') as f:
        pickle.dump(succinct, f)

    return org, len(aligned_map), len(succinct), flagged


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert aligned anchor pickles to succinct format. READ-ONLY for aligned pickles.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--aligned', required=True, type=str,
                        help='Path to aligned/ directory containing anchor pickles')
    parser.add_argument('--succinct', required=True, type=str,
                        help='Output directory for succinct pickles (must differ from --aligned)')
    parser.add_argument('--cores', type=int, default=1,
                        help='Number of parallel workers')
    parser.add_argument('--suffix', type=str, default='_with_syn_eval',
                        help='Filename suffix on aligned pickles (stripped to get org name). '
                             'Set to empty string if filenames are just the org name.')
    args = parser.parse_args()

    aligned_dir = os.path.abspath(args.aligned)
    succinct_dir = os.path.abspath(args.succinct)

    # Safety: output dir must not be the input dir
    if aligned_dir == succinct_dir:
        print("ERROR: --succinct must be a different directory than --aligned", file=sys.stderr)
        sys.exit(1)

    # Safety: output dir must not be inside input dir
    if succinct_dir.startswith(aligned_dir + os.sep):
        print("ERROR: --succinct must not be inside --aligned directory", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(aligned_dir):
        print(f"ERROR: --aligned directory does not exist: {aligned_dir}", file=sys.stderr)
        sys.exit(1)

    suffix = args.suffix

    # List all organism files in aligned dir (skip hidden files and directories)
    all_files = sorted([
        f for f in os.listdir(aligned_dir)
        if os.path.isfile(os.path.join(aligned_dir, f)) and not f.startswith('.')
    ])

    if not all_files:
        print(f"ERROR: No files found in {aligned_dir}", file=sys.stderr)
        sys.exit(1)

    # Strip suffix to get org names
    if suffix:
        all_orgs = []
        for f in all_files:
            if f.endswith(suffix):
                all_orgs.append(f[:-len(suffix)])
            else:
                print(f"WARNING: File '{f}' does not end with suffix '{suffix}', skipping", file=sys.stderr)
    else:
        all_orgs = all_files

    if not all_orgs:
        print(f"ERROR: No files matching suffix '{suffix}' found in {aligned_dir}", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(succinct_dir, exist_ok=True)

    print(f"Converting {len(all_orgs)} organisms: {aligned_dir} -> {succinct_dir} ({args.cores} cores)")
    if suffix:
        print(f"  File suffix: '{suffix}' (e.g. {all_orgs[0]}{suffix} -> {all_orgs[0]})")

    # Build task list
    tasks = [(org, aligned_dir, succinct_dir, all_orgs, suffix) for org in all_orgs]

    if args.cores == 1:
        # Sequential for easier debugging
        for task in tasks:
            org, n_aligned, n_succinct, flagged = convert_one(task)
            print(f"  {org}: {n_succinct}/{n_aligned} anchors kept, {flagged} relationships flagged")
    else:
        with Pool(args.cores) as pool:
            for org, n_aligned, n_succinct, flagged in pool.imap_unordered(convert_one, tasks):
                print(f"  {org}: {n_succinct}/{n_aligned} anchors kept, {flagged} relationships flagged")

    print("Done.")
