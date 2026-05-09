#!/usr/bin/env python3
"""
Generate statistics for aligned_succinct format.

Reports:
- Total anchors with valid matches
- Matches per organism pair
- Score distributions
- Orientation distributions

Usage: python statistics_succinct.py <organism> <work_dir> <anchor_dir>
"""

import sys
import pickle
from pathlib import Path
from statistics import mean, median, stdev
from collections import defaultdict


def compute_statistics(org, work_dir, anchor_dir):
    """
    Compute and print statistics for the succinct format.
    """
    succinct_path = Path(anchor_dir) / 'aligned_succinct' / org

    if not succinct_path.exists():
        print(f"ERROR: Succinct file not found: {succinct_path}", file=sys.stderr)
        sys.exit(1)

    # Load orgs list
    orgs_path = Path(work_dir) / 'orgs'
    orgs = []
    if orgs_path.exists():
        with open(orgs_path, 'r') as f:
            orgs = [line.strip() for line in f if line.strip()]

    # Load genome sizes for coverage calculation
    size = {}
    root_dir = Path(anchor_dir).parent
    for o in orgs:
        small_meta_path = root_dir / 'utils' / 'small_meta' / o
        if small_meta_path.exists():
            with open(small_meta_path, 'rb') as f:
                seqids, seqlen = pickle.load(f)
            size[o] = seqlen[-1]

    print(f"Loading succinct data for {org}...")
    with open(succinct_path, 'rb') as f:
        succinct = pickle.load(f)

    print('-----')
    print(f'Organism: {org}')
    print(f'Total anchors with valid matches: {len(succinct)}')

    # Collect per-org statistics
    matches_per_org = defaultdict(int)
    scores_per_org = defaultdict(list)
    aligned_bp_per_org = defaultdict(int)
    orientations = defaultdict(lambda: {'forward': 0, 'reverse': 0})

    for anchor_idx, anchor_data in succinct.items():
        for org2, org2_matches in anchor_data['matches'].items():
            for j, match_tuple in org2_matches.items():
                matches_per_org[org2] += 1

                score, coords, orientation = match_tuple
                scores_per_org[org2].append(int(score))
                aligned_bp_per_org[org2] += int(coords[1] - coords[0])
                orientations[org2][orientation] += 1

    print(f'\nMatches per organism:')
    for org2 in sorted(matches_per_org.keys()):
        count = matches_per_org[org2]
        scores = scores_per_org[org2]

        # Score statistics
        if len(scores) > 1:
            score_stats = f"avg={mean(scores):.1f}, med={median(scores):.1f}, std={stdev(scores):.1f}"
        elif len(scores) == 1:
            score_stats = f"score={scores[0]}"
        else:
            score_stats = "no scores"

        # Orientation distribution
        fwd = orientations[org2]['forward']
        rev = orientations[org2]['reverse']
        ori_stats = f"forward={fwd}, reverse={rev}"

        # Coverage if we have genome size
        if org in size:
            coverage = aligned_bp_per_org[org2] / size[org] * 100
            cov_stats = f", coverage={coverage:.2f}%"
        else:
            cov_stats = ""

        print(f"  {org2}: {count} matches ({score_stats}, {ori_stats}{cov_stats})")

    # Summary statistics
    total_matches = sum(matches_per_org.values())
    all_scores = [s for scores in scores_per_org.values() for s in scores]

    print(f'\nSummary:')
    print(f'  Total matches across all organisms: {total_matches}')
    if all_scores:
        print(f'  Overall score: avg={mean(all_scores):.1f}, med={median(all_scores):.1f}')

    # Anchor length statistics
    anchor_lengths = [int(a['end'] - a['start']) for a in succinct.values()]
    if anchor_lengths:
        print(f'  Anchor lengths: avg={mean(anchor_lengths):.1f}, med={median(anchor_lengths):.1f}')
        if len(anchor_lengths) > 1:
            print(f'                  std={stdev(anchor_lengths):.1f}')

    print('\nstatistics_succinct done')


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <organism> <work_dir> <anchor_dir>", file=sys.stderr)
        sys.exit(1)

    org = sys.argv[1]
    work_dir = sys.argv[2]
    anchor_dir = sys.argv[3]

    compute_statistics(org, work_dir, anchor_dir)
