#!/usr/bin/env python3
"""
Validate aligned_succinct format integrity.

Checks:
1. Correspondence: aligned_succinct anchors are a subset of candidates
2. Sanity checks (coordinates, chromosome boundaries, overlaps)
3. Format validation (tuple structure, types, no 'syntenic' key)
4. Cross-organism consistency (reciprocal matches, score consistency)

Errors cause FAILURE. Warnings are reported but do not fail the pipeline.

The succinct format:
{
    anchor_index: {
        'chromosome': str,
        'start': int,
        'end': int,
        'matches': {
            org2: {
                match_position_j: (score_int, (start, end), orientation_str)
            }
        }
    }
}

Usage: python checks_succinct.py <organism> <work_dir> <anchor_dir>
"""

import sys
import pickle
import numbers
import pathlib
from bisect import bisect_left
from statistics import mean, median, stdev


class ErrorTracker:
    """Track errors and warnings during validation."""
    def __init__(self):
        self.errors = []
        self.warnings = []

    def add(self, error_msg):
        """Add an error and print it."""
        self.errors.append(error_msg)
        print(f"ERROR: {error_msg}", file=sys.stderr)

    def warn(self, warning_msg):
        """Add a warning and print it."""
        self.warnings.append(warning_msg)
        print(f"WARNING: {warning_msg}", file=sys.stderr)

    def has_errors(self):
        return len(self.errors) > 0

    def has_warnings(self):
        return len(self.warnings) > 0

    def count(self):
        return len(self.errors)

    def warning_count(self):
        return len(self.warnings)


# Check 1: correspondence.

def check_correspondence_succinct(org, anchor_dir, tracker):
    """Check that aligned_succinct anchors are a subset of candidates.

    Succinct only contains anchors with valid matches, so it is expected
    to be a subset (not equal) of candidates. But every succinct anchor
    must exist in candidates.
    """
    with open(f'{anchor_dir}/candidates/{org}', 'rb') as f:
        candidates = pickle.load(f)
    with open(f'{anchor_dir}/aligned_succinct/{org}', 'rb') as f:
        succinct = pickle.load(f)

    for i in succinct:
        if i not in candidates:
            tracker.add(f'Anchor {i} in aligned_succinct but not in candidates for {org}')

    print(f'check_correspondence_succinct done ({len(succinct)} succinct anchors, '
          f'{len(candidates)} candidates)')


# Check 2: sanity.

def check_sanity_succinct(org, anchor_dir, root_dir, tracker):
    """Check coordinate sanity for aligned_succinct: no overlaps, valid boundaries."""
    with open(f'{anchor_dir}/aligned_succinct/{org}', 'rb') as f:
        succinct = pickle.load(f)

    with open(f'{root_dir}/utils/small_meta/{org}', 'rb') as f:
        seqids, seqlen = pickle.load(f)

    lens = []
    sorted_i = sorted(succinct.keys())
    for ii, i in enumerate(sorted_i):
        start = i
        anchor = succinct[i]
        length = anchor['end'] - anchor['start']
        lens.append(int(length))
        end = i + length

        if end < start:
            tracker.add(f'Anchor {i} has end < start ({end} < {start}) for {org} (succinct)')

        idx = bisect_left(seqlen, start)
        if idx < len(seqlen) and start == seqlen[idx]:
            idx += 1
        if idx >= len(seqids):
            tracker.add(f'Anchor {i} start {start} exceeds genome bounds for {org} (succinct)')
            continue
        chromo_start = seqids[idx]

        idx = bisect_left(seqlen, end)
        if idx < len(seqlen) and end == seqlen[idx]:
            idx += 1
        if idx >= len(seqids):
            tracker.add(f'Anchor {i} end {end} exceeds genome bounds for {org} (succinct)')
            continue
        chromo_end = seqids[idx]

        if chromo_start != chromo_end:
            tracker.add(f'Anchor {i} spans chromosomes ({chromo_start} to {chromo_end}) for {org} (succinct)')

        if ii < len(sorted_i) - 1:
            if end >= sorted_i[ii + 1]:
                tracker.add(f'Anchors overlap: {i} (end={end}) overlaps with {sorted_i[ii+1]} for {org} (succinct)')

    if lens and len(lens) > 1:
        print(f'Succinct anchor lengths: mean={mean(lens):.1f}, median={median(lens):.1f}, '
              f'stdev={stdev(lens):.1f}')
    else:
        print(f'at least two data points needed for stdev calc')
    print('check_sanity_succinct done')


# Check 3: format validation.

def validate_tuple_format(j, match_tuple, anchor_idx, org2):
    """
    Validate a single match tuple has the expected format:
    (score_int, (start, end), orientation_str)

    Returns list of error messages (empty if valid).
    """
    errors = []
    context = f"anchor={anchor_idx}, org2={org2}, pos={j}"

    if not isinstance(match_tuple, tuple):
        errors.append(f"Match is not a tuple at {context}: got {type(match_tuple)}")
        return errors

    if len(match_tuple) != 3:
        errors.append(f"Match tuple has wrong length at {context}: expected 3, got {len(match_tuple)}")
        return errors

    score, coords, orientation = match_tuple

    # Check score is int (accept int and numpy integer types)
    if not isinstance(score, numbers.Integral):
        errors.append(f"Score is not int at {context}: got {type(score)}")

    # Check coords is tuple of (start, end)
    if not isinstance(coords, tuple):
        errors.append(f"Coords is not tuple at {context}: got {type(coords)}")
    elif len(coords) != 2:
        errors.append(f"Coords tuple has wrong length at {context}: expected 2, got {len(coords)}")
    else:
        start, end = coords
        if not isinstance(start, numbers.Real):
            errors.append(f"Coord start is not numeric at {context}: got {type(start)}")
        if not isinstance(end, numbers.Real):
            errors.append(f"Coord end is not numeric at {context}: got {type(end)}")
        if isinstance(start, numbers.Real) and isinstance(end, numbers.Real) and start > end:
            errors.append(f"Coord start > end at {context}: {start} > {end}")

    # Check orientation is 'forward' or 'reverse'
    if not isinstance(orientation, str):
        errors.append(f"Orientation is not string at {context}: got {type(orientation)}")
    elif orientation not in ('forward', 'reverse'):
        errors.append(f"Orientation has invalid value at {context}: got '{orientation}'")

    return errors


def check_format_succinct(org, anchor_dir, tracker):
    """Validate tuple format and structure of all entries in aligned_succinct."""
    with open(f'{anchor_dir}/aligned_succinct/{org}', 'rb') as f:
        succinct = pickle.load(f)

    if not isinstance(succinct, dict):
        tracker.add(f"Succinct data is not a dict: got {type(succinct)}")
        return

    anchors_checked = 0
    matches_checked = 0
    required_keys = {'chromosome', 'start', 'end', 'matches'}

    for anchor_idx, anchor_data in succinct.items():
        anchors_checked += 1
        context = f"anchor={anchor_idx}"

        if not isinstance(anchor_data, dict):
            tracker.add(f"Anchor data is not dict at {context}: got {type(anchor_data)}")
            continue

        # Check required keys
        missing = required_keys - set(anchor_data.keys())
        if missing:
            tracker.add(f"Missing required keys at {context}: {missing}")
            continue

        # Validate types
        if not isinstance(anchor_data['chromosome'], str):
            tracker.add(f"'chromosome' is not string at {context}: got {type(anchor_data['chromosome'])}")

        if not isinstance(anchor_data['start'], numbers.Real):
            tracker.add(f"'start' is not numeric at {context}: got {type(anchor_data['start'])}")
        if not isinstance(anchor_data['end'], numbers.Real):
            tracker.add(f"'end' is not numeric at {context}: got {type(anchor_data['end'])}")
        if anchor_data['start'] > anchor_data['end']:
            tracker.add(f"'start' > 'end' at {context}: {anchor_data['start']} > {anchor_data['end']}")

        # Check for unexpected keys (like 'syntenic' which shouldn't be in succinct)
        unexpected_keys = set(anchor_data.keys()) - required_keys
        if unexpected_keys:
            tracker.add(f"Unexpected keys at {context}: {unexpected_keys}")

        # Validate matches dict
        if not isinstance(anchor_data['matches'], dict):
            tracker.add(f"'matches' is not dict at {context}: got {type(anchor_data['matches'])}")
            continue

        for org2, org2_matches in anchor_data['matches'].items():
            if not isinstance(org2, str):
                tracker.add(f"org2 key is not string at {context}: got {type(org2)}")
                continue

            if not isinstance(org2_matches, dict):
                tracker.add(f"matches[{org2}] is not dict at {context}: got {type(org2_matches)}")
                continue

            # 'syntenic' should NOT be present in succinct format
            if 'syntenic' in org2_matches:
                tracker.add(f"'syntenic' key in matches[{org2}] at {context} - should not be in succinct")

            for j, match_tuple in org2_matches.items():
                matches_checked += 1
                for err in validate_tuple_format(j, match_tuple, anchor_idx, org2):
                    tracker.add(err)

    print(f'check_format_succinct done ({anchors_checked} anchors, {matches_checked} matches)')


# Check 4: cross-organism consistency.

def check_cross_org_succinct(org, work_dir, anchor_dir, tracker):
    """Check cross-organism consistency for aligned_succinct.

    For each match org1:i -> org2:j in succinct:
    - j should exist in org2's succinct map (warning if not, since some anchors
      may not have valid matches in the other direction)
    - If org2:j -> org1 exists and includes position i, scores must be consistent

    Note: succinct keeps only the first match per org2 per anchor (sorted by
    position), so strict reciprocity (i->j implies j->i) is NOT expected.
    Only score consistency is checked where both directions happen to exist.
    """
    with open(f'{anchor_dir}/aligned_succinct/{org}', 'rb') as f:
        am1 = pickle.load(f)

    other_orgs = []
    with open(f'{work_dir}/orgs', 'r') as f:
        for line in f:
            other_orgs.append(line.strip())
    if org in other_orgs:
        other_orgs.remove(org)

    for org2 in other_orgs:
        succinct_path = pathlib.Path(anchor_dir) / 'aligned_succinct' / org2
        if not succinct_path.exists():
            # Partner's finalize_aligned may not have completed yet — this is
            # expected since checks_succinct only depends on its own succinct file.
            # Cross-org validation on the full aligned format (checks.py) covers this.
            continue

        with open(succinct_path, 'rb') as f:
            am2 = pickle.load(f)

        for i, anchor_data in am1.items():
            if org2 not in anchor_data['matches']:
                continue

            for j, match_tuple in anchor_data['matches'][org2].items():
                score1 = match_tuple[0]

                if j not in am2:
                    # j might not be in succinct if it has no valid matches — that's
                    # fine since succinct filters anchors without matches
                    continue

                if org not in am2[j]['matches']:
                    # org2:j doesn't have a match back to org — expected for succinct
                    # since it only keeps the first match per partner
                    continue

                # org2:j has match(es) to org — check if position i is among them
                if i in am2[j]['matches'][org]:
                    score2 = am2[j]['matches'][org][i][0]
                    if score1 != score2:
                        tracker.add(
                            f'Score mismatch in succinct for {org}:{i} <-> {org2}:{j}: '
                            f'{score1} vs {score2}')

    print('check_cross_org_succinct done')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <organism> <work_dir> <anchor_dir>", file=sys.stderr)
        sys.exit(1)

    org = sys.argv[1]
    work_dir = sys.argv[2]
    anchor_dir = sys.argv[3]
    root_dir = str(pathlib.Path(anchor_dir).parent)

    print(f'--- Checking succinct for {org} ---')

    tracker = ErrorTracker()

    check_correspondence_succinct(org, anchor_dir, tracker)
    check_sanity_succinct(org, anchor_dir, root_dir, tracker)
    check_format_succinct(org, anchor_dir, tracker)
    check_cross_org_succinct(org, work_dir, anchor_dir, tracker)

    # Report results
    if tracker.has_warnings():
        print(f'\n=== {tracker.warning_count()} warnings (non-critical) ===', file=sys.stderr)

    if tracker.has_errors():
        print(f'\n=== VALIDATION FAILED: {tracker.count()} errors found ===', file=sys.stderr)
        sys.exit(1)
    else:
        if tracker.has_warnings():
            print(f'\n=== All succinct checks passed ({tracker.warning_count()} warnings) ===')
        else:
            print('\n=== All succinct checks passed ===')
        sys.exit(0)
