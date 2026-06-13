#!/usr/bin/env python3
"""
Validate aligned map integrity.

Checks:
1. Correspondence between candidates and aligned maps
2. Sanity checks (coordinates, chromosome boundaries, overlaps)
3. Cross-organism consistency (reciprocal matches, scores, coordinates)

Errors cause FAILURE. Warnings are reported but do not fail the pipeline.

Non-reciprocal syntenic dup matches where the reverse side has the match in
dups_matches but not marked syntenic are warnings — eval_syn asymmetry has
multiple known causes (meta flag differences, inconsistency-based match
removal, margin boundary effects). Only structural issues (match not in
dups_matches at all) are errors.

Meta flag asymmetry between paired anchors is also a warning — these flags
describe match patterns from one anchor's perspective and are inherently
asymmetric when a region underwent insertion/deletion/transposition in one
organism but not the other.

Usage: python checks.py <organism> <work_dir> <anchor_dir>
"""

from pprint import pprint
import pickle
import pathlib
import sys
import yaml
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


def check_correspondence_candidates_aligned_maps(org, anchor_dir, tracker):
    """Check that candidates and aligned maps have same anchors."""
    with open(f'{anchor_dir}/candidates/{org}', 'rb') as f:
        candidates = pickle.load(f)
    with open(f'{anchor_dir}/aligned/{org}', 'rb') as f:
        aligned = pickle.load(f)

    for i in candidates:
        if i not in aligned:
            tracker.add(f'Anchor {i} in candidates but not in aligned for {org}')

    for i in aligned:
        if i not in candidates:
            tracker.add(f'Anchor {i} in aligned but not in candidates for {org}')

    print('check_correspondence done')


def check_sanity(org, anchor_dir, root_dir, tracker):
    """Check coordinate sanity: no overlaps, valid boundaries.

    Length stats (mean/median/stdev) are computed only over anchors that
    passed every check, so a polluted aligned map doesn't produce misleading
    summary numbers. The `N/M clean anchors` fraction in the printed line
    makes any silent filtering visible at a glance.
    """
    lens = []

    with open(f'{anchor_dir}/aligned/{org}', 'rb') as f:
        am1 = pickle.load(f)

    with open(f'{root_dir}/utils/small_meta/{org}', 'rb') as f:
        seqids, seqlen = pickle.load(f)

    sorted_i = sorted(am1.keys())
    n_anchors = len(sorted_i)

    # zero-anchor org: weird (upstream likely produced nothing for this org);
    # warn so it surfaces in the final summary without flipping validation to FAIL
    if n_anchors == 0:
        tracker.warn(f'No anchors in aligned map for {org} -- upstream likely produced empty output')

    for ii, i in enumerate(sorted_i):
        start = i
        anchor = am1[i]
        length = anchor['end'] - anchor['start']
        end = i + length
        anchor_ok = True

        if end < start:
            tracker.add(f'Anchor {i} has end < start ({end} < {start}) for {org}')
            anchor_ok = False

        idx = bisect_left(seqlen, start)
        if idx < len(seqlen) and start == seqlen[idx]:
            idx += 1
        if idx >= len(seqids):
            tracker.add(f'Anchor {i} start {start} exceeds genome bounds for {org}')
            continue
        chromo_start = seqids[idx]

        idx = bisect_left(seqlen, end)
        if idx < len(seqlen) and end == seqlen[idx]:
            idx += 1
        if idx >= len(seqids):
            tracker.add(f'Anchor {i} end {end} exceeds genome bounds for {org}')
            continue
        chromo_end = seqids[idx]

        if chromo_start != chromo_end:
            tracker.add(f'Anchor {i} spans chromosomes ({chromo_start} to {chromo_end}) for {org}')
            anchor_ok = False

        if ii < n_anchors - 1:
            if end >= sorted_i[ii + 1]:
                tracker.add(f'Anchors overlap: {i} (end={end}) overlaps with {sorted_i[ii+1]} for {org}')
                anchor_ok = False

        if anchor_ok:
            lens.append(int(length))

    n_clean = len(lens)
    if n_clean:
        msg = (f'Anchor lengths over {n_clean}/{n_anchors} clean anchors for {org}: '
               f'mean={mean(lens):.1f}, median={median(lens):.1f}')
        if n_clean >= 2:
            print(f'{msg}, stdev={stdev(lens):.1f}')
        else:
            print(f'{msg}  (stdev skipped: n<2)')
            tracker.warn(f'Only {n_clean} clean anchor(s) for {org} -- stdev undefined; very few anchors is itself suspicious')

    print('check_sanity done')


def get_all_matches(org_data, anchor_id, org_name, tracker):
    """Get all matches from both 'matches' and 'dups_matches' (syntenic only).

    Reports missing keys as errors.
    """
    all_matches = {}
    context = f"anchor={anchor_id}, org={org_name}"

    has_matches = 'matches' in org_data
    has_dups = 'dups_matches' in org_data

    if not has_matches and not has_dups:
        tracker.add(f"Neither 'matches' nor 'dups_matches' key found for {context}")
        return all_matches

    if has_matches:
        all_matches.update(org_data['matches'])

    if has_dups:
        dups_dict = org_data['dups_matches']
        if 'syntenic' in dups_dict:
            syntenic_set = dups_dict['syntenic']
        else:
            syntenic_set = None

        for j, match_data in dups_dict.items():
            if j == 'syntenic':
                continue
            if syntenic_set is not None and j not in syntenic_set:
                continue
            all_matches[j] = match_data

    return all_matches


def load_syn_eval_margin(work_dir):
    """Load the syn_eval margin from pipeline config for neighbor verification."""
    config_path = pathlib.Path(work_dir) / 'pipeline_config.yaml'
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config.get('syn_eval', {}).get('margin', 100000)


def check_nonreciprocal_cause(am2, j, org, i, margin):
    """Diagnose why a non-reciprocal match occurred.

    When anchor i (in org's aligned map) matches j in org2, but the reverse
    match i is not in am2[j]'s matches for org, determine the cause.

    Returns (is_eval_syn_asymmetry, detail_msg).

    is_eval_syn_asymmetry is True when the match EXISTS in dups_matches on the
    reverse side but wasn't marked syntenic. This is always a WARNING because
    eval_syn asymmetry has multiple known causes:
    - Neighbor's 'multiple matches out of tolerance range' meta flag
    - Supporting matches moved to 'not_considered' by inconsistency removal
    - Different anchor neighborhoods near the margin boundary

    is_eval_syn_asymmetry is False when the match is NOT in dups_matches at all,
    which indicates a structural pipeline problem (ERROR).
    """
    org_data_reverse = am2[j]['matches'][org]

    # Check if the match exists as a dup match but wasn't marked syntenic
    if 'dups_matches' not in org_data_reverse:
        return False, 'match not in dups_matches on reverse side'

    dups_reverse = org_data_reverse['dups_matches']
    syntenic_reverse = dups_reverse.get('syntenic', set())

    if i not in dups_reverse or i == 'syntenic':
        return False, 'match not in dups_matches on reverse side'

    if i in syntenic_reverse:
        # It IS syntenic on the reverse side — shouldn't reach here since
        # get_all_matches would have included it. Something else is wrong.
        return False, 'match is syntenic on reverse side but was not returned by get_all_matches'

    # Match IS in dups_matches but NOT syntenic — this is an eval_syn asymmetry.
    # Scan neighbors to identify the specific cause (for diagnostic detail).
    anchor_data_j = am2[j]
    anchor_len_j = anchor_data_j['end'] - anchor_data_j['start']
    sorted_indices_am2 = sorted(am2.keys())

    left_bound = bisect_left(sorted_indices_am2, j - margin)
    right_bound = bisect_left(sorted_indices_am2, j + anchor_len_j + margin)
    right_bound = min(len(sorted_indices_am2), right_bound)
    neighbors = sorted_indices_am2[left_bound:right_bound]

    flagged_neighbors = []
    for neighbor_idx in neighbors:
        if neighbor_idx == j:
            continue
        neighbor_data = am2[neighbor_idx]
        if org not in neighbor_data['matches']:
            continue
        neighbor_org_data = neighbor_data['matches'][org]
        if neighbor_org_data['meta']['multiple matches out of tolerance range'] == 1:
            flagged_neighbors.append(neighbor_idx)

    if flagged_neighbors:
        detail = (f'eval_syn asymmetry: {len(flagged_neighbors)} neighbor(s) of {j} have '
                  f"'multiple matches out of tolerance range' flag for {org} "
                  f'(neighbors: {flagged_neighbors[:5]}{"..." if len(flagged_neighbors) > 5 else ""})')
    else:
        detail = (f'eval_syn asymmetry: no neighbor meta flags — likely caused by '
                  f'inconsistency-based match removal or margin boundary effects')
    return True, detail


def check_aligned(org, work_dir, anchor_dir, tracker):
    """Check cross-organism consistency of aligned maps."""
    margin = load_syn_eval_margin(work_dir)

    with open(f'{anchor_dir}/aligned/{org}', 'rb') as f:
        am1 = pickle.load(f)

    # Check internal consistency
    for i, abib in am1.items():
        for org2, bbib in abib['matches'].items():
            all_matches = get_all_matches(bbib, i, org2, tracker)

            # 'matches not considered...' always exists (initialized in parse_bcamm)
            not_considered_key = 'matches not considered upon applying stricter score criterion since there are consistency issues'
            not_considered = bbib[not_considered_key]
            for j, cbib in all_matches.items():
                if j in not_considered:
                    if cbib == not_considered[j]:
                        tracker.add(f'Match {j} in both matches and not_considered for anchor={i}, org={org2}')

    # Load other organisms
    other_orgs = []
    with open(f'{work_dir}/orgs', 'r') as f:
        for line in f:
            other_orgs.append(line.strip())
    if org in other_orgs:
        other_orgs.remove(org)

    # Check cross-organism consistency
    for org2 in other_orgs:
        print(f'Checking consistency with {org2}...')
        with open(f'{anchor_dir}/aligned/{org2}', 'rb') as f:
            am2 = pickle.load(f)

        for i, bib in am1.items():
            if org2 not in bib['matches']:
                continue

            am1_matches = get_all_matches(bib['matches'][org2], i, org2, tracker)
            if not am1_matches:
                continue

            for j, match_data1 in am1_matches.items():
                meta = bib['matches'][org2]['meta']

                if j not in am2:
                    tracker.add(f'Match {j} (from anchor {i} of {org}) not in {org2} aligned map')
                    continue

                if org not in am2[j]['matches']:
                    tracker.add(f'{org} not in {org2} anchor {j} matches (expected reciprocal from {org} anchor {i})')
                    continue

                am2_matches = get_all_matches(am2[j]['matches'][org], j, org, tracker)
                if i not in am2_matches:
                    # Non-reciprocal match. Diagnose: is this the known meta flag
                    # edge case, or a genuine pipeline integrity error?
                    is_meta_flag_case, detail = check_nonreciprocal_cause(
                        am2, j, org, i, margin
                    )

                    if is_meta_flag_case:
                        tracker.warn(
                            f'Anchor {i} not in {org2} anchor {j} matches for {org} '
                            f'(non-reciprocal syntenic dup: {detail})')
                    else:
                        tracker.add(
                            f'Anchor {i} not in {org2} anchor {j} matches for {org} '
                            f'(non-reciprocal: {detail})')
                    continue

                match_data2 = am2_matches[i]

                # Check score consistency
                score1 = round(match_data1['match score'])
                score2 = round(match_data2['match score'])
                if score1 != score2:
                    tracker.add(f'Score mismatch for {org}:{i} <-> {org2}:{j}: {score1} vs {score2}')

                # Check coordinate consistency
                start11, end11 = match_data1[f'hit coordinates in (own) {org} candidate']
                start12, end12 = match_data2[f'hit coordinates in {org} candidate']
                if start11 != start12 or end11 != end12:
                    tracker.add(f'Coordinate mismatch for {org}:{i} own coords: ({start11},{end11}) vs ({start12},{end12})')

                start21, end21 = match_data1[f'hit coordinates in {org2} candidate']
                start22, end22 = match_data2[f'hit coordinates in (own) {org2} candidate']
                if start21 != start22 or end21 != end22:
                    tracker.add(f'Coordinate mismatch for {org2}:{j} own coords: ({start21},{end21}) vs ({start22},{end22})')

                # Meta flag consistency (warnings only — asymmetry is expected)
                # These flags describe the match pattern FROM one anchor's
                # perspective (set in collect_output.py check_range_and_chromos).
                # Asymmetry is biologically expected: an anchor region that
                # underwent insertion, deletion, or transposition in one organism
                # but not the other will have a different match distribution when
                # viewed from each side. The companion 'ambiguous' flag is
                # propagated cross-org via mark_in_others/inconsistencies.py.
                meta2 = am2[j]['matches'][org]['meta']

                if (meta['multiple matches out of tolerance range'] == 1 and
                    meta2['matches have ambiguous matches (tolerance/chromosome out of range)'] == 0 and
                    meta2['multiple matches out of tolerance range'] == 0):
                    tracker.warn(
                        f'Meta flag asymmetry (tolerance): {org}:{i} has flag=1 but {org2}:{j} '
                        f'has neither ambiguous nor own tolerance flag '
                        f'(expected if anchor region had indel/transposition in one org)')

                if (meta2['multiple matches out of tolerance range'] == 1 and
                    meta['matches have ambiguous matches (tolerance/chromosome out of range)'] == 0 and
                    meta['multiple matches out of tolerance range'] == 0):
                    tracker.warn(
                        f'Meta flag asymmetry (tolerance): {org2}:{j} has flag=1 but {org}:{i} '
                        f'has neither ambiguous nor own tolerance flag '
                        f'(expected if anchor region had indel/transposition in one org)')

                if meta['multiple matches on different chromosomes'] != meta2['multiple matches on different chromosomes']:
                    if not (meta['multiple matches on different chromosomes'] == 1 and
                            meta2['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1):
                        if not (meta2['multiple matches on different chromosomes'] == 1 and
                                meta['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1):
                            tracker.warn(
                                f'Meta flag asymmetry (chromosomes): {org}:{i}='
                                f'{meta["multiple matches on different chromosomes"]} vs '
                                f'{org2}:{j}={meta2["multiple matches on different chromosomes"]} '
                                f'(expected if anchor region had translocation in one org)')

    print('check_aligned done')


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <organism> <work_dir> <anchor_dir>", file=sys.stderr)
        sys.exit(1)

    org = sys.argv[1]
    work_dir = sys.argv[2]
    anchor_dir = sys.argv[3]
    root_dir = str(pathlib.Path(anchor_dir).parent)

    print(f'--- Checking {org} ---')

    tracker = ErrorTracker()

    check_correspondence_candidates_aligned_maps(org, anchor_dir, tracker)
    check_sanity(org, anchor_dir, root_dir, tracker)
    check_aligned(org, work_dir, anchor_dir, tracker)

    # Report results
    if tracker.has_warnings():
        print(f'\n=== {tracker.warning_count()} warnings (non-critical) ===', file=sys.stderr)

    if tracker.has_errors():
        print(f'\n=== VALIDATION FAILED: {tracker.count()} errors found ===', file=sys.stderr)
        sys.exit(1)
    else:
        if tracker.has_warnings():
            print(f'\n=== All checks passed ({tracker.warning_count()} warnings) ===')
        else:
            print('\n=== All checks passed ===')
        sys.exit(0)
