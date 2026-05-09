#!/usr/bin/env python3
"""
Verify aligned->succinct conversion and produce detailed filtering statistics.

This is an INDEPENDENT verification — it derives expected behavior from first
principles by examining the aligned data, NOT by copying the converter's code.

The specification for what should survive into succinct:
  1. An (anchor, org2) pair is kept if:
     a) It has regular matches OR reciprocal syntenic dups
     b) Neither ambiguity flag is set to 1 in meta
  2. Syntenic dups are reciprocal if for every position j in the syntenic set:
     - j exists as an anchor in org2's aligned map
     - org2's anchor j has a syntenic dup back to anchor_idx in org1
     If ANY j fails this check, the whole (anchor, org2) relationship is
     marked ambiguous (the converter does this in memory before filtering).
  3. An anchor is kept if at least one (anchor, org2) pair is kept.
  4. For kept pairs, only the first match (sorted by position) is stored.
  5. Regular matches take priority over syntenic dups.

Usage:
  python verify_succinct.py --aligned /path/to/aligned --succinct /path/to/succinct --cores 4
"""

import sys
import os
import argparse
import pickle
from pathlib import Path
from multiprocessing import Pool


def is_syntenic_dup_reciprocal(org, anchor_idx, org2, syntenic_set, partner_map):
    """
    Independently check whether syntenic dup markings for (anchor_idx -> org2)
    are reciprocal by examining the partner's aligned data.

    A syntenic relationship is reciprocal iff for EVERY position j in
    syntenic_set:
      - j is an anchor in partner_map
      - partner_map[j] has a match back to org
      - that match has dups_matches with a 'syntenic' set containing anchor_idx

    Returns True if reciprocal, False if any j fails.
    """
    for j in syntenic_set:
        # j must exist as an anchor in the partner organism
        if j not in partner_map:
            return False

        partner_anchor = partner_map[j]

        # partner must have matches back to our organism
        if 'matches' not in partner_anchor or org not in partner_anchor['matches']:
            return False

        partner_org_data = partner_anchor['matches'][org]

        # partner must have dups_matches with a syntenic set
        if 'dups_matches' not in partner_org_data:
            return False
        if 'syntenic' not in partner_org_data['dups_matches']:
            return False

        # our anchor_idx must be in partner's syntenic set
        if anchor_idx not in partner_org_data['dups_matches']['syntenic']:
            return False

    return True


def pair_should_be_kept(org, anchor_idx, org2, bib2, partner_maps):
    """
    Independently determine whether an (anchor, org2) pair should survive
    into the succinct output.

    The filtering pipeline is:
      1. No match data (no regular matches AND no syntenic dups) -> 'no_data'
      2. Ambiguity flags from meta (on-disk) -> 'ambiguous'
      3. Non-reciprocal syntenic dups (checked against partner data) -> 'ambiguous'
      4. Otherwise -> 'kept'

    Returns: ('kept', None) | ('no_data', None) | ('ambiguous', reason_str)
    """
    has_regular = 'matches' in bib2 and len(bib2['matches']) > 0
    has_syntenic = ('dups_matches' in bib2 and
                    'syntenic' in bib2['dups_matches'] and
                    len(bib2['dups_matches']['syntenic']) > 0)

    # No match data at all (neither regular matches nor syntenic dups exist)
    if not has_regular and not has_syntenic:
        return 'no_data', None

    # Start with on-disk ambiguity flags (set by earlier pipeline stages)
    meta = bib2.get('meta', {})
    tol = meta.get('multiple matches out of tolerance range', 0)
    ambig = meta.get('matches have ambiguous matches (tolerance/chromosome out of range)', 0)

    # Independently check syntenic dup reciprocity against partner data.
    # Non-reciprocal syntenic dups cause the ambiguity flag to be set.
    if has_syntenic and org2 in partner_maps:
        if not is_syntenic_dup_reciprocal(
                org, anchor_idx, org2,
                bib2['dups_matches']['syntenic'],
                partner_maps[org2]):
            ambig = 1

    if tol == 1 or ambig == 1:
        if tol == 1 and ambig == 1:
            return 'ambiguous', 'both'
        elif tol == 1:
            return 'ambiguous', 'tolerance'
        else:
            return 'ambiguous', 'chr'

    return 'kept', None


def compute_expected_match(org, bib2):
    """
    Independently compute what the succinct match entry should be for
    a kept (anchor, org2) pair.

    Returns dict of {j: (score, (start, end), orientation)} — should have
    exactly one entry (first match sorted by position).
    """
    expected = {}

    has_regular = 'matches' in bib2 and len(bib2['matches']) > 0
    has_syntenic = ('dups_matches' in bib2 and
                    'syntenic' in bib2['dups_matches'] and
                    len(bib2['dups_matches']['syntenic']) > 0)

    # Regular matches take priority — take only the first (sorted by position)
    if has_regular:
        for j, bib3 in sorted(bib2['matches'].items()):
            ori = 'reverse' if bib3['match is on other strand in other genome'] else 'forward'
            expected[j] = (
                int(bib3['match score']),
                tuple(bib3[f'hit coordinates in (own) {org} candidate']),
                ori
            )
            break

    # If no regular matches, try syntenic dups — first syntenic one sorted by position
    if not expected and has_syntenic:
        syn = bib2['dups_matches']['syntenic']
        dups_items = [(k, v) for k, v in bib2['dups_matches'].items() if k != 'syntenic']
        for j, bib3 in sorted(dups_items):
            if j not in syn:
                continue
            ori = 'reverse' if bib3['match is on other strand in other genome'] else 'forward'
            expected[j] = (
                int(bib3['match score']),
                tuple(bib3[f'hit coordinates in (own) {org} candidate']),
                ori
            )
            break

    return expected


def verify_one(args):
    """Verify and collect stats for a single organism."""
    org, aligned_dir, succinct_dir, suffix = args

    aligned_path = Path(aligned_dir) / (org + suffix)
    succinct_path = Path(succinct_dir) / org

    with open(aligned_path, 'rb') as f:
        aligned = pickle.load(f)

    with open(succinct_path, 'rb') as f:
        succinct = pickle.load(f)

    # Load partner aligned maps for reciprocity checking (READ-ONLY).
    # Only load partners that actually appear in our matches and have
    # syntenic dups, to avoid unnecessary I/O.
    needed_partners = set()
    for i, bib1 in aligned.items():
        if 'matches' not in bib1:
            continue
        for org2, bib2 in bib1['matches'].items():
            if ('dups_matches' in bib2 and
                    'syntenic' in bib2['dups_matches'] and
                    len(bib2['dups_matches']['syntenic']) > 0):
                needed_partners.add(org2)

    partner_maps = {}
    for org2 in needed_partners:
        partner_path = Path(aligned_dir) / (org2 + suffix)
        if partner_path.exists():
            with open(partner_path, 'rb') as f:
                partner_maps[org2] = pickle.load(f)

    errors = []

    # ── Per-(anchor, org2) pair stats ──
    pairs_total = 0
    pairs_no_regular_matches = 0
    pairs_no_syntenic_dups = 0
    pairs_has_nonsyn_dups_only = 0
    pairs_ambig_tolerance = 0
    pairs_ambig_chr = 0
    pairs_ambig_both = 0
    pairs_skipped_no_data = 0
    pairs_skipped_ambiguity = 0
    pairs_kept = 0
    pairs_nonreciprocal = 0

    # ── Per-anchor stats ──
    anchors_total = len(aligned)
    anchors_no_org2_entries = 0
    anchors_dropped_all_no_data = 0
    anchors_dropped_all_ambiguous = 0
    anchors_dropped_mixed = 0
    anchors_kept = 0
    anchors_should_be_kept = 0

    # ── Expected succinct: independently compute what should be in succinct ──
    expected_succinct = {}

    for i, bib1 in aligned.items():
        if 'matches' not in bib1 or len(bib1['matches']) == 0:
            anchors_no_org2_entries += 1
            continue

        pair_reasons = []
        expected_entry_matches = {}

        for org2, bib2 in bib1['matches'].items():
            pairs_total += 1

            # Diagnostic counters
            has_regular = 'matches' in bib2 and len(bib2['matches']) > 0
            has_dups = ('dups_matches' in bib2 and
                        len([k for k in bib2['dups_matches'] if k != 'syntenic']) > 0)
            has_syntenic = ('dups_matches' in bib2 and
                            'syntenic' in bib2['dups_matches'] and
                            len(bib2['dups_matches']['syntenic']) > 0)

            if not has_regular:
                pairs_no_regular_matches += 1
            if not has_syntenic:
                pairs_no_syntenic_dups += 1
            if has_dups and not has_syntenic:
                pairs_has_nonsyn_dups_only += 1

            # Check reciprocity independently
            if has_syntenic and org2 in partner_maps:
                if not is_syntenic_dup_reciprocal(
                        org, i, org2,
                        bib2['dups_matches']['syntenic'],
                        partner_maps[org2]):
                    pairs_nonreciprocal += 1

            # Determine expected outcome for this pair
            outcome, reason = pair_should_be_kept(org, i, org2, bib2, partner_maps)

            if outcome == 'no_data':
                pairs_skipped_no_data += 1
                pair_reasons.append('no_data')
            elif outcome == 'ambiguous':
                pairs_skipped_ambiguity += 1
                if reason == 'both':
                    pairs_ambig_both += 1
                elif reason == 'tolerance':
                    pairs_ambig_tolerance += 1
                else:
                    pairs_ambig_chr += 1
                pair_reasons.append('ambiguous')
            else:
                pairs_kept += 1
                pair_reasons.append('kept')
                expected_entry_matches[org2] = compute_expected_match(org, bib2)

        # Anchor-level outcome
        has_any_kept = 'kept' in pair_reasons
        if has_any_kept:
            anchors_should_be_kept += 1
            anchors_kept += 1
            expected_succinct[i] = {
                'chromosome': bib1['chromosome'],
                'start': bib1['start'],
                'end': bib1['end'],
                'matches': expected_entry_matches,
            }
        else:
            reasons_set = set(pair_reasons)
            if reasons_set == {'no_data'}:
                anchors_dropped_all_no_data += 1
            elif reasons_set == {'ambiguous'}:
                anchors_dropped_all_ambiguous += 1
            elif reasons_set <= {'no_data', 'ambiguous'}:
                anchors_dropped_mixed += 1
            else:
                errors.append(f"Anchor {i}: unexpected pair reasons {pair_reasons}")

    # Free partner maps
    del partner_maps

    # ── Cross-check: compare expected vs actual succinct ──
    succinct_missing = []   # in expected but not in actual succinct
    succinct_extra = []     # in actual succinct but not in expected
    value_mismatches = []

    for i in expected_succinct:
        if i not in succinct:
            succinct_missing.append(i)

    for i in succinct:
        if i not in expected_succinct:
            succinct_extra.append(i)

    # For anchors present in both, verify values match
    for i in succinct:
        if i not in expected_succinct:
            continue
        if i not in aligned:
            errors.append(f"Anchor {i}: in succinct but not in aligned at all")
            continue

        succ_entry = succinct[i]
        exp_entry = expected_succinct[i]

        # Check top-level fields
        if succ_entry['chromosome'] != exp_entry['chromosome']:
            value_mismatches.append(
                f"Anchor {i}: chromosome '{succ_entry['chromosome']}' "
                f"!= expected '{exp_entry['chromosome']}'")
        if succ_entry['start'] != exp_entry['start']:
            value_mismatches.append(
                f"Anchor {i}: start {succ_entry['start']} != expected {exp_entry['start']}")
        if succ_entry['end'] != exp_entry['end']:
            value_mismatches.append(
                f"Anchor {i}: end {succ_entry['end']} != expected {exp_entry['end']}")

        # Check matches: org2 set
        succ_orgs = set(succ_entry['matches'].keys())
        exp_orgs = set(exp_entry['matches'].keys())
        for org2 in exp_orgs - succ_orgs:
            value_mismatches.append(f"Anchor {i}: org2 {org2} expected but missing from succinct")
        for org2 in succ_orgs - exp_orgs:
            value_mismatches.append(f"Anchor {i}: org2 {org2} in succinct but not expected")

        # Check match values for shared org2s
        for org2 in succ_orgs & exp_orgs:
            succ_matches = succ_entry['matches'][org2]
            exp_matches = exp_entry['matches'][org2]

            succ_keys = set(succ_matches.keys())
            exp_keys = set(exp_matches.keys())

            for j in exp_keys - succ_keys:
                value_mismatches.append(
                    f"Anchor {i}, org2 {org2}, j={j}: expected but missing from succinct")
            for j in succ_keys - exp_keys:
                value_mismatches.append(
                    f"Anchor {i}, org2 {org2}, j={j}: in succinct but not expected")

            for j in succ_keys & exp_keys:
                s_score, s_coords, s_ori = succ_matches[j]
                e_score, e_coords, e_ori = exp_matches[j]
                if s_score != e_score:
                    value_mismatches.append(
                        f"Anchor {i}, org2 {org2}, j={j}: score {s_score} != expected {e_score}")
                if s_coords != e_coords:
                    value_mismatches.append(
                        f"Anchor {i}, org2 {org2}, j={j}: coords {s_coords} != expected {e_coords}")
                if s_ori != e_ori:
                    value_mismatches.append(
                        f"Anchor {i}, org2 {org2}, j={j}: ori '{s_ori}' != expected '{e_ori}'")

    return {
        'org': org,
        'anchors_total': anchors_total,
        'anchors_in_succinct': len(succinct),
        'anchors_no_org2_entries': anchors_no_org2_entries,
        'anchors_dropped_all_no_data': anchors_dropped_all_no_data,
        'anchors_dropped_all_ambiguous': anchors_dropped_all_ambiguous,
        'anchors_dropped_mixed': anchors_dropped_mixed,
        'anchors_kept': anchors_kept,
        'anchors_should_be_kept': anchors_should_be_kept,
        'pairs_total': pairs_total,
        'pairs_skipped_no_data': pairs_skipped_no_data,
        'pairs_skipped_ambiguity': pairs_skipped_ambiguity,
        'pairs_ambig_tolerance': pairs_ambig_tolerance,
        'pairs_ambig_chr': pairs_ambig_chr,
        'pairs_ambig_both': pairs_ambig_both,
        'pairs_no_regular_matches': pairs_no_regular_matches,
        'pairs_has_nonsyn_dups_only': pairs_has_nonsyn_dups_only,
        'pairs_kept': pairs_kept,
        'pairs_nonreciprocal': pairs_nonreciprocal,
        'succinct_missing': succinct_missing,
        'succinct_extra': succinct_extra,
        'value_mismatches': value_mismatches,
        'errors': errors,
    }


def print_org_report(s):
    """Print detailed report for one organism."""
    org = s['org']
    total = s['anchors_total']
    kept = s['anchors_in_succinct']
    dropped = total - kept
    pct_kept = 100 * kept / total if total > 0 else 0

    print(f"\n{'='*70}")
    print(f"  {org}")
    print(f"{'='*70}")

    print(f"\n  ANCHOR-LEVEL ({total} total):")
    print(f"    Kept in succinct:               {kept:>8}  ({pct_kept:.1f}%)")
    print(f"    Dropped:                        {dropped:>8}  ({100-pct_kept:.1f}%)")
    if dropped > 0:
        print(f"      No partner data at all:       {s['anchors_no_org2_entries']:>8}")
        print(f"      All pairs had no match data:  {s['anchors_dropped_all_no_data']:>8}")
        print(f"      All pairs were ambiguous:     {s['anchors_dropped_all_ambiguous']:>8}")
        print(f"      Mix (no data + ambiguous):    {s['anchors_dropped_mixed']:>8}")

    pt = s['pairs_total']
    print(f"\n  PAIR-LEVEL (anchor, org2) ({pt} total):")
    if pt > 0:
        print(f"    Kept:                           {s['pairs_kept']:>8}  ({100*s['pairs_kept']/pt:.1f}%)")
        print(f"    Skipped - no match data:        {s['pairs_skipped_no_data']:>8}  ({100*s['pairs_skipped_no_data']/pt:.1f}%)")
        print(f"      (no regular matches):         {s['pairs_no_regular_matches']:>8}")
        print(f"      (has non-syn dups only):      {s['pairs_has_nonsyn_dups_only']:>8}")
        print(f"    Skipped - ambiguity:            {s['pairs_skipped_ambiguity']:>8}  ({100*s['pairs_skipped_ambiguity']/pt:.1f}%)")
        print(f"      tolerance out of range:       {s['pairs_ambig_tolerance']:>8}")
        print(f"      chr/ambiguous out of range:   {s['pairs_ambig_chr']:>8}")
        print(f"      both flags:                   {s['pairs_ambig_both']:>8}")
        print(f"    Non-reciprocal syntenic dups:    {s['pairs_nonreciprocal']:>8}")

    # Verification
    ok = True
    if s['succinct_missing']:
        print(f"\n  VERIFICATION FAILURE: {len(s['succinct_missing'])} anchors should be in succinct but are missing!")
        for x in s['succinct_missing'][:5]:
            print(f"    - anchor {x}")
        if len(s['succinct_missing']) > 5:
            print(f"    ... and {len(s['succinct_missing'])-5} more")
        ok = False
    if s['succinct_extra']:
        print(f"\n  VERIFICATION FAILURE: {len(s['succinct_extra'])} anchors in succinct but should not be!")
        for x in s['succinct_extra'][:5]:
            print(f"    - anchor {x}")
        if len(s['succinct_extra']) > 5:
            print(f"    ... and {len(s['succinct_extra'])-5} more")
        ok = False
    if s['value_mismatches']:
        print(f"\n  VALUE MISMATCHES: {len(s['value_mismatches'])}")
        for x in s['value_mismatches'][:5]:
            print(f"    - {x}")
        if len(s['value_mismatches']) > 5:
            print(f"    ... and {len(s['value_mismatches'])-5} more")
        ok = False
    if s['errors']:
        print(f"\n  ERRORS: {len(s['errors'])}")
        for x in s['errors'][:5]:
            print(f"    - {x}")
        ok = False
    if s['anchors_in_succinct'] != s['anchors_should_be_kept']:
        print(f"\n  COUNT MISMATCH: succinct has {s['anchors_in_succinct']} but expected {s['anchors_should_be_kept']}")
        ok = False
    if ok:
        print(f"\n  VERIFICATION: OK")

    return ok


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Verify aligned->succinct conversion and produce filtering statistics.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--aligned', required=True, type=str,
                        help='Path to aligned/ directory')
    parser.add_argument('--succinct', required=True, type=str,
                        help='Path to succinct/ directory')
    parser.add_argument('--cores', type=int, default=1,
                        help='Number of parallel workers')
    parser.add_argument('--suffix', type=str, default='_with_syn_eval',
                        help='Filename suffix on aligned pickles')
    args = parser.parse_args()

    aligned_dir = os.path.abspath(args.aligned)
    succinct_dir = os.path.abspath(args.succinct)
    suffix = args.suffix

    if not os.path.isdir(aligned_dir):
        print(f"ERROR: {aligned_dir} not found", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(succinct_dir):
        print(f"ERROR: {succinct_dir} not found", file=sys.stderr)
        sys.exit(1)

    # List organisms from succinct dir (those are the ones we verify)
    orgs_to_verify = sorted([
        f for f in os.listdir(succinct_dir)
        if os.path.isfile(os.path.join(succinct_dir, f)) and not f.startswith('.')
    ])

    # Build full org list from aligned dir (for reciprocity context)
    all_aligned_files = sorted([
        f for f in os.listdir(aligned_dir)
        if os.path.isfile(os.path.join(aligned_dir, f)) and not f.startswith('.')
    ])
    if suffix:
        all_orgs = [f[:-len(suffix)] for f in all_aligned_files if f.endswith(suffix)]
    else:
        all_orgs = all_aligned_files

    # Verify each succinct file has a matching aligned file
    missing = [o for o in orgs_to_verify if not (Path(aligned_dir) / (o + suffix)).exists()]
    if missing:
        print(f"ERROR: {len(missing)} succinct files have no aligned counterpart:", file=sys.stderr)
        for m in missing[:5]:
            print(f"  {m} -> expected {m}{suffix}", file=sys.stderr)
        sys.exit(1)

    print(f"Verifying {len(orgs_to_verify)} organisms: {aligned_dir} vs {succinct_dir}")
    if len(all_orgs) != len(orgs_to_verify):
        print(f"  Note: aligned dir has {len(all_orgs)} organisms, succinct has {len(orgs_to_verify)}")

    tasks = [(org, aligned_dir, succinct_dir, suffix) for org in orgs_to_verify]
    all_stats = []

    if args.cores == 1:
        for task in tasks:
            all_stats.append(verify_one(task))
    else:
        with Pool(args.cores) as pool:
            for result in pool.imap_unordered(verify_one, tasks):
                all_stats.append(result)

    # Print per-org reports
    all_ok = True
    for s in sorted(all_stats, key=lambda x: x['org']):
        ok = print_org_report(s)
        if not ok:
            all_ok = False

    # Print summary
    tot_aligned = sum(s['anchors_total'] for s in all_stats)
    tot_succinct = sum(s['anchors_in_succinct'] for s in all_stats)
    tot_pairs = sum(s['pairs_total'] for s in all_stats)
    tot_pairs_kept = sum(s['pairs_kept'] for s in all_stats)
    tot_no_org2 = sum(s['anchors_no_org2_entries'] for s in all_stats)
    tot_no_data = sum(s['anchors_dropped_all_no_data'] for s in all_stats)
    tot_ambig = sum(s['anchors_dropped_all_ambiguous'] for s in all_stats)
    tot_mixed = sum(s['anchors_dropped_mixed'] for s in all_stats)

    print(f"\n{'='*70}")
    print(f"  SUMMARY ACROSS ALL {len(all_stats)} ORGANISMS")
    print(f"{'='*70}")
    print(f"  Total anchors in aligned:     {tot_aligned:>10}")
    print(f"  Total anchors in succinct:    {tot_succinct:>10}  ({100*tot_succinct/tot_aligned:.1f}%)")
    print(f"  Dropped anchors:              {tot_aligned-tot_succinct:>10}  ({100*(tot_aligned-tot_succinct)/tot_aligned:.1f}%)")
    print(f"    No partner data at all:     {tot_no_org2:>10}")
    print(f"    All pairs no match data:    {tot_no_data:>10}")
    print(f"    All pairs ambiguous:        {tot_ambig:>10}")
    print(f"    Mixed:                      {tot_mixed:>10}")
    print(f"  Total (anchor,org2) pairs:    {tot_pairs:>10}")
    print(f"  Pairs kept:                   {tot_pairs_kept:>10}  ({100*tot_pairs_kept/tot_pairs:.1f}%)" if tot_pairs > 0 else "")
    print(f"\n  Overall verification: {'PASSED' if all_ok else 'FAILED'}")
