#!/usr/bin/env python3
"""
Finalize aligned maps: run syntenic evaluation and create succinct output.

Two modes controlled by the 4th argument:

  eval_syn:
    Evaluate which dups_matches are syntenic (supported by flanking anchors).
    Updates aligned/{org} IN-PLACE with 'syntenic' sets in dups_matches.
    Per-organism, no barrier needed (only reads own aligned file).

  finalize_succinct:
    Check eval_syn reciprocity across ALL organisms (requires barrier),
    set 'matches have ambiguous matches' flag IN MEMORY for one-sided syntenic
    dups (not written to disk — avoids race condition with parallel partner reads),
    then create aligned_succinct/{org} with simplified tuple format.

The aligned_succinct format is identical to the former compressed_maps_multis_to_one:
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

Usage: python finalize_aligned.py <organism> <work_dir> <root_dir> <mode>
  mode = "eval_syn"          -> run eval_syntenic_dups + save aligned (no succinct)
  mode = "finalize_succinct" -> run reciprocity check + set flags + create succinct
"""

import sys
import pickle
import yaml
from pathlib import Path
from bisect import bisect_left


# Part 1: syntenic evaluation.

def load_config(work_dir):
    """Load pipeline configuration from work_dir/pipeline_config.yaml."""
    config_path = Path(work_dir) / 'pipeline_config.yaml'
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Get syn_eval parameters with defaults
    syn_eval = config.get('syn_eval', {})
    margin = syn_eval.get('margin', 100000)
    min_score = syn_eval.get('min_score', 500)

    return margin, min_score


def check_required_keys(data, keys, context):
    """
    Check that all required keys exist in data.
    Raises KeyError with descriptive message if any are missing.
    """
    missing = [k for k in keys if k not in data]
    if missing:
        msg = f"MISSING KEYS in {context}: {missing}. Available keys: {list(data.keys())}"
        print(f"ERROR: {msg}", file=sys.stderr)
        raise KeyError(msg)


def eval_syntenic_dups(org, aligned_map, reference_map, margin, min_score):
    """
    Evaluate which dups_matches are syntenic for a single organism.

    A dup match at position j is syntenic if neighboring anchors have regular
    matches that demonstrate two-sided flanking support in BOTH coordinate systems:

    1. Source two-sidedness: neighbors on both LEFT and RIGHT of anchor_idx
       (in this organism's coordinates) contribute qualifying matches
    2. Target two-sidedness: those matches land on both LEFT and RIGHT of dup_pos
       (in the partner organism's coordinates)
    3. All four accumulated scores >= min_score

    The same margin value is used for BOTH distance constraints:
    - Source side: |neighbor_idx - anchor_idx| < margin
    - Target side: |match_pos - dup_pos| < margin

    RECIPROCITY GUARANTEE: Using the same margin for both distance checks ensures
    the set of qualifying supporting matches is identical regardless of evaluation
    direction (org1->org2 vs org2->org1). The four scores are the same values
    (just relabeled), so the syntenic condition produces identical results from
    both directions. A rare edge case exists when the 'multiple matches out of
    tolerance range' meta flag differs between the two anchors in a supporting
    match pair, which can cause one direction to exclude a match the other includes.
    This is biologically correct (unreliable neighbors shouldn't be trusted as
    supporting evidence) but means perfect reciprocity is not 100% guaranteed.
    checks.py handles this case with a warning rather than an error.

    Modifies aligned_map IN-PLACE, adding 'syntenic' sets to dups_matches.
    Returns (syntenic_count, total_dups_matches) for statistics.
    """

    # Sorted list of anchor indices for binary search
    sorted_indices = sorted(list(reference_map.keys()))

    # Statistics
    total_dups_matches = 0
    syntenic_count = 0

    for anchor_idx, anchor_data in aligned_map.items():
        # Check required keys at anchor level
        check_required_keys(anchor_data, ['start', 'end', 'matches'], f"anchor {anchor_idx}")

        anchor_start = anchor_data['start']
        anchor_end = anchor_data['end']
        anchor_len = anchor_end - anchor_start

        for org2, org2_data in anchor_data['matches'].items():
            # Check required keys at org2 level
            check_required_keys(org2_data, ['meta'], f"anchor {anchor_idx} matches[{org2}]")

            # Skip if no dups_matches
            if 'dups_matches' not in org2_data:
                continue

            dups_matches = org2_data['dups_matches']

            # Find neighboring anchors within margin
            left_bound = bisect_left(sorted_indices, anchor_idx - margin)
            right_bound = bisect_left(sorted_indices, anchor_idx + anchor_len + margin)
            right_bound = min(len(sorted_indices), right_bound)
            local_indices = sorted_indices[left_bound:right_bound]

            # Initialize syntenic set for this org2 if not present
            if 'syntenic' not in dups_matches:
                dups_matches['syntenic'] = set()

            # Evaluate each dup match position
            for dup_pos, dup_match_data in list(dups_matches.items()):
                if dup_pos == 'syntenic':
                    continue  # Skip the metadata key

                total_dups_matches += 1

                # Accumulate support scores in BOTH coordinate systems.
                # Source side: position of neighbor relative to anchor_idx (this genome)
                # Target side: position of match relative to dup_pos (partner genome)
                # Both must show two-sided support for the match to be syntenic.
                source_left_score = 0
                source_right_score = 0
                target_left_score = 0
                target_right_score = 0

                for neighbor_idx in local_indices:
                    if neighbor_idx == anchor_idx:
                        continue

                    neighbor_data = reference_map[neighbor_idx]

                    # Check if neighbor has matches to org2
                    if org2 not in neighbor_data['matches']:
                        continue

                    neighbor_org2 = neighbor_data['matches'][org2]

                    # Skip if multiple matches out of tolerance
                    # 'meta' always exists when org2 entry exists (created in parse_bcamm)
                    if neighbor_org2['meta']['multiple matches out of tolerance range'] == 1:
                        continue

                    # Check regular matches (not dups_matches)
                    if 'matches' not in neighbor_org2:
                        continue

                    for match_pos, match_data in neighbor_org2['matches'].items():
                        # Exclude matches involving either partner being evaluated
                        if match_pos == dup_pos:
                            continue

                        # Check if match is within margin of our dup position
                        # (same margin used for source and target distance ensures reciprocity)
                        if abs(match_pos - dup_pos) < margin:
                            # 'match score' always exists in match entries (created in parse_bcamm)
                            score = match_data['match score']

                            # Source side: is neighbor left or right of our anchor?
                            if neighbor_idx < anchor_idx:
                                source_left_score += score
                            else:
                                source_right_score += score

                            # Target side: is match left or right of dup position?
                            if match_pos < dup_pos:
                                target_left_score += score
                            else:
                                target_right_score += score

                # Mark as syntenic if supported from both sides in BOTH coordinate systems
                if (source_left_score >= min_score and source_right_score >= min_score and
                        target_left_score >= min_score and target_right_score >= min_score):
                    dups_matches['syntenic'].add(dup_pos)
                    syntenic_count += 1

    return syntenic_count, total_dups_matches


# Part 2: create succinct output.

def check_forward_ambiguity(match_bib):
    """Check if match has forward ambiguity (should be filtered out)."""
    if match_bib['meta']['multiple matches out of tolerance range'] == 1 or \
       match_bib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
        return True
    return False


def create_succinct(org, aligned_map):
    """
    Create succinct format from aligned map.

    This produces the EXACT format as compressed_maps_multis_to_one:
    - Only includes anchors with at least one valid match
    - Only includes first match per org2 (sorted by position)
    - For dups_matches, only includes syntenic ones and only if no regular matches
    - Each match becomes tuple: (score_int, (start, end), orientation_str)
    - No 'syntenic' key in output

    Returns the succinct dictionary.
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
            # Skip entries without matches or with forward ambiguity
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

            # Try regular matches first - take only the first one (sorted by position)
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
                # Need to iterate without modifying original
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


# Part 3: reciprocity check (cross-organism, requires barrier).

def check_eval_syn_reciprocity(org, aligned_map, anchor_dir, work_dir):
    """
    Check that syntenic dup markings are reciprocal across organisms.

    For each anchor with syntenic dups to org2, verify that the partner's
    aligned map also marks the reverse as syntenic. If ANY syntenic dup in
    a (anchor, org2) relationship is NOT reciprocal, set the ambiguous flag
    on the org2_data['meta'] for that relationship. This causes
    check_forward_ambiguity() to skip the entire pairwise relationship
    in succinct creation.

    Args:
        org: This organism's name
        aligned_map: This organism's aligned map (modified in-place)
        anchor_dir: Path to anchors/ directory
        work_dir: Working directory (for reading orgs file)

    Returns:
        Number of (anchor, org2) relationships flagged as ambiguous.
    """
    # Read orgs file to get all partner organisms
    orgs_path = Path(work_dir) / 'orgs'
    with open(orgs_path, 'r') as f:
        all_orgs = [line.strip() for line in f if line.strip()]

    other_orgs = [o for o in all_orgs if o != org]
    flagged_count = 0
    aligned_dir = Path(anchor_dir) / 'aligned'

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

        # Load partner's aligned map
        partner_path = aligned_dir / org2
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
                # Check: does partner_map[j] exist?
                if j not in partner_map:
                    relationship_nonreciprocal = True
                    break

                # Check: does partner have matches to org?
                if org not in partner_map[j]['matches']:
                    relationship_nonreciprocal = True
                    break

                partner_org_data = partner_map[j]['matches'][org]

                # Check: does partner have dups_matches?
                if 'dups_matches' not in partner_org_data:
                    relationship_nonreciprocal = True
                    break

                partner_dups = partner_org_data['dups_matches']

                # Check: does partner have syntenic set?
                if 'syntenic' not in partner_dups:
                    relationship_nonreciprocal = True
                    break

                # Check: is anchor_idx in partner's syntenic set?
                if anchor_idx not in partner_dups['syntenic']:
                    relationship_nonreciprocal = True
                    break

            if relationship_nonreciprocal:
                org2_data['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] = 1
                flagged_count += 1

        # Free partner map
        del partner_map

    return flagged_count


def run_eval_syn(org, work_dir, root_dir):
    """
    Mode 1: Run syntenic evaluation only.

    1. Load aligned/{org}
    2. Run eval_syn logic to mark syntenic dups (updates in-place)
    3. Save updated aligned/{org}

    No succinct creation — that happens in finalize_succinct after the barrier.
    """
    margin, min_score = load_config(work_dir)

    aligned_path = Path(root_dir) / 'anchors' / 'aligned' / org

    if not aligned_path.exists():
        print(f"ERROR: Aligned file not found: {aligned_path}", file=sys.stderr)
        sys.exit(1)

    print(f"[eval_syn] Processing {org} with margin={margin}, min_score={min_score}")

    # Load aligned map
    with open(aligned_path, 'rb') as f:
        aligned_map = pickle.load(f)

    # Keep a reference copy for looking up neighboring anchors
    # (we modify aligned_map in place, so need unmodified reference)
    with open(aligned_path, 'rb') as f:
        reference_map = pickle.load(f)

    # Run syntenic evaluation
    syntenic_count, total_dups = eval_syntenic_dups(
        org, aligned_map, reference_map, margin, min_score
    )

    # Save updated aligned map back to the same file
    with open(aligned_path, 'wb') as f:
        pickle.dump(aligned_map, f)

    print(f"[eval_syn] {syntenic_count}/{total_dups} dups_matches marked as syntenic")


def run_finalize_succinct(org, work_dir, root_dir):
    """
    Mode 2: Reciprocity check + create succinct output.

    Requires ALL organisms to have completed eval_syn (barrier in Snakemake).

    1. Load aligned/{org} (already has syntenic markings from eval_syn)
    2. Run reciprocity check — set ambiguous flag IN MEMORY for one-sided syntenic dups
    3. Create aligned_succinct/{org} with compress format (uses in-memory flags)

    NOTE: Does NOT write back to aligned/{org}. The ambiguous flags are applied in
    memory only, consumed by create_succinct() in the same process. Writing back would
    cause a race condition: all finalize_succinct rules run in parallel, and one
    instance writing aligned/{org1} while another reads it as a partner would corrupt
    the pickle file. The aligned/{org} on disk retains the eval_syn state (syntenic
    markings) which is the last safe write point.
    """
    aligned_path = Path(root_dir) / 'anchors' / 'aligned' / org
    anchor_dir = Path(root_dir) / 'anchors'
    succinct_dir = Path(root_dir) / 'anchors' / 'aligned_succinct'
    succinct_path = succinct_dir / org

    if not aligned_path.exists():
        print(f"ERROR: Aligned file not found: {aligned_path}", file=sys.stderr)
        sys.exit(1)

    print(f"[finalize_succinct] Processing {org}")

    # Load aligned map (already has syntenic markings)
    with open(aligned_path, 'rb') as f:
        aligned_map = pickle.load(f)

    # Run reciprocity check — sets ambiguous flag IN MEMORY for one-sided syntenic dups
    # (not written back to disk to avoid race condition with parallel partner reads)
    flagged = check_eval_syn_reciprocity(org, aligned_map, anchor_dir, work_dir)
    print(f"[finalize_succinct] Reciprocity check: {flagged} (anchor, org2) relationships flagged as ambiguous")

    # Create succinct output (uses in-memory flags from reciprocity check)
    succinct_dir.mkdir(parents=True, exist_ok=True)
    succinct = create_succinct(org, aligned_map)

    with open(succinct_path, 'wb') as f:
        pickle.dump(succinct, f)

    print(f"[finalize_succinct] Succinct: {len(succinct)} anchors with valid matches written to {succinct_path}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <organism> <work_dir> <root_dir> <mode>", file=sys.stderr)
        print(f"  mode = eval_syn | finalize_succinct", file=sys.stderr)
        sys.exit(1)

    org = sys.argv[1]
    work_dir = sys.argv[2]
    root_dir = sys.argv[3]
    mode = sys.argv[4]

    if mode == 'eval_syn':
        run_eval_syn(org, work_dir, root_dir)
    elif mode == 'finalize_succinct':
        run_finalize_succinct(org, work_dir, root_dir)
    else:
        print(f"ERROR: Unknown mode '{mode}'. Use 'eval_syn' or 'finalize_succinct'.", file=sys.stderr)
        sys.exit(1)
