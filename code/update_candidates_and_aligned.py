#!/usr/bin/env python3
"""
Unified update_candidates and update_aligned_stage1 script.
Combines candidate update with aligned map update and cross-genome notifications.

Flow:
1. Load existing to_del from filter_windows
2. Load candidates and aligned (if exist)
3. Delete old candidates while tracking matches
4. Add new candidates from indices
5. Integrate dups while tracking matches
6. Write: candidates, aligned, from_* notifications
"""

import pickle
import os
import sys
import pathlib
from bisect import bisect_left, bisect_right
from copy import deepcopy
from subprocesses import get_min_anchor_size, get_score_threshold, get_overlap_threshold, get_gap_size


def which_chromo(seqlen, seqids, i):
    """
    Determine which chromosome a position belongs to.
    seqlen is cumulative lengths, seqids is list of chromosome IDs.
    """
    pos = bisect_left(seqlen, i)
    if pos < len(seqlen) and i == seqlen[pos]:
        pos += 1
    return seqids[pos]


def load_indices(org, work_dir):
    """
    Load all *_indices files and merge into single list.

    Index format from self-blast-clasp.py:
    [[(abs_start, abs_end), [[[region_start, region_end], score], ...]], ...]

    Returns: sorted list of indices
    """
    indices_dir = f'{work_dir}/indices/{org}'
    global_idx = []

    # Find all indices files (both new and old)
    for name in os.listdir(indices_dir):
        if '_indices' in name and 'global' not in name:
            try:
                with open(f'{indices_dir}/{name}', 'rb') as f:
                    chunk_indices = pickle.load(f)
                # Indices files contain lists, not dicts
                global_idx += chunk_indices
                print(f'Loaded {len(chunk_indices)} indices from {name}')
            except Exception as e:
                print(f'WARNING: Could not load {name}: {e}')
                raise

    # Also check old_indices if exists
    old_indices_dir = f'{work_dir}/old_indices/{org}'
    if os.path.exists(old_indices_dir):
        for name in os.listdir(old_indices_dir):
            if '_indices' in name and 'global' not in name:
                try:
                    with open(f'{old_indices_dir}/{name}', 'rb') as f:
                        chunk_indices = pickle.load(f)
                    global_idx += chunk_indices
                    print(f'Loaded {len(chunk_indices)} indices from old_indices/{name}')
                except Exception as e:
                    print(f'WARNING: Could not load old_indices/{name}: {e}')

    # Save merged indices
    with open(f'{indices_dir}/indices_global', 'wb') as f:
        pickle.dump(global_idx, f)

    print(f'Merged indices: {len(global_idx)} total')
    return sorted(global_idx)


def add_new_candidates(bib, bib_aligned, org, work_dir, root, ws):
    """
    Add new candidates from indices.

    Index format: [(abs_start, abs_end), [[[region_start, region_end], score], ...]]

    Builds candidate entries with:
    - regions: [1, ...] - region boundaries (relative, starts at 1)
    - scores_regions: [...] - scores between regions
    - chromosome: chromo name
    - start: chromosome-relative start
    - end: chromosome-relative end
    - blast_word_size: ws
    - matches: {}
    """
    # Load metadata
    with open(f'{root}/utils/metadata_genomes/{org}', 'rb') as f:
        seqids, seqlen, s = pickle.load(f)

    with open(f'{work_dir}/indices/{org}/indices_global', 'rb') as f:
        idx = pickle.load(f)

    idx = sorted(idx)
    score_threshold = get_score_threshold()
    added_count = 0

    for i in idx:
        abs_start = i[0][0]
        abs_end = i[0][1]

        if abs_start in bib:
            print(f'Index {abs_start} already in bib (should have been deleted), overwriting')

        # Determine chromosome
        chromo = which_chromo(seqlen, seqids, abs_start)
        chromo_idx = seqids.index(chromo)
        if chromo_idx > 0:
            chromo_offset = seqlen[chromo_idx - 1]
        else:
            chromo_offset = 0
        chromo_start = abs_start - chromo_offset
        chromo_end = abs_end - chromo_offset

        # Build regions and scores_regions
        # regions starts with 1 (relative position within window)
        # regions contains the boundaries, scores_regions contains scores between them
        regions = [1]
        scores_regions = []

        # Sort region/score pairs by start position
        region_scores = sorted(i[1], key=lambda x: x[0][0])

        for region_pair in region_scores:
            region_bounds, score = region_pair
            reg_start, reg_end = region_bounds

            if reg_start < regions[-1]:
                # Overlapping region, skip
                pass
            elif reg_start == regions[-1]:
                # Adjacent region
                pass
            else:
                # Gap before this region - fill with default score
                scores_regions.append(score_threshold)
                regions.append(reg_start)

            scores_regions.append(score)
            regions.append(reg_end)

        # Fill any gap at the end
        window_len = abs_end - abs_start
        if regions[-1] != window_len:
            scores_regions.append(score_threshold)
            regions.append(window_len)

        # Create candidate entry
        entry = {
            'regions': regions,
            'scores_regions': scores_regions,
            'chromosome': chromo,
            'start': chromo_start,
            'end': chromo_end,
            'blast_word_size': ws,
            'matches': {}
        }

        bib[abs_start] = entry
        bib_aligned[abs_start] = deepcopy(entry)
        added_count += 1

    print(f'Added {added_count} new candidates from indices')
    return bib, bib_aligned


def overlaps(x, y):
    """Check if two intervals overlap."""
    if x[0] >= y[0] and x[0] <= y[1]:
        return True
    if x[1] >= y[0] and x[1] <= y[1]:
        return True
    if y[0] >= x[0] and y[0] <= x[1]:
        return True
    if y[1] >= x[0] and y[1] <= x[1]:
        return True
    return False


def integrate_dups_and_track_matches(bib, aligned, org, work_dir, matches_to_del):
    """
    Integrate dups into candidates map.
    Track deletions for cross-genome cleanup.

    Logic (with 30% threshold):
    - Small overlap (< 30% of candidate): Shrink dup, KEEP candidate
    - Large overlap (>= 30%): Don't shrink, DELETE candidate
    - If shrunk dup < min_size (300bp), don't add it
    """
    # Check if dups file exists
    dups_file = f'{work_dir}/dups/{org}'
    if not os.path.exists(dups_file):
        print(f'No dups file found for {org}, skipping dups processing')
        return bib, aligned

    with open(dups_file, 'rb') as f:
        dups = pickle.load(f)

    print(f'Integrating {len(dups)} dups for {org}')

    # Get config values
    overlap_threshold = get_overlap_threshold()  # 0.3 = 30%
    gap_size = get_gap_size()  # 100 bp
    min_size = get_min_anchor_size()  # 300 bp

    # Get sorted candidates (include ALL - dups and non-dups) for overlap checking
    # CRITICAL: Must include existing dups to detect old-dup vs new-dup overlaps in incremental runs
    ses = sorted([(i, i + x['end'] - x['start']) for i, x in bib.items()])
    starts = [x[0] for x in ses]
    ends = [x[1] for x in ses]

    candidates_to_delete = []
    dups_added = 0

    for i, d in dups.items():
        end = i + d['end'] - d['start']

        # Find potentially overlapping candidates using binary search
        end_bigger_start_of_dup = bisect_left(ends, i)
        start_smaller_end_of_dup = bisect_left(starts, end)

        # Create new entry for the dup (may be modified)
        new_entry = {
            'i': i,
            'bib': deepcopy(d)
        }

        cur_start = i
        cur_end = end

        # Check overlaps with nearby candidates
        for idx in range(max(0, end_bigger_start_of_dup - 3), min(len(ends), start_smaller_end_of_dup + 4)):
            if not overlaps((starts[idx], ends[idx]), (cur_start, cur_end)):
                continue

            # Calculate overlap
            overlap_amount = min(cur_end, ends[idx]) - max(cur_start, starts[idx])
            cand_len = ends[idx] - starts[idx]

            if overlap_amount < cand_len * overlap_threshold:
                # Small overlap (< 30%): Shrink dup, KEEP candidate
                # Determine which side of dup overlaps the candidate

                if cur_end >= starts[idx] and cur_end <= ends[idx]:
                    # Dup's RIGHT side overlaps candidate's LEFT side
                    # Shrink dup's END (right side)
                    shrink_amount = cur_end - starts[idx] + gap_size
                    new_entry['bib']['end'] -= shrink_amount
                    cur_end -= shrink_amount
                elif cur_start >= starts[idx] and cur_start <= ends[idx]:
                    # Dup's LEFT side overlaps candidate's RIGHT side
                    # Shrink dup's START (left side)
                    shrink_amount = ends[idx] - cur_start + gap_size
                    new_entry['i'] += shrink_amount
                    new_entry['bib']['start'] += shrink_amount
                    cur_start += shrink_amount
                # Don't delete candidate for small overlaps
            else:
                # Large overlap (>= 30%): Don't shrink, DELETE candidate
                candidates_to_delete.append(starts[idx])

        # Validate and add adjusted dup
        dup_len = cur_end - cur_start
        dup_bib_len = new_entry['bib']['end'] - new_entry['bib']['start']

        # Check various invalid conditions (chromosome boundaries, too small, etc.)
        if (dup_len < min_size or cur_end < cur_start or cur_end < 0 or
            new_entry['bib']['end'] < 0 or new_entry['bib']['start'] < 0 or
            dup_bib_len < 0):
            continue

        # Mark as dup and add to both bib and aligned (CRITICAL: must mirror exactly)
        new_entry['bib']['dup'] = 1
        bib[new_entry['i']] = new_entry['bib']
        # CRITICAL FIX: Also add to aligned to maintain alignment with candidates
        # aligned must mirror bib at all times (same candidates, same basic structure)
        aligned[new_entry['i']] = deepcopy(new_entry['bib'])
        dups_added += 1

    # Delete overlapping candidates and track matches for cross-genome cleanup
    deleted_count = 0
    for cand_idx in set(candidates_to_delete):
        if cand_idx in bib:
            # Check if this candidate had matches before deleting
            if aligned and cand_idx in aligned:
                track_matches_for_deletion(cand_idx, aligned, org, matches_to_del)
                del aligned[cand_idx]

            del bib[cand_idx]
            deleted_count += 1

    print(f'Integrated {dups_added} dups, deleted {deleted_count} overlapping candidates')
    # CRITICAL FIX: Return both bib and aligned to maintain synchronization
    # aligned must mirror bib candidates exactly at all times
    return bib, aligned


def track_matches_for_deletion(cand_idx, aligned, org, matches_to_del):
    """
    Track matches for a candidate being deleted.
    Extracted helper function to avoid duplication.
    """
    if cand_idx not in aligned:
        return

    # 'matches' at anchor level always exists
    for org2, match_data in aligned[cand_idx]['matches'].items():
        # Collect all matches (regular and not considered)
        matched_indices = set()

        # 'matches not considered...' always exists (initialized in parse_bcamm)
        for j in match_data['matches not considered upon applying stricter score criterion since there are consistency issues']:
                matched_indices.add(j)

        if 'matches' in match_data:
            for j in match_data['matches']:
                matched_indices.add(j)

        # Handle dups_matches if present - track ALL of them
        if 'dups_matches' in match_data:
            for j, match_info in match_data['dups_matches'].items():
                if j == 'syntenic':
                    continue  # Skip the metadata key (it's a set, not a match)
                if isinstance(match_info, dict):
                    matched_indices.add(j)  # Track ALL dups_matches for deletion notification

        # Add to notifications
        for j in matched_indices:
            if org2 not in matches_to_del:
                matches_to_del[org2] = {}
            if j not in matches_to_del[org2]:
                matches_to_del[org2][j] = {org: set([cand_idx])}
            else:
                if org not in matches_to_del[org2][j]:
                    matches_to_del[org2][j][org] = set([cand_idx])
                else:
                    matches_to_del[org2][j][org].add(cand_idx)


def process_deletions_and_track_matches(to_del_list, bib, aligned, org, matches_to_del):
    """
    Delete candidates from to_del while tracking their matches for cross-genome cleanup.
    Updates both bib and aligned dictionaries.
    """
    deleted_from_bib = 0
    deleted_from_aligned = 0

    for cand_idx in to_del_list:
        # Delete from candidates
        if cand_idx in bib:
            del bib[cand_idx]
            deleted_from_bib += 1

        # Check aligned map for matches before deleting
        if aligned and cand_idx in aligned:
            # Track matches for notification
            track_matches_for_deletion(cand_idx, aligned, org, matches_to_del)

            # Delete from aligned
            del aligned[cand_idx]
            deleted_from_aligned += 1

    print(f'Deleted {deleted_from_bib} from candidates, {deleted_from_aligned} from aligned')
    return bib, aligned


def main():
    if len(sys.argv) < 3:
        print("Usage: update_candidates_and_aligned.py <organism> <work_dir>")
        sys.exit(1)

    org = sys.argv[1]
    work_dir = sys.argv[2]

    root = str(pathlib.Path(__file__).parents[1])
    anchor_dir = root + '/anchors'

    print(f'=' * 80)
    print(f'Unified update_candidates_and_aligned for {org}')
    print(f'Work dir: {work_dir}')
    print(f'Anchor dir: {anchor_dir}')
    print(f'=' * 80)

    # Get blast word size from params file
    try:
        with open(f'{work_dir}/params', 'r') as f:
            params = f.readline().strip()
            params_list = params.split()
            ws = params_list[5]
    except FileNotFoundError:
        print('WARNING: No params file found, using default word_size 11')
        ws = '11'

    # Initialize tracking for cross-genome notifications
    matches_to_del = {}

    # 1. Load to_del from filter_windows
    to_del_list = []
    try:
        with open(f'{work_dir}/to_del/{org}/to_del', 'rb') as f:
            to_del_list = pickle.load(f)
        print(f'Loaded {len(to_del_list)} candidates to delete from filter_windows')
    except FileNotFoundError:
        print(f'No to_del file from filter_windows (first run)')
    except Exception as e:
        print(f'Error loading to_del: {e}')
        raise

    # 2. Load existing candidates (if exist)
    bib = {}
    try:
        with open(f'{anchor_dir}/candidates/{org}', 'rb') as f:
            bib = pickle.load(f)
        print(f'Loaded {len(bib)} existing candidates')
    except FileNotFoundError:
        print(f'No existing candidates (first run)')

    # 3. Load existing aligned map (if exists)
    aligned = {}
    aligned_exists = False
    try:
        with open(f'{anchor_dir}/aligned/{org}', 'rb') as f:
            aligned = pickle.load(f)
        print(f'Loaded existing aligned map with {len(aligned)} candidates')
        aligned_exists = True
    except FileNotFoundError:
        # First run - will create aligned from candidates after update
        print(f'No aligned map found (will create after update)')

    # 4. Delete old candidates while tracking matches
    if to_del_list:
        bib, aligned = process_deletions_and_track_matches(
            to_del_list, bib, aligned, org, matches_to_del
        )

    # 5. Load indices
    indices = load_indices(org, work_dir)

    if len(indices) == 0:
        print(f'No indices found for {org}, exiting')
        return

    # 6. Add new candidates from indices
    bib, aligned = add_new_candidates(bib, aligned, org, work_dir, root, ws)

    # 7. Integrate dups (if exist) - also tracks matches for deletions
    # CRITICAL FIX: Capture both bib and aligned to maintain synchronization
    bib, aligned = integrate_dups_and_track_matches(bib, aligned, org, work_dir, matches_to_del)

    # 8. Verify no overlaps - FAIL if any found (data corruption)
    ses = sorted([(i, i + x['end'] - x['start']) for i, x in bib.items()])
    overlaps_found = []
    for enu, se in enumerate(ses[:-1]):
        if se[1] >= ses[enu + 1][0]:
            overlaps_found.append((se, ses[enu+1]))
            print(f'ERROR: Overlapping candidates after processing: {se} and {ses[enu+1]}', file=sys.stderr)

    if overlaps_found:
        print(f'FATAL: Found {len(overlaps_found)} overlapping candidate pairs. This indicates data corruption.', file=sys.stderr)
        print(f'Cannot proceed with corrupted data. Exiting.', file=sys.stderr)
        sys.exit(1)

    # 9. Write final candidates
    print(f'Writing {len(bib)} candidates')
    with open(f'{anchor_dir}/candidates/{org}', 'wb') as f:
        pickle.dump(bib, f)

    # 10. Write updated aligned map
    with open(f'{anchor_dir}/aligned/{org}', 'wb') as f:
        pickle.dump(aligned, f)
        f.flush()
        os.fsync(f.fileno())
    print(f'Wrote aligned map with {len(aligned)} candidates')

    # 11. Write cross-genome notifications
    for org2, notification_data in matches_to_del.items():
        try:
            os.makedirs(f'{work_dir}/to_del/{org2}', exist_ok=True)
            with open(f'{work_dir}/to_del/{org2}/from_{org}', 'wb') as f:
                pickle.dump(notification_data, f)
            print(f'Wrote notification to {org2} ({len(notification_data)} candidates affected)')
        except Exception as e:
            print(f'ERROR: Could not write notification to {org2}: {e}')
            raise

    print(f'Wrote notifications to {len(matches_to_del)} organisms')
    print(f'Unified update complete for {org}')


if __name__ == "__main__":
    main()
