"""
Phase 3a: Update Aligned Maps
Cross-genome synchronization to clean up references to deleted candidates.

Two-stage process with barrier:
1. Stage 1: NOW INTEGRATED INTO update_candidates
2. BARRIER: Wait for all orgs to write notifications
3. Stage 2: Read notifications and clean aligned maps (per genome)

Runs in parallel with Phase 3b (bcamm).
"""

import os

# Rule order is defined in main Snakefile to avoid conflicts
# (removed duplicate ruleorder that conflicted with main Snakefile)

# Per-pair checkpoint scheduling: each pair (org1, org2) gets its own
# aggregate_pair rule whose input function calls checkpoint.get() for only its
# 2 organisms. When both complete, that pair's bcamm jobs are immediately
# discovered and scheduled. collect_output depends on static pair_done touch
# files (no checkpoints). This replaces an earlier global-barrier approach
# where collect_output expand()ed all update_compared checkpoints at once and
# a single IncompleteCheckpointException aborted the whole function.

def get_pairs_to_compute_static():
    """
    Get list of (org1, org2) pairs to compute STATICALLY (no filesystem access).
    Uses the same redundancy elimination logic as get_pairwise_tuples().
    """
    if PAIRWISE_FILE is not None and os.path.exists(PAIRWISE_FILE):
        return load_and_validate_custom_pairs(PAIRWISE_FILE, ALL_ORGS)
    else:
        pairs = []
        for org1 in COMPUTE_ANCHORS_FOR:
            for org2 in ALL_ORGS:
                if org2 == org1:
                    continue
                # OPTIMIZATION: Skip symmetric pair if both in compute_for
                if org2 in COMPUTE_ANCHORS_FOR and org2 < org1:
                    continue
                pairs.append((org1, org2))
        return pairs

# Compute pairs at Snakefile parse time (all variables available from common.smk)
_PAIRS_TO_COMPUTE = get_pairs_to_compute_static()


def get_bcamm_outputs_for_pair(wildcards):
    """
    Enumerate parse_bcamm outputs for a SINGLE (org1, org2) pair.

    Only triggers checkpoint.get() for the two organisms in this pair.
    This means bcamm jobs are discovered per-pair: as soon as BOTH organisms'
    update_compared checkpoints complete, their bcamm jobs are scheduled.

    KEY ADVANTAGE over the old global function:
    - Old: ONE function loops over ALL pairs, hits first incomplete checkpoint,
      aborts. Only after ALL checkpoints complete can it enumerate chunks.
    - New: Each pair has its OWN function with at most 2 checkpoint.get() calls.
      Pairs whose organisms are ready get scheduled immediately.

    NOTE: If org1's checkpoint is incomplete, the exception aborts before
    reaching org2's get(). On re-evaluation after org1 completes, org2 is
    checked. This is at most a 2-cycle delay per pair (vs N-cycle globally).
    """
    org1 = wildcards.org1
    org2 = wildcards.org2

    # Trigger checkpoints for both organisms
    # checkpoint.get() raises IncompleteCheckpointException if not done,
    # causing Snakemake to schedule the checkpoint and defer this function
    checkpoints.update_compared.get(organism=org1)
    checkpoints.update_compared.get(organism=org2)

    # Both checkpoints complete - enumerate chunks from org2's forward_split
    chunk_dir = f"{WORK_DIR}/sequences_to_compare/{org2}/forward_split"
    chunks = sorted([
        name.split('forward.')[1]
        for name in os.listdir(chunk_dir)
        if name.startswith('forward.')
    ])

    parse_outputs = []
    orgs_tuple = f"{org1}---{org2}"
    for chunk in chunks:
        # parse_bcamm creates TWO files per comparison
        parse_outputs.append(f"{WORK_DIR}/parse_bcamm/{org1}/{orgs_tuple}_{chunk}")
        parse_outputs.append(f"{WORK_DIR}/parse_bcamm/{org2}/{orgs_tuple}_{chunk}")

    return parse_outputs


rule aggregate_pair:
    """
    Per-pair aggregation: wait for all parse_bcamm outputs for one (org1, org2) pair.

    This rule breaks the global checkpoint barrier. Instead of one function
    needing ALL organisms' checkpoints to complete, each pair only needs its
    TWO organisms' checkpoints. bcamm jobs start per-pair as soon as both
    organisms complete update_compared, regardless of other organisms' progress.

    The touch file serves as a static dependency for collect_output.
    """
    input:
        get_bcamm_outputs_for_pair
    output:
        done = touch(f"{WORK_DIR}/touch/pair_done_{{org1}}---{{org2}}")
    resources:
        mem_mb = 100  # Trivial synchronization rule


def get_pair_done_files_for_organism(wildcards):
    """
    Get all pair_done touch files for pairs involving this organism.
    STATIC enumeration - no checkpoint calls, no filesystem scanning.
    Uses _PAIRS_TO_COMPUTE computed at Snakefile parse time.
    """
    org = wildcards.organism
    return [
        f"{WORK_DIR}/touch/pair_done_{org1}---{org2}"
        for org1, org2 in _PAIRS_TO_COMPUTE
        if org in (org1, org2)
    ]

# Phase 3a: update aligned maps.
# Stage 1 (creating aligned/{org}) is handled by update_candidates rule in
# 03_update_candidates.smk. Only stage 2 remains here.

rule update_aligned_stage2:
    """
    Stage 2: Read notifications and clean aligned map.

    BARRIER: Waits for ALL organisms to finish update_candidates (which now writes notifications).
    Reads all from_* notifications and removes stale match references.
    Updates aligned map IN-PLACE.
    """
    input:
        # CRITICAL: Wait for ALL orgs to complete update_candidates (writes from_* notifications)
        [ancient(f"{WORK_DIR}/touch/update_candidates_done_{organism}")
         for organism in ALL_ORGS]
    output:
        # Aligned file updated in-place (explicit about modification)
        aligned = update(f"{ROOT_DIR}/anchors/aligned/{{organism}}"),
        done = touch(f"{WORK_DIR}/touch/update_aligned_stage2_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/update_aligned_stage2/{{organism}}.err",
        out = f"{WORK_DIR}/logs/update_aligned_stage2/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['update_aligned_stage2'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./update_aligned.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            stage2 \
            > {{log.out}} 2> {{log.err}}"""

# Phase 4: match collection and validation.

rule collect_output:
    """
    Collect all pairwise comparison results for each organism.

    MUST wait for:
    1. update_aligned_stage2 (clean aligned map)
    2. All per-pair aggregate_pair rules (which ensure bcamm/parse_bcamm are done)

    NOTE: The global expand() barrier is REMOVED. Per-pair aggregate_pair rules
    handle checkpoint triggering individually via get_bcamm_outputs_for_pair(),
    enabling bcamm jobs to start as soon as each pair's two organisms complete
    update_compared (not waiting for ALL organisms).
    """
    input:
        # Wait for stage2 to complete (barrier)
        [ancient(f"{WORK_DIR}/touch/update_aligned_stage2_done_{organism}")
         for organism in ALL_ORGS],
        # Wait for all pairs involving this organism (static file paths, no checkpoints)
        pair_done = get_pair_done_files_for_organism
    output:
        # Aligned file updated in-place with collected matches (explicit about modification)
        aligned = update(f"{ROOT_DIR}/anchors/aligned/{{organism}}"),
        done = touch(f"{WORK_DIR}/touch/collect_output_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/collect_output/{{organism}}.err",
        out = f"{WORK_DIR}/logs/collect_output/{{organism}}.out"
    resources:
        # partner-load rule: accumulates candidates/{org2} for every partner in a
        # `bibs` dict; partner candidates are small (~genome/50) so max_genome with
        # MULT=5 comfortably covers focal_aligned + sum(partner_candidates)
        mem_mb = lambda wildcards, attempt: get_max_genome_size_mb() * MEM_MULT['collect_output'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./collect_output.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}"""

rule resolve_inconsistencies:
    """
    Process inconsistency flags across all organisms.
    SYNCHRONIZATION BARRIER: Requires ALL organisms to finish collect_output.
    Updates anchors/aligned/{org} IN-PLACE to resolve cross-genome inconsistencies.

    NOTE: This is an in-place update chain using Snakemake's update() wrapper:
      1. update_candidates CREATES aligned/{org}
      2. update_aligned_stage2 UPDATES aligned/{org} in-place
      3. collect_output UPDATES aligned/{org} in-place
      4. resolve_inconsistencies UPDATES aligned/{org} in-place
      5. eval_syn UPDATES aligned/{org} in-place (syntenic markers, LAST WRITE)
      6. finalize_succinct READS aligned/{org} + partners, creates aligned_succinct/{org}
         (does NOT write aligned/{org} — flags applied in memory only to avoid race condition)
    """
    input:
        # CRITICAL: Global dependency - all orgs must finish collecting
        [ancient(f"{WORK_DIR}/touch/collect_output_done_{organism}")
         for organism in ALL_ORGS]
    output:
        # Aligned file finalized with resolved inconsistencies (explicit about modification)
        aligned = update(f"{ROOT_DIR}/anchors/aligned/{{organism}}"),
        done = touch(f"{WORK_DIR}/touch/inconsistencies_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/inconsistencies/{{organism}}.err",
        out = f"{WORK_DIR}/logs/inconsistencies/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['resolve_inconsistencies'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./inconsistencies.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}"""

# Phase 4.5: eval_syn + finalize_succinct.

rule eval_syn:
    """
    Evaluate which dups_matches are syntenic (supported by flanking anchors).
    Updates aligned/{org} IN-PLACE with 'syntenic' sets in dups_matches.
    Per-organism, no barrier needed (only reads own aligned file).
    """
    input:
        inconsistencies_done = ancient(f"{WORK_DIR}/touch/inconsistencies_done_{{organism}}")
    output:
        aligned = update(f"{ROOT_DIR}/anchors/aligned/{{organism}}"),
        done = touch(f"{WORK_DIR}/touch/eval_syn_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/eval_syn/{{organism}}.err",
        out = f"{WORK_DIR}/logs/eval_syn/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['eval_syn'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./finalize_aligned.py \
            {{wildcards.organism}} {WORK_DIR} {ROOT_DIR} eval_syn \
            > {{log.out}} 2> {{log.err}}"""

rule finalize_succinct:
    """
    Check eval_syn reciprocity and create succinct output.

    BARRIER: Requires ALL organisms' eval_syn to complete, because cross-org
    reciprocity check loads partner aligned maps that eval_syn modifies.

    Sets 'matches have ambiguous matches' flag IN MEMORY for one-sided syntenic
    dups (not written back to aligned/{org} to avoid race condition with parallel
    partner reads), then creates aligned_succinct/{org} with simplified tuple format.

    Memory: Loads aligned/{org} + partner aligned maps one at a time (~2x map).
    """
    input:
        eval_syn_done = ancient(f"{WORK_DIR}/touch/eval_syn_done_{{organism}}"),
        all_eval_syn = [ancient(f"{WORK_DIR}/touch/eval_syn_done_{org}") for org in ALL_ORGS]
    output:
        succinct = f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}"
    log:
        err = f"{WORK_DIR}/logs/finalize_succinct/{{organism}}.err",
        out = f"{WORK_DIR}/logs/finalize_succinct/{{organism}}.out"
    resources:
        # partner-load rule: check_eval_syn_reciprocity loads aligned/{org2} for
        # every partner one at a time; mem_mb sized from the largest possible partner
        mem_mb = lambda wildcards, attempt: get_max_genome_size_mb() * MEM_MULT['finalize_succinct'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./finalize_aligned.py \
            {{wildcards.organism}} {WORK_DIR} {ROOT_DIR} finalize_succinct \
            > {{log.out}} 2> {{log.err}}"""
