"""
Phase 2: Update Candidates
Combines update_candidates + update_aligned_stage1 + add_dups logic into single per-genome step.

Flow:
1. Delete old candidates from to_del (filter_windows)
2. Add new candidates from indices (blast-clasp)
3. Integrate dups (if exists)
4. Update aligned map
5. Track all deletions and write cross-genome notifications

This runs AFTER Phase 1 (filter + self-blast + get_dups) completes.
"""

def has_dups_params():
    """Check if dups_params.txt exists and is non-empty."""
    try:
        return os.path.exists(f"{WORK_DIR}/dups_params.txt") and os.path.getsize(f"{WORK_DIR}/dups_params.txt") > 0
    except:
        return False

def get_to_del_input(wildcards):
    """
    Get to_del file for update_candidates.
    OPTIMIZATION: If candidates already exist, return empty list to skip filter_windows.
    This allows pairwise_downstream target to skip anchor generation stages.
    """
    candidates_file = f"{ROOT_DIR}/anchors/candidates/{wildcards.organism}"
    touch_file = f"{WORK_DIR}/touch/update_candidates_done_{wildcards.organism}"
    if os.path.exists(candidates_file) and os.path.exists(touch_file):
        return []
    return f"{WORK_DIR}/to_del/{wildcards.organism}/to_del"

def get_dups_input(wildcards):
    """
    Get dups file for update_candidates.
    OPTIMIZATION: If candidates already exist, return empty list to skip get_dups.
    Only returns dups file if:
    1. Outputs don't exist yet (need to run update_candidates)
    2. dups_params.txt exists and is non-empty (dups processing enabled)
    """
    candidates_file = f"{ROOT_DIR}/anchors/candidates/{wildcards.organism}"
    touch_file = f"{WORK_DIR}/touch/update_candidates_done_{wildcards.organism}"
    if os.path.exists(candidates_file) and os.path.exists(touch_file):
        return []
    if not has_dups_params():
        return []
    return ancient(f"{WORK_DIR}/dups/{wildcards.organism}")

def get_metadata_input(wildcards):
    """
    Get metadata file for update_candidates.
    OPTIMIZATION: If candidates already exist, return empty list to skip dependency.

    This is CRITICAL for pairwise_downstream target: when anchors exist but
    metadata_genomes is regenerated (because update_compared needs it), we must
    NOT trigger update_candidates. Without this, ancient() alone doesn't help
    because Snakemake still re-runs rules when inputs are produced in the same execution.
    """
    candidates_file = f"{ROOT_DIR}/anchors/candidates/{wildcards.organism}"
    touch_file = f"{WORK_DIR}/touch/update_candidates_done_{wildcards.organism}"
    if os.path.exists(candidates_file) and os.path.exists(touch_file):
        return []  # Skip metadata dependency - we won't be running anyway
    return [ancient(f"{ROOT_DIR}/utils/metadata_genomes/{wildcards.organism}")]

rule update_candidates:
    """
    Unified update of candidates AND aligned map with cross-genome notifications.
    Combines former update_candidates + update_aligned_stage1 logic:
    - Merges all .fasta_indices from BLAST chunks
    - Deletes old candidates while tracking matches
    - Adds new candidates from indices
    - Integrates dups
    - Updates aligned map
    - Writes cross-genome notifications (from_* files)

    Can run per-organism as soon as its BLAST completes.
    """
    input:
        # Wait for ALL self-BLAST chunks to complete (skipped if outputs exist)
        indices_files = lambda wildcards: get_all_blast_outputs(wildcards),
        # Metadata file (skipped if outputs exist - prevents pairwise_downstream from triggering this rule)
        metadata = lambda wildcards: get_metadata_input(wildcards),
        # to_del created by filter_windows WITHIN THIS RUN (skipped if outputs exist)
        to_del_in = lambda wildcards: get_to_del_input(wildcards),
        # Optional: wait for dups if dups_params.txt exists and is non-empty (skipped if outputs exist)
        dups = lambda wildcards: get_dups_input(wildcards)
    output:
        # Primary outputs (use update() for incremental updates)
        candidates = update(f"{ROOT_DIR}/anchors/candidates/{{organism}}"),
        # Aligned file created or updated (use update() for re-runs)
        aligned = update(f"{ROOT_DIR}/anchors/aligned/{{organism}}"),
        # NOTE: indices_global is created internally by the script but NOT declared as output
        # because it's only used within the same script and not by any other rule.
        # Declaring it as output would cause Snakemake to re-run the rule if it's missing
        # even when candidates/aligned already exist from a previous run.
        # Touch file for sequencing
        done = touch(f"{WORK_DIR}/touch/update_candidates_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/update_candidates/{{organism}}.err",
        out = f"{WORK_DIR}/logs/update_candidates/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['update_candidates'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./update_candidates_and_aligned.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}"""
