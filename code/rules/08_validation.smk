"""
Stage 8: Validation and Statistics
Final quality control and statistics generation
These can run in parallel per organism
"""

rule checks:
    """
    Validate anchor integrity:
    - Check for overlapping candidates
    - Verify chromosome boundaries
    - Validate match metadata
    - Cross-organism consistency (reciprocal matches, scores, coordinates)

    Depends on ALL organisms' finalize_aligned completing, because cross-org checks
    read aligned/{org2} which finalize_aligned modifies in-place. Without this barrier,
    checks could read a partially-modified aligned file.

    Memory: Loads aligned/{org} (kept in memory) + aligned/{org2} one at a time for each
    partner organism. Peak memory is ~2x a single aligned map. Uses MEM_MULT['checks']
    which is set higher than finalize_aligned since two maps are held simultaneously.
    """
    input:
        succinct = f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}",
        all_finalized = expand(f"{ROOT_DIR}/anchors/aligned_succinct/{{org}}", org=ALL_ORGS)
    output:
        touch(f"{WORK_DIR}/touch/checks_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/checks/{{organism}}.err",
        out = f"{WORK_DIR}/logs/checks/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['checks'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./checks.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            {ROOT_DIR}/anchors \
            > {{log.out}} 2> {{log.err}}"""

rule statistics:
    """
    Generate anchor statistics and coverage metrics.

    Depends on finalize_aligned to ensure aligned file is complete (not being modified).
    By depending on aligned_succinct output, we ensure finalize_aligned has fully completed.

    Memory: Loads aligned/{org} + small_meta for all orgs. Single aligned map.
    """
    input:
        succinct = f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}"
    output:
        touch(f"{WORK_DIR}/touch/statistics_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/statistics/{{organism}}.err",
        out = f"{WORK_DIR}/logs/statistics/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['statistics'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./aligned_statistics.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            {ROOT_DIR}/anchors \
            > {{log.out}} 2> {{log.err}}"""

# Succinct format validation (runs after finalize_aligned).

rule checks_succinct:
    """
    Validate aligned_succinct format, sanity, correspondence, and cross-org consistency.

    Checks:
    - Correspondence: succinct anchors are a subset of candidates
    - Sanity: coordinates, chromosome boundaries, overlaps
    - Format: tuple structure (score_int, (start, end), orientation_str), no 'syntenic' key
    - Cross-organism: score consistency for reciprocal matches

    Depends on ALL organisms' finalize_aligned completing, because cross-org checks
    load partner succinct files that only exist after finalize_aligned creates them.

    Memory: Loads succinct/{org} + candidates/{org} + small_meta/{org} + succinct/{org2}
    one at a time. Similar to checks but succinct maps are smaller than aligned maps.
    """
    input:
        succinct = f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}",
        all_finalized = expand(f"{ROOT_DIR}/anchors/aligned_succinct/{{org}}", org=ALL_ORGS)
    output:
        touch(f"{WORK_DIR}/touch/checks_succinct_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/checks_succinct/{{organism}}.err",
        out = f"{WORK_DIR}/logs/checks_succinct/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['checks_succinct'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./checks_succinct.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            {ROOT_DIR}/anchors \
            > {{log.out}} 2> {{log.err}}"""

rule statistics_succinct:
    """
    Generate statistics for aligned_succinct format.

    Reports match counts, score distributions, and coverage metrics.

    Depends on finalize_aligned to ensure succinct file exists.

    Memory: Loads succinct/{org} + small_meta. Lighter than aligned statistics.
    """
    input:
        succinct = f"{ROOT_DIR}/anchors/aligned_succinct/{{organism}}"
    output:
        touch(f"{WORK_DIR}/touch/statistics_succinct_done_{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/statistics_succinct/{{organism}}.err",
        out = f"{WORK_DIR}/logs/statistics_succinct/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['statistics_succinct'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./statistics_succinct.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            {ROOT_DIR}/anchors \
            > {{log.out}} 2> {{log.err}}"""
