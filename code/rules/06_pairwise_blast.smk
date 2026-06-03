"""
Stage 6: Pairwise BLAST + CLASP + AWK + Match (BCAMM)
Find synteny anchors between genome pairs
This is the most computationally intensive stage
"""

rule bcamm:
    """
    BLAST-CLASP-AWK-MergeMatches for pairwise genome comparison.
    For each (org1, org2, chunk) combination:
    - BLAST org1 candidates against org2 chunk
    - Run CLASP to chain hits
    - Filter by score threshold (>40)

    The {file} wildcard represents split chunk names (e.g., part-1.fasta).
    bcamm.py reads from forward_split/forward.{file} and reverse_split/reverse.{file}.
    """
    input:
        # Depend on org1's BLAST database (to query from)
        # BLAST requires ALL 3 files (.nin, .nhr, .nsq) - declare all explicitly
        # NOTE: No checkpoint.get() needed because bcamm jobs are only discovered via
        # get_bcamm_outputs_for_pair() which triggers checkpoint.get() for both organisms.
        # By the time bcamm is scheduled, org1's checkpoint is guaranteed complete.
        org1_blastdb_nsq = ancient(f"{WORK_DIR}/blastdbs/anchor_candidates_{{org1}}_forward.nsq"),
        org1_blastdb_nin = ancient(f"{WORK_DIR}/blastdbs/anchor_candidates_{{org1}}_forward.nin"),
        org1_blastdb_nhr = ancient(f"{WORK_DIR}/blastdbs/anchor_candidates_{{org1}}_forward.nhr"),
        # NOTE: No explicit org2 checkpoint dependency needed either because:
        # 1. bcamm jobs are only discovered via get_bcamm_outputs_for_pair() in aggregate_pair
        # 2. That function triggers checkpoint.get() for both organisms before enumerating chunks
        # 3. The file dependencies below (org2_fwd_split, org2_rev_split) implicitly wait for org2's files
        # Depend on org2's split files (to query against)
        # The {file} wildcard is filled by downstream checkpoint discovery
        org2_fwd_split = ancient(f"{WORK_DIR}/sequences_to_compare/{{org2}}/forward_split/forward.{{file}}"),
        org2_rev_split = ancient(f"{WORK_DIR}/sequences_to_compare/{{org2}}/reverse_split/reverse.{{file}}")
    output:
        fwd_awk = f"{WORK_DIR}/clasp_out_forward/{{org1}}---{{org2}}/{{file}}_awk",
        rev_awk = f"{WORK_DIR}/clasp_out_reverse/{{org1}}---{{org2}}/{{file}}_awk"
    log:
        err = f"{WORK_DIR}/logs/bcamm/{{org1}}---{{org2}}/{{file}}.err",
        out = f"{WORK_DIR}/logs/bcamm/{{org1}}---{{org2}}/{{file}}.out"
    params:
        ws = WS_ANCHOR_MAKING
    shell:
        # Reconstruct orgs_tuple for Python script (which doesn't change)
        f"""cd {CODE_DIR} && {PYTHON} ./bcamm.py \
            {{wildcards.org1}}---{{wildcards.org2}} \
            {{wildcards.file}} \
            {{params.ws}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}} \
        && awk '{{{{ if ($1 != "#" && ($8+0)>40) print }}}}' \
            {WORK_DIR}/clasp_out_forward/{{wildcards.org1}}---{{wildcards.org2}}/{{wildcards.file}} \
            > {{output.fwd_awk}} \
        && awk '{{{{ if ($1 != "#" && ($8+0)>40) print }}}}' \
            {WORK_DIR}/clasp_out_reverse/{{wildcards.org1}}---{{wildcards.org2}}/{{wildcards.file}} \
            > {{output.rev_awk}}"""

rule parse_bcamm:
    """
    Parse CLASP output and identify synteny anchors.
    Applies score thresholds based on self-BLAST results.
    Creates parse_bcamm outputs for both organisms in the pair.

    Uses {org1} and {org2} wildcards instead of {orgs_tuple} to avoid lambda in output.
    Both output files use same {org1}---{org2} pattern in filename, different directories.
    """
    input:
        fwd_awk = ancient(rules.bcamm.output.fwd_awk),
        rev_awk = ancient(rules.bcamm.output.rev_awk),
        # parse_bcamm.py reads candidate dictionaries for both organisms
        candidates_org1 = ancient(f"{ROOT_DIR}/anchors/candidates/{{org1}}"),
        candidates_org2 = ancient(f"{ROOT_DIR}/anchors/candidates/{{org2}}")
    output:
        # NO LAMBDA! Static wildcards only
        # parse_bcamm.py creates TWO files: one in org1's directory, one in org2's directory
        # Both files use {org1}---{org2}_{file} naming pattern
        org1_file = f"{WORK_DIR}/parse_bcamm/{{org1}}/{{org1}}---{{org2}}_{{file}}",
        org2_file = f"{WORK_DIR}/parse_bcamm/{{org2}}/{{org1}}---{{org2}}_{{file}}"
    params:
        ws = WS_ANCHOR_MAKING
    shell:
        f"""
        # Defensive directory creation (normally handled by make_directories.py)
        # Protects against manual deletion or --continue-run with incomplete setup
        mkdir -p {WORK_DIR}/parse_bcamm/{{wildcards.org1}}
        mkdir -p {WORK_DIR}/parse_bcamm/{{wildcards.org2}}

        cd {CODE_DIR} && {PYTHON} ./parse_bcamm.py \
            {{wildcards.org1}}---{{wildcards.org2}} \
            {{wildcards.file}} \
            {{params.ws}} \
            {WORK_DIR}"""
