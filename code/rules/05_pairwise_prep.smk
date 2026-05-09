"""
Stage 5: Prepare Pairwise Comparisons
Create split sequences for efficient pairwise BLAST
"""

checkpoint update_compared:
    """
    Prepare sequences for pairwise comparison.
    Depends on updated anchor candidates (Phase 2).
    CHECKPOINT: Creates split files that bcamm rules will process.
    """
    input:
        candidates = ancient(f"{ROOT_DIR}/anchors/candidates/{{organism}}"),
        metadata = ancient(f"{ROOT_DIR}/utils/metadata_genomes/{{organism}}")
    output:
        # Actual files created by update_compared.py
        fwd_fasta = f"{WORK_DIR}/sequences_to_compare/{{organism}}/forward.fasta",
        rev_fasta = f"{WORK_DIR}/sequences_to_compare/{{organism}}/reverse.fasta",
        blastdb = f"{WORK_DIR}/blastdbs/anchor_candidates_{{organism}}_forward.nsq"
        # Note: forward_split/* and reverse_split/* are created by fasta-splitter
        # These are discovered dynamically after checkpoint completion
    log:
        err = f"{WORK_DIR}/logs/update_compared/{{organism}}.err",
        out = f"{WORK_DIR}/logs/update_compared/{{organism}}.out"
    params:
        ws = WS_ANCHOR_MAKING
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['update_compared'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./update_compared.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}"""
