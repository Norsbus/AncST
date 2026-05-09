"""
Stage 1: Filter Candidate Windows
Using GenMap k-mer counts and macle complexity scores
"""

rule filter_windows:
    """
    Filter genomic windows using GenMap frequency data and macle complexity scores.
    Uses organism-specific K, E, L, I, P parameters from genmap_params.txt and/or macle_params.txt.
    At least one parameter file must be present.

    On subsequent runs: If candidates/{org} exists, marks deleted candidates in to_del/{org}/to_del
    """
    input:
        unpack(get_filter_inputs),
        # Parameter files (global - not per organism)
        # NOTE: ancient() because timestamp doesn't matter - only content within this run
        genmap_params = ancient(f"{WORK_DIR}/genmap_params.txt") if os.path.exists(f"{WORK_DIR}/genmap_params.txt") else [],
        macle_params = ancient(f"{WORK_DIR}/macle_params.txt") if os.path.exists(f"{WORK_DIR}/macle_params.txt") else []
    output:
        indices = f"{WORK_DIR}/windows_low_kmers_indices/{{organism}}",
        fasta_fwd = f"{WORK_DIR}/windows_low_kmers_fastas/{{organism}}.fasta",
        fasta_rev = f"{WORK_DIR}/windows_low_kmers_fastas_rev/{{organism}}.fasta",
        # to_del created if candidates/{org} exists (subsequent runs)
        to_del = f"{WORK_DIR}/to_del/{{organism}}/to_del"
    params:
        # Pass first file for backwards compatibility (script ignores and reads from params files)
        genmap_file = lambda w, input: input.genmap[0] if isinstance(input.genmap, list) else input.genmap,
        macle_file = lambda w, input: input.macle[0] if isinstance(input.macle, list) else input.macle,
        k = get_genmap_k,
        e = get_genmap_e,
        L = get_genmap_L,
        I = get_genmap_I,
        P = get_genmap_percentile
    log:
        err = f"{WORK_DIR}/logs/filter/{{organism}}.err",
        out = f"{WORK_DIR}/logs/filter/{{organism}}.out"
    threads: 1
    shell:
        f"""cd {WORK_DIR} && \
            {PYTHON} {ROOT_DIR}/code/filter_windows_macle_and_genmap_manual.py \
            {{wildcards.organism}} \
            {{params.genmap_file}} \
            {{params.macle_file}} \
            {{output.indices}} \
            {{threads}} \
            {{params.k}} {{params.e}} {{params.L}} {{params.I}} {{params.P}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}"""

checkpoint split_fastas:
    """
    Split candidate FASTA files into chunks for parallel BLAST processing.
    Part size is configurable via pipeline_config.yaml (splitting.self_blast_part_size).
    This is a CHECKPOINT to enable dynamic file discovery for BLAST stage.
    """
    input:
        fwd = ancient(rules.filter_windows.output.fasta_fwd),
        rev = ancient(rules.filter_windows.output.fasta_rev)
    output:
        touch(f"{WORK_DIR}/touch/split_done_{{organism}}")
    params:
        split_dir_fwd = f"{WORK_DIR}/split_out_forward/{{organism}}",
        split_dir_rev = f"{WORK_DIR}/split_out_reverse/{{organism}}",
        part_size = SELF_BLAST_PART_SIZE
    group: "initial_split"
    shell:
        """
        fasta-splitter --part-size {params.part_size} --out-dir {params.split_dir_fwd} {input.fwd}
        fasta-splitter --part-size {params.part_size} --out-dir {params.split_dir_rev} {input.rev}
        """
