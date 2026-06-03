"""
Phase 0: Utils Preprocessing
Creates indices and metadata required for candidate filtering.
All rules are per-genome and fully parallel.

Memory allocation is dynamic based on genome file size.
Multipliers are defined in common.smk (MEM_MULT dict).
Use --low-mem flag for reduced multipliers, --high-mem for 2x multipliers.
"""

import os

rule macle_index:
    """
    Create macle index for organism.
    Required for macle complexity scoring.

    Memory: MEM_MULT['macle_index'] * genome size (suffix array construction)
    """
    input:
        genome = ancient(f"{ROOT_DIR}/utils/genomes/{{organism}}.fasta")
    output:
        index = f"{ROOT_DIR}/utils/macle_indices/{{organism}}"
    log:
        err = f"{WORK_DIR}/logs/macle_index/{{organism}}.err",
        out = f"{WORK_DIR}/logs/macle_index/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['macle_index'] * attempt
    shell:
        """
        macle -s {input.genome} > {output.index} \
            2> {log.err}
        echo "Created macle index for {wildcards.organism}" > {log.out}
        """

rule macle_Cm:
    """
    Compute macle complexity scores for organism with specific parameters.
    Creates one output per parameter set defined in macle_params.txt.

    Memory: MEM_MULT['macle_Cm'] * genome size (complexity calculation)
    """
    input:
        index = ancient(rules.macle_index.output.index)
    output:
        scores = f"{ROOT_DIR}/utils/macle_out/{{organism}}/{{w}}_{{p}}.txt"
    params:
        w = lambda w: w.w,
        p = lambda w: w.p
    log:
        err = f"{WORK_DIR}/logs/macle_Cm/{{organism}}_{{w}}_{{p}}.err",
        out = f"{WORK_DIR}/logs/macle_Cm/{{organism}}_{{w}}_{{p}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['macle_Cm'] * attempt
    shell:
        """
        mkdir -p $(dirname {output.scores})
        macle -i {input.index} -w {params.w} -k {params.p} > {output.scores} \
            2> {log.err}
        echo "Created macle Cm for {wildcards.organism} w={params.w} p={params.p}" > {log.out}
        """

rule genmap_index:
    """
    Create GenMap index for organism.
    Required for k-mer frequency counting.

    Memory: MEM_MULT['genmap_index'] * genome size (FM-index construction)
    """
    input:
        genome = ancient(f"{ROOT_DIR}/utils/genomes/{{organism}}.fasta")
    output:
        index = directory(f"{ROOT_DIR}/utils/genmap_indices/{{organism}}")
    log:
        err = f"{WORK_DIR}/logs/genmap_index/{{organism}}.err",
        out = f"{WORK_DIR}/logs/genmap_index/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['genmap_index'] * attempt
    shell:
        """
        genmap index -F {input.genome} -I {output.index} \
            > {log.out} 2> {log.err}
        """

rule genmap_map:
    """
    Compute k-mer frequencies using GenMap.
    Creates one .freq16 file per (k, e) parameter set.

    Memory: MEM_MULT['genmap_map'] * genome size (k-mer frequency mapping)
    Threads: adaptive (cores / num_species, allows parallel jobs)
    """
    input:
        index = ancient(rules.genmap_index.output.index)
    output:
        freq = f"{ROOT_DIR}/utils/genmap_out/{{organism}}/{{k}}_{{e}}.freq16"
    params:
        k = lambda w: w.k,
        e = lambda w: w.e,
        # GenMap with -fl -r adds .freq16 automatically, so pass prefix without extension
        prefix = lambda w, output: output.freq.replace('.freq16', '')
    log:
        err = f"{WORK_DIR}/logs/genmap_map/{{organism}}_{{k}}_{{e}}.err",
        out = f"{WORK_DIR}/logs/genmap_map/{{organism}}_{{k}}_{{e}}.out"
    threads: lambda wildcards: max(1, workflow.cores // sum(len(params) for params in GENMAP_PARAMS.values())) if hasattr(workflow, 'cores') and workflow.cores and GENMAP_PARAMS else 16
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['genmap_map'] * attempt
    shell:
        """
        mkdir -p $(dirname {output.freq})
        genmap map -K {params.k} -E {params.e} \
            -I {input.index} -O {params.prefix} -T {threads} -fl -r \
            > {log.out} 2> {log.err}
        """

rule make_blastdb:
    """
    Create BLAST database for organism.
    makeblastdb creates .nin, .nhr, and .nsq files atomically.
    Using .nsq as output marker (always present, never split for large DBs).
    """
    input:
        genome = ancient(f"{ROOT_DIR}/utils/genomes/{{organism}}.fasta")
    output:
        blastdb = f"{ROOT_DIR}/utils/blastdbs/{{organism}}.nsq"
    log:
        err = f"{WORK_DIR}/logs/blastdb/{{organism}}.err",
        out = f"{WORK_DIR}/logs/blastdb/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['make_blastdb'] * attempt
    shell:
        """
        makeblastdb -in {input.genome} \
            -out {ROOT_DIR}/utils/blastdbs/{wildcards.organism} \
            -dbtype nucl \
            >> {log.out} 2>> {log.err}
        """

rule make_metadata_genomes:
    """
    Create full metadata pickle with chromosome IDs, lengths, and full sequence.
    """
    input:
        genome = ancient(f"{ROOT_DIR}/utils/genomes/{{organism}}.fasta")
    output:
        metadata = f"{ROOT_DIR}/utils/metadata_genomes/{{organism}}"
    log:
        err = f"{WORK_DIR}/logs/metadata_genomes/{{organism}}.err",
        out = f"{WORK_DIR}/logs/metadata_genomes/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['make_metadata_genomes'] * attempt
    shell:
        """
        cd {CODE_DIR} && {PYTHON} ./create_metadata.py {wildcards.organism} {ROOT_DIR} --task metadata \
            >> {log.out} 2>> {log.err}
        """

rule make_small_meta:
    """
    Create small metadata pickle with chromosome IDs and cumulative lengths only (no sequence).
    """
    input:
        genome = ancient(f"{ROOT_DIR}/utils/genomes/{{organism}}.fasta")
    output:
        small_meta = f"{ROOT_DIR}/utils/small_meta/{{organism}}"
    log:
        err = f"{WORK_DIR}/logs/small_meta/{{organism}}.err",
        out = f"{WORK_DIR}/logs/small_meta/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['make_small_meta'] * attempt
    shell:
        """
        cd {CODE_DIR} && {PYTHON} ./create_metadata.py {wildcards.organism} {ROOT_DIR} --task small_meta \
            >> {log.out} 2>> {log.err}
        """

def get_all_preprocessing_outputs(wildcards):
    """
    Enumerate all preprocessing outputs for all organisms.
    Includes macle, genmap, and metadata for each organism.

    Note on organism usage:
    - BLASTdbs: ALL_ORGS (needed for pairwise BLAST against all organisms)
    - small_meta: ALL_ORGS (needed by downstream processing)
    - metadata_genomes: COMPUTE_ANCHORS_FOR only (only for anchor generation)
    - GenMap/macle: COMPUTE_ANCHORS_FOR only (filtered by GENMAP_PARAMS/MACLE_PARAMS)
    """
    outputs = []

    # BLASTdbs and small_meta needed for ALL organisms
    for org in ALL_ORGS:
        outputs.append(f"{ROOT_DIR}/utils/blastdbs/{org}.nsq")
        outputs.append(f"{ROOT_DIR}/utils/small_meta/{org}")

    # metadata_genomes only needed for anchor-generating organisms
    for org in COMPUTE_ANCHORS_FOR:
        outputs.append(f"{ROOT_DIR}/utils/metadata_genomes/{org}")

    # GenMap outputs (only for organisms in GENMAP_PARAMS, which is subset of COMPUTE_ANCHORS_FOR)
    for org in GENMAP_PARAMS:
        for params in GENMAP_PARAMS[org]:
            k, e = params['k'], params['e']
            outputs.append(f"{ROOT_DIR}/utils/genmap_out/{org}/{k}_{e}.freq16")

    # Macle outputs (only for organisms in MACLE_PARAMS, which is subset of COMPUTE_ANCHORS_FOR)
    for org in MACLE_PARAMS:
        for params in MACLE_PARAMS[org]:
            w, p = params['w'], params['p']
            outputs.append(f"{ROOT_DIR}/utils/macle_out/{org}/{w}_{p}.txt")

    for org in DUPS_PARAMS:
        for params in DUPS_PARAMS[org]:
            k1, e1, k2, e2 = params[0], params[1], params[2], params[3]
            outputs.append(f"{ROOT_DIR}/utils/genmap_out/{org}/{k1}_{e1}.freq16")
            outputs.append(f"{ROOT_DIR}/utils/genmap_out/{org}/{k2}_{e2}.freq16")

    return outputs

# NOTE: preprocessing_complete touch file has been removed.
# Rules now depend directly on actual preprocessing outputs via get_all_preprocessing_outputs()
