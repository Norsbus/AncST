"""
Stage 4: Find Duplications (Optional)
Identify duplicate/paralogous regions using different parameters
This stage is INDEPENDENT and could run in parallel with earlier stages!
"""

import os

def get_dups_genmap_inputs(wildcards):
    """
    Dynamically discover genmap files needed for dups based on dups_params.txt.
    Returns list of genmap frequency files for the organism.
    Uses set to avoid duplicate file declarations.
    """
    org = wildcards.organism
    genmap_files = set()

    if os.path.exists(f"{WORK_DIR}/dups_params.txt"):
        with open(f"{WORK_DIR}/dups_params.txt", "r") as f:
            for line in f:
                if org in line:
                    parts = line.strip().split()[1:]
                    if len(parts) >= 4:
                        k1 = parts[0]
                        e1 = parts[1]
                        k2 = parts[2]
                        e2 = parts[3]

                        # Add both k1/e1 and k2/e2 genmap files
                        genmap_files.add(ancient(f"{ROOT_DIR}/utils/genmap_out/{org}/{k1}_{e1}.freq16"))
                        genmap_files.add(ancient(f"{ROOT_DIR}/utils/genmap_out/{org}/{k2}_{e2}.freq16"))

    return list(genmap_files)

rule get_dups:
    """
    Find potential duplications using GenMap with different parameters.
    Only runs if dups_params.txt exists and is non-empty.

    Dynamic genmap inputs are discovered at runtime from dups_params.txt using
    get_dups_genmap_inputs() function (similar to checkpoint + lambda pattern).
    """
    input:
        genome = ancient(f"{ROOT_DIR}/utils/genomes/{{organism}}.fasta"),
        metadata = ancient(f"{ROOT_DIR}/utils/metadata_genomes/{{organism}}"),
        # Parameter file (ancient - timestamp doesn't matter within a run)
        dups_params = ancient(f"{WORK_DIR}/dups_params.txt") if os.path.exists(f"{WORK_DIR}/dups_params.txt") else [],
        # Dynamic genmap files discovered from dups_params.txt
        genmap_files = lambda w: get_dups_genmap_inputs(w)
    output:
        f"{WORK_DIR}/dups/{{organism}}"
    log:
        err = f"{WORK_DIR}/logs/dups/{{organism}}.err",
        out = f"{WORK_DIR}/logs/dups/{{organism}}.out"
    resources:
        mem_mb = lambda wildcards, attempt: get_genome_size_mb(wildcards) * MEM_MULT['get_dups'] * attempt
    shell:
        f"""cd {CODE_DIR} && {PYTHON} ./find_dups.py \
            {{wildcards.organism}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}"""
