"""
Stage 2: Self-BLAST and CLASP
BLAST each candidate against its own genome to establish edit distance thresholds
"""

import os

# Dynamic file discovery after checkpoint
def get_blast_chunk_files(wildcards):
    """
    Discover all FASTA chunks for an organism after split_fastas checkpoint completes.
    Returns list of chunk filenames for use in expand().

    NOTE: Callers should check if outputs exist before calling this to avoid
    unnecessary checkpoint triggers.
    """
    # Wait for checkpoint to complete
    checkpoint_output = checkpoints.split_fastas.get(**wildcards).output[0]

    # Discover forward chunks
    forward_dir = f"{WORK_DIR}/split_out_forward/{wildcards.organism}"
    chunks = []
    if os.path.exists(forward_dir):
        chunks = [f for f in os.listdir(forward_dir) if f.endswith('.fasta')]

    return chunks

def get_all_blast_outputs(wildcards):
    """
    Get all expected BLAST indices files for an organism.
    Used by update_candidates to wait for all BLAST jobs.

    OPTIMIZATION: If candidates already exist (e.g., from previous run),
    skip checkpoint evaluation to avoid triggering the full preprocessing chain.
    This allows pairwise_downstream target to skip anchor generation stages.
    """
    # Check if outputs already exist - if so, skip checkpoint to avoid DAG explosion
    candidates_file = f"{ROOT_DIR}/anchors/candidates/{wildcards.organism}"
    touch_file = f"{WORK_DIR}/touch/update_candidates_done_{wildcards.organism}"
    if os.path.exists(candidates_file) and os.path.exists(touch_file):
        # Outputs exist, no need for BLAST indices
        return []

    chunks = get_blast_chunk_files(wildcards)

    # Extract just the filename without extension
    # Files are like: GCF_000001215.4.part-001.fasta
    # We want: GCF_000001215.4.part-001
    chunk_names = [os.path.splitext(f)[0] for f in chunks]

    return expand(
        f"{WORK_DIR}/indices/{{organism}}/{{file}}.fasta_indices",
        organism=wildcards.organism,
        file=chunk_names
    )

rule blast_and_clasp_self:
    """
    BLAST candidate chunks against own genome, run CLASP to chain hits.
    Determines edit distance thresholds for pairwise comparison.

    Input files are discovered dynamically from split_out_forward/reverse directories.
    File wildcard matches the complete chunk filename (e.g., GCF_000001215.4.part-001).
    """
    input:
        # Split FASTA files that will be processed
        fwd_split = lambda wildcards: ancient(f"{WORK_DIR}/split_out_forward/{wildcards.organism}/{wildcards.file}.fasta"),
        rev_split = lambda wildcards: ancient(f"{WORK_DIR}/split_out_reverse/{wildcards.organism}/{wildcards.file}.fasta"),
        # Files actually read by self-blast-clasp.py
        small_meta = lambda wildcards: ancient(f"{ROOT_DIR}/utils/small_meta/{wildcards.organism}"),
        blastdb = lambda wildcards: ancient(f"{ROOT_DIR}/utils/blastdbs/{wildcards.organism}.nsq")
    output:
        # Actual indices file created by self-blast-clasp.py
        # Note: wildcard {file} does NOT include .fasta, but actual file DOES
        indices = f"{WORK_DIR}/indices/{{organism}}/{{file}}.fasta_indices"
    wildcard_constraints:
        file="[^/]+\\.part-\\d+"  # Match full chunk name like: organism.part-###
    log:
        err = f"{WORK_DIR}/logs/first_self/{{organism}}_{{file}}.err",
        out = f"{WORK_DIR}/logs/first_self/{{organism}}_{{file}}.out"
    params:
        ws = WS_ANCHOR_MAKING
    shell:
        f"""
        # Defensive directory creation (normally handled by make_directories.py)
        # Protects against manual deletion or --continue-run with incomplete setup
        mkdir -p {WORK_DIR}/blast_out_forward/{{wildcards.organism}}
        mkdir -p {WORK_DIR}/blast_out_reverse/{{wildcards.organism}}
        mkdir -p {WORK_DIR}/clasp_out_forward/{{wildcards.organism}}
        mkdir -p {WORK_DIR}/clasp_out_reverse/{{wildcards.organism}}
        mkdir -p {WORK_DIR}/indices/{{wildcards.organism}}

        # Run BLAST and CLASP
        cd {CODE_DIR} && {PYTHON} ./self-blast-clasp.py \
            {{wildcards.file}}.fasta \
            {{params.ws}} \
            {WORK_DIR} \
            > {{log.out}} 2> {{log.err}}
        """

# NOTE: No aggregate rule needed!
# update_candidates depends directly on .fasta_indices files and creates indices_global itself
