"""
Common configuration and setup for CAncST Snakefile.
"""

import os
import pathlib
import sys
import yaml

# Pipeline config loading (split sizes, score thresholds, etc.)

def load_pipeline_config():
    """Load pipeline configuration from YAML file."""
    try:
        code_dir = pathlib.Path(__file__).parent.parent.resolve()  # code/rules -> code
        root_dir = code_dir.parent  # code -> CAncST
        config_file = root_dir / 'template' / 'pipeline_config.yaml'

        with open(config_file, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        # Return defaults if config not found
        return {
            'splitting': {'self_blast_part_size': 100000, 'pairwise_part_size': 10000000},
        }

PIPELINE_CONFIG = load_pipeline_config()

# Extract split sizes for fasta-splitter
# NOTE: SELF_BLAST_PART_SIZE is defined here via PIPELINE_CONFIG
# Default: 100000 if template/pipeline_config.yaml missing or doesn't specify
# This value controls how sequences are split for self-BLAST operations
SELF_BLAST_PART_SIZE = PIPELINE_CONFIG.get('splitting', {}).get('self_blast_part_size', 100000)
PAIRWISE_PART_SIZE = PIPELINE_CONFIG.get('splitting', {}).get('pairwise_part_size', 10000000)

# Directories from config
CODE_DIR = config['code_dir']
WORK_DIR = config['work_dir']
ROOT_DIR = config['root_dir']

# Environment setup
os.environ['TMPDIR'] = f'{ROOT_DIR}/utils/tmp'

# Per-organism parameter loading.

# Data structures to hold per-organism parameters
GENMAP_PARAMS = {}  # org -> [{'k': k, 'e': e, 'L': L, 'I': I, 'percentile': P, 'flag': flag}, ...]
MACLE_PARAMS = {}   # org -> [{'w': w, 'p': p, 'percentile': perc, 'flag': flag}, ...]
DUPS_PARAMS = {}    # org -> [[param_list], [param_list], ...]

# Load GenMap parameters (OPTIONAL)
genmap_file = f"{WORK_DIR}/genmap_params.txt"
if os.path.exists(genmap_file):
    with open(genmap_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 7:
                org, k, e, L, I, percentile, param_flag = parts[:7]
                if org not in GENMAP_PARAMS:
                    GENMAP_PARAMS[org] = []
                GENMAP_PARAMS[org].append({
                    'k': k,
                    'e': e,
                    'L': L,
                    'I': I,
                    'percentile': percentile,
                    'flag': param_flag
                })
    total_param_sets = sum(len(params) for params in GENMAP_PARAMS.values())
    print(f"Loaded GenMap parameters: {len(GENMAP_PARAMS)} organisms, {total_param_sets} parameter sets")
else:
    print(f"WARNING: {genmap_file} not found - GenMap filtering will be disabled")

# Load macle parameters (OPTIONAL)
macle_file = f"{WORK_DIR}/macle_params.txt"
if os.path.exists(macle_file):
    with open(macle_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 5:
                org, w, p, percentile, param_flag = parts[:5]
                if org not in MACLE_PARAMS:
                    MACLE_PARAMS[org] = []
                MACLE_PARAMS[org].append({
                    'w': w,
                    'p': p,
                    'percentile': percentile,
                    'flag': param_flag
                })
    total_param_sets = sum(len(params) for params in MACLE_PARAMS.values())
    print(f"Loaded macle parameters: {len(MACLE_PARAMS)} organisms, {total_param_sets} parameter sets")
else:
    print(f"WARNING: {macle_file} not found - macle filtering will be disabled")

# Require at least one parameter file
if not GENMAP_PARAMS and not MACLE_PARAMS:
    raise FileNotFoundError(
        f"At least one parameter file required:\n"
        f"  - {genmap_file} (for GenMap filtering)\n"
        f"  - {macle_file} (for macle filtering)\n"
        f"Both are missing!"
    )

# Load dups parameters (OPTIONAL)
dups_file = f"{WORK_DIR}/dups_params.txt"
if os.path.exists(dups_file):
    with open(dups_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 11:
                org = parts[0]
                if org not in DUPS_PARAMS:
                    DUPS_PARAMS[org] = []
                DUPS_PARAMS[org].append(parts[1:])  # Append each parameter line

# Global parameters (not organism-specific)
if os.path.exists(f"{WORK_DIR}/params"):
    with open(f"{WORK_DIR}/params", "r") as f:
        line = f.read()
        params = line.split()
        if len(params) >= 7:
            WS_ANCHOR_MAKING = params[6]  # BLAST word size for anchor making
        else:
            WS_ANCHOR_MAKING = "11"  # Default
else:
    WS_ANCHOR_MAKING = "11"

# Load all organisms
ALL_ORGS = []
with open(f"{WORK_DIR}/orgs", "r") as f:
    for line in f:
        ALL_ORGS.append(line.strip())

# Load organisms to compute anchors for
COMPUTE_ANCHORS_FOR = []
with open(f"{WORK_DIR}/compute_anchors_for", "r") as f:
    for line in f:
        COMPUTE_ANCHORS_FOR.append(line.strip())

# Note: ORGS_TO_PROCESS variable removed as redundant
# Use COMPUTE_ANCHORS_FOR directly for organisms needing anchor generation

# Optional: Custom pairwise comparisons file
PAIRWISE_FILE = config.get('pairwise_file', None)
if PAIRWISE_FILE:
    print(f"Custom pairwise comparisons file specified: {PAIRWISE_FILE}")

def get_genmap_k(wildcards):
    """Get GenMap K value for organism. Returns 0 if no params (dummy).
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in GENMAP_PARAMS:
        return "0"
    return GENMAP_PARAMS[org][0]['k']  # First parameter set

def get_genmap_e(wildcards):
    """Get GenMap E value for organism. Returns 0 if no params (dummy).
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in GENMAP_PARAMS:
        return "0"
    return GENMAP_PARAMS[org][0]['e']  # First parameter set

def get_genmap_L(wildcards):
    """Get window length L for organism. Returns 0 if no params (dummy).
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in GENMAP_PARAMS:
        return "0"
    return GENMAP_PARAMS[org][0]['L']  # First parameter set

def get_genmap_I(wildcards):
    """Get window interval I for organism. Returns 0 if no params (dummy).
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in GENMAP_PARAMS:
        return "0"
    return GENMAP_PARAMS[org][0]['I']  # First parameter set

def get_genmap_percentile(wildcards):
    """Get percentile threshold for organism. Returns 0 if no params (dummy).
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in GENMAP_PARAMS:
        return "0"
    return GENMAP_PARAMS[org][0]['percentile']  # First parameter set

def get_macle_w(wildcards):
    """Get macle window size W for organism.
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in MACLE_PARAMS:
        raise ValueError(f"No macle parameters found for {org}")
    return MACLE_PARAMS[org][0]['w']  # First parameter set

def get_macle_p(wildcards):
    """Get macle pitch P for organism.
    NOTE: Returns FIRST parameter set (filter_windows reads all sets from file)."""
    org = wildcards.organism
    if org not in MACLE_PARAMS:
        raise ValueError(f"No macle parameters found for {org}")
    return MACLE_PARAMS[org][0]['p']  # First parameter set

def get_genmap_file(wildcards):
    """Get GenMap output file path for organism. Returns None if no params.
    NOTE: Returns FIRST parameter set file (filter_windows reads all from params file)."""
    org = wildcards.organism
    if org not in GENMAP_PARAMS:
        return None
    k = GENMAP_PARAMS[org][0]['k']  # First parameter set
    e = GENMAP_PARAMS[org][0]['e']
    return f"{ROOT_DIR}/utils/genmap_out/{org}/{k}_{e}.freq16"

def get_macle_file(wildcards):
    """Get macle output file path for organism. Returns None if no params.
    NOTE: Returns FIRST parameter set file (filter_windows reads all from params file)."""
    org = wildcards.organism
    if org not in MACLE_PARAMS:
        return None
    w = MACLE_PARAMS[org][0]['w']  # First parameter set
    p = MACLE_PARAMS[org][0]['p']
    return f"{ROOT_DIR}/utils/macle_out/{org}/{w}_{p}.txt"

def get_filter_inputs(wildcards):
    """
    Get ALL preprocessing inputs for a single organism.
    Waits for ALL genmap/macle parameter sets, not just first.

    NOTE: The filter_windows script IGNORES these input arguments and reads
    directly from genmap_params.txt and macle_params.txt. These inputs are
    ONLY for Snakemake dependency tracking to ensure all preprocessing is complete.

    All returned paths are wrapped with ancient() to disable timestamp checking
    when using --force-recompute with copied utils.
    """
    org = wildcards.organism

    inputs = {
        "genome": ancient(f"{ROOT_DIR}/utils/genomes/{org}.fasta"),
        "metadata": ancient(f"{ROOT_DIR}/utils/metadata_genomes/{org}"),
        "small_meta": ancient(f"{ROOT_DIR}/utils/small_meta/{org}"),
        "blastdb": ancient(f"{ROOT_DIR}/utils/blastdbs/{org}.nsq")
    }

    # Collect ALL GenMap files for this organism (for dependency tracking)
    genmap_files = []
    if org in GENMAP_PARAMS:
        for params in GENMAP_PARAMS[org]:  # ALL parameter sets
            k, e = params['k'], params['e']
            genmap_files.append(ancient(f"{ROOT_DIR}/utils/genmap_out/{org}/{k}_{e}.freq16"))

    if genmap_files:
        inputs["genmap"] = genmap_files  # List of all files
    else:
        # No genmap params - use genome as dummy for compatibility
        inputs["genmap"] = [ancient(f"{ROOT_DIR}/utils/genomes/{org}.fasta")]

    # Collect ALL macle files for this organism (for dependency tracking)
    macle_files = []
    if org in MACLE_PARAMS:
        for params in MACLE_PARAMS[org]:  # ALL parameter sets
            w, p = params['w'], params['p']
            macle_files.append(ancient(f"{ROOT_DIR}/utils/macle_out/{org}/{w}_{p}.txt"))

    if macle_files:
        inputs["macle"] = macle_files  # List of all files
    else:
        # No macle params - use genome as dummy for compatibility
        inputs["macle"] = [ancient(f"{ROOT_DIR}/utils/genomes/{org}.fasta")]

    return inputs

# Memory multipliers for resource allocation (genome_size_mb * multiplier * attempt)
# Default: Conservative values for large genomes / memory-constrained systems
# --low-mem: Reduced values (~3x) to allow more parallel jobs when memory isn't bottleneck
# --high-mem: Doubled values (2x) for very large genomes or OOM-prone environments

LOW_MEM = config.get("low_mem", False)
HIGH_MEM = config.get("high_mem", False)

if LOW_MEM:
    # Reduced multipliers for lenient mode
    MEM_MULT = {
        'genmap_index': 10,       # down from 30
        'macle_index': 5,         # down from 20
        'genmap_map': 5,          # down from 15
        'macle_Cm': 3,            # down from 10
        'get_dups': 3,            # down from 10
        'update_candidates': 2,   # down from 5
        'update_compared': 2,     # down from 5
        'update_aligned_stage2': 2,  # down from 5
        'collect_output': 2,      # down from 5
        'resolve_inconsistencies': 2,  # down from 5
        'eval_syn': 2,            # down from 5 (single aligned map + reference copy)
        'finalize_succinct': 4,   # loads own map + partner maps one at a time
        'checks': 4,             # down from 10 (holds 2 aligned maps simultaneously)
        'checks_succinct': 3,    # down from 5 (succinct maps smaller than aligned)
        'statistics': 2,         # down from 5 (single aligned map + small_meta)
        'statistics_succinct': 1, # down from 3 (single succinct map + small_meta)
        'make_blastdb': 1,        # down from 3
        'make_metadata_genomes': 1,  # down from 3
        'make_small_meta': 1,     # down from 2
    }
    print("Using LOW memory multipliers (--low-mem: reduced for more parallelism)")
elif HIGH_MEM:
    # 2x default multipliers for large genomes / OOM-prone jobs
    MEM_MULT = {
        'genmap_index': 60,
        'macle_index': 40,
        'genmap_map': 30,
        'macle_Cm': 20,
        'get_dups': 20,
        'update_candidates': 10,
        'update_compared': 10,
        'update_aligned_stage2': 10,
        'collect_output': 10,
        'resolve_inconsistencies': 10,
        'eval_syn': 10,
        'finalize_succinct': 20,
        'checks': 20,
        'checks_succinct': 10,
        'statistics': 10,
        'statistics_succinct': 6,
        'make_blastdb': 6,
        'make_metadata_genomes': 6,
        'make_small_meta': 4,
    }
    print("Using HIGH memory multipliers (2x default for large genomes)")
else:
    # Default (conservative) multipliers
    MEM_MULT = {
        'genmap_index': 30,
        'macle_index': 20,
        'genmap_map': 15,
        'macle_Cm': 10,
        'get_dups': 10,
        'update_candidates': 5,
        'update_compared': 5,
        'update_aligned_stage2': 5,
        'collect_output': 5,
        'resolve_inconsistencies': 5,
        'eval_syn': 5,            # single aligned map + reference copy
        'finalize_succinct': 10,  # loads own map + partner maps one at a time
        'checks': 10,            # holds 2 aligned maps simultaneously (am1 + am2)
        'checks_succinct': 5,   # holds 2 succinct maps + candidates + small_meta
        'statistics': 5,        # single aligned map + small_meta
        'statistics_succinct': 3, # single succinct map + small_meta
        'make_blastdb': 3,
        'make_metadata_genomes': 3,
        'make_small_meta': 2,
    }

def get_genome_size_mb(wildcards):
    """
    Get genome file size in MB for dynamic resource allocation.
    FAILS LOUDLY if genome file doesn't exist - genomes must be available before running pipeline!
    """
    genome_path = f"{ROOT_DIR}/utils/genomes/{wildcards.organism}.fasta"
    if not os.path.exists(genome_path):
        raise FileNotFoundError(
            f"Genome file not found: {genome_path}\n"
            f"Genomes must be downloaded before running the pipeline!\n"
            f"Organism: {wildcards.organism}"
        )
    size_bytes = os.path.getsize(genome_path)
    size_mb = int(size_bytes / (1024 * 1024))
    return size_mb

def get_organism_from_file(filename):
    """Extract organism name from split file path."""
    # Split files are in split_out_forward/{org}/filename
    # This will be called from rules that need to determine org from file wildcard
    # For now, return None - will be implemented when needed
    return None

def get_split_files_for_org(org, work_dir=WORK_DIR):
    """Get all split files for a given organism."""
    forward_dir = f"{work_dir}/split_out_forward/{org}"
    if os.path.isdir(forward_dir):
        return [name for name in os.listdir(forward_dir)]
    return []

def get_all_split_files(orgs, work_dir=WORK_DIR):
    """Get all split files across all organisms."""
    files = []
    for org in orgs:
        files += get_split_files_for_org(org, work_dir)
    return files

def load_and_validate_custom_pairs(filepath, all_orgs):
    """
    Load custom pairwise comparisons from TSV file.

    File format: org1\torg2 (one pair per line, # for comments)

    Applies:
    - Strict validation (all organisms must exist in all_orgs)
    - Lexicographic normalization (ensures org1 < org2)
    - Automatic deduplication

    Returns: List of (org1, org2) tuples, sorted and unique
    """
    pairs = []
    errors = []

    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) != 2:
                errors.append(f"Line {line_num}: Expected 2 tab-separated columns, got {len(parts)}")
                continue

            org1, org2 = parts[0].strip(), parts[1].strip()

            # Strict validation
            if org1 not in all_orgs:
                errors.append(f"Line {line_num}: Organism '{org1}' not found in orgs file")
            if org2 not in all_orgs:
                errors.append(f"Line {line_num}: Organism '{org2}' not found in orgs file")
            if org1 == org2:
                errors.append(f"Line {line_num}: Self-comparison not allowed: {org1}")

            if errors:
                continue

            # Normalize to lexicographic order (org1 < org2)
            if org1 < org2:
                pairs.append((org1, org2))
            else:
                pairs.append((org2, org1))

    # Fail immediately if any errors
    if errors:
        error_msg = f"Custom pairwise comparisons validation failed ({filepath}):\n" + "\n".join(errors)
        raise ValueError(error_msg)

    # Deduplicate (handles user specifying both A->B and B->A)
    pairs = sorted(list(set(pairs)))

    print(f"Loaded {len(pairs)} unique pairwise comparisons from {filepath}")
    return pairs

def get_pairwise_tuples(compute_for, all_orgs, work_dir=WORK_DIR, custom_pairs_file=None):
    """
    Generate (org1, org2, chunk) tuples for pairwise comparison.

    Two modes:
    1. DEFAULT: Automatic pairs from compute_for vs all_orgs
       - Species in compute_for compared against all_orgs
       - Redundancy elimination: only one direction when both in compute_for
    2. CUSTOM: Explicit pairs from file (overrides default completely)
       - Read pairs from custom_pairs_file
       - Automatic redundancy elimination via lexicographic ordering

    Returns: (org1_list, org2_list, files_list) - three parallel lists
    """
    org1_list = []
    org2_list = []
    files = []

    # Determine which pairs to compute
    if custom_pairs_file is not None and os.path.exists(custom_pairs_file):
        # CUSTOM MODE: Read and validate explicit pairs
        print(f"Custom pairwise mode: using pairs from {custom_pairs_file}")
        pairs_to_compute = load_and_validate_custom_pairs(custom_pairs_file, all_orgs)
    else:
        # DEFAULT MODE: compute_for vs all_orgs with optimization
        print(f"Default pairwise mode: compute_anchors_for ({len(compute_for)}) vs all_orgs ({len(all_orgs)})")
        pairs_to_compute = []
        for org1 in compute_for:
            for org2 in all_orgs:
                if org2 == org1:
                    continue

                # OPTIMIZATION: Skip symmetric pair if both in compute_for
                # parse_bcamm.py writes to BOTH directories, so only need one direction
                if org2 in compute_for and org2 < org1:
                    continue

                pairs_to_compute.append((org1, org2))

        print(f"  Generated {len(pairs_to_compute)} unique pairs (redundancy eliminated)")

    # Enumerate chunks for each pair
    for org1, org2 in pairs_to_compute:
        chunk_dir = f"{work_dir}/sequences_to_compare/{org2}/forward_split"
        if os.path.isdir(chunk_dir):
            chunks = [name.split('forward.')[1]
                     for name in os.listdir(chunk_dir)
                     if name.startswith('forward.')]
            for chunk in chunks:
                org1_list.append(org1)
                org2_list.append(org2)
                files.append(chunk)

    return org1_list, org2_list, files

# All rules use the same conda environment
CONDA_ENV = "snakemake"

# Python interpreter - use the same interpreter running Snakemake
# sys.executable is more robust than PATH resolution (handles venvs, conda, etc.)
PYTHON = sys.executable
