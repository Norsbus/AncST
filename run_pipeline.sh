#!/usr/bin/env bash
#
# AncST Pipeline Runner
# Complete workflow wrapper for synteny anchor detection
#
# Usage: ./run_pipeline.sh [OPTIONS]
#
# This script:
# 1. Optionally downloads genomes from NCBI
# 2. Validates genome files and configurations
# 3. Generates or validates parameters
# 4. Runs the complete Snakemake pipeline
#

set -euo pipefail

# Increase file descriptor limit. Snakemake with parallel jobs needs many file
# handles open simultaneously; raising the soft limit only affects this process.

current_soft=$(ulimit -n)
current_hard=$(ulimit -H -n)

# Target: Maximum safe for 64-bit systems
# - Linux/Mac workstations: 1048576 (1M, typical for modern systems)
# - Compute servers: Often higher (4M-16M)
# - Strategy: Use 1048576 or hard limit, whichever is lower
if [[ "$current_hard" == "unlimited" ]]; then
    target_limit=1048576  # 1M for unlimited systems
elif [[ "$current_hard" -lt 1048576 ]]; then
    target_limit="$current_hard"  # Use hard limit if lower
else
    target_limit=1048576  # 1M is safe for all 64-bit systems
fi

if [[ "$current_soft" != "unlimited" && "$current_soft" -lt "$target_limit" ]]; then
    if ulimit -n "$target_limit" 2>/dev/null; then
        echo "[INFO] Increased file descriptor limit: $current_soft -> $target_limit"
    else
        echo "[WARNING] Could not increase file descriptor limit (current: $current_soft, hard limit: $current_hard)"
        echo "[WARNING] Pipeline may encounter 'Too many open files' errors with many parallel jobs"
        echo "[WARNING] Consider running: ulimit -n $target_limit"
    fi
fi

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$SCRIPT_DIR"

# Default values
WORK_DIR="${PROJECT_ROOT}/template"
CODE_DIR="${PROJECT_ROOT}/code"
GENOMES_DIR="${PROJECT_ROOT}/utils/genomes"  # Fixed location, not configurable
UTILS_DIR="${PROJECT_ROOT}/utils"

SPECIES_FILE=""  # User-provided species file (any name)
ORGS_FILE="${WORK_DIR}/orgs"  # Internal copy (always named "orgs")
COMPUTE_ANCHORS_FILE="${WORK_DIR}/compute_anchors_for"
USER_PROVIDED_COMPUTE_ANCHORS=false

# Environment configuration
ENV_NAME="AncST"  # Default conda environment name
SKIP_CONDA=false  # Whether to skip conda activation
CONTINUE_RUN=false  # Whether to preserve template/ files (resume interrupted run)

# Parameter handling (NEW: cleaner interface)
GENMAP_PARAMS_FILE=""         # Path to manual genmap params (empty = auto-generate)
MACLE_PARAMS_FILE=""          # Path to manual macle params (empty = auto-generate)
DUPS_PARAMS_FILE=""           # Path to manual dups params (empty = auto-generate)
SKIP_GENMAP=false             # Skip genmap auto-generation
SKIP_MACLE=false              # Skip macle auto-generation
SKIP_DUPS=false               # Skip dups auto-generation

CORES=1
MAX_MEM_MB=""     # Maximum total memory in MB (passed to Snakemake --resources)
PROC_MEM_LIMIT=false   # --process-mem-limit: OPT-IN per-process RLIMIT_AS + OOM splitting (off by default)
PROC_TIME_LIMIT=false  # --process-time-limit: OPT-IN per-process timeout + timeout splitting (off by default)
SLURM_MODE=false  # Whether to use SLURM profile instead of --cores
SLURM_CONFIG=""   # Optional path to custom SLURM config.yaml
# Pipeline stage range -- replaces the old --target flag with a true stage-range model.
# FROM_STAGE = where the run starts; TO_STAGE = where it ends (inclusive on both sides).
# Default is the full pipeline. Set via --from-target STAGE or --target-from-to FROM,TO.
FROM_STAGE="preprocessing"
TO_STAGE="downstream"
# Canonical ordered stages. Index matters for validation (FROM <= TO).
# (Plain indexed array -- avoid `declare -A` so the script stays bash-3.2 compatible,
#  a portability target for the macOS dev environment.)
STAGES=(preprocessing anchors pairwise validation downstream)
# Return the index of a stage name in STAGES (echoed on stdout); rc=1 if unknown.
stage_idx() {
    local s=$1 i
    for i in "${!STAGES[@]}"; do
        [[ "${STAGES[$i]}" == "$s" ]] && { echo "$i"; return 0; }
    done
    return 1
}
# Return the snakemake top-level rule that produces the endpoint for a given TO stage.
# (Case statement instead of an associative array -- bash 3.2 compatibility.)
to_rule_for_stage() {
    case "$1" in
        preprocessing) echo "preprocessing_only" ;;
        anchors)       echo "anchors_only" ;;
        pairwise)      echo "pairwise_endpoint" ;;
        validation)    echo "validation_only" ;;
        downstream)    echo "all" ;;
        *)             return 1 ;;
    esac
}
# Return rc 0 if STAGE ($1) lies within the requested [FROM_STAGE, TO_STAGE] range.
# Used to make pre-flight tool validation stage-aware (a downstream-only run does
# not need genmap/macle/blast, etc.). Explicit if/return -- NOT a bare trailing
# `(( ))` -- because under `set -e` a bare arithmetic command that evaluates to 0
# returns rc 1 and could abort the script if this were ever called outside a
# condition context. This form is safe in any context.
stage_in_range() {
    local s_idx from_idx to_idx
    s_idx=$(stage_idx "$1") || return 1
    from_idx=$(stage_idx "$FROM_STAGE")
    to_idx=$(stage_idx "$TO_STAGE")
    if (( s_idx >= from_idx && s_idx <= to_idx )); then
        return 0
    fi
    return 1
}
# Did the user explicitly pass --from-target or --target-from-to ? Used to
# detect "both flags together" error.
FROM_TARGET_SET=false
TARGET_FROM_TO_SET=false

PAIRWISE_FILE=""  # Optional custom pairwise comparisons file

# Snakemake pass-through
SNAKEMAKE_EXTRA_FLAGS=""  # Additional flags passed directly to snakemake command

# Memory configuration
LOW_MEM=false         # Whether to use reduced memory multipliers for rules
HIGH_MEM=false        # Whether to use 2x memory multipliers for rules
AUTO_FIX_CORRUPT=false  # Whether to auto-remove corrupt utils/ files (requires --auto-fix-corrupt)
ABORT_HARD=false  # Whether the filter raises on degenerate genmap/macle budgets instead of degrading (--abort-hard)

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_phase() {
    echo -e "\n${GREEN}=========================================="
    echo -e "$1"
    echo -e "==========================================${NC}"
}

error_exit() {
    log_error "$1"
    exit 1
}

# Portable realpath function (works on macOS and Linux)
# macOS may not have realpath, and readlink -f doesn't work on BSD
get_realpath() {
    local path="$1"
    if [[ -z "$path" ]]; then
        echo ""
        return 1
    fi
    # Try realpath first (Linux, newer macOS)
    if command -v realpath &>/dev/null; then
        realpath "$path" 2>/dev/null && return 0
    fi
    # Fallback: use Python (always available)
    python3 -c "import os; print(os.path.realpath('$path'))" 2>/dev/null && return 0
    # Last resort: just echo the path
    echo "$path"
}

# Usage message
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required:
  --species FILE              Path to species file (one genome per line)
                              Can contain NCBI accessions (GCF_*/GCA_*) or custom genome names
                              Custom genomes must exist in utils/genomes/{name}.fasta

Optional:
  --compute-anchors FILE      Path to compute_anchors_for file
                              (default: same as species file)
  --env-name NAME             Name of conda environment to use (default: AncST)
  --skip-conda                Skip conda activation (if managing env manually)

  Parameters (auto-generate by default):
  --genmap-params FILE       Use manual GenMap parameters from FILE
  --macle-params FILE        Use manual macle parameters from FILE
  --dups-params FILE         Use manual dups parameters from FILE
  --skip-genmap              Skip GenMap parameter generation
  --skip-macle               Skip macle parameter generation
  --skip-dups                Skip dups parameter generation

  NOTE: At least one of genmap or macle must be provided
        (either via manual file or auto-generation)

  Pipeline options:
  --cores N                  Number of cores to use (default: 1, local execution only)
  --max-mem MB               Maximum TOTAL memory in MB across all concurrent jobs,
                             passed to Snakemake job scheduling (--resources mem_mb=MB).
                             Does NOT by itself limit blast/clasp processes; it only
                             becomes the per-process RLIMIT_AS base (max-mem/cores) when
                             you ALSO pass --process-mem-limit.
                             Example: --max-mem 250000 -> Snakemake total budget 250000 MB
  --process-mem-limit        OPT-IN (OFF by default). Cap each blast/clasp process with
                             RLIMIT_AS = max-mem/cores and enable OOM splitting/retry.
                             Base defaults to 100000 MB if --max-mem omitted (warns).
                             Off by default because RLIMIT_AS caps VIRTUAL address space,
                             which blast's mmap'd DB blows through on many-core nodes;
                             prefer letting SLURM cgroups bound real memory instead.
  --process-time-limit       OPT-IN (OFF by default). Enable per-process timeouts and
                             timeout splitting/retry for blast/clasp. Durations come from
                             template/pipeline_config.yaml (timeout_self_minutes=360,
                             timeout_bcamm_minutes=180); edit there to change them.
  --pairwise-comparisons FILE  Custom pairwise comparisons (TSV: org1<TAB>org2)
                               Overrides default pairing logic (default: auto-generate from compute_anchors_for)
  --slurm PROFILE_DIR        Use SLURM cluster execution with profile directory
                             PROFILE_DIR must contain a config.yaml without deprecated flags
                             (Required for SLURM - no default provided)
  Stage range (NEW model -- replaces the old --target flag):
  --from-target STAGE        Start at STAGE, run through end of pipeline.
                             STAGE one of: preprocessing | anchors | pairwise | validation | downstream
                             Default: preprocessing (= full pipeline).
                             Earlier stages must already have produced their outputs
                             (the script fabricates the touch files for them).
  --target-from-to FROM,TO   Run from stage FROM through stage TO (inclusive).
                             Mutually exclusive with --from-target.
                             Examples:
                               --target-from-to pairwise,downstream     # stages 5-9 incl validation
                               --target-from-to anchors,validation       # stages 1-8
                               --target-from-to downstream,downstream    # downstream only
  --continue-run             Resume interrupted run without cleaning template/.
                             AUTO-BEHAVIOUR: injects --rerun-triggers mtime so
                             code-hash changes don't cascade unwanted reruns
                             (overridable by passing --rerun-triggers ... via
                             --snakemake-flags). NO automatic handling of
                             IncompleteFilesException -- snakemake's own error
                             is shown verbatim; for a thorough integrity check
                             of pipeline outputs run:
                               python utils/check_pipeline_state.py <PROJECT_ROOT>
                             WARNING: Assumes template/ is in consistent state.
                             Skips setup_directories.py call entirely.
  --snakemake-flags "FLAGS"  Pass additional flags directly to Snakemake.
                             Enclose in quotes. Common examples:
                               --snakemake-flags "--dry-run"
                               --snakemake-flags "--forcerun rule1 rule2"
                               --snakemake-flags "--forceall"
                               --snakemake-flags "--dry-run --quiet"
  --low-mem                  Use reduced memory multipliers for all rules (~3x lower).
                             Useful when memory isn't the bottleneck and you want
                             more parallel jobs.
  --high-mem                 Use 2x memory multipliers for all rules.
                             Use when jobs are getting OOM-killed on large genomes.
                             Mutually exclusive with --low-mem.
  --auto-fix-corrupt         Automatically remove corrupt/incomplete files in utils/.
                             Without this flag, corrupt files are reported and the
                             pipeline exits with an error (safe default).
                             Use this when you want automatic cleanup after interrupted runs.
  --abort-hard               Strict mode. At several points the pipeline makes a
                             reasonable assumption to keep going when a parameter
                             setting and the input combine into a degenerate case;
                             by default it logs and continues. With this flag those
                             cases raise instead. Use it when you don't fully trust
                             your data/parameters.

  Directories:
  --work-dir DIR             Work directory (default: template/)

  Help:
  -h, --help                 Show this help message

Species File Format:
  One genome per line. Lines can be:
  - NCBI accessions (e.g., GCF_000001215.4, GCA_028878055.3) - will be downloaded
  - Custom genome names (e.g., my_genome) - must exist in utils/genomes/my_genome.fasta
  - Comments starting with # (ignored)
  - Blank lines (ignored)

Examples:
  # Download 3 NCBI genomes and run with auto-generated parameters
  echo -e "GCF_000001215.4\\nGCF_016746365.2\\nGCF_016746395.2" > my_species.txt
  $0 --species my_species.txt

  # Mix NCBI and custom genomes
  cat > mixed_species.txt << 'END'
  # NCBI genomes (will download)
  GCF_000001215.4
  GCA_028878055.3
  # Custom genomes (must be in utils/genomes/)
  my_custom_fly
  another_genome
  END
  $0 --species mixed_species.txt --cores 4

  # Use manual genmap params, auto-generate macle and dups
  $0 --species my_species.txt --genmap-params my_genmap.txt

  # Skip dups entirely (only genmap and macle)
  $0 --species my_species.txt --skip-dups

  # Use custom pairwise comparisons
  cat > my_pairs.txt << 'END'
  GCF_000001215.4	GCF_016746365.2
  GCF_000001215.4	GCA_028878055.3
  END
  $0 --species my_species.txt --pairwise-comparisons my_pairs.txt --cores 4

EOF
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --species)
            SPECIES_FILE="$2"
            shift 2
            ;;
        --compute-anchors)
            COMPUTE_ANCHORS_FILE="$2"
            USER_PROVIDED_COMPUTE_ANCHORS=true
            shift 2
            ;;
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --skip-conda)
            SKIP_CONDA=true
            shift
            ;;
        --genmap-params)
            GENMAP_PARAMS_FILE="$2"
            shift 2
            ;;
        --macle-params)
            MACLE_PARAMS_FILE="$2"
            shift 2
            ;;
        --dups-params)
            DUPS_PARAMS_FILE="$2"
            shift 2
            ;;
        --skip-genmap)
            SKIP_GENMAP=true
            shift
            ;;
        --skip-macle)
            SKIP_MACLE=true
            shift
            ;;
        --skip-dups)
            SKIP_DUPS=true
            shift
            ;;
        --pairwise-comparisons)
            PAIRWISE_FILE="$2"
            shift 2
            ;;
        --cores)
            CORES="$2"
            shift 2
            ;;
        --max-mem)
            MAX_MEM_MB="$2"
            # Validate it's a number
            if ! [[ "$MAX_MEM_MB" =~ ^[0-9]+$ ]]; then
                error_exit "--max-mem requires a positive integer (memory in MB)"
            fi
            shift 2
            ;;
        --process-mem-limit)
            PROC_MEM_LIMIT=true
            shift
            ;;
        --process-time-limit)
            PROC_TIME_LIMIT=true
            shift
            ;;
        --slurm)
            SLURM_MODE=true
            # --slurm now REQUIRES a profile path
            if [[ $# -lt 2 ]] || [[ "$2" =~ ^-- ]]; then
                error_exit "--slurm requires a profile directory path. Usage: --slurm /path/to/profile"
            fi
            # Convert to absolute path immediately to avoid resolution issues
            # If path is already absolute, use as-is; otherwise prepend PWD
            if [[ "$2" = /* ]]; then
                SLURM_CONFIG="$2"
            else
                SLURM_CONFIG="$(pwd)/$2"
            fi
            shift 2
            ;;
        --from-target)
            if [[ "$TARGET_FROM_TO_SET" == "true" ]]; then
                log_error "--from-target and --target-from-to are mutually exclusive"
                exit 1
            fi
            if ! stage_idx "$2" >/dev/null; then
                log_error "Invalid stage: $2"
                log_error "Valid stages: ${STAGES[*]}"
                exit 1
            fi
            FROM_STAGE="$2"
            FROM_TARGET_SET=true
            shift 2
            ;;
        --target-from-to)
            if [[ "$FROM_TARGET_SET" == "true" ]]; then
                log_error "--from-target and --target-from-to are mutually exclusive"
                exit 1
            fi
            # Split on comma; reject if no comma was present
            if [[ "$2" != *","* ]]; then
                log_error "--target-from-to expects FROM,TO (comma-separated), got: $2"
                exit 1
            fi
            _TFT_FROM="${2%%,*}"
            _TFT_TO="${2##*,}"
            if ! stage_idx "$_TFT_FROM" >/dev/null; then
                log_error "Invalid FROM stage: $_TFT_FROM"
                log_error "Valid stages: ${STAGES[*]}"
                exit 1
            fi
            if ! stage_idx "$_TFT_TO" >/dev/null; then
                log_error "Invalid TO stage: $_TFT_TO"
                log_error "Valid stages: ${STAGES[*]}"
                exit 1
            fi
            if (( $(stage_idx "$_TFT_FROM") > $(stage_idx "$_TFT_TO") )); then
                log_error "--target-from-to: FROM ($_TFT_FROM, idx $(stage_idx "$_TFT_FROM")) must come at or before TO ($_TFT_TO, idx $(stage_idx "$_TFT_TO"))"
                exit 1
            fi
            FROM_STAGE="$_TFT_FROM"
            TO_STAGE="$_TFT_TO"
            TARGET_FROM_TO_SET=true
            unset _TFT_FROM _TFT_TO
            shift 2
            ;;
        --continue-run)
            CONTINUE_RUN=true
            shift
            ;;
        --snakemake-flags)
            # Pass any flags directly to Snakemake
            SNAKEMAKE_EXTRA_FLAGS="$2"
            shift 2
            ;;
        --low-mem)
            LOW_MEM=true
            shift
            ;;
        --high-mem)
            HIGH_MEM=true
            shift
            ;;
        --auto-fix-corrupt)
            AUTO_FIX_CORRUPT=true
            shift
            ;;
        --abort-hard)
            ABORT_HARD=true
            shift
            ;;
        --work-dir)
            WORK_DIR="$(cd "$2" 2>/dev/null && pwd)" || error_exit "--work-dir not accessible: $2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Environment activation and tool validation.

if [ "$SKIP_CONDA" = false ]; then
    log_info "Activating conda environment: $ENV_NAME"

    # Initialize conda
    if [[ -f "$HOME/opt/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/opt/anaconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)"
    else
        error_exit "Cannot find conda installation. Please ensure conda is installed."
    fi

    # Activate the environment
    # Temporarily disable errexit and nounset for conda activate
    # (conda scripts aren't compatible with strict mode)
    set +eu
    conda activate "$ENV_NAME" 2>&1
    ACTIVATE_STATUS=$?
    set -eu

    if [ $ACTIVATE_STATUS -ne 0 ]; then
        error_exit "Failed to activate conda environment '$ENV_NAME'. Run './install.sh --env-name $ENV_NAME' first."
    fi

    log_success "Environment '$ENV_NAME' activated"

    # Debug PATH and binary locations (helpful for troubleshooting)
    log_info "Debugging PATH setup:"
    log_info "  PATH=$PATH"
    log_info "  macle location: $(command -v macle || echo 'NOT FOUND')"
    log_info "  clasp.x location: $(command -v clasp.x || echo 'NOT FOUND')"
    log_info "  CONDA_PREFIX=$CONDA_PREFIX"
else
    log_info "Skipping conda activation (--skip-conda specified)"
fi

# Validate all required tools
log_info "Validating required tools..."

check_tool() {
    local tool="$1"
    local check_cmd="$2"

    if ! eval "$check_cmd" &>/dev/null; then
        log_error "Tool '$tool' not found or not working"
        return 1
    fi
    log_success "$tool"
    return 0
}

# Track validation failures
VALIDATION_FAILED=false

# Determine which param families will be available, from the parsed flags alone
# (manual params file given, OR not skipped -> auto-generated). Single source of
# truth: these booleans gate BOTH the stage-aware tool validation below AND the
# "at least one of genmap or macle" check further down. They depend only on CLI
# flags, so computing them here (before validation) is safe.
GENMAP_AVAILABLE=false
if [[ -n "$GENMAP_PARAMS_FILE" ]] || [[ "$SKIP_GENMAP" != "true" ]]; then
    GENMAP_AVAILABLE=true
fi
MACLE_AVAILABLE=false
if [[ -n "$MACLE_PARAMS_FILE" ]] || [[ "$SKIP_MACLE" != "true" ]]; then
    MACLE_AVAILABLE=true
fi
# dups also drives the genmap binary (dups_params -> genmap freq16 outputs), so it
# must be considered when deciding whether genmap needs validating.
DUPS_AVAILABLE=false
if [[ -n "$DUPS_PARAMS_FILE" ]] || [[ "$SKIP_DUPS" != "true" ]]; then
    DUPS_AVAILABLE=true
fi

# Core Python packages -- needed by every stage, including downstream.
check_tool "Python" "python --version" || VALIDATION_FAILED=true
check_tool "NumPy" "python -c 'import numpy'" || VALIDATION_FAILED=true
check_tool "BioPython" "python -c 'import Bio'" || VALIDATION_FAILED=true
check_tool "PyYAML" "python -c 'import yaml'" || VALIDATION_FAILED=true
check_tool "Snakemake" "snakemake --version" || VALIDATION_FAILED=true

# Bioinformatics tools -- validated ONLY for the stages that actually invoke them,
# so e.g. a downstream-only or validation-only run (which use neither) is not blocked
# by a missing genmap/macle/blast. Stage->binary map verified against the rule files.

# datasets: only used to download NCBI genomes during preprocessing. Resume runs
# (--from-target > preprocessing) require pre-existing genomes, so it isn't needed.
if stage_in_range preprocessing; then
    check_tool "datasets (NCBI)" "datasets --version" || VALIDATION_FAILED=true
fi

# genmap/macle build indices in preprocessing, but only run when their params exist,
# which requires the anchors (filter) stage to be in range AND the family enabled.
# genmap is additionally driven by dups params (dups_params -> genmap freq16 outputs).
if stage_in_range preprocessing && stage_in_range anchors; then
    if [[ "$GENMAP_AVAILABLE" == "true" || "$DUPS_AVAILABLE" == "true" ]]; then
        check_tool "genmap" "genmap --version" || VALIDATION_FAILED=true
    fi
    if [[ "$MACLE_AVAILABLE" == "true" ]]; then
        # Robust check: macle's help flag is -h (--help is undefined and returns
        # non-zero on some builds); grep the output and ignore the exit code, like
        # the clasp.x check below. Catches a missing/broken binary, build-portable.
        check_tool "macle" "( macle -h 2>&1 || true ) | grep -qi 'macle\|usage'" || VALIDATION_FAILED=true
    fi
fi

# makeblastdb: genome BLASTdbs in preprocessing (unconditional for ALL_ORGS) and
# candidate BLASTdbs in the pairwise stage.
if stage_in_range preprocessing || stage_in_range pairwise; then
    check_tool "makeblastdb" "makeblastdb -version" || VALIDATION_FAILED=true
fi

# blastn/clasp.x/fasta-splitter: used by self-blast (anchors) and pairwise matching.
if stage_in_range anchors || stage_in_range pairwise; then
    check_tool "BLAST" "blastn -version" || VALIDATION_FAILED=true
    check_tool "clasp.x" "( clasp.x 2>&1 || true ) | grep -qi 'usage\|clasp'" || VALIDATION_FAILED=true
    check_tool "fasta-splitter" "( fasta-splitter --help 2>&1 || true ) | grep -qi 'usage\|fasta'" || VALIDATION_FAILED=true
fi

if [[ "$VALIDATION_FAILED" == "true" ]]; then
    log_error "Some tools are missing or not working"
    log_info "Please run './install.sh --env-name $ENV_NAME' to set up the environment"
    exit 1
fi

log_success "All required tools validated"

# Validate required arguments
if [[ -z "$SPECIES_FILE" ]]; then
    error_exit "--species FILE is required"
fi

if [[ "$LOW_MEM" == "true" && "$HIGH_MEM" == "true" ]]; then
    error_exit "--low-mem and --high-mem are mutually exclusive"
fi

if [[ ! -f "$SPECIES_FILE" ]]; then
    error_exit "Species file not found: $SPECIES_FILE"
fi

# Validate parameter configuration
# GENMAP_AVAILABLE / MACLE_AVAILABLE / DUPS_AVAILABLE are computed once above
# (before tool validation), from the parsed flags. Reuse them here.

# At least one of genmap or macle must be available
if [[ "$GENMAP_AVAILABLE" != "true" ]] && [[ "$MACLE_AVAILABLE" != "true" ]]; then
    error_exit "At least one of genmap or macle parameters must be provided (either via manual file or auto-generation). You cannot skip both."
fi

# Species file validation and processing.

validate_and_process_species() {
    log_info "Validating species file: $SPECIES_FILE"

    # Arrays to track genomes
    local -a ncbi_accessions=()
    local -a custom_genomes=()
    local -a errors=()
    local -a all_genomes=()

    # NCBI accession pattern: GC[AF]_[0-9]{9}.[0-9]+
    local ncbi_pattern='^GC[AF]_[0-9]{9}\.[0-9]+$'

    # Parse species file
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Skip empty lines and comments
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue

        # Trim whitespace
        genome=$(echo "$line" | tr -d '[:space:]')
        [[ -z "$genome" ]] && continue

        all_genomes+=("$genome")

        # Check if it's an NCBI accession
        if [[ "$genome" =~ $ncbi_pattern ]]; then
            ncbi_accessions+=("$genome")
        else
            custom_genomes+=("$genome")
        fi
    done < "$SPECIES_FILE"

    # Check if file is empty
    if [[ ${#all_genomes[@]} -eq 0 ]]; then
        error_exit "Species file is empty or contains only comments: $SPECIES_FILE"
    fi

    log_info "Found ${#all_genomes[@]} genomes: ${#ncbi_accessions[@]} NCBI, ${#custom_genomes[@]} custom"

    # Validate NCBI accessions - skip if genome already exists locally
    if [[ ${#ncbi_accessions[@]} -gt 0 ]]; then
        log_info "Validating ${#ncbi_accessions[@]} NCBI accessions..."
        for acc in "${ncbi_accessions[@]}"; do
            # First check if genome already exists locally (any extension)
            local exists_locally=false
            for ext in fasta fna fa faa; do
                if [[ -f "${GENOMES_DIR}/${acc}.${ext}" ]]; then
                    exists_locally=true
                    log_info "  $acc - found locally, skipping NCBI check"
                    break
                fi
            done

            # Only query NCBI if not found locally
            if [[ "$exists_locally" == "false" ]]; then
                if ! datasets summary genome accession "$acc" 2>/dev/null | grep -q '"total_count":[[:space:]]*[1-9]'; then
                    errors+=("NCBI accession not found or invalid: $acc")
                fi
            fi
        done
    fi

    # Validate custom genomes exist in utils/genomes/
    # Check for multiple possible extensions: .fasta, .fna, .fa, .faa
    if [[ ${#custom_genomes[@]} -gt 0 ]]; then
        log_info "Validating ${#custom_genomes[@]} custom genomes..."
        for genome in "${custom_genomes[@]}"; do
            local found=false
            for ext in fasta fna fa faa; do
                if [[ -f "${GENOMES_DIR}/${genome}.${ext}" ]]; then
                    found=true
                    break
                fi
            done
            if [[ "$found" == "false" ]]; then
                errors+=("Custom genome not found: ${GENOMES_DIR}/${genome}.{fasta,fna,fa,faa}")
            fi
        done
    fi

    # Report all errors at once
    if [[ ${#errors[@]} -gt 0 ]]; then
        log_error "===== Species Validation Errors ====="
        for err in "${errors[@]}"; do
            log_error "  - $err"
        done
        log_error "======================================"
        error_exit "Species file validation failed with ${#errors[@]} error(s)"
    fi

    log_success "All ${#all_genomes[@]} genomes validated successfully"

    # Copy species file to orgs (internal name), sorted + deduped.
    # LC_ALL=C gives byte-wise order identical on Linux/Mac/Intel/ARM.
    # sort -o handles the case where $SPECIES_FILE == $ORGS_FILE (reads fully before writing).
    if [[ "$(get_realpath "$SPECIES_FILE")" != "$(get_realpath "$ORGS_FILE")" ]]; then
        LC_ALL=C sort -u "$SPECIES_FILE" -o "$ORGS_FILE"
        log_info "Sorted and copied species file to: $ORGS_FILE"
    else
        LC_ALL=C sort -u "$ORGS_FILE" -o "$ORGS_FILE"
        log_info "Sorted species file in place: $ORGS_FILE"
    fi

    # Copy orgs to utils/ locations (real copies, not symlinks)
    # Required by utils/make_directories.py and utils/util_code/make_anchor_directories.py
    # Skip for --continue-run to preserve existing state
    if [[ "$CONTINUE_RUN" != "true" ]]; then
        mkdir -p "${UTILS_DIR}" "${UTILS_DIR}/util_code"
        rm -f "${UTILS_DIR}/orgs" "${UTILS_DIR}/util_code/orgs"
        cp "$ORGS_FILE" "${UTILS_DIR}/orgs"
        cp "$ORGS_FILE" "${UTILS_DIR}/util_code/orgs"
        log_info "Copied orgs to utils/orgs and utils/util_code/orgs"
    fi

    # Download NCBI genomes if any
    if [[ ${#ncbi_accessions[@]} -gt 0 ]]; then
        download_ncbi_genomes "${ncbi_accessions[@]}"
    fi
}

download_ncbi_genomes() {
    local accessions=("$@")
    log_info "Downloading ${#accessions[@]} NCBI genomes..."

    mkdir -p "$GENOMES_DIR"

    # Create temporary accessions file
    local temp_acc_file="${WORK_DIR}/temp_ncbi_accessions.txt"
    printf "%s\n" "${accessions[@]}" > "$temp_acc_file"

    # Download using get_genomes_mp.py
    local download_script="${UTILS_DIR}/util_code/get_genomes_mp.py"
    if [[ ! -f "$download_script" ]]; then
        error_exit "Download script not found: $download_script"
    fi

    python "$download_script" "$temp_acc_file" "$GENOMES_DIR" || error_exit "Failed to download NCBI genomes"

    # Extract the downloaded genomes
    local extract_script="${UTILS_DIR}/util_code/extract_fastas_mp.py"
    if [[ ! -f "$extract_script" ]]; then
        error_exit "Extract script not found: $extract_script"
    fi

    log_info "Extracting NCBI genomes..."
    python "$extract_script" "$temp_acc_file" "$GENOMES_DIR" || error_exit "Failed to extract NCBI genomes"

    rm -f "$temp_acc_file"
    log_success "NCBI genomes downloaded and extracted successfully"
}

# Header
log_info "=========================================="
log_info "AncST Pipeline Runner"
log_info "=========================================="
log_info "Work directory: $WORK_DIR"
log_info "Genomes directory: $GENOMES_DIR"
log_info "Cores: $CORES"

# Step 1: Validate and standardize genome files
validate_genomes() {
    log_info "Validating genome files for this job..."

    # Check if genomes directory exists (follows the symlink in the slot layout).
    if [[ ! -d "$GENOMES_DIR" ]]; then
        error_exit "Genomes directory not found: $GENOMES_DIR"
    fi
    # ORGS_FILE is written by validate_and_process_species, which runs first.
    if [[ ! -f "$ORGS_FILE" ]]; then
        error_exit "Orgs file not found (expected from validate_and_process_species): $ORGS_FILE"
    fi

    # Validate ONLY the genomes this job needs (the orgs in ORGS_FILE), not the
    # whole utils/genomes directory. In the comp-server slot layout that directory
    # is a symlink into a cache shared across slots, so scanning/renaming all of it
    # every job was both wasteful and a cross-slot mutation hazard. Each org name
    # equals its on-disk basename ${org}.fasta (same contract used by
    # validate_and_process_species when it locates ${org}.{fasta,fna,fa,faa}).
    local count=0
    local -a missing=()
    while IFS= read -r org || [[ -n "$org" ]]; do
        # Mirror validate_and_process_species parsing: skip comments/blanks, trim.
        [[ "$org" =~ ^[[:space:]]*# ]] && continue
        org=$(echo "$org" | tr -d '[:space:]')
        [[ -z "$org" ]] && continue

        local target="${GENOMES_DIR}/${org}.fasta"
        if [[ ! -f "$target" ]]; then
            # Standardize a non-.fasta extension for THIS org only. Server-side
            # extraction already yields .fasta (so this is a no-op there); this
            # keeps standalone use working when a user supplies .fna/.fa/.faa.
            local found=""
            for ext in fna fa faa; do
                if [[ -f "${GENOMES_DIR}/${org}.${ext}" ]]; then
                    found="${GENOMES_DIR}/${org}.${ext}"
                    break
                fi
            done
            if [[ -n "$found" ]]; then
                log_info "Renaming: $(basename "$found") -> ${org}.fasta"
                mv "$found" "$target"
            else
                missing+=("$org")
                continue
            fi
        fi

        # Basic FASTA format check (warning only).
        if ! head -1 "$target" | grep -q "^>"; then
            log_warning "File may not be valid FASTA: ${org}.fasta"
        fi
        ((count++)) || true
    done < "$ORGS_FILE"

    if [[ ${#missing[@]} -gt 0 ]]; then
        error_exit "Missing genome file(s) for this job: ${missing[*]} (looked in $GENOMES_DIR for .fasta/.fna/.fa/.faa)"
    fi
    if [[ $count -eq 0 ]]; then
        error_exit "No genomes required by this job (orgs file empty?): $ORGS_FILE"
    fi

    log_success "Validated $count genome file(s) for this job"
}

# Step 3: Create or validate orgs file

# Step 4: Setup compute_anchors_for file
setup_compute_anchors() {
    log_info "Setting up compute_anchors_for file..."

    # ALWAYS write to the work dir ($WORK_DIR), whether user-provided or default
    # This ensures Snakemake reads the correct file (not a stale placeholder)
    # Sorted + deduped for deterministic pair enumeration (LC_ALL=C for cross-platform byte-wise order).
    if [[ "$USER_PROVIDED_COMPUTE_ANCHORS" == "true" ]]; then
        log_info "Using user-provided compute_anchors_for file"
        LC_ALL=C sort -u "$COMPUTE_ANCHORS_FILE" -o "${WORK_DIR}/compute_anchors_for"
    else
        log_info "Using orgs file as compute_anchors_for"
        cp "$ORGS_FILE" "${WORK_DIR}/compute_anchors_for"
    fi

    # Validate compute_anchors_for is subset of orgs
    log_info "Validating that compute_anchors_for is a subset of species..."

    # Collect all valid organisms from species file
    local all_species=()
    while IFS= read -r org || [[ -n "$org" ]]; do
        [[ -z "$org" ]] && continue
        [[ "$org" =~ ^# ]] && continue
        all_species+=("$org")
    done < "$ORGS_FILE"

    # Check each organism in compute_anchors_for
    local compute_list=()
    local invalid_orgs=()
    while IFS= read -r org || [[ -n "$org" ]]; do
        [[ -z "$org" ]] && continue
        [[ "$org" =~ ^# ]] && continue
        compute_list+=("$org")

        if ! grep -qx "$org" "$ORGS_FILE"; then
            invalid_orgs+=("$org")
        fi
    done < "$COMPUTE_ANCHORS_FILE"

    # If any invalid organisms found, show verbose error
    if [ ${#invalid_orgs[@]} -gt 0 ]; then
        log_error "ERROR: compute_anchors_for contains organisms not present in species file!"
        echo ""
        echo "  Invalid organisms in compute_anchors_for:"
        for org in "${invalid_orgs[@]}"; do
            echo "    - $org"
        done
        echo ""
        echo "  Available organisms in species file (${#all_species[@]} total):"
        for org in "${all_species[@]}"; do
            echo "    - $org"
        done
        echo ""
        echo "  Organisms to compute anchors for must be a subset of the species file."
        echo "  Either:"
        echo "    1. Add missing organisms to your species file"
        echo "    2. Remove invalid organisms from compute_anchors_for file"
        echo "    3. Don't use --compute-anchors flag (defaults to all species)"
        echo ""
        error_exit "Validation failed: compute_anchors_for is not a subset of species"
    fi

    log_success "Subset validation passed: ${#compute_list[@]} organisms to compute anchors for (out of ${#all_species[@]} total species)"

    log_success "compute_anchors_for validated: $COMPUTE_ANCHORS_FILE"
}

# Validate parameter file format
validate_params_file() {
    local file=$1
    local type=$2  # genmap, macle, or dups
    local expected_fields=0

    case $type in
        genmap) expected_fields=7 ;;
        macle) expected_fields=5 ;;
        dups) expected_fields=11 ;;
    esac

    log_info "Validating $type parameters: $file"

    # Check file exists
    if [[ ! -f "$file" ]]; then
        error_exit "$type parameters file not found: $file"
    fi

    # Load valid organism names
    local valid_orgs=()
    while IFS= read -r org || [[ -n "$org" ]]; do
        [[ -z "$org" ]] && continue
        [[ "$org" =~ ^# ]] && continue
        valid_orgs+=("$org")
    done < "$ORGS_FILE"

    # Validate each line
    local line_num=0
    while IFS= read -r line || [[ -n "$line" ]]; do
        ((++line_num))

        # Skip empty lines and comments
        [[ -z "$line" ]] && continue
        [[ "$line" =~ ^# ]] && continue

        # Split into fields
        read -ra fields <<< "$line"
        local field_count=${#fields[@]}

        # Check field count
        if [[ $field_count -ne $expected_fields ]]; then
            error_exit "$type params line $line_num: Expected $expected_fields fields, got $field_count"
        fi

        local org="${fields[0]}"

        # Check organism name is valid
        local org_valid=false
        for valid_org in "${valid_orgs[@]}"; do
            if [[ "$org" == "$valid_org" ]]; then
                org_valid=true
                break
            fi
        done
        if [[ "$org_valid" == "false" ]]; then
            log_warning "$type params line $line_num: Organism '$org' not in orgs file"
        fi

        # Type-specific validation
        case $type in
            genmap)
                # genmap: org k e L I percentile flag
                local window_l="${fields[3]}"
                local pitch="${fields[4]}"
                if ! [[ "$window_l" =~ ^[0-9]+$ ]] || ! [[ "$pitch" =~ ^[0-9]+$ ]]; then
                    error_exit "$type params line $line_num: L and I must be integers"
                fi
                if [[ $pitch -ge $window_l ]]; then
                    error_exit "$type params line $line_num: pitch ($pitch) must be < window_length ($window_l)"
                fi
                ;;

            macle)
                # macle: org w p percentile flag
                local window_w="${fields[1]}"
                local pitch="${fields[2]}"
                if ! [[ "$window_w" =~ ^[0-9]+$ ]] || ! [[ "$pitch" =~ ^[0-9]+$ ]]; then
                    error_exit "$type params line $line_num: w and p must be integers"
                fi
                if [[ $pitch -ge $window_w ]]; then
                    error_exit "$type params line $line_num: pitch ($pitch) must be < window_size ($window_w)"
                fi
                ;;

            dups)
                # dups: org k1 e1 k2 e2 w p x1 y1 x2 y2
                # Check y1, y2 are percentages (0-100)
                local y1="${fields[8]}"
                local y2="${fields[10]}"
                if ! [[ "$y1" =~ ^[0-9]+$ ]] || [[ $y1 -lt 0 ]] || [[ $y1 -gt 100 ]]; then
                    error_exit "$type params line $line_num: y1 must be 0-100, got: $y1"
                fi
                if ! [[ "$y2" =~ ^[0-9]+$ ]] || [[ $y2 -lt 0 ]] || [[ $y2 -gt 100 ]]; then
                    error_exit "$type params line $line_num: y2 must be 0-100, got: $y2"
                fi

                local window_w="${fields[5]}"
                local pitch="${fields[6]}"
                if ! [[ "$window_w" =~ ^[0-9]+$ ]] || ! [[ "$pitch" =~ ^[0-9]+$ ]]; then
                    error_exit "$type params line $line_num: w and p must be integers"
                fi
                if [[ $pitch -ge $window_w ]]; then
                    error_exit "$type params line $line_num: pitch ($pitch) must be < window_size ($window_w)"
                fi
                ;;
        esac
    done < "$file"

    log_success "$type parameters validated"
}

# Step 5: Generate or validate parameters
setup_parameters() {
    log_info "Setting up parameters..."

    # Check if generate_params.py exists
    local GENERATE_SCRIPT="${WORK_DIR}/generate_params.py"
    if [[ ! -f "$GENERATE_SCRIPT" ]]; then
        GENERATE_SCRIPT="${CODE_DIR}/generate_params.py"
    fi

    # --- GenMap Parameters ---
    local GENMAP_FINAL="${WORK_DIR}/genmap_params.txt"
    if [[ -n "$GENMAP_PARAMS_FILE" ]]; then
        # Manual file provided
        log_info "Using manual GenMap parameters: $GENMAP_PARAMS_FILE"
        if [[ ! -f "$GENMAP_PARAMS_FILE" ]]; then
            error_exit "GenMap parameters file not found: $GENMAP_PARAMS_FILE"
        fi
        validate_params_file "$GENMAP_PARAMS_FILE" "genmap"
        if [[ "$(get_realpath "$GENMAP_PARAMS_FILE")" != "$(get_realpath "$GENMAP_FINAL" 2>/dev/null)" ]]; then
            cp "$GENMAP_PARAMS_FILE" "$GENMAP_FINAL"
        fi
    elif [[ "$SKIP_GENMAP" != "true" ]]; then
        # Auto-generate
        log_info "Auto-generating GenMap parameters..."
        if [[ ! -f "$GENERATE_SCRIPT" ]]; then
            error_exit "generate_params.py not found"
        fi
        python "$GENERATE_SCRIPT" --genmap --orgs "$COMPUTE_ANCHORS_FILE" --output "$GENMAP_FINAL" --work-dir "$WORK_DIR" --genomes-dir "$GENOMES_DIR" || error_exit "Failed to generate GenMap parameters"
        log_success "GenMap parameters generated: $GENMAP_FINAL"
    else
        log_info "Skipping GenMap parameters"
        # Remove any stale genmap_params.txt so common.smk doesn't load it and
        # 00_preprocessing enumerate genmap targets for orgs in the stale file.
        if [[ -f "$GENMAP_FINAL" ]]; then
            log_warning "Removing stale $GENMAP_FINAL to honor --skip-genmap"
            rm -f "$GENMAP_FINAL"
        fi
    fi

    # --- macle Parameters ---
    local MACLE_FINAL="${WORK_DIR}/macle_params.txt"
    if [[ -n "$MACLE_PARAMS_FILE" ]]; then
        # Manual file provided
        log_info "Using manual macle parameters: $MACLE_PARAMS_FILE"
        if [[ ! -f "$MACLE_PARAMS_FILE" ]]; then
            error_exit "macle parameters file not found: $MACLE_PARAMS_FILE"
        fi
        validate_params_file "$MACLE_PARAMS_FILE" "macle"
        if [[ "$(get_realpath "$MACLE_PARAMS_FILE")" != "$(get_realpath "$MACLE_FINAL" 2>/dev/null)" ]]; then
            cp "$MACLE_PARAMS_FILE" "$MACLE_FINAL"
        fi
    elif [[ "$SKIP_MACLE" != "true" ]]; then
        # Auto-generate
        log_info "Auto-generating macle parameters..."
        if [[ ! -f "$GENERATE_SCRIPT" ]]; then
            error_exit "generate_params.py not found"
        fi
        python "$GENERATE_SCRIPT" --macle --orgs "$COMPUTE_ANCHORS_FILE" --output "$MACLE_FINAL" --work-dir "$WORK_DIR" --genomes-dir "$GENOMES_DIR" || error_exit "Failed to generate macle parameters"
        log_success "macle parameters generated: $MACLE_FINAL"
    else
        log_info "Skipping macle parameters"
        # Remove any stale macle_params.txt so common.smk doesn't load it and
        # 00_preprocessing enumerate macle targets for orgs in the stale file.
        if [[ -f "$MACLE_FINAL" ]]; then
            log_warning "Removing stale $MACLE_FINAL to honor --skip-macle"
            rm -f "$MACLE_FINAL"
        fi
    fi

    # --- Dups Parameters ---
    local DUPS_FINAL="${WORK_DIR}/dups_params.txt"
    if [[ -n "$DUPS_PARAMS_FILE" ]]; then
        # Manual file provided
        log_info "Using manual dups parameters: $DUPS_PARAMS_FILE"
        if [[ ! -f "$DUPS_PARAMS_FILE" ]]; then
            error_exit "dups parameters file not found: $DUPS_PARAMS_FILE"
        fi
        validate_params_file "$DUPS_PARAMS_FILE" "dups"
        if [[ "$(get_realpath "$DUPS_PARAMS_FILE")" != "$(get_realpath "$DUPS_FINAL" 2>/dev/null)" ]]; then
            cp "$DUPS_PARAMS_FILE" "$DUPS_FINAL"
        fi
    elif [[ "$SKIP_DUPS" != "true" ]]; then
        # Auto-generate
        log_info "Auto-generating dups parameters..."
        if [[ ! -f "$GENERATE_SCRIPT" ]]; then
            error_exit "generate_params.py not found"
        fi
        python "$GENERATE_SCRIPT" --dups --orgs "$COMPUTE_ANCHORS_FILE" --output "$DUPS_FINAL" --work-dir "$WORK_DIR" --genomes-dir "$GENOMES_DIR" || error_exit "Failed to generate dups parameters"
        log_success "dups parameters generated: $DUPS_FINAL"
    else
        log_info "Skipping dups parameters"
        # Remove any existing dups_params.txt from a previous run so Snakemake
        # rules (get_dups_genmap_inputs, has_dups_params) don't see it and
        # request dups-related genmap files
        if [[ -f "$DUPS_FINAL" ]]; then
            log_info "Removing existing dups_params.txt (--skip-dups active)"
            rm "$DUPS_FINAL"
        fi
    fi

    # config lives in the work dir; write the runtime cores/mem/timeout into it
    PIPELINE_CONFIG_DEST="${WORK_DIR}/pipeline_config.yaml"
    sed -i.bak "s/^  cores:.*/  cores: $CORES/" "$PIPELINE_CONFIG_DEST"
    sed -i.bak "s/^  abort_hard:.*/  abort_hard: $ABORT_HARD/" "$PIPELINE_CONFIG_DEST"

    # Per-process limits for blast/clasp are OPT-IN (off by default). RLIMIT_AS caps VIRTUAL
    # address space, which blast's mmap'd DB + glibc thread arenas blow through on many-core
    # nodes -> spurious failures; SLURM cgroups bound real RSS instead. Enable explicitly.
    if [[ "$PROC_MEM_LIMIT" == "true" ]]; then
        EFFECTIVE_MEM=${MAX_MEM_MB:-100000}
        if [[ -z "$MAX_MEM_MB" ]]; then
            log_warning "--process-mem-limit set without --max-mem; using default base ${EFFECTIVE_MEM}MB (pass --max-mem to control)"
        fi
        log_info "Per-process memory limit ENABLED: ${EFFECTIVE_MEM}MB / ${CORES} cores = $((EFFECTIVE_MEM / CORES))MB (RLIMIT_AS)"
        sed -i.bak "s/^  max_mem_mb:.*/  max_mem_mb: $EFFECTIVE_MEM/" "$PIPELINE_CONFIG_DEST"
    else
        log_info "Per-process memory limit OFF (default; enable with --process-mem-limit)"
        sed -i.bak "s/^  max_mem_mb:.*/  max_mem_mb: null/" "$PIPELINE_CONFIG_DEST"
    fi

    if [[ "$PROC_TIME_LIMIT" == "true" ]]; then
        # keep the template's timeout_self_minutes / timeout_bcamm_minutes (edit there to change)
        log_info "Per-process timeouts ENABLED (durations from pipeline_config.yaml)"
    else
        log_info "Per-process time limit OFF (default; enable with --process-time-limit)"
        sed -i.bak "s/^  timeout_self_minutes:.*/  timeout_self_minutes: null/" "$PIPELINE_CONFIG_DEST"
        sed -i.bak "s/^  timeout_bcamm_minutes:.*/  timeout_bcamm_minutes: null/" "$PIPELINE_CONFIG_DEST"
    fi

    rm -f "${PIPELINE_CONFIG_DEST}.bak"
}

# Step 5b: Setup custom pairwise comparisons (if provided)
setup_pairwise_comparisons() {
    if [[ -n "$PAIRWISE_FILE" ]]; then
        log_info "Setting up custom pairwise comparisons..."

        # Validate file exists
        if [[ ! -f "$PAIRWISE_FILE" ]]; then
            error_exit "Pairwise comparisons file not found: $PAIRWISE_FILE"
        fi

        # Validate file format (basic check)
        if ! grep -qE $'^\t|^[^#\t]+\t[^#\t]+' "$PAIRWISE_FILE"; then
            log_warning "Pairwise file may have formatting issues (expected: org1<TAB>org2)"
        fi

        # Copy to work directory
        local PAIRWISE_DEST="${WORK_DIR}/pairwise_comparisons.txt"
        if [[ "$(get_realpath "$PAIRWISE_FILE" 2>/dev/null)" != "$(get_realpath "$PAIRWISE_DEST" 2>/dev/null)" ]]; then
            cp "$PAIRWISE_FILE" "$PAIRWISE_DEST"
            log_success "Custom pairwise comparisons copied to: $PAIRWISE_DEST"
        else
            log_info "Pairwise comparisons file already in work directory"
        fi
    fi
}

# Step 5d: Validate integrity of existing utils/ preprocessing files
# Detects corrupt/partial files left by interrupted runs.
# Behavior depends on --auto-fix-corrupt flag:
#   WITH flag:    removes corrupt files so Snakemake regenerates them
#   WITHOUT flag: reports corrupt files and exits with error (safe default)
# Works on both Linux and Mac (POSIX-compatible checks only).
validate_utils_integrity() {
    log_phase "Validating utils/ file integrity"

    local utils_dir="${PROJECT_ROOT}/utils"
    local corrupt_count=0
    # Collect corrupt file paths for reporting/removal
    local corrupt_files=()
    # Collect genmap index dirs separately (need rm -rf, not rm -f)
    local corrupt_dirs=()
    # Collect blastdb base paths separately (need wildcard removal)
    local corrupt_blastdb_bases=()

    # --- Check macle indices for 0-byte files ---
    if [[ -d "${utils_dir}/macle_indices" ]]; then
        for f in "${utils_dir}/macle_indices/"*; do
            [[ -e "$f" ]] || continue  # skip if glob matched nothing
            [[ -d "$f" ]] && continue  # skip directories
            if [[ ! -s "$f" ]]; then
                log_warning "Corrupt macle index (0 bytes): $(basename "$f")"
                corrupt_files+=("$f")
                corrupt_count=$((corrupt_count + 1))
            fi
        done
    fi

    # --- Check macle outputs for 0-byte files ---
    if [[ -d "${utils_dir}/macle_out" ]]; then
        for org_dir in "${utils_dir}/macle_out/"*/; do
            [[ -d "$org_dir" ]] || continue
            for f in "${org_dir}"*.txt; do
                [[ -e "$f" ]] || continue
                if [[ ! -s "$f" ]]; then
                    log_warning "Corrupt macle output (0 bytes): $(basename "$(dirname "$f")")/$(basename "$f")"
                    corrupt_files+=("$f")
                    corrupt_count=$((corrupt_count + 1))
                fi
            done
        done
    fi

    # --- Check genmap map outputs for 0-byte files ---
    if [[ -d "${utils_dir}/genmap_out" ]]; then
        for org_dir in "${utils_dir}/genmap_out/"*/; do
            [[ -d "$org_dir" ]] || continue
            for f in "${org_dir}"*.freq16; do
                [[ -e "$f" ]] || continue
                if [[ ! -s "$f" ]]; then
                    log_warning "Corrupt genmap map output (0 bytes): $(basename "$(dirname "$f")")/$(basename "$f")"
                    corrupt_files+=("$f")
                    corrupt_count=$((corrupt_count + 1))
                fi
            done
        done
    fi

    # --- Check genmap index directories for required core files ---
    # A complete genmap index must have: index.sa.val, index.lf.drv, index.rev.lf.drv
    # These are written last during indexing — missing means interrupted build
    if [[ -d "${utils_dir}/genmap_indices" ]]; then
        local required_files=("index.sa.val" "index.lf.drv" "index.rev.lf.drv")
        for org_dir in "${utils_dir}/genmap_indices/"*/; do
            [[ -d "$org_dir" ]] || continue
            local org_name
            org_name=$(basename "$org_dir")
            for req_file in "${required_files[@]}"; do
                if [[ ! -s "${org_dir}${req_file}" ]]; then
                    log_warning "Incomplete genmap index for ${org_name} (missing/empty ${req_file})"
                    corrupt_dirs+=("$org_dir")
                    corrupt_count=$((corrupt_count + 1))
                    break  # one missing file is enough to flag the directory
                fi
            done
        done
    fi

    # --- Check blastdb marker files for 0-byte ---
    if [[ -d "${utils_dir}/blastdbs" ]]; then
        for f in "${utils_dir}/blastdbs/"*.nsq; do
            [[ -e "$f" ]] || continue
            if [[ ! -s "$f" ]]; then
                local db_base="${f%.nsq}"
                log_warning "Corrupt blastdb (0-byte .nsq): $(basename "$db_base")"
                corrupt_blastdb_bases+=("$db_base")
                corrupt_count=$((corrupt_count + 1))
            fi
        done
    fi

    # --- Check metadata pickles for 0-byte ---
    for meta_dir in "metadata_genomes" "small_meta"; do
        if [[ -d "${utils_dir}/${meta_dir}" ]]; then
            for f in "${utils_dir}/${meta_dir}/"*; do
                [[ -e "$f" ]] || continue
                [[ -d "$f" ]] && continue
                if [[ ! -s "$f" ]]; then
                    log_warning "Corrupt ${meta_dir} file (0 bytes): $(basename "$f")"
                    corrupt_files+=("$f")
                    corrupt_count=$((corrupt_count + 1))
                fi
            done
        fi
    done

    # --- Act on findings ---
    if [[ $corrupt_count -gt 0 ]]; then
        if [[ "$AUTO_FIX_CORRUPT" == "true" ]]; then
            # Remove corrupt files
            for f in "${corrupt_files[@]}"; do
                rm -f "$f"
            done
            for d in "${corrupt_dirs[@]}"; do
                rm -rf "$d"
            done
            for base in "${corrupt_blastdb_bases[@]}"; do
                rm -f "${base}".n*
            done
            log_warning "Removed ${corrupt_count} corrupt/incomplete file(s) from utils/"
            log_info "Snakemake will regenerate them automatically"
        else
            log_error "Found ${corrupt_count} corrupt/incomplete file(s) in utils/"
            log_error "Pipeline cannot proceed with corrupt preprocessing files."
            log_error ""
            log_error "Options:"
            log_error "  1. Re-run with --auto-fix-corrupt to auto-remove and regenerate"
            log_error "  2. Manually inspect and remove the listed files"
            log_error ""
            error_exit "Corrupt utils/ files detected. Use --auto-fix-corrupt to enable automatic cleanup."
        fi
    else
        log_success "All utils/ files passed integrity checks"
    fi
}

# Step 6: Setup all directories
setup_all_directories() {
    log_phase "STEP 6: Directory Setup"

    # Mode 2: Continue run (explicit preservation)
    if [[ "$CONTINUE_RUN" == "true" ]]; then
        log_info "Mode: Continue Run (--continue-run)"
        log_info "Preserving template/ directories for resume"
        log_info "Skipping setup_directories.py - assuming consistent state"
        return 0
    fi

    # Mode 1: Clean run (DEFAULT - always cleans template/)
    log_info "Mode: Clean Run (default)"
    log_info "Cleaning template/ directories"
    python3 "${PROJECT_ROOT}/setup_directories.py" "$WORK_DIR" || error_exit "Directory setup failed"
    log_success "All directories ready"
}

# Step 5c: Skip parameter generation if the filter stage isn't in the requested range.
# Replaces the old skip_params_for_pairwise_targets which hard-coded three TARGET names.
# Filter (= the "anchors" stage in our 5-stage model) is the only consumer of
# genmap/macle/dups params, so we skip generating them whenever anchors is not
# in [FROM_STAGE, TO_STAGE]. The new common.smk:103 gate uses the same logic
# at parse time to skip the FileNotFoundError check.
skip_params_if_filter_not_in_range() {
    local from_idx to_idx anchors_idx
    from_idx=$(stage_idx "$FROM_STAGE")
    to_idx=$(stage_idx "$TO_STAGE")
    anchors_idx=$(stage_idx "anchors")
    if (( from_idx > anchors_idx )) || (( to_idx < anchors_idx )); then
        log_info "Anchors stage not in requested range ($FROM_STAGE..$TO_STAGE)"
        log_info "Skipping parameter file generation (filter rule won't run)"
        SKIP_GENMAP=true
        SKIP_MACLE=true
        SKIP_DUPS=true
    fi
}

# Step 6b: Create touch files for existing anchors NOT in compute_anchors_for
# This handles the incremental scenario: some orgs have pre-existing anchors,
# while others (in compute_anchors_for) are being computed fresh.
# Without these touch files, Snakemake would try to run update_candidates
# for ALL_ORGS, triggering unwanted anchor computation for existing orgs.
setup_existing_anchor_touch_files() {
    log_info "Checking for existing anchors outside compute_anchors_for..."

    local touch_dir="${WORK_DIR}/touch"
    mkdir -p "$touch_dir"

    local orgs_with_existing=0

    # Verify required files exist
    if [[ ! -f "${WORK_DIR}/compute_anchors_for" ]]; then
        error_exit "compute_anchors_for not found at ${WORK_DIR}/compute_anchors_for - run setup_compute_anchors first"
    fi
    if [[ ! -f "${WORK_DIR}/orgs" ]]; then
        error_exit "orgs not found at ${WORK_DIR}/orgs - run validate_and_process_species first"
    fi

    # Check each organism in orgs file
    while IFS= read -r org || [[ -n "$org" ]]; do
        [[ -z "$org" || "$org" =~ ^# ]] && continue

        # Skip if this org is in compute_anchors_for (will create its own touch file)
        # Using grep instead of associative array for bash 3.x compatibility
        if grep -qx "$org" "${WORK_DIR}/compute_anchors_for"; then
            continue
        fi

        local candidates_file="${PROJECT_ROOT}/anchors/candidates/${org}"
        local aligned_file="${PROJECT_ROOT}/anchors/aligned/${org}"

        if [[ -f "$candidates_file" && -f "$aligned_file" ]]; then
            # Org has existing anchors but is NOT in compute_anchors_for
            # Create touch file so Snakemake doesn't try to recompute
            touch "${touch_dir}/update_candidates_done_${org}"
            ((orgs_with_existing++)) || true
            log_info "  Found existing anchors for: $org (touch file created)"
        else
            # Org is not in compute_anchors_for AND has no anchors
            # This is an error for any run that includes the pairwise stage
            # (pairwise needs anchors for every org in ALL_ORGS, not just
            # compute_anchors_for). Stages preprocessing/anchors don't need it.
            local _from_idx _to_idx _pairwise_idx
            _from_idx=$(stage_idx "$FROM_STAGE")
            _to_idx=$(stage_idx "$TO_STAGE")
            _pairwise_idx=$(stage_idx "pairwise")
            if (( _from_idx <= _pairwise_idx )) && (( _to_idx >= _pairwise_idx )); then
                log_error "Organism '$org' is in orgs but not in compute_anchors_for, and has no existing anchors!"
                log_error "  Either add it to compute_anchors_for, or provide existing anchors"
                log_error "  Expected: $candidates_file"
                log_error "  Expected: $aligned_file"
                exit 1
            fi
        fi
    done < "${WORK_DIR}/orgs"

    if [[ $orgs_with_existing -gt 0 ]]; then
        log_success "Created touch files for $orgs_with_existing organisms with existing anchors"
        log_info "These will be used for pairwise comparisons but not recomputed"
    fi
}

# Step 6c (part 1 of 2): READ-ONLY validation of --from-target prerequisites.
# Split out of the old setup_touch_files_for_from_stage so it can run BEFORE the
# destructive directory clean (setup_all_directories). Previously the
# --continue-run guard and the prerequisite check lived together with the touch
# fabrication and ran AFTER the clean -- so a `--from-target validation|downstream`
# WITHOUT --continue-run would wipe template/ intermediates (touch/, parse_bcamm/,
# sequences_to_compare/) and only THEN report "requires --continue-run", destroying
# the very state the guard exists to protect. This function performs NO writes
# (no mkdir, no touch); see fabricate_touch_files_for_from_stage for the writes.
#
# Cases by FROM:
#   preprocessing | anchors  -> nothing to validate (anchors_only rule creates its
#                               own update_candidates_done inputs).
#   pairwise                 -> anchors candidates/aligned must exist per ALL_ORGS.
#   validation | downstream  -> REQUIRES --continue-run because the pairwise stage
#                               produces template/ intermediates (parse_bcamm/,
#                               sequences_to_compare/, pair_done_*) we cannot
#                               fabricate; --continue-run preserves them. Also
#                               requires aligned_succinct/{org} to exist.
validate_from_stage_prerequisites() {
    local from_idx pairwise_idx validation_idx
    from_idx=$(stage_idx "$FROM_STAGE")
    pairwise_idx=$(stage_idx "pairwise")
    validation_idx=$(stage_idx "validation")

    # Nothing to validate for FROM=preprocessing or FROM=anchors
    if (( from_idx < pairwise_idx )); then
        return 0
    fi

    # FROM >= validation requires --continue-run so template/touch and template's
    # per-pair intermediate state from a prior pairwise run are preserved.
    # Without --continue-run, setup_directories.py would clean template/ and
    # pair_done_* / parse_bcamm/ / sequences_to_compare/ would be gone -> snakemake
    # would cascade back into the pairwise chain. This check MUST run before the
    # clean (it does now: called pre-setup_all_directories in main()).
    if (( from_idx >= validation_idx )) && [[ "$CONTINUE_RUN" != "true" ]]; then
        log_error "--from-target $FROM_STAGE requires --continue-run"
        log_error "  Skipping the pairwise stage means snakemake must see prior pair_done_*"
        log_error "  touch files and other template/ intermediates from an earlier complete"
        log_error "  pairwise run. Without --continue-run, template/ is cleaned and that"
        log_error "  state is lost -> snakemake would cascade back into pairwise."
        log_error "  Re-run with: --from-target $FROM_STAGE --continue-run"
        exit 1
    fi

    log_phase "STEP 6c (validate): Checking --from-target prerequisites for FROM=$FROM_STAGE"

    if [[ ! -f "${WORK_DIR}/orgs" ]]; then
        error_exit "orgs not found at ${WORK_DIR}/orgs"
    fi

    # aligned_succinct/{org} is only a prerequisite for FROM >= validation.
    local need_aligned_succinct_file=false
    if (( from_idx >= validation_idx )); then
        need_aligned_succinct_file=true
    fi

    # Per-org existence check on anchor outputs. These live OUTSIDE the work dir
    # (in ${PROJECT_ROOT}/anchors), so reading them before the clean is safe.
    local orgs_ok=0 orgs_missing=0
    while IFS= read -r org || [[ -n "$org" ]]; do
        [[ -z "$org" || "$org" =~ ^# ]] && continue

        local candidates_file="${PROJECT_ROOT}/anchors/candidates/${org}"
        local aligned_file="${PROJECT_ROOT}/anchors/aligned/${org}"
        local succinct_file="${PROJECT_ROOT}/anchors/aligned_succinct/${org}"
        local missing=""

        if [[ ! -f "$candidates_file" ]]; then missing="$missing  $candidates_file"$'\n'; fi
        if [[ ! -f "$aligned_file" ]];    then missing="$missing  $aligned_file"$'\n';    fi
        if [[ "$need_aligned_succinct_file" == "true" && ! -f "$succinct_file" ]]; then
            missing="$missing  $succinct_file"$'\n'
        fi

        if [[ -n "$missing" ]]; then
            log_warning "Missing prerequisites for $org:"
            printf "%b" "$missing" | sed 's/^/   /' >&2
            ((orgs_missing++)) || true
            continue
        fi
        ((orgs_ok++)) || true
    done < "${WORK_DIR}/orgs"

    if (( orgs_missing > 0 )); then
        log_error "Cannot run with --from-target $FROM_STAGE: $orgs_missing organism(s) missing prerequisite files"
        log_error "Run earlier stages first, or provide the missing files"
        exit 1
    fi

    log_success "Prerequisite state verified for $orgs_ok organisms (FROM=$FROM_STAGE)"
    if (( from_idx >= validation_idx )); then
        log_info "  (assumes prior pairwise run left pair_done_* / parse_bcamm/ / sequences_to_compare/ intact in template/)"
    fi
}

# Step 6c (part 2 of 2): WRITE the fabricated "stage X done" touch files so
# snakemake's DAG sees stages before FROM_STAGE as completed and short-circuits
# them. Runs AFTER the directory clean, because touch/ lives in the work dir and
# is wiped/recreated by the clean -- fabricating earlier would just be deleted.
# Prerequisites are already verified by validate_from_stage_prerequisites (called
# pre-clean), so this assumes the anchor files exist; it fabricates touch files
# with mtime matching anchors/aligned/{org} (via touch -r) to avoid a spurious
# mtime-cascade, and never clobbers an existing touch file.
fabricate_touch_files_for_from_stage() {
    local from_idx pairwise_idx validation_idx downstream_idx
    from_idx=$(stage_idx "$FROM_STAGE")
    pairwise_idx=$(stage_idx "pairwise")
    validation_idx=$(stage_idx "validation")
    downstream_idx=$(stage_idx "downstream")

    # Nothing to do for FROM=preprocessing or FROM=anchors
    if (( from_idx < pairwise_idx )); then
        return 0
    fi

    log_phase "STEP 6c (fabricate): Writing touch files for stages before FROM=$FROM_STAGE"

    local touch_dir="${WORK_DIR}/touch"
    mkdir -p "$touch_dir"

    if [[ ! -f "${WORK_DIR}/orgs" ]]; then
        error_exit "orgs not found at ${WORK_DIR}/orgs"
    fi

    while IFS= read -r org || [[ -n "$org" ]]; do
        [[ -z "$org" || "$org" =~ ^# ]] && continue

        local aligned_file="${PROJECT_ROOT}/anchors/aligned/${org}"
        # Defensive: prerequisites were validated pre-clean; skip silently if the
        # anchor file is somehow absent (validation owns the error reporting).
        [[ -f "$aligned_file" ]] || continue

        # Build the list of touch files this FROM needs (per-org)
        local touch_basenames=()
        touch_basenames+=("update_candidates_done_${org}")
        if (( from_idx >= validation_idx )); then
            touch_basenames+=("update_aligned_stage2_done_${org}")
            touch_basenames+=("collect_output_done_${org}")
            touch_basenames+=("inconsistencies_done_${org}")
            touch_basenames+=("eval_syn_done_${org}")
        fi
        if (( from_idx >= downstream_idx )); then
            touch_basenames+=("checks_done_${org}")
            touch_basenames+=("statistics_done_${org}")
            touch_basenames+=("checks_succinct_done_${org}")
            touch_basenames+=("statistics_succinct_done_${org}")
        fi

        # Only fabricate touch files that DON'T already exist (don't clobber
        # existing valid timestamps from a prior --continue-run state).
        # Use touch -r aligned_file so the fabricated touch's mtime matches the
        # anchor file -> snakemake's mtime check sees outputs >= inputs -> no
        # spurious rerun cascade.
        local n_fabricated=0
        local tf
        for tf in "${touch_basenames[@]}"; do
            if [[ ! -f "$touch_dir/$tf" ]]; then
                touch -r "$aligned_file" "$touch_dir/$tf"
                ((n_fabricated++)) || true
            fi
        done
        if (( n_fabricated > 0 )); then
            log_info "  $org: fabricated $n_fabricated missing touch file(s) (mtime via touch -r aligned)"
        fi
    done < "${WORK_DIR}/orgs"
}

# Step 7: Run Snakemake pipeline
run_pipeline() {
    log_info "Running Snakemake pipeline..."
    log_info "This may take 1-2 hours depending on genome count and size..."

    # Verify installation before running pipeline
    log_info "Verifying installation dependencies..."
    if ! python "${PROJECT_ROOT}/verify_installation.py"; then
        log_error "Installation verification failed! Please check dependencies."
        log_error "Run './install.sh' to install missing dependencies."
        exit 1
    fi
    log_success "Installation verified successfully"

    # Validate utils/ files before running pipeline
    validate_utils_integrity

    cd "$CODE_DIR"

    # Build snakemake command
    # Note: Using default rerun triggers (mtime, params, code, software_env) for better change detection
    # --keep-incomplete: Preserve failed outputs for debugging
    # --allow-ambiguity: Required because Snakemake 9's ruleorder doesn't properly handle
    #                    multiple rules with update() on the same output file
    # -p: Print shell commands for debugging
    # NOTE: We do NOT use --scheduler greedy, we want the ILP solver (via pulp) for optimal scheduling

    # Build config string
    CONFIG_STR="work_dir=$WORK_DIR code_dir=$CODE_DIR root_dir=$PROJECT_ROOT"

    # Add pairwise file if provided
    if [[ -n "$PAIRWISE_FILE" ]]; then
        CONFIG_STR="$CONFIG_STR pairwise_file=${WORK_DIR}/pairwise_comparisons.txt"
    fi

    # Add memory mode flags if specified
    if [ "$LOW_MEM" = true ]; then
        CONFIG_STR="$CONFIG_STR low_mem=true"
        log_info "Using low memory specifications (reduced multipliers)"
    elif [ "$HIGH_MEM" = true ]; then
        CONFIG_STR="$CONFIG_STR high_mem=true"
        log_info "Using high memory specifications (2x multipliers)"
    fi

    if [ "$SLURM_MODE" = true ]; then
        # SLURM mode: use profile instead of --cores
        if [ -n "$SLURM_CONFIG" ]; then
            # User provided a custom config path - auto-detect file vs directory
            if [ -f "$SLURM_CONFIG" ]; then
                # It's a file - use parent directory as profile (--profile expects directory)
                PROFILE_PATH="$(dirname "$SLURM_CONFIG")"
                log_warning "Note: --slurm expects directory, not file. Using parent: $PROFILE_PATH"
            elif [ -d "$SLURM_CONFIG" ]; then
                # It's a directory - use as-is
                PROFILE_PATH="$SLURM_CONFIG"
            else
                error_exit "--slurm path does not exist: $SLURM_CONFIG"
            fi
        else
            # --slurm requires a profile path
            error_exit "--slurm requires a profile path. Usage: --slurm /path/to/profile"
        fi
        log_info "Using SLURM profile: $PROFILE_PATH"
        SNAKEMAKE_CMD="snakemake --snakefile Snakefile --profile $PROFILE_PATH --keep-incomplete --latency-wait 60 --allow-ambiguity -p"
    else
        # Local mode: use --cores
        SNAKEMAKE_CMD="snakemake --snakefile Snakefile --cores $CORES --keep-incomplete --latency-wait 60 --allow-ambiguity -p"

    fi

    # Pick the snakemake top-level rule from the TO stage (Snakemake 8+ requires
    # targets before options). For (FROM,TO) == (preprocessing,downstream) i.e.
    # the full pipeline, that's "all" which is the default first rule.
    local _TARGET_RULE
    _TARGET_RULE=$(to_rule_for_stage "$TO_STAGE") \
        || error_exit "internal error: no rule mapped for TO_STAGE=$TO_STAGE"
    log_info "Stage range: $FROM_STAGE..$TO_STAGE (snakemake target rule: $_TARGET_RULE)"
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD $_TARGET_RULE"

    # Add config AFTER target -- includes from_stage/to_stage so common.smk
    # can gate the param-file requirement check on whether anchors is in range.
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --config $CONFIG_STR from_stage=$FROM_STAGE to_stage=$TO_STAGE"

    # Add memory limit AFTER config (must not precede positional target args,
    # because --resources uses nargs='+' and would consume the target name)
    if [[ -n "$MAX_MEM_MB" ]]; then
        log_info "Setting maximum total memory: ${MAX_MEM_MB} MB"
        SNAKEMAKE_CMD="$SNAKEMAKE_CMD --resources mem_mb=$MAX_MEM_MB"
    fi

    # --- --continue-run auto-behaviours (see plan: stage-range CLI + auto-behaviours)
    # 1) inject --rerun-triggers mtime BEFORE the user's --snakemake-flags so a
    #    user-passed --rerun-triggers can still override ours (snakemake CLI:
    #    later occurrences win).
    # We do NOT try to handle IncompleteFilesException ourselves -- that would
    # require parsing snakemake's stderr format and poking at .snakemake/
    # internals (path schemes, backup layout) that are not API. Instead we
    # surface snakemake's own error verbatim and point at utils/check_pipeline_state.py
    # for a thorough integrity check the user can run independently.
    if [[ "$CONTINUE_RUN" == "true" ]]; then
        SNAKEMAKE_CMD="$SNAKEMAKE_CMD --rerun-triggers mtime"
        log_info "Continue-run: adding --rerun-triggers mtime (overridable via --snakemake-flags)"
        log_info "  If snakemake aborts with IncompleteFilesException, run:"
        log_info "    python ${PROJECT_ROOT}/utils/check_pipeline_state.py ${PROJECT_ROOT} --work-dir $WORK_DIR --cores $CORES"
        log_info "  for a structural integrity report on pipeline outputs."

        # Code-edit-during-resume detection: --rerun-triggers mtime makes snakemake
        # ignore code-hash changes between runs. If the user edited any code/*.py
        # AFTER an output was last produced, that output won't be regenerated even
        # though the rule's logic changed. We warn (not block) here so the user
        # can decide. Walks done in Python for mac/linux portability.
        _CODE_EDIT_WARN=$(PROJECT_ROOT="$PROJECT_ROOT" python3 - <<'PY'
import os, sys
root = os.environ['PROJECT_ROOT']
code_dir = os.path.join(root, 'code')
out_dir  = os.path.join(root, 'anchors', 'aligned')
if not (os.path.isdir(code_dir) and os.path.isdir(out_dir)):
    sys.exit(0)
newest_code, newest_code_file = -1.0, None
for dp, _, files in os.walk(code_dir):
    if os.sep + '.' in dp:        # skip dot dirs (e.g. .snakemake/)
        continue
    for f in files:
        if not f.endswith('.py'):
            continue
        p = os.path.join(dp, f)
        try:
            m = os.path.getmtime(p)
        except OSError:
            continue
        if m > newest_code:
            newest_code, newest_code_file = m, p
if newest_code_file is None:
    sys.exit(0)
oldest_out, oldest_out_file = float('inf'), None
for f in os.listdir(out_dir):
    p = os.path.join(out_dir, f)
    if not os.path.isfile(p):
        continue
    try:
        m = os.path.getmtime(p)
    except OSError:
        continue
    if m < oldest_out:
        oldest_out, oldest_out_file = m, p
if oldest_out_file is None:
    sys.exit(0)
if newest_code > oldest_out:
    # tab-separated: code_file, output_file, delta_seconds
    print(f"{newest_code_file}\t{oldest_out_file}\t{int(newest_code - oldest_out)}")
PY
)
        if [[ -n "$_CODE_EDIT_WARN" ]]; then
            local _cf _of _dt
            IFS=$'\t' read -r _cf _of _dt <<<"$_CODE_EDIT_WARN"
            log_warning "Code edit detected AFTER an existing pipeline output:"
            log_warning "  newest code:    $_cf"
            log_warning "  oldest output:  $_of"
            log_warning "  delta:          ${_dt}s (code is newer)"
            log_warning "Under --rerun-triggers mtime, snakemake will NOT rerun rules whose"
            log_warning "code body changed. If your edit changes pipeline outputs you should:"
            log_warning "  a) override the trigger set:  --snakemake-flags '--rerun-triggers mtime code'"
            log_warning "  b) delete the specific output(s) you want regenerated"
            log_warning "  c) only ignore if the edit was a non-semantic refactor"
        fi
    fi

    # Add user-provided Snakemake flags if specified (LAST so they can override ours)
    if [ -n "$SNAKEMAKE_EXTRA_FLAGS" ]; then
        log_info "Additional Snakemake flags: $SNAKEMAKE_EXTRA_FLAGS"
        SNAKEMAKE_CMD="$SNAKEMAKE_CMD $SNAKEMAKE_EXTRA_FLAGS"
    fi

    # Run snakemake; let its stdout/stderr go to the user verbatim so they see
    # IncompleteFilesException, MissingInputException, etc. as snakemake reports them.
    local _SM_RC=0
    eval "$SNAKEMAKE_CMD" || _SM_RC=$?

    if (( _SM_RC != 0 )); then
        log_error "Snakemake exited with status $_SM_RC (output above)."
        if [[ "$CONTINUE_RUN" == "true" ]]; then
            log_error "  If the cause is incomplete files from a prior interrupted run:"
            log_error "    1) inspect integrity:  python ${PROJECT_ROOT}/utils/check_pipeline_state.py ${PROJECT_ROOT} --work-dir $WORK_DIR --cores $CORES"
            log_error "    2) clear markers for files you trust:  snakemake --cleanup-metadata <FILE1> <FILE2> ..."
            log_error "    3) delete files you want regenerated, then re-run with --continue-run"
        fi
        error_exit "Snakemake pipeline failed"
    fi

    cd "$PROJECT_ROOT"

    log_success "Pipeline completed successfully!"
}

# Step 8: Report results
report_results() {
    log_info "=========================================="
    log_info "Pipeline Results"
    log_info "=========================================="

    # Count anchors. NOTE: anchors live at the project root (anchors/), NOT under
    # the work dir -- scripts use root + '/anchors' (root-fixed, shared across runs).
    if [[ -d "${PROJECT_ROOT}/anchors" ]]; then
        local anchor_count=$(find "${PROJECT_ROOT}/anchors" -name "*.fasta" 2>/dev/null | wc -l)
        log_success "Anchor files generated: $anchor_count"
    fi

    # --- Validation Summary ---
    # Aggregate check results per organism
    local checks_dir="${WORK_DIR}/logs/checks"
    local stats_dir="${WORK_DIR}/logs/statistics"
    local stats_succinct_dir="${WORK_DIR}/logs/statistics_succinct"

    if [[ -d "$checks_dir" ]]; then
        log_info ""
        log_info "=========================================="
        log_info "Validation Summary"
        log_info "=========================================="

        local passed=0
        local passed_with_warnings=0
        local failed=0
        local failed_orgs=""
        local warned_orgs=""

        # Read orgs list
        while IFS= read -r org || [[ -n "$org" ]]; do
            [[ -z "$org" || "$org" =~ ^# ]] && continue

            local err_file="${checks_dir}/${org}.err"
            if [[ ! -f "$err_file" ]]; then
                continue  # Check didn't run for this org
            fi

            local has_errors=false
            local has_warnings=false

            # Check for ERRORs (pipeline failures)
            if grep -q "^ERROR:" "$err_file" 2>/dev/null; then
                has_errors=true
            fi

            # Check for WARNINGs (non-critical)
            if grep -q "^WARNING:" "$err_file" 2>/dev/null; then
                has_warnings=true
            fi

            if [[ "$has_errors" == "true" ]]; then
                ((failed++)) || true
                if [[ -n "$failed_orgs" ]]; then
                    failed_orgs="${failed_orgs}, ${org}"
                else
                    failed_orgs="${org}"
                fi
            elif [[ "$has_warnings" == "true" ]]; then
                ((passed_with_warnings++)) || true
                if [[ -n "$warned_orgs" ]]; then
                    warned_orgs="${warned_orgs}, ${org}"
                else
                    warned_orgs="${org}"
                fi
            else
                ((passed++)) || true
            fi
        done < "${WORK_DIR}/orgs"

        local total=$((passed + passed_with_warnings + failed))

        if [[ $failed -eq 0 && $passed_with_warnings -eq 0 ]]; then
            log_success "All $total species passed all checks"
        else
            if [[ $passed -gt 0 ]]; then
                log_success "$passed species passed all checks"
            fi
            if [[ $passed_with_warnings -gt 0 ]]; then
                log_warning "$passed_with_warnings species passed with warnings: $warned_orgs"
                log_info "  (Non-reciprocal syntenic dups are warnings, not errors."
                log_info "   This is rare and biologically justified. See logs for details.)"
            fi
            if [[ $failed -gt 0 ]]; then
                log_error "$failed species FAILED validation: $failed_orgs"
            fi
        fi

        log_info ""
        log_info "Check logs:      ${checks_dir}/"
        if [[ -d "$stats_dir" ]]; then
            log_info "Statistics logs:  ${stats_dir}/"
        fi
        if [[ -d "$stats_succinct_dir" ]]; then
            log_info "Succinct stats:   ${stats_succinct_dir}/"
        fi
    fi

    # List output directories
    log_info ""
    log_info "Output locations:"
    log_info "  Anchor candidates: ${PROJECT_ROOT}/anchors/candidates/"
    log_info "  Aligned maps:      ${PROJECT_ROOT}/anchors/aligned/"
    log_info "  Succinct maps:     ${PROJECT_ROOT}/anchors/aligned_succinct/"
    log_info "  BLAST results:     ${WORK_DIR}/clasp_out_forward/, clasp_out_reverse/"
    log_info "  Parse results:     ${WORK_DIR}/parse_bcamm/"

    # Report downstream outputs if they exist
    if [[ -f "${PROJECT_ROOT}/downstream/full.tar.gz" ]]; then
        log_info ""
        log_info "Downstream outputs:"
        log_info "  Full archive: ${PROJECT_ROOT}/downstream/full.tar.gz"
        log_info "  Main archive: ${PROJECT_ROOT}/downstream/main.tar.gz"
        log_info "  Packaged output: ${PROJECT_ROOT}/downstream/out/"
        log_info "  Phylogenetic trees: ${PROJECT_ROOT}/downstream/*.nwk"
        log_info "  GFF files: ${PROJECT_ROOT}/downstream/gff/"
        log_info "  MCScanX files: ${PROJECT_ROOT}/downstream/MCScanX.*"
    fi
}

# Main execution
validate_work_dir() {
    log_info "Validating work dir: $WORK_DIR"
    [[ -d "$WORK_DIR" ]] || error_exit "work dir not found: $WORK_DIR"
    local report rc=0
    report=$(WORK_DIR="$WORK_DIR" python3 - <<'PY'
import os, re, json, sys
wd = os.environ['WORK_DIR']
problems = []
cfg = os.path.join(wd, 'pipeline_config.yaml')
if not os.path.isfile(cfg):
    problems.append("missing file: pipeline_config.yaml")
else:
    text = open(cfg, encoding='utf-8', errors='replace').read()
    for key in ('cores', 'max_mem_mb', 'timeout_self_minutes', 'timeout_bcamm_minutes'):
        if not re.search(r'^  ' + key + r':', text, re.M):
            problems.append("pipeline_config.yaml: missing '  " + key + ":' line (2-space indent)")
dc = os.path.join(wd, 'directories_config.json')
if not os.path.isfile(dc):
    problems.append("missing file: directories_config.json")
else:
    try:
        d = json.load(open(dc, encoding='utf-8'))
        if not isinstance(d, dict) or not isinstance(d.get('directories'), list) or not d['directories']:
            problems.append("directories_config.json: needs a non-empty 'directories' list")
    except Exception as e:
        problems.append("directories_config.json: invalid JSON (" + str(e) + ")")
for p in problems:
    print(p)
sys.exit(1 if problems else 0)
PY
) || rc=$?
    if [[ $rc -ne 0 ]]; then
        log_error "Work dir not usable: $WORK_DIR"
        while IFS= read -r l; do [[ -n "$l" ]] && log_error "  - $l"; done <<< "$report"
        error_exit "every required file must be present and correctly formatted in the work dir"
    fi
    log_success "Work dir OK"
}

main() {
    validate_work_dir             # work dir must have the required files (fail loud)
    validate_and_process_species  # Validates species file, downloads NCBI genomes, creates orgs file
    validate_genomes              # Normalizes extensions (.fna/.fa/.faa -> .fasta) and validates FASTA format
    validate_from_stage_prerequisites  # READ-ONLY: --from-target guard + anchor prereq check. MUST run
                                       # before the clean so a missing --continue-run errors BEFORE
                                       # template/ intermediates are wiped (see function header).
    setup_all_directories         # Destructive clean of template/ (skipped under --continue-run).
                                  # Kept after all read-only checks that can fail.
    skip_params_if_filter_not_in_range  # Skip param generation when anchors stage not in [FROM,TO]
    setup_compute_anchors         # Now safe - directory cleanup already done
    setup_parameters
    setup_pairwise_comparisons    # Setup custom pairwise comparisons if provided
    setup_existing_anchor_touch_files  # Create touch files for pre-existing anchors (incremental mode)
    fabricate_touch_files_for_from_stage  # WRITE touch files for stages before FROM_STAGE (post-clean:
                                          # touch/ is wiped by the clean, so fabricate after it)
    run_pipeline
    report_results

    log_success "=========================================="
    log_success "AncST Pipeline Complete!"
    log_success "=========================================="
}

# Run main function
main
