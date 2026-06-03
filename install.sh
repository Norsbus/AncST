#!/usr/bin/env bash

############################################################################
# AncST Installation Script - Cross-Platform Edition
#
# This script:
# 1. Detects OS (Linux/macOS) and architecture (x86_64/ARM64)
# 2. Creates conda environment from environment.yml (all dependencies)
# 3. Compiles macle and clasp.x from source with proper arch flags
#
# Usage: ./install.sh [--env-name NAME] [--force-recreate]
############################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get absolute path to project root
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"

# Default environment name
ENV_NAME="AncST"
FORCE_RECREATE=false
CLEAN_CACHE=false
CONDA_PATH=""
FORCE_OLD_CONDA=false

# Build log file
BUILD_LOG="$PROJECT_ROOT/tmp/install_build.log"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --force-recreate)
            FORCE_RECREATE=true
            shift
            ;;
        --clean-cache)
            CLEAN_CACHE=true
            shift
            ;;
        --conda-path)
            CONDA_PATH="$2"
            shift 2
            ;;
        --force-old-conda)
            FORCE_OLD_CONDA=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --env-name NAME       Name for conda environment (default: AncST)"
            echo "  --force-recreate      Remove and recreate environment from scratch"
            echo "                        (also cleans conda package cache)"
            echo "  --clean-cache         Clean conda package cache before installation"
            echo "                        Use if you see 'corrupted cache' or '404' errors"
            echo "  --conda-path PATH     Explicitly specify conda executable path"
            echo "                        (use if auto-detection fails or picks wrong version)"
            echo "  --force-old-conda     Allow using conda < 23.1.0 (bypass version check)"
            echo "                        WARNING: May fail during environment creation"
            echo "  --help, -h            Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Logging functions
log_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }
error_exit() {
    log_error "$1"
    if [ -f "$BUILD_LOG" ]; then
        echo ""
        log_info "Build log saved at: $BUILD_LOG"
        log_info "Check the log for detailed error information"
    fi
    exit 1
}

############################################################################
# STEP 0: PRE-FLIGHT CHECKS AND PLATFORM DETECTION
############################################################################

log_info "AncST Installation Starting..."
log_info "Project root: $PROJECT_ROOT"
echo ""

# Check required build tools in CURRENT environment (before PATH modifications)
log_info "Checking required build tools..."

# Detect OS and architecture
log_info "Detecting platform..."
OS_TYPE=$(uname -s)
case "$OS_TYPE" in
    Linux*)     OS_NAME="Linux";;
    Darwin*)    OS_NAME="macOS";;
    *)          error_exit "Unsupported OS: $OS_TYPE (only Linux and macOS are supported)";;
esac

ARCH_TYPE=$(uname -m)
case "$ARCH_TYPE" in
    x86_64)     ARCH_NAME="x86_64";;
    aarch64)    ARCH_NAME="aarch64";;  # Linux ARM64
    arm64)      ARCH_NAME="arm64";;    # macOS ARM64
    i686|i386|armv7*)
        error_exit "Unsupported architecture: $ARCH_TYPE

  AncST requires a 64-bit system (x86_64 or ARM64/aarch64)
  Your system appears to be 32-bit, which is not supported."
        ;;
    *)
        error_exit "Unknown architecture: $ARCH_TYPE

  Please report this issue at: https://github.com/Norsbus/AncST/issues"
        ;;
esac

log_success "Platform detected: $OS_NAME on $ARCH_NAME"

# Check for gcc FIRST (needs special verbose error)
if [ "$OS_NAME" != "macOS" ] && ! command -v gcc &>/dev/null; then
    error_exit "gcc (C compiler) not found!

  gcc is required to compile macle and clasp from source.

  Installation options:

  1. System installation (RECOMMENDED for HPC/clusters):
     - On RHEL/CentOS/Rocky Linux: sudo yum groupinstall 'Development Tools'
     - On Ubuntu/Debian: sudo apt-get install build-essential
     - On macOS: xcode-select --install

  2. HPC module system (if available on your cluster):
     - module load gcc
     - module load GCC/11.3.0  (or similar - check: module avail gcc)
     - Add to ~/.bashrc: module load gcc

  3. Conda installation (NOT recommended on HPC - may cause conflicts):
     - conda install -n $ENV_NAME gcc_linux-64 gxx_linux-64

  After installing or loading gcc, re-run this script."
fi

log_success "Found gcc: $(gcc --version 2>/dev/null | head -1 || echo gcc)"

# Check other critical build tools
MISSING_CRITICAL=()
for cmd in git make curl; do
    if ! command -v $cmd &>/dev/null; then
        MISSING_CRITICAL+=("$cmd")
    fi
done

if [ ${#MISSING_CRITICAL[@]} -gt 0 ]; then
    error_exit "Missing required tools: ${MISSING_CRITICAL[*]}

  Please install them first:
    - On Ubuntu/Debian: sudo apt-get install ${MISSING_CRITICAL[*]}
    - On RHEL/CentOS: sudo yum install ${MISSING_CRITICAL[*]}
    - On macOS: brew install ${MISSING_CRITICAL[*]}"
fi

log_success "All required build tools found (git, make, curl)"

# Check for cmake in system (BEFORE modifying PATH or creating conda env)
SYSTEM_CMAKE=""
if command -v cmake &>/dev/null; then
    SYSTEM_CMAKE=$(which cmake)
    log_info "Found system cmake: $SYSTEM_CMAKE"
else
    log_info "No system cmake found (will install to conda environment)"
fi

# Detect compiler strategy (conda vs system)
USE_SYSTEM_COMPILERS=false
CONDA_ARCH=""

if [ "$OS_NAME" = "macOS" ]; then
    # On macOS: Check for conda/system architecture mismatch
    if command -v conda &> /dev/null; then
        CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
        if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/bin/python" ]; then
            CONDA_ARCH=$(file "$CONDA_BASE/bin/python" | grep -o "x86_64\|arm64" | head -1)
            if [ -z "$CONDA_ARCH" ]; then
                log_warning "Could not determine conda Python architecture; assuming native"
                CONDA_ARCH="$SYSTEM_ARCH"
            fi

            SYSTEM_ARCH=$ARCH_NAME

            if [ "$CONDA_ARCH" != "$SYSTEM_ARCH" ]; then
                log_warning "Architecture Mismatch Detected!"
                log_warning "    System: $SYSTEM_ARCH"
                log_warning "    Conda:  $CONDA_ARCH (running via Rosetta)"
                log_warning ""
                log_warning "    SOLUTION: Using system compilers (Xcode) for native $SYSTEM_ARCH builds"
                log_warning "    This avoids conflicts with conda's $CONDA_ARCH compiler wrappers"
                log_warning ""
                USE_SYSTEM_COMPILERS=true
            else
                log_info "macOS: Native conda detected ($CONDA_ARCH), can use conda compilers"
            fi
        fi
    else
        log_info "macOS: No conda found, will use system compilers"
        USE_SYSTEM_COMPILERS=true
    fi
else
    # Linux: No Rosetta issues, conda compilers work fine
    log_info "Linux: Will use available compilers (conda or system)"
fi

# Determine macOS build architecture (critical for CMake + static libs)
# This must be set AFTER USE_SYSTEM_COMPILERS detection above
if [ "$OS_NAME" = "macOS" ]; then
    if [ "$USE_SYSTEM_COMPILERS" = true ]; then
        # Using system compilers: build for native system architecture
        MACOS_BUILD_ARCH="$ARCH_NAME"
        log_info "macOS: System compilers -> building for native architecture: $MACOS_BUILD_ARCH"
    else
        # Using conda compilers: build for conda's architecture
        # (This handles case where conda and system match)
        MACOS_BUILD_ARCH="${CONDA_ARCH:-$ARCH_NAME}"
        log_info "macOS: Conda compilers -> building for conda architecture: $MACOS_BUILD_ARCH"
    fi
fi

# Normalize architecture names for display
if [ "$ARCH_NAME" = "arm64" ]; then
    ARCH_DISPLAY="ARM64 (Apple Silicon)"
elif [ "$ARCH_NAME" = "aarch64" ]; then
    ARCH_DISPLAY="ARM64"
else
    ARCH_DISPLAY="$ARCH_NAME"
fi

# Set compilation flags based on platform

NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
[ -z "$NPROC" ] && NPROC=4

export MAKEFLAGS="-j$NPROC"
if [ "$OS_NAME" = "macOS" ]; then
    # Use MACOS_BUILD_ARCH (determined above based on compiler strategy)
    BUILD_ARCH="$MACOS_BUILD_ARCH"

    if [[ "$BUILD_ARCH" == "x86_64" ]]; then
        # Building for x86_64: Use -march=core2 optimizations if supported
        # Note: This only happens on Intel Mac OR when building x86_64 with conda compilers
        if command -v ${CXX:-clang++} &>/dev/null; then
            TEST_COMPILER="${CXX:-clang++}"
        else
            TEST_COMPILER="clang++"
        fi

        if $TEST_COMPILER -march=core2 -x c++ -c /dev/null -o /dev/null 2>/dev/null; then
            CFLAGS_NATIVE="-arch x86_64 -O2 -pipe -march=core2 -mtune=haswell -mssse3"
            CXXFLAGS_NATIVE="-arch x86_64 -O2 -pipe -march=core2 -mtune=haswell -mssse3"
            LDFLAGS_NATIVE="-arch x86_64"
            log_info "macOS x86_64: Using optimized flags (-march=core2)"
        else
            # Fallback for x86_64 without core2 support
            CFLAGS_NATIVE="-arch x86_64 -O2 -pipe"
            CXXFLAGS_NATIVE="-arch x86_64 -O2 -pipe"
            LDFLAGS_NATIVE="-arch x86_64"
            log_info "macOS x86_64: Using generic flags"
        fi
    else
        # Building for arm64: NEVER use -march= flags (not supported on ARM)
        CFLAGS_NATIVE="-arch $BUILD_ARCH -O2 -pipe"
        CXXFLAGS_NATIVE="-arch $BUILD_ARCH -O2 -pipe"
        LDFLAGS_NATIVE="-arch $BUILD_ARCH"
        log_info "macOS arm64: Using native ARM flags"
    fi
else
    # Linux: native compilation with optimizations
    # WARNING: -march=native compiles for THIS machine's CPU. On SLURM clusters,
    # login nodes often have newer CPUs than compute nodes. If you compile on a
    # login node, binaries may crash with SIGILL (illegal instruction) on compute
    # nodes. Always run install.sh on a compute node:
    #   srun --mem=8G --cpus-per-task=4 ./install.sh
    if command -v sinfo &>/dev/null && [ -z "${SLURM_JOB_ID:-}" ]; then
        log_warning "============================================================"
        log_warning "SLURM cluster detected but running on LOGIN NODE."
        log_warning "-march=native will compile for this node's CPU."
        log_warning "If compute nodes have older CPUs, binaries will crash (SIGILL)."
        log_warning "Recommended: srun --mem=8G --cpus-per-task=4 ./install.sh"
        log_warning "============================================================"
        log_warning "Continuing in 10 seconds (Ctrl+C to cancel)..."
        sleep 10
    fi
    CFLAGS_NATIVE="-march=native -O3"
    CXXFLAGS_NATIVE="-march=native -O3"
    LDFLAGS_NATIVE=""
fi

# Fallback flags (if -march=native not supported)
CFLAGS_FALLBACK="-O3"
CXXFLAGS_FALLBACK="-O3"
LDFLAGS_FALLBACK=""

log_info "Compilation will use $NPROC parallel jobs"

############################################################################
# CLEANUP: Remove old build artifacts (clean on start)
############################################################################

log_info "Cleaning previous build artifacts..."
rm -rf "$PROJECT_ROOT/tmp/macle_build"
rm -rf "$PROJECT_ROOT/tmp/clasp_build"
rm -f "$BUILD_LOG"
mkdir -p "$PROJECT_ROOT/tmp"
log_success "Build environment cleaned"
echo ""

############################################################################
# STEP 0.5: FIND AND INITIALIZE BEST CONDA INSTALLATION
############################################################################

find_and_initialize_best_conda() {
    local min_version="23.1.0"
    local best_conda=""
    local best_version=""
    local version_check_failed=false

    log_info "Searching for conda installations..."

    # STRATEGY 0: Use explicitly specified conda path (highest priority)
    if [ -n "$CONDA_PATH" ]; then
        log_info "Using explicitly specified conda: $CONDA_PATH"

        if [ ! -x "$CONDA_PATH" ]; then
            error_exit "Specified conda path is not executable: $CONDA_PATH

  Please check the path and try again."
        fi

        local version=$("$CONDA_PATH" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "")

        if [ -z "$version" ]; then
            error_exit "Could not determine version of specified conda: $CONDA_PATH

  Please ensure it's a valid conda executable."
        fi

        log_success "Using conda $version (user-specified)"
        best_conda="$CONDA_PATH"
        best_version="$version"
        version_check_failed=false

        # Skip to initialization
    else
        # STRATEGY 1: Check HPC module system first
        if command -v module &>/dev/null; then
            log_info "Detected HPC module system - checking for conda modules..."

            # Get available conda modules
            local conda_modules=$(module avail conda 2>&1 | grep -oE 'conda/[0-9.]+' || true)

            if [ -n "$conda_modules" ]; then
                log_info "Found conda modules:"
                echo "$conda_modules" | while IFS= read -r mod; do
                    log_info "  - $mod"
                done

                # Find highest version module
                local highest_module=$(echo "$conda_modules" | sort -Vr | head -1)
                log_info "Attempting to load: $highest_module"

                if module load "$highest_module" &>/dev/null; then
                    log_success "Loaded HPC module: $highest_module"

                    # Check if conda is now available
                    if command -v conda &>/dev/null; then
                        best_conda=$(command -v conda)
                        best_version=$(conda --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "")
                        if [ -n "$best_version" ]; then
                            log_success "Module provided conda $best_version at $best_conda"
                            # Skip to version check
                        fi
                    fi
                else
                    log_warning "Failed to load $highest_module - falling back to PATH search"
                fi
            else
                log_info "No conda modules found in HPC system"
            fi
        fi

        # STRATEGY 2: Check if conda is already initialized as shell function
        if [ -z "$best_conda" ]; then
            if type conda &>/dev/null 2>&1; then
                local conda_type=$(type -t conda 2>/dev/null || echo "")
                if [ "$conda_type" = "function" ]; then
                    log_info "Detected conda as shell function (already initialized)"

                    # Try to get actual conda executable from function variables
                    local func_conda="${CONDA_EXE:-${__conda_exe:-}}"

                    if [ -n "$func_conda" ] && [ -x "$func_conda" ]; then
                        local version=$("$func_conda" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "")
                        if [ -n "$version" ]; then
                            log_info "Found conda $version at $func_conda (from shell function)"
                            best_conda="$func_conda"
                            best_version="$version"
                        fi
                    fi
                fi
            fi
        fi

        # STRATEGY 3: Search PATH and common locations
        if [ -z "$best_conda" ]; then
            log_info "Searching for conda in PATH and common locations..."

            local conda_candidates=()

            # Check all conda executables in PATH (using which -a if available)
            if command -v which &>/dev/null; then
                while IFS= read -r conda_path; do
                    [ -n "$conda_path" ] && [ -x "$conda_path" ] && conda_candidates+=("$conda_path")
                done < <(which -a conda 2>/dev/null || true)
            fi

            # Check common installation locations
            local common_locations=(
                "$HOME/miniconda3/bin/conda"
                "$HOME/anaconda3/bin/conda"
                "$HOME/mambaforge/bin/conda"
                "$HOME/.conda/bin/conda"
                "/opt/conda/bin/conda"
                "/opt/miniconda3/bin/conda"
                "/opt/anaconda3/bin/conda"
                "/usr/local/miniconda3/bin/conda"
                "/usr/local/anaconda3/bin/conda"
            )

            for loc in "${common_locations[@]}"; do
                [ -x "$loc" ] && conda_candidates+=("$loc")
            done

            # Remove duplicates
            local unique_candidates=($(printf '%s\n' "${conda_candidates[@]}" | sort -u))

            if [ ${#unique_candidates[@]} -eq 0 ]; then
                error_exit "No conda installations found

  Please install Miniconda or Anaconda:
    - Download: https://docs.conda.io/en/latest/miniconda.html
    - Or use: curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-${OS_NAME}-${ARCH_NAME}.sh"
            fi

            # Check version of each candidate and pick most recent
            log_info "Evaluating ${#unique_candidates[@]} conda installation(s):"

            for conda_path in "${unique_candidates[@]}"; do
                local version=$("$conda_path" --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo "")

                if [ -z "$version" ]; then
                    log_warning "  $conda_path - could not determine version (skipping)"
                    continue
                fi

                log_info "  $conda_path - version $version"

                # Update best if this is newer (using sort -V for version comparison)
                if [ -z "$best_version" ] || [ "$(printf '%s\n' "$version" "$best_version" | sort -V | tail -1)" = "$version" ]; then
                    best_conda="$conda_path"
                    best_version="$version"
                fi
            done
        fi
    fi

    # Validate we found something
    if [ -z "$best_conda" ] || [ -z "$best_version" ]; then
        error_exit "Could not find any usable conda installation

  Please ensure conda is installed and accessible.

  Or specify conda path explicitly:
    $0 --conda-path /path/to/conda"
    fi

    # Check if best version meets minimum requirement
    if [ "$(printf '%s\n' "$best_version" "$min_version" | sort -V | head -1)" != "$min_version" ]; then
        if [ "$FORCE_OLD_CONDA" = true ]; then
            log_warning "WARNING: Using old conda version $best_version "
            log_warning ""
            log_warning "  Minimum recommended: $min_version"
            log_warning "  You are using: $best_version"
            log_warning "  Location: $best_conda"
            log_warning ""
            log_warning "  This may cause environment creation to FAIL with errors like:"
            log_warning "    - ResolvePackageNotFound"
            log_warning "    - Dependency conflicts"
            log_warning "    - Package version incompatibilities"
            log_warning ""
            log_warning "  If installation fails, update conda or use a newer version."
            log_warning ""
            log_info "Continuing anyway due to --force-old-conda flag..."
        else
            error_exit "Best conda found is too old: $best_version (minimum required: $min_version)

  Location: $best_conda

  SOLUTIONS:
    1. Update conda: conda update -n base conda
    2. Install newer miniconda: https://docs.conda.io/en/latest/miniconda.html
    3. On HPC: Ask admin to provide conda >= $min_version module
    4. Or specify a different conda manually:
       $0 --conda-path /path/to/newer/conda
    5. Or force use of this old conda (may fail):
       $0 --force-old-conda"
        fi
    fi

    log_success "Selected conda $best_version at $best_conda"

    # Initialize conda shell integration
    local conda_base=$("$best_conda" info --base 2>/dev/null || echo "")

    if [ -n "$conda_base" ] && [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
        # Source conda shell integration (standard conda activation pattern)
        source "$conda_base/etc/profile.d/conda.sh"
        log_success "Initialized conda shell integration from $conda_base"
    else
        log_warning "Could not initialize conda shell integration (may still work)"
    fi

    # Ensure selected conda is first in PATH
    local conda_bin_dir=$(dirname "$best_conda")
    export PATH="$conda_bin_dir:$PATH"

    return 0
}

# Run conda detection (unless explicitly skipped for testing)
log_info "=== Finding best conda installation ==="
if ! find_and_initialize_best_conda; then
    error_exit "Failed to find suitable conda installation

  Please ensure conda >= 23.1.0 is available, or specify one manually:
    $0 --conda-path /path/to/conda"
fi
echo ""

############################################################################
# STEP 1: CREATE CONDA ENVIRONMENT
############################################################################

# Check for conda/mamba
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    log_success "Using mamba (faster): $(mamba --version 2>&1 | head -1)"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    CONDA_VERSION=$(conda --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
    log_success "Using conda: conda $CONDA_VERSION"

    # Check conda version (need 23.10+ for libmamba solver)
    CONDA_MAJOR=$(echo "$CONDA_VERSION" | cut -d. -f1)
    CONDA_MINOR=$(echo "$CONDA_VERSION" | cut -d. -f2)

    if [ "$CONDA_MAJOR" -lt 23 ] || ([ "$CONDA_MAJOR" -eq 23 ] && [ "$CONDA_MINOR" -lt 10 ]); then
        log_warning "Conda version $CONDA_VERSION is older than recommended (23.10+)"
        log_warning ""
        log_warning "  Older conda versions may fail to resolve environment.yml dependencies."
        log_warning ""
        log_warning "  RECOMMENDED: Update conda first:"
        log_warning "    conda update -n base conda"
        log_warning ""
        log_warning "  Or use mamba (faster solver):"
        log_warning "    conda install -n base conda-libmamba-solver"
        log_warning "    conda config --set solver libmamba"
        log_warning ""
        log_info "Continuing anyway - if environment creation fails, update conda..."
    else
        log_success "Conda version OK (>= 23.10)"
    fi
else
    error_exit "Neither conda nor mamba found after detection!

  This is unexpected after STEP 0.5 succeeded. Possible causes:
    - Conda initialization failed
    - PATH was modified after detection
    - Conda executable is broken

  SOLUTIONS:
    1. Specify conda path explicitly:
       $0 --conda-path /path/to/conda
    2. Install fresh miniconda:
       curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-${OS_NAME}-${ARCH_NAME}.sh
    3. Check if conda module was unloaded (HPC): module load conda"
fi

# Check for environment.yml
ENV_FILE="$PROJECT_ROOT/environment.yml"
if [ ! -f "$ENV_FILE" ]; then
    error_exit "environment.yml not found at $ENV_FILE

  This file should be in the project root directory.
  If you're in the wrong directory, cd to the AncST root first."
fi

# Create or update the environment
if [ "$FORCE_RECREATE" = true ] && conda env list | grep -q "^$ENV_NAME "; then
    log_info "Force recreating environment '$ENV_NAME'..."

    # Clean conda cache when force recreating (prevents corrupted cache issues)
    log_info "Cleaning conda package cache..."
    conda clean --all -y 2>&1 | tee -a "$BUILD_LOG" || log_warning "Cache cleaning had some issues but continuing..."
    log_success "Cache cleaned"

    conda env remove -n "$ENV_NAME" -y || error_exit "Failed to remove environment"
    log_success "Existing environment removed"
    CLEAN_CACHE=false  # Already cleaned, don't clean again
elif [ "$CLEAN_CACHE" = true ]; then
    # Clean cache if explicitly requested (without force recreate)
    log_info "Cleaning conda package cache..."
    conda clean --all -y 2>&1 | tee -a "$BUILD_LOG" || log_warning "Cache cleaning had some issues but continuing..."
    log_success "Cache cleaned"
fi

log_info "Creating/updating conda environment from environment.yml..."
if conda env list | grep -q "^$ENV_NAME "; then
    log_info "Environment '$ENV_NAME' exists. Updating..."

    # Try normal update first
    if ! $CONDA_CMD env update -n "$ENV_NAME" -f "$ENV_FILE" --prune 2>&1 | tee -a "$BUILD_LOG"; then
        log_warning "Update failed - checking if it's a plugin issue..."

        # Check if error mentions plugins
        if grep -qi "plugin\|CONDA_NO_PLUGINS" "$BUILD_LOG" 2>/dev/null; then
            log_warning "Detected conda plugin error - retrying with --no-plugins..."

            # Retry with --no-plugins
            if $CONDA_CMD --no-plugins env update -n "$ENV_NAME" -f "$ENV_FILE" --prune 2>&1 | tee -a "$BUILD_LOG"; then
                log_success "Environment updated (with --no-plugins workaround)"
            else
                error_exit "Failed to update conda environment even with --no-plugins

  SOLUTIONS:
    1. Clean cache and recreate: $0 --env-name $ENV_NAME --force-recreate
    2. Just clean cache: conda clean --all && $0 --env-name $ENV_NAME
    3. Manually remove: conda env remove -n $ENV_NAME
    4. Check conda configuration: conda config --show

  Error details in: $BUILD_LOG"
            fi
        else
            error_exit "Failed to update conda environment

  SOLUTIONS:
    1. Remove and recreate: $0 --force-recreate
    2. Clean corrupted cache: $0 --clean-cache
    3. Or manually: conda clean --all && $0

  Error details in: $BUILD_LOG"
        fi
    else
        log_success "Environment updated!"
    fi
else
    if ! $CONDA_CMD env create -n "$ENV_NAME" -f "$ENV_FILE" 2>&1 | tee -a "$BUILD_LOG"; then
        log_error "Failed to create conda environment"
        echo ""

        # Check if error was due to environment already existing
        if grep -q "prefix already exists" "$BUILD_LOG" 2>/dev/null || \
           grep -q "CondaValueError.*already exists" "$BUILD_LOG" 2>/dev/null; then
            log_error "Environment already exists!"
            log_error ""
            log_error "The environment name '$ENV_NAME' (or default 'AncST') already exists."
            log_error ""
            log_error "SOLUTIONS:"
            log_error "  1. Remove existing environment first:"
            log_error "     conda env remove -n $ENV_NAME"
            log_error "     Then re-run: $0"
            log_error ""
            log_error "  2. Or use --force-recreate to remove and recreate:"
            log_error "     $0 --force-recreate"
            log_error ""
            log_error "  3. Or use a different environment name:"
            log_error "     $0 --env-name MyNewEnvName"
            exit 1
        fi

        # Check if error was due to conda plugins
        if grep -qi "plugin\|CONDA_NO_PLUGINS" "$BUILD_LOG" 2>/dev/null; then
            log_warning "Detected conda plugin error - retrying with --no-plugins..."

            # Retry with --no-plugins
            if $CONDA_CMD --no-plugins env create -n "$ENV_NAME" -f "$ENV_FILE" 2>&1 | tee -a "$BUILD_LOG"; then
                log_success "Environment created (with --no-plugins workaround)"
                # Skip to next section (don't re-execute the error handling below)
            else
                log_error "Failed even with --no-plugins"
                log_info ""
                log_info "Error details saved to: $BUILD_LOG"
                exit 1
            fi
        else
            # Not a plugin issue - show generic troubleshooting
            log_info "Common solutions:"
            log_info ""
            log_info "1. Clean corrupted package cache (RECOMMENDED FIRST):"
            log_info "     $0 --clean-cache"
            log_info "     OR: conda clean --all && ./install.sh"
            log_info ""
            log_info "2. Update conda to latest version (24.x):"
            log_info "     conda update -n base conda"
            log_info ""
            log_info "3. Install and use libmamba solver (faster, better dependency resolution):"
            log_info "     conda install -n base conda-libmamba-solver"
            log_info "     conda config --set solver libmamba"
            log_info "     Then re-run: ./install.sh"
            log_info ""
        log_info "4. If specific packages can't be found:"
        log_info "     - Update package index: conda update --all -n base -c conda-forge"
        log_info "     - Verify channels: conda config --show channels"
        log_info "       Should include: conda-forge, bioconda, defaults"
        log_info ""
        log_info "5. Check internet connection and proxy/firewall settings"
        log_info ""
            log_info "6. If conda version issues persist, try different conda installation:"
            log_info "     $0 --conda-path /path/to/different/conda"
            log_info ""
            log_info "Error details saved to: $BUILD_LOG"
            exit 1
        fi  # End plugin check
    fi  # End create failure check
    log_success "Environment created!"
fi  # End environment exists check

# Get conda environment path
CONDA_ENV_PATH=$(conda env list | grep "^$ENV_NAME " | awk '{print $NF}')
if [ -z "$CONDA_ENV_PATH" ]; then
    # Try common locations
    for base_path in "$HOME/opt/anaconda3" "$HOME/anaconda3" "$HOME/miniconda3" "$HOME/mambaforge"; do
        if [ -d "$base_path/envs/$ENV_NAME" ]; then
            CONDA_ENV_PATH="$base_path/envs/$ENV_NAME"
            break
        fi
    done
fi
log_info "Environment path: $CONDA_ENV_PATH"
if [ -z "$CONDA_ENV_PATH" ] || [ ! -d "$CONDA_ENV_PATH" ]; then
    error_exit "Could not determine conda environment path for '$ENV_NAME'"
fi

############################################################################
# STEP 1.5: CMAKE SETUP AND PATH CONFIGURATION
############################################################################

# If no system cmake was found, install it to conda environment
if [ -z "$SYSTEM_CMAKE" ]; then
    log_info "Installing cmake (>=3.5,<4.0) from conda-forge into $ENV_NAME environment..."
    conda install -n "$ENV_NAME" -y -c conda-forge "cmake>=3.5,<4.0" >> "$BUILD_LOG" 2>&1 || \
        error_exit "Failed to install cmake to conda environment

  Try manually:
    conda activate $ENV_NAME
    conda install -c conda-forge 'cmake>=3.5,<4.0'"
    log_success "cmake installed to conda environment (modern but backward-compatible)"
fi

# NOW export PATH to use conda environment binaries
# This ensures we use conda's genmap, blastn, etc. for the pipeline
# and conda's cmake if we just installed it
export PATH="$CONDA_ENV_PATH/bin:$PATH"

# Architecture-agnostic note: Compiler detection logic removed from here
# Instead, compilers are checked/unset right before building macle/clasp
# This is more robust and works across all architectures and conda versions

# macOS: If conda architecture mismatches system, force system compilers
if [ "$OS_NAME" = "macOS" ] && [ "$USE_SYSTEM_COMPILERS" = true ]; then
    log_info "macOS: Switching to system compilers (bypassing conda toolchain)"

    # Unset ALL conda compiler variables
    unset CC CXX FC
    unset CFLAGS CXXFLAGS CPPFLAGS LDFLAGS
    unset CONDA_BUILD_SYSROOT

    # Set to system compilers
    export CC="/usr/bin/clang"
    export CXX="/usr/bin/clang++"

    # Ensure system toolchain precedes conda wrappers
    export PATH="/usr/bin:/bin:/usr/sbin:/sbin:$PATH"

    log_info "Cleared conda compiler flags to use clean system toolchain"
fi

log_info "CC in use: $(command -v ${CC:-cc})"
log_info "CXX in use: $(command -v ${CXX:-c++})"

log_info "Using conda environment binaries for pipeline and build tools"

# Verify cmake is now available (either system or from conda env)
if ! command -v cmake &>/dev/null; then
    error_exit "cmake not found after installation attempt

  This is unexpected. Please report this issue at:
  https://github.com/Norsbus/AncST/issues"
fi

CMAKE_LOCATION=$(which cmake)
if [[ "$CMAKE_LOCATION" == "$CONDA_ENV_PATH"* ]]; then
    log_info "Using cmake from conda environment: $CMAKE_LOCATION"
else
    log_info "Using system cmake: $CMAKE_LOCATION"
fi

echo ""

############################################################################
# HELPER: Try compilation with fallback
############################################################################

try_compile_with_fallback() {
    local component_name="$1"
    local make_target="$2"

    log_info "Building $component_name with optimized flags..."

    # Try with native optimization first
    if make $make_target CFLAGS="$CFLAGS_NATIVE" CXXFLAGS="$CXXFLAGS_NATIVE" LDFLAGS="$LDFLAGS_NATIVE" >> "$BUILD_LOG" 2>&1; then
        log_success "$component_name built with native optimizations"
        return 0
    else
        log_warning "Native optimization flags failed for $component_name"
        log_warning "Falling back to generic -O3 optimization..."
        log_warning "Performance may be 5-15% slower than with -march=native"
        log_warning "    Consider upgrading GCC or filing an issue if this persists"

        # Clean and try with fallback flags
        make clean >> "$BUILD_LOG" 2>&1 || true

        if make $make_target CFLAGS="$CFLAGS_FALLBACK" CXXFLAGS="$CXXFLAGS_FALLBACK" LDFLAGS="$LDFLAGS_FALLBACK" >> "$BUILD_LOG" 2>&1; then
            log_success "$component_name built with fallback optimization"
            return 0
        else
            error_exit "Failed to build $component_name even with fallback flags

  See build log for details: $BUILD_LOG"
            return 1
        fi
    fi
}

############################################################################
# STEP 2: COMPILE MACLE FROM SOURCE
############################################################################

log_info "=== Compiling macle from source for $OS_NAME $ARCH_DISPLAY ==="

# Remove old binary from conda env (if exists from previous install)
rm -f "$CONDA_ENV_PATH/bin/macle"

MACLE_DIR="$PROJECT_ROOT/tmp/macle_build"
mkdir -p "$MACLE_DIR"

# Clone if needed
if [ ! -f "$MACLE_DIR/macle/Makefile" ]; then
    log_info "Cloning macle repository..."
    rm -rf "$MACLE_DIR/macle"
    git clone https://github.com/EvolBioInf/macle.git "$MACLE_DIR/macle" >> "$BUILD_LOG" 2>&1 || \
        error_exit "Failed to clone macle repository

  Repository: https://github.com/EvolBioInf/macle.git
  Possible causes:
    - No internet connection
    - GitHub is down or repository moved
    - Git authentication issues

  Try manually: git clone https://github.com/EvolBioInf/macle.git"
    log_success "macle repository cloned"
fi

cd "$MACLE_DIR/macle"

# Clean any previous builds
log_info "Cleaning previous macle builds..."
make clean >> "$BUILD_LOG" 2>&1 || true
rm -rf build 2>/dev/null || true

# CRITICAL for macOS: Set CMAKE_OSX_ARCHITECTURES for all CMake-based builds
# This prevents universal binary (fat file) issues with static libraries
# added 22.01.
if [ "$OS_NAME" = "macOS" ]; then
    export CMAKE_OSX_ARCHITECTURES="$MACOS_BUILD_ARCH"
    log_info "macOS: Setting CMAKE_OSX_ARCHITECTURES=$CMAKE_OSX_ARCHITECTURES for divsufsort"
fi

# CRITICAL: Clean up broken compiler variables before building
# When using custom --env-name, conda's activation scripts may set CC/CXX/FC to paths
# that reference the wrong environment name (from environment.yml's name: field)
# This is architecture-agnostic: we don't try to "fix" paths, we just unset broken ones
# and let the build system find compilers naturally via PATH or system defaults
log_info "Checking compiler environment variables..."
COMPILER_CLEANED=false

for var in CC CXX FC; do
    # Use ${!var:-} to safely handle unset variables (returns empty string if unset)
    compiler_path="${!var:-}"

    if [ -n "$compiler_path" ]; then
        # Check if compiler path exists and is executable
        if [ ! -x "$compiler_path" ]; then
            log_warning "$var points to non-existent or non-executable file: $compiler_path"
            unset $var
            COMPILER_CLEANED=true
            log_info "Unset $var - build system will find compiler naturally"
        # Additional check: if path contains /envs/ but not our actual environment name
        elif [[ "$compiler_path" == *"/envs/"* ]] && [[ "$compiler_path" != *"/envs/$ENV_NAME/"* ]]; then
            log_warning "$var points to wrong environment: $compiler_path"
            log_warning "    Expected environment: $ENV_NAME"
            unset $var
            COMPILER_CLEANED=true
            log_info "Unset $var - build system will use compilers from correct environment"
        else
            log_info "$var is valid: $compiler_path"
        fi
    fi
done

if [ "$COMPILER_CLEANED" = true ]; then
    log_success "Compiler variables cleaned - builds will use correct compilers from PATH"
    log_info "PATH includes: $CONDA_ENV_PATH/bin (correct environment)"
else
    log_success "All compiler variables are valid"
fi

# Build divsufsort dependency (uses CMake)
try_compile_with_fallback "divsufsort" "divsufsort" || log_warning "divsufsort build issues (may be non-critical)"

# Prepare and build SDSL
log_info "Preparing SDSL source code..."

# Clone SDSL if needed (don't build yet)
if [ ! -d "sdsl-lite" ]; then
    log_info "Cloning SDSL repository with submodules..."
    # Clone recursively to get libdivsufsort submodule
    git clone --recursive https://github.com/simongog/sdsl-lite.git >> "$BUILD_LOG" 2>&1 || \
        error_exit "Failed to clone SDSL repository"
else
    # If already cloned, make sure submodules are initialized
    cd sdsl-lite
    if [ ! -f "external/libdivsufsort/CMakeLists.txt" ]; then
        log_info "Initializing SDSL git submodules (libdivsufsort)..."
        git submodule update --init --recursive >> "$BUILD_LOG" 2>&1 || \
            log_warning "Failed to initialize submodules - SDSL build may fail"
    fi
    cd ..
fi

# Patch SDSL BEFORE building (known compatibility issue with modern compilers)
# Bug: louds_tree.hpp references m_select1/m_select0 but should be m_bv_select1/m_bv_select0
# See: https://www.mail-archive.com/freebsd-pkg-fallout@freebsd.org/msg2487620.html
if [ -f "sdsl-lite/include/sdsl/louds_tree.hpp" ]; then
    log_info "Applying SDSL compatibility patches for louds_tree.hpp bug..."
    if [ "$OS_NAME" = "macOS" ]; then
        # macOS uses BSD sed (requires '' for in-place edit)
        sed -i '' 's/tree\.m_select1/tree.m_bv_select1/g' sdsl-lite/include/sdsl/louds_tree.hpp
        sed -i '' 's/tree\.m_select0/tree.m_bv_select0/g' sdsl-lite/include/sdsl/louds_tree.hpp
        log_success "SDSL patches applied for macOS compatibility"
    else
        # Linux uses GNU sed
        sed -i 's/tree\.m_select1/tree.m_bv_select1/g' sdsl-lite/include/sdsl/louds_tree.hpp
        sed -i 's/tree\.m_select0/tree.m_bv_select0/g' sdsl-lite/include/sdsl/louds_tree.hpp
        log_success "SDSL patches applied for Linux compatibility"
    fi
fi

# NOW build SDSL (patches applied, submodules initialized)
# SDSL will use its own libdivsufsort from external/libdivsufsort submodule
log_info "Building SDSL library..."

# Check if SDSL already built successfully (make sdsl is not idempotent)
if [ -f "sdsl/lib/libsdsl.a" ]; then
    log_success "SDSL already built successfully - skipping rebuild"

    # On macOS, verify architecture matches current system
    if [ "$OS_NAME" = "macOS" ]; then
        CURRENT_ARCH=$(uname -m)
        LIBRARY_ARCH=$(lipo -archs sdsl/lib/libsdsl.a 2>/dev/null | xargs)

        if [ "$LIBRARY_ARCH" != "$CURRENT_ARCH" ]; then
            log_warning "Architecture mismatch detected!"
            log_warning "    Current system: $CURRENT_ARCH"
            log_warning "    SDSL library built for: $LIBRARY_ARCH"
            log_warning "    This may cause linking errors. Rebuilding SDSL..."

            # Remove old library and rebuild
            rm -rf sdsl sdsl-lite/build
            log_info "Cleaned old SDSL build - will rebuild for $CURRENT_ARCH"
        fi
    fi
fi

# Build SDSL if not already built (or if architecture mismatch on macOS)
if [ ! -f "sdsl/lib/libsdsl.a" ]; then
    # CRITICAL for macOS: Force single-architecture build
    # Conda's compiler on macOS defaults to universal binaries (arm64+x86_64)
    # This creates fat object files which CANNOT be put in static archives (.a)
    # See: https://cmake.org/pipermail/cmake/2013-September/055916.html
    # added 22.01.
    if [ "$OS_NAME" = "macOS" ]; then
        export CMAKE_OSX_ARCHITECTURES="$MACOS_BUILD_ARCH"
        log_info "macOS: Building for single architecture: $CMAKE_OSX_ARCHITECTURES"
        log_info "    (Prevents 'fat file not allowed in archive' error)"
    fi


    # CRITICAL: Call sdsl-lite/install.sh DIRECTLY (not through macle Makefile)
    # Reason: macle's Makefile has LDFLAGS with -ldivsufsort which pollutes SDSL's CMake build
    # This causes CMake compiler test to fail on macOS (strict linker rejects missing libs)
    # See: https://github.com/Norsbus/AncST/issues/... (LDFLAGS pollution bug)
    log_info "Building SDSL via install.sh (bypassing macle Makefile to avoid LDFLAGS pollution)..."

    if (cd sdsl-lite && ./install.sh "$MACLE_DIR/macle/sdsl" >> "$BUILD_LOG" 2>&1); then
        log_success "SDSL built and installed successfully"

        # Verify architecture on macOS
        if [ "$OS_NAME" = "macOS" ]; then
            BUILT_ARCH=$(lipo -archs sdsl/lib/libsdsl.a 2>/dev/null | xargs)
            log_info "SDSL library architecture: $BUILT_ARCH"
        fi
    else
        log_warning "SDSL build had issues - check if library exists anyway"
        # Sometimes build fails on second attempt but library is already there
        if [ -f "sdsl/lib/libsdsl.a" ]; then
            log_success "SDSL library found despite build warnings"
        else
            error_exit "SDSL build failed - see log: $BUILD_LOG"
        fi
    fi

    # Clean up environment variable after SDSL build
    if [ "$OS_NAME" = "macOS" ]; then
        unset CMAKE_OSX_ARCHITECTURES
    fi
fi

# Build macle - CRITICAL: Set up library paths for linking
log_info "Configuring macle build environment..."
export SDSL_ROOT="$MACLE_DIR/macle/sdsl"
export CPLUS_INCLUDE_PATH="$SDSL_ROOT/include:$MACLE_DIR/macle/src:${CPLUS_INCLUDE_PATH:-}"
export C_INCLUDE_PATH="$SDSL_ROOT/include:$MACLE_DIR/macle/src:${C_INCLUDE_PATH:-}"
export LIBRARY_PATH="$SDSL_ROOT/lib:${LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$SDSL_ROOT/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

# CRITICAL: Append divsufsort + sdsl library linking flags to LDFLAGS_NATIVE and LDFLAGS_FALLBACK
# These must be appended to the flags that try_compile_with_fallback passes to make
MACLE_LIB_FLAGS="-L$SDSL_ROOT/lib -lsdsl -ldivsufsort64 -ldivsufsort"
LDFLAGS_NATIVE="$LDFLAGS_NATIVE $MACLE_LIB_FLAGS"
LDFLAGS_FALLBACK="$LDFLAGS_FALLBACK $MACLE_LIB_FLAGS"

try_compile_with_fallback "macle" "" || error_exit "Failed to build macle - see log: $BUILD_LOG"

# Test and install binary into conda environment
if [ -f "build/macle" ]; then
    log_info "Testing macle binary..."
    # Test by calling --help (may exit with non-zero, use || true for pipefail)
    if ( ./build/macle --help 2>&1 || true ) | grep -qi "usage\|macle"; then
        log_info "Installing macle to conda environment bin directory..."
        cp build/macle "$CONDA_ENV_PATH/bin/"
        chmod +x "$CONDA_ENV_PATH/bin/macle"
        log_success "macle installed to $CONDA_ENV_PATH/bin/macle"
    else
        error_exit "macle compiled but failed execution test

  This suggests an architecture mismatch or missing runtime libraries.
  See build log: $BUILD_LOG"
    fi
else
    error_exit "macle binary not found after build

  Build may have failed silently. See log: $BUILD_LOG"
fi

cd "$PROJECT_ROOT"
echo ""

############################################################################
# STEP 3: COMPILE CLASP.X FROM SOURCE
############################################################################

log_info "=== Compiling clasp.x from source for $OS_NAME $ARCH_DISPLAY ==="

# Remove old binary from conda env (if exists from previous install)
rm -f "$CONDA_ENV_PATH/bin/clasp.x"

CLASP_DIR="$PROJECT_ROOT/tmp/clasp_build"
CLASP_URL="http://legacy.bioinf.uni-leipzig.de/Software/clasp/clasp_v1_1.tar.gz"

mkdir -p "$CLASP_DIR"
cd "$CLASP_DIR"

# Download if needed
if [ ! -f "clasp_v1_1.tar.gz" ]; then
    log_info "Downloading clasp..."
    curl -L -o clasp_v1_1.tar.gz "$CLASP_URL" >> "$BUILD_LOG" 2>&1 || \
        error_exit "Failed to download clasp

  URL: $CLASP_URL
  Possible causes:
    - No internet connection
    - Server is down
    - Firewall blocking access

  Try manually: curl -L -O $CLASP_URL"
    log_success "clasp downloaded"
fi

# Extract
log_info "Extracting clasp..."
tar -xzf clasp_v1_1.tar.gz >> "$BUILD_LOG" 2>&1 || error_exit "Failed to extract clasp archive"

# Find directory
CLASP_SRC=""
for dir in clasp clasp_v1_1 clasp_v1.1; do
    if [ -d "$dir" ]; then
        CLASP_SRC="$dir"
        break
    fi
done

if [ -z "$CLASP_SRC" ]; then
    error_exit "Could not find clasp source directory after extraction

  Expected one of: clasp, clasp_v1_1, clasp_v1.1
  See: $BUILD_LOG"
fi

cd "$CLASP_SRC"

# Clean previous builds
make clean >> "$BUILD_LOG" 2>&1 || true

# Patch Makefile for conda SSL and architecture flags
if [ -f "Makefile" ]; then
    log_info "Patching clasp Makefile..."
    cp Makefile Makefile.orig

    # Try native flags first
    if [ -d "$CONDA_ENV_PATH/lib" ]; then
        # Add conda paths and our CFLAGS
        if [ "$OS_NAME" = "macOS" ]; then
            sed -i '' "s|^CFLAGS=|CFLAGS=$CFLAGS_NATIVE -L$CONDA_ENV_PATH/lib -I$CONDA_ENV_PATH/include |" Makefile
        else
            sed -i "s|^CFLAGS=|CFLAGS=$CFLAGS_NATIVE -L$CONDA_ENV_PATH/lib -I$CONDA_ENV_PATH/include |" Makefile
        fi
    else
        # Just add our CFLAGS
        if [ "$OS_NAME" = "macOS" ]; then
            sed -i '' "s|^CFLAGS=|CFLAGS=$CFLAGS_NATIVE |" Makefile
        else
            sed -i "s|^CFLAGS=|CFLAGS=$CFLAGS_NATIVE |" Makefile
        fi
    fi
fi

# Build clasp
log_info "Building clasp.x..."
if ! make >> "$BUILD_LOG" 2>&1; then
    log_warning "Build with native flags failed, trying fallback..."

    # Restore Makefile and try with fallback
    if [ -f "Makefile.orig" ]; then
        cp Makefile.orig Makefile

        if [ -d "$CONDA_ENV_PATH/lib" ]; then
            if [ "$OS_NAME" = "macOS" ]; then
                sed -i '' "s|^CFLAGS=|CFLAGS=$CFLAGS_FALLBACK -L$CONDA_ENV_PATH/lib -I$CONDA_ENV_PATH/include |" Makefile
            else
                sed -i "s|^CFLAGS=|CFLAGS=$CFLAGS_FALLBACK -L$CONDA_ENV_PATH/lib -I$CONDA_ENV_PATH/include |" Makefile
            fi
        else
            if [ "$OS_NAME" = "macOS" ]; then
                sed -i '' "s|^CFLAGS=|CFLAGS=$CFLAGS_FALLBACK |" Makefile
            else
                sed -i "s|^CFLAGS=|CFLAGS=$CFLAGS_FALLBACK |" Makefile
            fi
        fi

        make clean >> "$BUILD_LOG" 2>&1 || true
        make >> "$BUILD_LOG" 2>&1 || error_exit "Failed to build clasp.x - see log: $BUILD_LOG"

        log_warning "clasp.x built with fallback optimization (may be 5-15% slower)"
    else
        error_exit "Failed to build clasp.x - see log: $BUILD_LOG"
    fi
else
    log_success "clasp.x built with native optimizations"
fi

# Test and install binary into conda environment
if [ -f "clasp.x" ]; then
    log_info "Testing clasp.x binary..."
    # clasp.x exits with error when run without required args, but outputs usage
    # Use || true to prevent pipefail from killing the test
    if ( ./clasp.x 2>&1 || true ) | grep -qi "usage\|clasp"; then
        log_info "Installing clasp.x to conda environment bin directory..."
        cp clasp.x "$CONDA_ENV_PATH/bin/"
        chmod +x "$CONDA_ENV_PATH/bin/clasp.x"
        log_success "clasp.x installed to $CONDA_ENV_PATH/bin/clasp.x"
    else
        error_exit "clasp.x compiled but failed execution test

  This suggests an architecture mismatch or missing runtime libraries.
  See build log: $BUILD_LOG"
    fi
else
    error_exit "clasp.x binary not found after build

  Build may have failed silently. See log: $BUILD_LOG"
fi

cd "$PROJECT_ROOT"
echo ""

############################################################################
# STEP 4: CREATE ACTIVATION SCRIPT
############################################################################

log_info "Creating activation script..."

cat > "$PROJECT_ROOT/activate_$ENV_NAME.sh" << EOF
#!/usr/bin/env bash
# AncST Environment Activation Script

PROJECT_ROOT="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"

# Activate conda environment
eval "\$(conda shell.bash hook)"
conda activate $ENV_NAME

# Note: macle and clasp.x are installed in conda env bin (automatically in PATH)

echo ""
echo "AncST environment activated!"
echo ""
echo "  Platform:    $OS_NAME $ARCH_DISPLAY"
echo "  Conda env:   $ENV_NAME (\$CONDA_PREFIX)"
echo "  Project:     \$PROJECT_ROOT"
echo ""
echo "Available tools (from conda environment):"
command -v snakemake &>/dev/null && echo "  snakemake \$(snakemake --version 2>&1 | head -1)"
command -v genmap &>/dev/null && echo "  genmap \$(genmap --version 2>&1 | head -1)"
command -v blastn &>/dev/null && echo "  blastn \$(blastn -version 2>&1 | head -1 | cut -d: -f2 | xargs)"
command -v macle &>/dev/null && echo "  macle (compiled for $OS_NAME $ARCH_DISPLAY)"
command -v clasp.x &>/dev/null && echo "  clasp.x (compiled for $OS_NAME $ARCH_DISPLAY)"
echo ""
echo "To run the pipeline:"
echo "  cd template && python run_complete_snakemake_manual.py"
echo ""
EOF

chmod +x "$PROJECT_ROOT/activate_$ENV_NAME.sh"
log_success "Activation script created"

############################################################################
# SUCCESS: CLEANUP AND FINAL VERIFICATION
############################################################################

# Clean up build artifacts on success
log_info "Cleaning build artifacts..."
rm -rf "$PROJECT_ROOT/tmp/macle_build"
rm -rf "$PROJECT_ROOT/tmp/clasp_build"
log_success "Build artifacts cleaned"

echo ""
log_info "================================================================"
log_success "                 Installation Complete!                        "
log_info "================================================================"
echo ""
log_info "Built for:         $OS_NAME $ARCH_DISPLAY"
log_info "Conda environment: $ENV_NAME"
log_info "Build log:         $BUILD_LOG (preserved for reference)"
echo ""
log_info "To activate the environment:"
echo ""
echo "  source $PROJECT_ROOT/activate_$ENV_NAME.sh"
echo ""
log_info "Or manually:"
echo ""
echo "  conda activate $ENV_NAME"
echo ""
log_info "To verify installation:"
echo ""
log_info "To run the pipeline:"
echo ""
echo "  cd template && python run_complete_snakemake_manual.py"
echo ""
