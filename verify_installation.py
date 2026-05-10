#!/usr/bin/env python3
"""
CAncST installation verification.
Checks Python version, packages, PuLP solver, external tools, and directory layout.
"""

import sys
import subprocess
from pathlib import Path

def check_python_version():
    print("=" * 70)
    print("PYTHON VERSION")
    print("=" * 70)
    version = sys.version_info
    print(f"Python {version.major}.{version.minor}.{version.micro}")
    if version.major == 3 and version.minor >= 8:
        print("OK: Python version is compatible (3.8+)")
        return True
    else:
        print("FAIL: Python version too old (need 3.8+)")
        return False

def check_python_packages():
    print("\n" + "=" * 70)
    print("PYTHON PACKAGES (Core Pipeline)")
    print("=" * 70)

    packages = [
        ('snakemake', '9.14.5'),
        ('Bio', '1.86'),
        ('numpy', '2.4.0'),
        ('pulp', '2.8.0'),
        ('yaml', '6.0.3'),
        ('pandas', '2.3.3'),
        ('psutil', '7.2.1'),
        ('requests', '2.32.5'),
    ]

    all_ok = True
    for pkg, min_version in packages:
        try:
            module = __import__(pkg)
            version = getattr(module, '__version__', 'unknown')
            print(f"OK   {pkg:20} {version}")
        except ImportError:
            print(f"FAIL {pkg:20} NOT INSTALLED (need >={min_version})")
            all_ok = False

    return all_ok

def check_pulp_solver():
    print("\n" + "=" * 70)
    print("PULP SOLVER")
    print("=" * 70)

    try:
        import pulp
        solvers = pulp.listSolvers(onlyAvailable=True)
        print(f"Available solvers: {', '.join(solvers)}")

        if 'COIN_CMD' in solvers:
            print("OK: COIN_CMD solver available (required)")
            coin_ok = True
        else:
            print("FAIL: COIN_CMD solver NOT available (required)")
            coin_ok = False

        if 'GUROBI_CMD' in solvers:
            print("OK: GUROBI_CMD solver available (optional, faster)")
        else:
            print("     GUROBI_CMD solver not available (optional)")

        return coin_ok
    except ImportError:
        print("FAIL: PuLP not installed")
        return False

def check_external_tools():
    print("\n" + "=" * 70)
    print("EXTERNAL TOOLS")
    print("=" * 70)

    tools = [
        ('genmap', 'required'),
        ('macle', 'required'),
        ('blastn', 'required'),
        ('makeblastdb', 'required'),
        ('fasta-splitter', 'required'),
    ]

    all_ok = True
    for tool, status in tools:
        result = subprocess.run(['which', tool], capture_output=True, text=True)
        if result.returncode == 0:
            path = result.stdout.strip()
            print(f"OK   {tool:20} {path}")
        else:
            print(f"FAIL {tool:20} NOT FOUND ({status})")
            if status == 'required':
                all_ok = False

    return all_ok

def check_downstream_packages():
    print("\n" + "=" * 70)
    print("PYTHON PACKAGES (Downstream - Optional)")
    print("=" * 70)

    packages = [
        'matplotlib',
        'scipy',
        'networkx',
        'ete3',
        'pygenomeviz',
        'dbscan1d',
        'markov_clustering',
        'tralda',
    ]

    installed = []
    missing = []

    for pkg in packages:
        try:
            module = __import__(pkg)
            version = getattr(module, '__version__', 'unknown')
            installed.append(f"OK   {pkg:20} {version}")
        except ImportError:
            missing.append(f"     {pkg:20} not installed")

    for item in installed:
        print(item)
    for item in missing:
        print(item)

    if missing:
        print("\nNote: These packages are optional for visualization/analysis.")
        print("      Install with the downstream/install_downstream_requirements.sh script.")

def check_directory_structure():
    print("\n" + "=" * 70)
    print("DIRECTORY STRUCTURE")
    print("=" * 70)

    base = Path(".")
    dirs = ['code', 'utils', 'template', 'downstream']

    all_ok = True
    for d in dirs:
        path = base / d
        if path.exists():
            print(f"OK   {d}/")
        else:
            print(f"FAIL {d}/ NOT FOUND")
            all_ok = False

    return all_ok

def main():
    print("\n" + "=" * 70)
    print("CAncST INSTALLATION VERIFICATION")
    print("=" * 70)

    results = {
        'Python version': check_python_version(),
        'Python packages': check_python_packages(),
        'PuLP solver': check_pulp_solver(),
        'External tools': check_external_tools(),
        'Directory structure': check_directory_structure(),
    }

    check_downstream_packages()

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    for category, status in results.items():
        marker = "OK  " if status else "FAIL"
        print(f"{marker} {category}")

    all_ok = all(results.values())

    if all_ok:
        print("\nAll core dependencies are installed. You can run the CAncST pipeline.")
        return 0
    else:
        print("\nSome dependencies are missing. See above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
