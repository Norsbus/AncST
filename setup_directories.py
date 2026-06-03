#!/usr/bin/env python3
"""
Unified Directory Setup for AncST Pipeline
Orchestrates directory creation across template/, utils/, and anchors/

Behavior:
- template/ directories: CLEANED and recreated each run (working dirs)
- utils/ directories: Created if missing, NOT deleted (preprocessing outputs)
- anchors/ directories: Created if missing, NOT deleted (results)
"""

import os
import sys
import pathlib
from subprocess import run

def main():
    # Get project root
    script_dir = pathlib.Path(__file__).parent.resolve()

    print("=" * 60)
    print("AncST Pipeline - Directory Setup")
    print("=" * 60)

    # Step 1: Setup template/ directories (CLEAN + RECREATE)
    print("\n[1/3] Setting up template/ working directories (clean run)...")
    template_script = script_dir / 'template' / 'make_directories.py'
    if template_script.exists():
        result = run(['python', str(template_script)], cwd=script_dir / 'template')
        if result.returncode != 0:
            print(f"ERROR: Failed to setup template directories")
            sys.exit(1)
        print("Template directories ready")
    else:
        print(f"WARNING: {template_script} not found, skipping")

    # Step 2: Setup utils/ directories (CREATE IF MISSING, don't delete)
    # Note: orgs files in utils/ and utils/util_code/ are copied by run_pipeline.sh
    # before this script runs. They must already exist.
    print("\n[2/3] Setting up utils/ preprocessing directories...")
    utils_script = script_dir / 'utils' / 'make_directories.py'
    if utils_script.exists():
        result = run(['python', str(utils_script)], cwd=script_dir / 'utils')
        if result.returncode != 0:
            print(f"ERROR: Failed to setup utils directories")
            sys.exit(1)
        print("Utils directories ready")
    else:
        print(f"WARNING: {utils_script} not found, skipping")

    # Step 3: Setup anchors/ directories (CREATE IF MISSING, don't delete)
    print("\n[3/3] Setting up anchors/ result directories...")
    anchors_script = script_dir / 'utils' / 'util_code' / 'make_anchor_directories.py'
    if anchors_script.exists():
        result = run(['python', str(anchors_script)], cwd=script_dir / 'utils' / 'util_code')
        if result.returncode != 0:
            print(f"ERROR: Failed to setup anchors directories")
            sys.exit(1)
        print("Anchors directories ready")
    else:
        print(f"WARNING: {anchors_script} not found, skipping")

    print("\n" + "=" * 60)
    print("Directory setup complete!")
    print("=" * 60)

    return 0

if __name__ == "__main__":
    sys.exit(main())
