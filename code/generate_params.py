#!/usr/bin/env python3
"""
Auto-generate parameters for GenMap, macle, and dups based on genome characteristics.

Usage:
    python generate_params.py --genmap --orgs orgs.txt --output genmap_params.txt
    python generate_params.py --macle --orgs orgs.txt --output macle_params.txt
    python generate_params.py --dups --orgs orgs.txt --output dups_params.txt
"""

import argparse
import os
import sys
import pathlib
import yaml
from math import ceil, log

def load_pipeline_config(work_dir):
    config_file = pathlib.Path(work_dir) / 'pipeline_config.yaml'
    if not config_file.is_file():
        print(f"ERROR: Config file not found: {config_file}", file=sys.stderr)
        sys.exit(1)
    with open(config_file,'r') as f:
        return yaml.safe_load(f)

DEFAULTS = {}

def get_genome_size(genome_path):
    """Calculate total genome size from FASTA file."""
    total_size = os.path.getsize(f'{genome_path}')
    return total_size

def calculate_k_value(genome_size):
    """
    Calculate optimal K value for GenMap based on genome size.
    Formula: K = ceil(log₄(genome_size))
    """
    k = ceil(log(genome_size, 4))
    return k

def generate_genmap_params(orgs_file, genomes_dir, output_file):
    """Generate GenMap parameters for each genome."""
    print(f"Generating GenMap parameters...")

    params = []
    with open(orgs_file, 'r') as f:
        for line in f:
            org = line.strip()
            if not org or org.startswith('#'):
                continue

            genome_file = os.path.join(genomes_dir, f"{org}.fasta")
            if not os.path.exists(genome_file):
                print(f"Warning: Genome file not found: {genome_file}", file=sys.stderr)
                continue

            genome_size = get_genome_size(genome_file)
            k_value = calculate_k_value(genome_size)

            # GenMap parameters: organism k e L I percentile flag
            # Defaults loaded from pipeline_config.yaml
            gm = DEFAULTS.get('genmap', {})
            e = gm.get('e', 0)
            L = gm.get('L', 300)
            I = gm.get('I', 50)
            percentile = gm.get('percentile', 42)
            params.append(f"{org} {k_value} {e} {L} {I} {percentile} 0")
            print(f"  {org}: genome_size={genome_size:,}, K={k_value}")

    with open(output_file, 'w') as f:
        f.write('\n'.join(params) + '\n')

    print(f"GenMap parameters written to: {output_file}")
    return len(params)

def generate_macle_params(orgs_file, genomes_dir, output_file):
    """Generate macle parameters for each genome."""
    print(f"Generating macle parameters...")

    params = []
    with open(orgs_file, 'r') as f:
        for line in f:
            org = line.strip()
            if not org or org.startswith('#'):
                continue

            genome_file = os.path.join(genomes_dir, f"{org}.fasta")
            if not os.path.exists(genome_file):
                print(f"Warning: Genome file not found: {genome_file}", file=sys.stderr)
                continue

            # macle parameters: organism w p percentile flag
            # Defaults loaded from pipeline_config.yaml
            mc = DEFAULTS.get('macle', {})
            w = mc.get('w', 1000)
            p = mc.get('p', 50)
            percentile = mc.get('percentile', 42)
            params.append(f"{org} {w} {p} {percentile} 0")
            print(f"  {org}: w={w}, p={p}, percentile={percentile}")

    with open(output_file, 'w') as f:
        f.write('\n'.join(params) + '\n')

    print(f"macle parameters written to: {output_file}")
    return len(params)

def generate_dups_params(orgs_file, genomes_dir, output_file):
    """Generate duplication detection parameters."""
    print(f"Generating dups parameters...")

    params = []
    with open(orgs_file, 'r') as f:
        for line in f:
            org = line.strip()
            if not org or org.startswith('#'):
                continue

            genome_file = os.path.join(genomes_dir, f"{org}.fasta")
            if not os.path.exists(genome_file):
                print(f"Warning: Genome file not found: {genome_file}", file=sys.stderr)
                continue

            # dups parameters: organism k1 e1 k2 e2 w p x1 y1 x2 y2
            # Defaults loaded from pipeline_config.yaml
            dp = DEFAULTS.get('dups', {})
            k1 = dp.get('k1', 50)
            e1 = dp.get('e1', 4)
            k2 = dp.get('k2', 21)
            e2 = dp.get('e2', 2)
            w = dp.get('w', 300)
            p = dp.get('p', 50)
            x1 = dp.get('x1', 16)
            y1 = dp.get('y1', 90)
            x2 = dp.get('x2', 10)
            y2 = dp.get('y2', 10)
            params.append(f"{org} {k1} {e1} {k2} {e2} {w} {p} {x1} {y1} {x2} {y2}")
            print(f"  {org}: k1={k1}, e1={e1}, k2={k2}, e2={e2}, w={w}, p={p}, x1={x1}, y1={y1}, x2={x2}, y2={y2}")

    with open(output_file, 'w') as f:
        f.write('\n'.join(params) + '\n')

    print(f"dups parameters written to: {output_file}")
    return len(params)

def main():
    parser = argparse.ArgumentParser(
        description='Auto-generate pipeline parameters based on genome characteristics'
    )

    # Parameter type selection (mutually exclusive)
    param_group = parser.add_mutually_exclusive_group(required=True)
    param_group.add_argument('--genmap', action='store_true',
                            help='Generate GenMap parameters')
    param_group.add_argument('--macle', action='store_true',
                            help='Generate macle parameters')
    param_group.add_argument('--dups', action='store_true',
                            help='Generate dups parameters')

    # Required arguments
    parser.add_argument('--orgs', required=True,
                       help='Path to orgs file (list of genome IDs)')
    parser.add_argument('--output', required=True,
                       help='Output file path for parameters')
    parser.add_argument('--work-dir', required=True,
                       help='Work dir containing pipeline_config.yaml')

    # Optional arguments
    parser.add_argument('--genomes-dir', default='../utils/genomes',
                       help='Directory containing genome FASTA files')

    args = parser.parse_args()

    global DEFAULTS
    DEFAULTS = load_pipeline_config(args.work_dir).get('defaults', {})

    # Validate inputs
    if not os.path.exists(args.orgs):
        print(f"Error: Orgs file not found: {args.orgs}", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.genomes_dir):
        print(f"Error: Genomes directory not found: {args.genomes_dir}", file=sys.stderr)
        sys.exit(1)

    # Generate parameters
    try:
        if args.genmap:
            count = generate_genmap_params(args.orgs, args.genomes_dir, args.output)
        elif args.macle:
            count = generate_macle_params(args.orgs, args.genomes_dir, args.output)
        elif args.dups:
            count = generate_dups_params(args.orgs, args.genomes_dir, args.output)

        print(f"\nSuccess! Generated parameters for {count} genomes")
        return 0

    except Exception as e:
        print(f"Error generating parameters: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    sys.exit(main())
