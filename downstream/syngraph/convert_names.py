#!/usr/bin/env python3
"""
Convert between internal (Xorg/XorgchrN) and NCBI accession names in syngraph
output files, using the AncST mapping file.

Supports single files or syngraph output prefixes (converts all matching files).

Usage:
  # Single file
  python convert_names.py \\
    --mapping utils/mapping \\
    --input syngraph_output.tsv \\
    --output converted_output.tsv \\
    --direction internal_to_ncbi

  # Prefix mode: converts all files matching prefix.*
  python convert_names.py \\
    --mapping utils/mapping \\
    --prefix infer_out_syngraph_build_foo \\
    --direction internal_to_ncbi
"""

import argparse
import glob
import os
import re
import sys


def parse_mapping(filepath):
    """
    Parse the AncST utils/mapping file.

    Returns:
        org_mapping:  dict  accession <-> alias  (bidirectional)
        chr_mapping:  dict  alias <-> accession  (bidirectional, flat)
    """
    org_mapping = {}
    chr_flat = {}

    with open(filepath) as f:
        reading_species = False
        current_accession = None

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if 'species' in line:
                    reading_species = True
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            if reading_species:
                accession, alias = parts[0], parts[1]
                org_mapping[accession] = alias
                org_mapping[alias] = accession
                current_accession = accession
                reading_species = False
            else:
                if current_accession is None:
                    continue
                acc, alias = parts[0], parts[1]
                chr_flat[acc] = alias
                chr_flat[alias] = acc

    return org_mapping, chr_flat


def build_replacement_map(org_mapping, chr_flat, direction):
    """
    Build a replacement map.

    direction: 'internal_to_ncbi' or 'ncbi_to_internal'
    """
    replacements = {}

    if direction == 'internal_to_ncbi':
        for key, val in chr_flat.items():
            if re.match(r'^\d+orgchr\d+$', key):
                replacements[key] = val
        for key, val in org_mapping.items():
            if re.match(r'^\d+org$', key):
                replacements[key] = val
    else:
        for key, val in chr_flat.items():
            if not re.match(r'^\d+orgchr\d+$', key):
                replacements[key] = val
        for key, val in org_mapping.items():
            if not re.match(r'^\d+org$', key):
                replacements[key] = val

    return replacements


def convert_text(text, replacements):
    """Replace all occurrences, longest keys first to avoid partial matches."""
    sorted_keys = sorted(replacements.keys(), key=len, reverse=True)
    if not sorted_keys:
        return text
    pattern = '|'.join(re.escape(k) for k in sorted_keys)
    return re.sub(pattern, lambda m: replacements[m.group(0)], text)


def convert_file(inpath, outpath, replacements):
    """Convert a single file, writing to outpath."""
    with open(inpath) as fin, open(outpath, 'w') as fout:
        for line in fin:
            fout.write(convert_text(line, replacements))
    print(f"  {inpath} -> {outpath}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert names between internal (Xorg) and NCBI '
                    'accessions in syngraph files.',
    )
    parser.add_argument('--mapping', required=True,
                        help='Path to AncST utils/mapping file')
    parser.add_argument('--input',
                        help='Single input file to convert')
    parser.add_argument('--output',
                        help='Single output file path (required with --input)')
    parser.add_argument('--prefix',
                        help='Syngraph output prefix. Converts all files '
                             'matching prefix.* and writes to '
                             'prefix.converted.* (or use --output-prefix)')
    parser.add_argument('--output-prefix',
                        help='Output prefix for prefix mode '
                             '(default: {prefix}.converted)')
    parser.add_argument('--direction', required=True,
                        choices=['internal_to_ncbi', 'ncbi_to_internal'],
                        help='Conversion direction')

    args = parser.parse_args()

    if not args.input and not args.prefix:
        parser.error("Either --input or --prefix is required")
    if args.input and not args.output:
        parser.error("--output is required with --input")
    if args.input and args.prefix:
        parser.error("Use either --input or --prefix, not both")

    org_mapping, chr_flat = parse_mapping(args.mapping)
    replacements = build_replacement_map(org_mapping, chr_flat, args.direction)
    print(f"Loaded {len(replacements)} name mappings ({args.direction})")

    if args.input:
        convert_file(args.input, args.output, replacements)
    else:
        files = sorted(glob.glob(args.prefix + '.*'))
        if not files:
            print(f"No files found matching {args.prefix}.*", file=sys.stderr)
            sys.exit(1)

        files = [f for f in files if not f.endswith('.pickle')]

        out_prefix = args.output_prefix or args.prefix + '.converted'
        print(f"Converting {len(files)} files:")

        for inpath in files:
            suffix = inpath[len(args.prefix):]
            outpath = out_prefix + suffix
            convert_file(inpath, outpath, replacements)

        print(f"Done. Output prefix: {out_prefix}")


if __name__ == '__main__':
    main()
