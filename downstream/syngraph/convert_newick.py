#!/usr/bin/env python3
"""
Convert leaf names in Newick tree files between NCBI accessions and internal
(Xorg) aliases using the AncST mapping file.

Usage:
  python convert_newick.py \\
    --mapping utils/mapping \\
    --input tree.nwk \\
    --output converted_tree.nwk \\
    --direction ncbi_to_internal

  python convert_newick.py \\
    --mapping utils/mapping \\
    --input tree.nwk \\
    --direction ncbi_to_internal   # prints to stdout if --output not given
"""

import argparse
import re
import sys


def parse_mapping(filepath):
    """Parse the AncST mapping file. Returns bidirectional accession <-> alias dict."""
    org_mapping = {}
    with open(filepath) as f:
        reading_species = False
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
                reading_species = False
    return org_mapping


def tokenize_newick(nwk_string):
    """
    Tokenize a Newick string into meaningful tokens.

    Returns list of (token_type, value) tuples:
        'LPAREN', 'RPAREN', 'COMMA', 'COLON', 'SEMICOLON',
        'LABEL', 'NUMBER', 'COMMENT'
    """
    tokens = []
    i = 0
    n = len(nwk_string)

    while i < n:
        c = nwk_string[i]

        if c in ' \t\n\r':
            i += 1
            continue

        if c == '(':
            tokens.append(('LPAREN', c))
            i += 1
        elif c == ')':
            tokens.append(('RPAREN', c))
            i += 1
        elif c == ',':
            tokens.append(('COMMA', c))
            i += 1
        elif c == ':':
            tokens.append(('COLON', c))
            i += 1
        elif c == ';':
            tokens.append(('SEMICOLON', c))
            i += 1
        elif c == '[':
            j = nwk_string.index(']', i)
            tokens.append(('COMMENT', nwk_string[i:j + 1]))
            i = j + 1
        elif c == "'":
            j = i + 1
            while j < n and nwk_string[j] != "'":
                if nwk_string[j] == '\\':
                    j += 1
                j += 1
            label = nwk_string[i + 1:j]
            tokens.append(('LABEL', label))
            i = j + 1
        elif c == '"':
            j = i + 1
            while j < n and nwk_string[j] != '"':
                if nwk_string[j] == '\\':
                    j += 1
                j += 1
            label = nwk_string[i + 1:j]
            tokens.append(('LABEL', label))
            i = j + 1
        else:
            j = i
            while j < n and nwk_string[j] not in '(),:;[] \t\n\r':
                j += 1
            text = nwk_string[i:j]
            if tokens and tokens[-1][0] == 'COLON':
                tokens.append(('NUMBER', text))
            else:
                tokens.append(('LABEL', text))
            i = j

    return tokens


def convert_newick(nwk_string, org_mapping, direction, strip_internal=False):
    """
    Convert leaf labels in a Newick string.

    direction: 'ncbi_to_internal' or 'internal_to_ncbi'
    strip_internal: if True, remove internal node labels
    """
    conversion = {}
    for key, val in org_mapping.items():
        if direction == 'ncbi_to_internal':
            if not re.match(r'^\d+org$', key):
                conversion[key] = val
        else:
            if re.match(r'^\d+org$', key):
                conversion[key] = val

    tokens = tokenize_newick(nwk_string)

    result_parts = []
    for idx, (ttype, tval) in enumerate(tokens):
        if ttype == 'LABEL':
            prev_meaningful = None
            for pidx in range(idx - 1, -1, -1):
                if tokens[pidx][0] not in ('COMMENT',):
                    prev_meaningful = tokens[pidx][0]
                    break

            if prev_meaningful == 'RPAREN':
                if not strip_internal:
                    result_parts.append(tval)
            else:
                converted = conversion.get(tval, tval)
                result_parts.append(converted)
        elif ttype == 'NUMBER':
            result_parts.append(tval)
        elif ttype == 'COMMENT':
            result_parts.append(tval)
        else:
            result_parts.append(tval)

    return ''.join(result_parts)


def direction_matches(key, direction):
    """Check if a key belongs to the source side of the conversion."""
    is_internal = bool(re.match(r'^\d+org$', key))
    if direction == 'ncbi_to_internal':
        return not is_internal
    else:
        return is_internal


def main():
    parser = argparse.ArgumentParser(
        description='Convert leaf names in Newick trees between NCBI '
                    'accessions and internal (Xorg) aliases.',
    )
    parser.add_argument('--mapping', required=True,
                        help='Path to AncST utils/mapping file')
    parser.add_argument('--input', required=True,
                        help='Input Newick file')
    parser.add_argument('--output',
                        help='Output Newick file (stdout if not given)')
    parser.add_argument('--direction', required=True,
                        choices=['ncbi_to_internal', 'internal_to_ncbi'],
                        help='Conversion direction')
    parser.add_argument('--strip-internal', action='store_true',
                        help='Remove internal node labels (e.g. Inner14)')

    args = parser.parse_args()

    org_mapping = parse_mapping(args.mapping)

    with open(args.input) as f:
        nwk_content = f.read()

    output_lines = []
    for line in nwk_content.splitlines():
        stripped = line.strip()
        if stripped and ';' in stripped:
            converted = convert_newick(stripped, org_mapping, args.direction,
                                       strip_internal=args.strip_internal)
            output_lines.append(converted)
        else:
            output_lines.append(line)

    result = '\n'.join(output_lines)
    if not result.endswith('\n'):
        result += '\n'

    if args.output:
        with open(args.output, 'w') as f:
            f.write(result)
        print(f"Converted: {args.input} -> {args.output} ({args.direction})",
              file=sys.stderr)
    else:
        sys.stdout.write(result)


if __name__ == '__main__':
    main()
