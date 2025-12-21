#! /usr/bin/env python3

"""
Draw verbose synteny visualizations with pygenomeviz.

This script creates detailed synteny plots showing genomic regions of interest
with alignment links between organisms. It reads coordinates from a coords file,
pairwise alignment data, and determines chromosome orientations automatically.

Usage:
    ./draw_verbose.py MARGIN [REF_ORG REF_CHROMO REF_START] [OPTIONS]

    # Basic usage with margin
    ./draw_verbose.py 50000

    # With reference element
    ./draw_verbose.py 50000 GCF_016746395.2 NC_052520.2 80300

    # With custom paths
    ./draw_verbose.py 50000 --coords-file ../coords --alignments-file pairwise_alignments_table

    # Auto-generate plotting order
    ./draw_verbose.py 50000 --auto-order

    # Use GFF files instead of coords
    ./draw_verbose.py 50000 --gff-dir ~/test_out_dir/out/stable/synthology/halos/gff

    # Draw all elements from GFF (not just connected component)
    ./draw_verbose.py 50000 --gff-dir ~/test_out_dir/out/stable/synthology/halos/gff --draw-all-elements

Arguments:
    MARGIN: Integer margin around elements (required)
    REF_ORG, REF_CHROMO, REF_START: Reference element to focus on (optional)

Options:
    --coords-file: Path to coords file (default: ../coords)
    --alignments-file: Path to pairwise alignments table (default: pairwise_alignments_table)
    --output-dir: Output directory for images (default: images_draw_verbose)
    --auto-order: Auto-generate plotting_order file before drawing

Reads from:
    - coords file: genomic coordinates of elements of interest
    - pairwise_alignments_table: alignment data between organisms
    - plotting_order (optional): order of tracks in figure
    - try_to_draw_all (optional): flag to draw all chromosomes

Writes to:
    - images_draw_verbose/{margin}/{ref}_*.svg: synteny figure
    - tables_drawing/{org}.tsv: gene feature tables
"""

from sys import argv
import argparse
from matplotlib.pyplot import cm
from pygenomeviz import GenomeViz
from random import randint
import numpy as np
from subprocess import run
from pprint import pprint
import pickle
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from os.path import isfile

### IMPORTANT NOTES $$$

# 1. since the chromsome/contig lengths is not explicitely known here, the drawing can be longer than a chromosome/contig

def parse_gff_directory(gff_dir):
    """
    Parse all GFF files in a directory and convert to coords format.

    Returns same structure as coords file parsing:
    - orgs: list of organism IDs
    - coords: dict of {org: {chromo: [(start, end, strand, gene_name), ...]}}
    - to_ori: set of (org, chromo) tuples
    """
    import os
    from pathlib import Path

    orgs = []
    coords = {}
    to_ori = set()

    gff_files = sorted(Path(gff_dir).glob('*.gff'))
    if not gff_files:
        gff_files = sorted(Path(gff_dir).glob('*.gff3'))

    for gff_file in gff_files:
        # Extract organism ID from filename (e.g., GCF_000001215.4.gff -> GCF_000001215.4)
        org = gff_file.stem
        orgs.append(org)
        coords[org] = {}

        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                chromo = fields[0]
                feature_type = fields[2]
                start = fields[3]
                end = fields[4]
                strand_symbol = fields[6]
                attributes = fields[8]

                # Only process CDS features
                if feature_type != 'CDS':
                    continue

                # Convert strand to orientation
                if strand_symbol == '+':
                    strand = 1
                    ori_gene = 'forward'
                elif strand_symbol == '-':
                    strand = -1
                    ori_gene = 'reverse'
                else:
                    continue

                # Extract gene name from attributes
                gene_name = None
                for attr in attributes.split(';'):
                    if attr.startswith('gene='):
                        gene_name = attr.split('=')[1]
                        break
                    elif attr.startswith('Name='):
                        gene_name = attr.split('=')[1]
                        break

                if not gene_name:
                    # Fall back to protein_id
                    for attr in attributes.split(';'):
                        if attr.startswith('protein_id='):
                            gene_name = attr.split('=')[1]
                            break

                if not gene_name:
                    gene_name = f"{chromo}_{start}_{end}"

                # Add to coords structure
                if chromo not in coords[org]:
                    coords[org][chromo] = []
                    to_ori.add((org, chromo))

                coords[org][chromo].append((int(start), int(end), strand, gene_name))

    return orgs, coords, to_ori


def turn(start,end,segments):
    max_e = 0
    min_s = 1e15
    for segment in segments:
        max_e = max(segment[1],max_e)
        min_s = min(segment[0],min_s)
    len_line = max_e - min_s
    new_end = len_line - start
    new_start = new_end - (end - start)
    return new_start, new_end

def overlaps(ranges):
    """Merge overlapping ranges and combine their color sets and strand info."""
    ranges = sorted(ranges)
    it = iter(ranges)
    try:
        curr_start, curr_stop, curr_color, curr_strands = next(it)
    except StopIteration:
        return
    for start, stop, color, strands in it:
        if curr_start <= start <= curr_stop:
            curr_stop = max(curr_stop, stop)
            curr_color.update(color)
            curr_strands.update(strands)
        else:
            yield curr_start, curr_stop, curr_color, curr_strands
            curr_start, curr_stop, curr_color, curr_strands = start, stop, color, strands
    yield curr_start, curr_stop, curr_color, curr_strands

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description='Draw verbose synteny visualizations with pygenomeviz',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  %(prog)s 50000
  %(prog)s 50000 GCF_016746395.2 NC_052520.2 80300
  %(prog)s 50000 --auto-order --coords-file ../coords
    """
)
parser.add_argument('margin', type=int, help='Margin around elements (bp)')
parser.add_argument('ref_org', nargs='?', default=None, help='Reference organism')
parser.add_argument('ref_chromo', nargs='?', default=None, help='Reference chromosome')
parser.add_argument('ref_start', nargs='?', default=None, help='Reference start position')
parser.add_argument('--coords-file', default='../coords', help='Path to coords file (default: ../coords)')
parser.add_argument('--gff-dir', default=None, help='Path to directory containing GFF files (use instead of coords-file)')
parser.add_argument('--alignments-file', default='pairwise_alignments_table', help='Path to pairwise alignments table (default: pairwise_alignments_table)')
parser.add_argument('--output-dir', default='images_draw_verbose', help='Output directory for images (default: images_draw_verbose)')
parser.add_argument('--auto-order', action='store_true', help='Auto-generate plotting_order file before drawing')
parser.add_argument('--draw-all-elements', action='store_true', help='Draw all elements instead of filtering by connectivity to reference')

args = parser.parse_args()

# Auto-generate plotting order if requested
if args.auto_order:
    print('Auto-generating plotting_order...')
    result = run('./generate_plotting_order.py', shell=True)
    if result.returncode != 0:
        print('Warning: generate_plotting_order.py failed, continuing anyway')

margin = args.margin

# Set reference element if provided
if args.ref_org and args.ref_chromo and args.ref_start:
    ref = (args.ref_org, args.ref_chromo, args.ref_start)
else:
    ref = 0

# Parse coords from either GFF directory or coords file
orgs = []
coords = {}
to_ori = set()
ori = {}

if args.gff_dir:
    import os
    gff_dir_expanded = os.path.expanduser(args.gff_dir)
    print(f'Parsing GFF files from {gff_dir_expanded}...')
    orgs, coords, to_ori = parse_gff_directory(gff_dir_expanded)
    print(f'Found {len(orgs)} organisms with {sum(len(coords[org]) for org in orgs)} chromosomes')
    # Set default reference if not specified
    if ref == 0 and orgs:
        first_org = orgs[0]
        first_chromo = list(coords[first_org].keys())[0]
        first_start = coords[first_org][first_chromo][0][0]
        ref = (first_org, first_chromo, str(first_start))
        print(f'Using default reference: {ref}')
else:
    # Parse coords file (original format)
    with open(args.coords_file) as f:
        for line in f:
            line = line.strip().split()
            if len(line) == 5:
                org,chromo,start,end,ori_gene = line
                hit_name = f'{org}_{chromo}_{start}_{end}'
            elif len(line) == 6:
                org,chromo,start,end,ori_gene,hit_name = line
            else:
                continue

            if (chromo,org,start) == ref:
                ori[ref] = ori_gene
            if org not in coords:
                coords[org] = {}
                orgs.append(org)
            if chromo not in coords[org]:
                coords[org][chromo] = []
                to_ori.add((org,chromo))
            if ori_gene == 'forward':
                strand = 1
            else:
                strand = -1
            if ref == 0:
                ref = (org,chromo,start)
            coords[org][chromo].append((int(start),int(end),strand,hit_name))

run(f'mkdir -p tables_drawing',shell=True)
tables = {}
for org in orgs:
    tables[org] = open(f'tables_drawing/{org}.tsv','w')
    tables[org].write(f'species\tchromosome\tstart\tend\tstrand\tlabel\toriginal protein name (in ncbi annotation/proteome)\n')

tot_lens = {}
no_anchors = {}
for org in orgs:
    tot_lens[org] = 0
    no_anchors[org] = 0

bib = {}

no_colors = 0
with open(args.alignments_file) as f:
    for line in f:
        if line[0] == '#':
            continue
        no_colors += 1

colors = cm.rainbow(np.linspace(0, 1, no_colors))
used = {}
orientation = {}
alignments = {}

pairs = {}

with open(args.alignments_file) as f:
    for line in f:
        if line[0] == '#':
            continue
        no,org1,chromo1,start1,end1,org2,chromo2,start2,end2,ori_alignment,score,snd1,snd2 = line.split()
        if org1 not in coords or org2 not in coords or chromo1 not in coords[org1] or chromo2 not in coords[org2]:
            continue
        start1 = int(start1)
        end1 = int(end1)
        start2 = int(start2)
        end2 = int(end2)
        if start1 > end1:
            start1,end1 = end1,start1
        if start2 > end2:
            start2,end2 = end2,start2
        relevant1 = 0
        rel1 = []
        relevant2 = 0
        rel2 = []
        for t in coords[org1][chromo1]:
            start,end,strand,hit_name = t
            if start1 > start - margin and end1 < end + margin:
                relevant1 = 1
                rel1.append((org1,chromo1,t[0]))
        for t in coords[org2][chromo2]:
            start,end,strand,hit_name = t
            if start2 > start - margin and end2 < end + margin:
                relevant2 = 1
                rel2.append((org2,chromo2,t[0]))
        if relevant1 == 0 or relevant2 == 0:
            continue
        for ele1 in rel1:
            if ele1 not in pairs:
                pairs[ele1] = set()
            for ele2 in rel2:
                if ele2 not in pairs:
                    pairs[ele2] = set()
                pairs[ele1].add(ele2)
                pairs[ele2].add(ele1)

        if (org1,chromo1,start1,end1) not in alignments:
            alignments[(org1,chromo1,start1,end1)] = []
        alignments[(org1,chromo1,start1,end1)].append([org2,chromo2,start2,end2,ori_alignment])
        if (org2,chromo2,start2,end2) not in alignments:
            alignments[(org2,chromo2,start2,end2)] = []
        alignments[(org2,chromo2,start2,end2)].append([org1,chromo1,start1,end1,ori_alignment])
        if org1 not in bib:
            bib[org1] = {}
            orientation[org1] = {}
        if chromo1 not in bib[org1]:
            bib[org1][chromo1] = []
            orientation[org1][chromo1] = {}
        if org2 not in bib:
            bib[org2] = {}
            orientation[org2] = {}
        if chromo2 not in bib[org2]:
            bib[org2][chromo2] = []
            orientation[org2][chromo2] = {}
        if org1 not in orientation[org2][chromo2]:
            orientation[org2][chromo2][org1] = {}
        if chromo1 not in orientation[org2][chromo2][org1]:
            orientation[org2][chromo2][org1][chromo1] = []
        if org2 not in orientation[org1][chromo1]:
            orientation[org1][chromo1][org2] = {}
        if chromo2 not in orientation[org1][chromo1][org2]:
            orientation[org1][chromo1][org2][chromo2] = []
        ri = randint(0,no_colors-1)
        while ri in used:
            ri = randint(0,no_colors-1)
        used[ri] = 1
        color = ri
        bib[org1][chromo1].append((start1,end1,set([color])))
        bib[org2][chromo2].append((start2,end2,set([color])))
        orientation[org1][chromo1][org2][chromo2].append(ori_alignment)
        orientation[org2][chromo2][org1][chromo1].append(ori_alignment)

# Find connected components using BFS
def get_connected_component(start_ele, pairs_graph, visited):
    """BFS to find one connected component starting from start_ele."""
    component = set()
    component.add(start_ele)
    new_eles = pairs_graph[start_ele] if start_ele in pairs_graph else set()
    visited.add(start_ele)

    while len(new_eles) > 0:
        ne = set()
        for cur_ele in new_eles:
            if cur_ele in visited:
                continue
            visited.add(cur_ele)
            component.add(cur_ele)
            if cur_ele in pairs_graph:
                for ele2 in pairs_graph[cur_ele]:
                    if ele2 not in visited:
                        ne.add(ele2)
        new_eles = ne

    return component

# Determine which elements to draw
if args.draw_all_elements:
    # Find all connected components
    print('Finding all connected components (--draw-all-elements flag set)')
    all_elements = set(pairs.keys())
    visited = set()
    components = []

    for start_ele in all_elements:
        if start_ele not in visited:
            component = get_connected_component(start_ele, pairs, visited)
            components.append(component)

    print(f'Found {len(components)} connected components')
    for i, comp in enumerate(components):
        chromo_counts = {}
        for org, chromo, start in comp:
            key = f'{org}:{chromo}'
            chromo_counts[key] = chromo_counts.get(key, 0) + 1
        print(f'  Component {i+1}: {len(comp)} elements across {len(chromo_counts)} chromosomes')
        print(f'    Chromosomes: {dict(list(chromo_counts.items())[:5])}')

    # Will draw each component as a separate image
    if len(components) == 0:
        print('No elements with alignments found')
        exit(0)
else:
    # Filter by connectivity to reference element (default behavior)
    # Start from reference, expand to all elements connected through alignment pairs
    if (ref[0],ref[1],int(ref[2])) not in pairs or len(pairs[(ref[0],ref[1],int(ref[2]))]) == 0:
        print('exiting because focal element does not have alignments')
        exit(0)
    else:
        visited = set()
        component = get_connected_component((ref[0],ref[1],int(ref[2])), pairs, visited)
        components = [component]

# Save original bib before component loop (will be filtered per component)
import copy
original_bib = copy.deepcopy(bib)

# Process each connected component separately
for comp_idx, to_draw in enumerate(components):
    if args.draw_all_elements:
            print(f'\n=== Processing Component {comp_idx+1}/{len(components)} ({len(to_draw)} elements) ===')

    # Restore original bib for this component
    bib = copy.deepcopy(original_bib)

    new_bib = {}
    tdbib = {}
    for org,chromo,start in to_draw:
        if org not in tdbib:
            tdbib[org] = {}
        if chromo not in tdbib[org]:
            tdbib[org][chromo] = []
        for startele,end,strand,hit_name in coords[org][chromo]:
            if start == startele:
                tdbib[org][chromo].append((start,end,strand))

    for org,x in bib.items():
        if org not in tdbib:
            continue
        for chromo,l in x.items():
            if chromo not in tdbib[org]:
                continue
            for ele in l:
                for gene in tdbib[org][chromo]:
                    if ele[0] >= gene[0]-margin and ele[1] <= gene[1]+margin:
                        if org not in new_bib:
                            new_bib[org] = {}
                        if chromo not in new_bib[org]:
                            new_bib[org][chromo] = []
                        new_bib[org][chromo].append(ele)

    bib = new_bib

    # Debug: show what's in bib after filtering
    if args.draw_all_elements:
        print(f'After filtering, bib contains {sum(len(v) for v in bib.values())} chromosome entries')
        for org in bib:
            print(f'  {org}: {list(bib[org].keys())}')

    # Determine chromosome orientations per connected component
    # Goal: minimize inversions across tracks in the figure
    # Pick first element in component as reference for this component
    comp_ref_ele = list(to_draw)[0]
    comp_ref = (comp_ref_ele[0], comp_ref_ele[1])  # (org, chromo)

    if comp_ref not in ori:
        ori[comp_ref] = 'forward'

    # Orient other chromosomes in this component relative to component reference
    for org,x in orientation.items():
        for chromo,y in x.items():
            if (org,chromo) in ori:
                continue
            # Check if this chromosome is in current component
            in_component = any((org, chromo, start) in to_draw for start in [t[2] for t in to_draw if t[0]==org and t[1]==chromo])
            if not in_component:
                continue

            if org not in orientation[comp_ref[0]][comp_ref[1]]:
                continue
            if chromo not in orientation[comp_ref[0]][comp_ref[1]][org]:
                continue
            oris = orientation[comp_ref[0]][comp_ref[1]][org][chromo]
            # Orient based on majority of alignment orientations to reference
            if oris.count('forward') > oris.count('reverse'):
                ori[(org,chromo)] = 'forward'
            else:
                ori[(org,chromo)] = 'reverse'

    # Orient remaining chromosomes based on already-oriented chromosomes
    # Propagate orientations through the alignment graph
    for org,x in orientation.items():
        for chromo,y in x.items():
            if (org,chromo) in ori:
                continue
            if org not in orientation[ref[0]][ref[1]]:
                continue
            if chromo not in orientation[ref[0]][ref[1]][org]:
                continue
            oris = orientation[ref[0]][ref[1]][org][chromo]
            # Orient based on majority of alignment orientations to reference
            if oris.count('forward') > oris.count('reverse'):
                ori[(org,chromo)] = 'forward'
            else:
                ori[(org,chromo)] = 'reverse'


    # Orient remaining chromosomes based on already-oriented chromosomes
    # Propagate orientations through the alignment graph
    for org,x in orientation.items():
        for chromo,y in x.items():
            if (org,chromo) in ori:
                continue
            for org2,z in y.items():
                for chromo2,oris in z.items():
                    if (org2,chromo2) not in ori:
                        continue
                    if org2 not in orientation[org][chromo]:
                        continue
                    if chromo2 not in orientation[org][chromo][org2]:
                        continue
                    oris = orientation[org][chromo][org2][chromo2]
                    # Determine orientation relative to already-oriented chromosome
                    if oris.count('forward') > oris.count('reverse'):
                        if ori[(org2,chromo2)] == 'forward':
                            ori[(org,chromo)] = 'forward'
                        else:
                            ori[(org,chromo)] = 'reverse'
                    else:
                        if ori[(org2,chromo2)] == 'forward':
                            ori[(org,chromo)] = 'reverse'
                        else:
                            ori[(org,chromo)] = 'forward'

    first = 1
    if isfile('try_to_draw_all'):
        while len(ori) != len(to_ori):
            if first == 1:
                for org,chromo in to_ori:
                    if (org,chromo) not in ori:
                        ori[(org,chromo)] = 'forward'
                        first = 0
                        break
            else:
                first = 1
                for org,x in orientation.items():
                    for chromo,y in x.items():
                        if (org,chromo) in ori:
                            continue
                        for org2,z in y.items():
                            for chromo2,oris in z.items():
                                if (org2,chromo2) not in ori:
                                    continue
                                if org2 not in orientation[org][chromo]:
                                    continue
                                if chromo2 not in orientation[org][chromo][org2]:
                                    continue
                                oris = orientation[org][chromo][org2][chromo2]
                                if oris.count('forward') > oris.count('reverse'):
                                    if ori[(org2,chromo2)] == 'forward':
                                        ori[(org,chromo)] = 'forward'
                                    else:
                                        ori[(org,chromo)] = 'reverse'
                                else:
                                    if ori[(org2,chromo2)] == 'forward':
                                        ori[(org,chromo)] = 'reverse'
                                    else:
                                        ori[(org,chromo)] = 'forward'

    for org,x in bib.items():
        for chromo,l in x.items():
            new_l = []
            for start,end,color_set in l:
                strands = set()
                alis = alignments[(org,chromo,start,end)]
                for ali in alis:
                    org2,chromo2,start2,end2,orient = ali
                    if orient == '*':
                        strands.add('reverse')
                        continue
                    if (org,chromo) not in ori:
                        continue
                    if (org2,chromo2) not in ori:
                        continue
                    ori1 = ori[(org,chromo)]
                    ori2 = ori[(org2,chromo2)]
                    if (ori1 == ori2 and ali[-1] == 'forward') or (ori1 != ori2 and ali[-1] == 'reverse'):
                        strands.add('forward')
                    else:
                        if ali[-1] == 'forward':
                            if ori1 == 'forward':
                                strands.add('forward')
                            else:
                                strands.add('reverse')
                        else:
                            strands.add('forward')
                            for enu,t in enumerate(alignments[(org2,chromo2,start2,end2)]):
                                a,b,c,d,e = t
                                if a == org and b == chromo and c == start and d == end:
                                    alignments[(org2,chromo2,start2,end2)][enu][-1] = '*'
                        
                new_l.append((start,end,color_set,strands))
            bib[org][chromo] = new_l

    # Merge overlapping alignment fragments and assign consistent colors
    # Fragments that overlap get the same color in the final figure
    color_sets = {}
    members = {}
    counter = 1

    for org,x in bib.items():
        for chromo,l in x.items():
            l = list(overlaps(l))
            for start,end,color_set,strands in l:
                # Find which color groups are already taken by overlapping fragments
                colors_taken = set()
                for color in color_set:
                    if color in members:
                        colors_taken.add(members[color])
                # Assign new color group or merge with existing groups
                if len(colors_taken) == 0:
                    # New color group
                    color_sets[counter] = color_set
                    for color in color_set:
                        members[color] = counter
                    counter += 1
                elif len(colors_taken) == 1:
                    # Merge into existing color group
                    taken = list(colors_taken)[0]
                    color_sets[taken].update(color_set)
                    for color in color_set:
                        members[color] = taken
                else:
                    # Multiple color groups need to be merged
                    color_sets[counter] = set()
                    for color in colors_taken:
                        color_sets[counter].update(color_sets[color])
                        for color2 in color_sets[color]:
                            members[color2] = counter
                        del color_sets[color]
                    color_sets[counter].update(color_set)
                    for color in color_set:
                        members[color] = counter
                    counter += 1
            bib[org][chromo] = l

    new_colors = iter(cm.rainbow(np.linspace(0, 1, len(color_sets))))
    colors = {}

    for org,x in bib.items():
        for chromo,l in x.items():
            new_l = []
            for start,end,color_set,strands in l:
                cont = 0
                for coord in coords[org][chromo]:
                    if (start >= coord[0] and start <= coord[1]) or (end >= coord[0] and end <= coord[1]) or (start <= coord[0] and end >= coord[1]) or (start >= coord[0] and end <= coord[1]):
                        cont = 1
                        break
                if cont == 1:
                    continue
                if len(set([members[color] for color in color_set])) != 1:
                    input('ERROR collecting all overlapping fragments and their colors')
                else:
                    color_no = list(set([members[color] for color in color_set]))[0]
                if color_no not in colors:
                    colors[color_no] = next(new_colors)
                color = colors[color_no]
                new_l.append((start,end,color,strands))
            bib[org][chromo] = new_l

    color_counter = 1
    used = {}
    gv = GenomeViz()
    shifts = {}
    genome_sizes = {}

    org_chromo_order = []
    if isfile('plotting_order'):
        with open('plotting_order') as f:
            for line in f:
                org,chromo = line.strip().split()
                org_chromo_order.append((org,chromo))
    else:
        for org,d in bib.items():
            for chromo in d:
                org_chromo_order.append((org,chromo))

    new_order = []
    for org,chromo in org_chromo_order:
        if org in bib and chromo in bib[org]:
            new_order.append((org,chromo))
    org_chromo_order = new_order

    # Assign default orientation 'forward' to any chromosomes that still don't have one
    for org, chromo in org_chromo_order:
        if (org, chromo) not in ori:
            ori[(org, chromo)] = 'forward'

    # Debug: show final plotting order
    if args.draw_all_elements:
        print(f'Final plotting order has {len(org_chromo_order)} tracks')
        unoriented_count = sum(1 for org, chromo in org_chromo_order if ori.get((org, chromo)) == 'forward')
        print(f'  ({unoriented_count} with default forward orientation)')

    for org,chromo in org_chromo_order:
        l = bib[org][chromo]
        if (org,chromo) not in ori or len(l) == 0:
            continue
        # first find elements and their size to make figures and get shift
        start = min(y[0] for y in l)
        end = max(y[1] for y in l)
        # i take anchors to actually check which genes are in the cluster
        if org not in coords or chromo not in coords[org]:
            continue
        start_min_ele = 1e15
        end_max_ele = 0
        for start_l,end_l,strand,hit_name in coords[org][chromo]:
            start_l,end_l = int(start_l),int(end_l)
            if (end_l < start and start-end_l > margin) or (start_l > end and start_l - end > margin):
                continue
            if start_l < start_min_ele:
                start_min_ele = start_l
            if end_l > end_max_ele:
                end_max_ele = end_l
        if start_min_ele == 1e15 and end_max_ele == 0:
            continue
        start = start_min_ele - margin
        local_margin_up = margin
        local_margin_down = margin
        if start < 0:
            local_margin_up = start_min_ele
            start = 0
        end = end_max_ele + margin
        if ori[(org,chromo)] == 'reverse':
            local_margin_up,local_margin_down = local_margin_down,local_margin_up
        shift = start
        genome_size = end - start
        name = org+' '+chromo + ' ' + ori[(org,chromo)]
        name = name.replace('_',' ')
        shifts[name] = shift

        gene_drawings = []
        for start_l,end_l,strand,hit_name in coords[org][chromo]:
            start_l,end_l = int(start_l),int(end_l)
            if end_l < start or start_l > end:
                continue
            start_l -= shift
            end_l -= shift
            if ori[(org,chromo)] == 'reverse':
                start_l,end_l = turn(start_l,end_l,[(0,genome_size)])
                if strand == -1:
                    strand = 1
                else:
                    strand = -1
            if start_l > genome_size or end_l > genome_size or start_l < 0 or end_l < 0:
                continue
            gene_drawings.append((max(0,start_l-local_margin_up),min(end_l+local_margin_down,genome_size)))

        gene_drawings = sorted(gene_drawings)
        target_ranges = []
        i = 1
        s1 = gene_drawings[0][0]
        e1 = gene_drawings[0][1]
        while i < len(gene_drawings):
            s2 = gene_drawings[i][0]
            e2 = gene_drawings[i][1]
            if s2 - e1 > 500000:
                target_ranges.append((s1,e1))
                s1 = s2
            e1 = e2
            i+=1
        if len(target_ranges) == 0:
            target_ranges = (gene_drawings[0][0],gene_drawings[-1][1])
        else:
            if gene_drawings[-1][1] >  target_ranges[-1][1]:
                if target_ranges[-1][0] != s1:
                    target_ranges.append((s1,gene_drawings[-1][1]))
                else:
                    target_ranges.append((gene_drawings[-1][0],gene_drawings[-1][1]))

        track = gv.add_feature_track(name,segments=target_ranges)
        seg_labs = []
        genome_sizes[name] = []
        for segment in track.segments:
            if ori[(org,chromo)] == 'reverse':
                seg_labs.append((segment.start,segment.end))
            else:
                segment.add_sublabel(f'{segment.start+shift}-{segment.end+shift}')
            genome_sizes[name].append((int(segment.start),int(segment.end)))
            tot_lens[org] += int(segment.end) - int(segment.start)
        if len(seg_labs) > 0:
            # For reverse tracks, convert visualization coords back to genomic coords
            seg_labels = []
            seg_labs = sorted(seg_labs)
            for s,e in seg_labs:
                # Convert from visualization space back to genomic space
                # turn() did: new_end = genome_size - start; new_start = new_end - (end-start)
                # So reverse: genomic_start = shift + (genome_size - viz_end)
                #            genomic_end = shift + (genome_size - viz_start)
                genomic_start = shift + (genome_size - e)
                genomic_end = shift + (genome_size - s)
                # Show in decreasing order for reverse strand
                seg_labels.append(f'{int(genomic_end)} - {int(genomic_start)}')
            for i,s in enumerate(track.segments):
                s.add_sublabel(seg_labels[i])
        track.set_segment_sep(symbol="//")

        for start_l,end_l,strand,hit_name in coords[org][chromo]:
            orig_start = start_l
            orig_end = end_l
            start_l,end_l = int(start_l),int(end_l)
            if end_l < start or start_l > end:
                continue
            start_l -= shift
            end_l -= shift
            if strand == -1:
                vpos = 'bottom'
                ymargin = 2.5
            else:
                vpos = 'top'
                ymargin = 0.8
            if ori[(org,chromo)] == 'reverse':
                start_l,end_l = turn(start_l,end_l,genome_sizes[name])
                if strand == -1:
                    strand = 1
                    vpos = 'top'
                    ymargin = 0.8
                else:
                    strand = -1
                    vpos = 'bottom'
                    ymargin = 2.5
            
            label = hit_name
            facecolor = 'black'
            for segment in track.segments:
                try:
                    segment.add_feature(start_l,end_l,facecolor='black',label = label, text_kws=dict(color=facecolor,vpos=vpos,ymargin=ymargin,hpos='center'),strand = strand,plotstyle='arrow')
                    tables[org].write(f'{org}\t{chromo}\t{orig_start}\t{orig_end}\t{strand}\t{label}\t{hit_name}\n')
                except:
                    pass
    
        for start_l,end_l,color,strands in l:
            if start_l > end_l:
                start_l,end_l = end_l,start_l
            cont = 0
            for coord in coords[org][chromo]:
                if (start_l >= coord[0] and start_l <= coord[1]) or (end_l >= coord[0] and end_l <= coord[1]) or (start_l <= coord[0] and end_l >= coord[1]) or (start_l >= coord[0] and end_l <= coord[1]):
                    cont = 1
                    break
            if cont == 1:
                continue
            start_l = start_l - shift
            end_l = end_l - shift
            if ori[(org,chromo)] == 'reverse':
                start_l,end_l = turn(start_l,end_l,genome_sizes[name])
            if 'reverse' in strands:
                strand_l = -1
            else:
                strand_l = 1
            if tuple(color) not in used:
                used[tuple(color)] = color_counter
                color_counter += 1
            label = ''
            strand_l = 1
            added = 0
            for segment in track.segments:
                try:
                    segment.add_feature(start_l,end_l,facecolor=color,strand = strand_l,plotstyle='arrow',label=label)
                    added = 1
                    no_anchors[org] += 1
                except:
                    pass
            if added == 0:
                print(f'CANNOT ADD anchor of {org} {chromo} {start_l} {end_l} {strands}...probably out of range: genome_size = {genome_size}')

    done = {}
    for t1,l in alignments.items():
        org,chromo,start1,end1 = t1
        if (org,chromo) not in ori:
            continue
        orig_start1 = start1
        orig_end1 = end1
        orig_start1 = int(orig_start1)
        orig_end1 = int(orig_end1)
        cont = 0
        for coord in coords[org][chromo]:
            if (start1 >= coord[0] and start1 <= coord[1]) or (end1 >= coord[0] and end1 <= coord[1]) or (start1 <= coord[0] and end1 >= coord[1]) or (start1 >= coord[0] and end1 <= coord[1]):
                cont = 1
                break
        if cont == 1:
            continue
        name1 = org+' '+chromo+' '+ori[(org,chromo)]
        name1 = name1.replace('_',' ')
        if name1 not in shifts:
            continue
        shift1 = shifts[name1]
        start1 -= shift1
        end1 -= shift1
        if ori[(org,chromo)] == 'reverse':
            start1,end1 = turn(start1,end1,genome_sizes[name1])
        for org2,chromo2,start2,end2,ori_alignment in l:
            if (org2,chromo2) not in ori:
                continue
            orig_start2 = start2
            orig_end2 = end2
            orig_start2 = int(orig_start2)
            orig_end2 = int(orig_end2)
            if ((org,chromo,orig_start1,orig_end1),(org2,chromo2,orig_start2,orig_end2)) in done or ((org2,chromo2,orig_start2,orig_end2),(org,chromo,orig_start1,orig_end1)) in done:
                continue
            else:
                done[((org,chromo,orig_start1,orig_end1),(org2,chromo2,orig_start2,orig_end2))] = 1
                done[((org2,chromo2,orig_start2,orig_end2),(org,chromo,orig_start1,orig_end1))] = 1
            cont = 0
            for coord in coords[org2][chromo2]:
                if (start2 >= coord[0] and start2 <= coord[1]) or (end2 >= coord[0] and end2 <= coord[1]) or (start2 <= coord[0] and end2 >= coord[1]) or (start2 >= coord[0] and end2 <= coord[1]):
                    cont = 1 
                    break
            if cont == 1:
                continue
            name2 = org2+' '+chromo2+ ' ' + ori[(org2,chromo2)]
            name2 = name2.replace('_',' ')
            if name2 not in shifts:
                continue
            shift2 = shifts[name2]
            start2 -= shift2
            end2 -= shift2
            if ori_alignment == 'reverse' and ori[(org,chromo)] == ori[(org2,chromo2)] or ori_alignment == 'forward' and ori[(org,chromo)] != ori[(org2,chromo2)]:
                start2,end2 = end2,start2
            if ori[(org2,chromo2)] == 'reverse':
                start2,end2 = turn(start2,end2,genome_sizes[name2])

            try:
                for track in gv.feature_tracks:
                    if track.name != name1:
                        continue
                    for segment in track.segments:
                        segstart = int(segment.start)
                        segend = int(segment.end)
                        if start1 < segstart or end1 > segend:
                            continue
                        else:
                            seg1 = segment.name
                for track in gv.feature_tracks:
                    if track.name != name2:
                        continue
                    for segment in track.segments:
                        segstart = int(segment.start)
                        segend = int(segment.end)
                        if start2 < segstart or end2 > segend:
                            continue
                        else:
                            seg2 = segment.name
                gv.add_link((name1,seg1,start1,end1),(name2,seg2,start2,end2), color="skyblue", inverted_color="lime", curve=True)
            except:
                print(f'link FROM {name1} {start1} {end1} TO {name2} {start2} {end2} not adjacent')

    run(f'mkdir -p {args.output_dir} && mkdir -p {args.output_dir}/{margin}',shell=True)
    fig = gv.plotfig()
    if args.draw_all_elements:
        # Find a representative element for naming
        first_ele = list(to_draw)[0]
        output_filename = f'{args.output_dir}/{margin}/component_{comp_idx+1}_of_{len(components)}_ref_{first_ele[0]}_{first_ele[1]}_{first_ele[2]}.svg'
    else:
        output_filename = f'{args.output_dir}/{margin}/{ref}_verbose_with_alignments_with_prot_labels.svg'
    fig.savefig(output_filename)
    print(f'Saved figure to {output_filename}')
    plt.close(fig)

    for org in orgs:
        tables[org].close()
