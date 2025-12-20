#!/usr/bin/env python3

"""
Draw SynOrthogroup synteny visualizations with pygenomeviz.

This script creates synteny plots for each SynOrthogroup, showing genes
across species with alignment links. It reads SynOrthogroups.tsv output
from export_synorthogroups.py, looks up coordinates, and draws figures.

Usage:
    ./draw_synorthogroups.py MARGIN SYNORTHO_DIR [OPTIONS]

    # Basic usage
    ./draw_synorthogroups.py 50000 synorthogroups/no_new_edges

    # With custom paths
    ./draw_synorthogroups.py 50000 synorthogroups/no_new_edges --parsed-gff parsed_faa_gff

    # Draw specific SOGs only
    ./draw_synorthogroups.py 50000 synorthogroups/no_new_edges --sogs SOG0000001,SOG0000003

Arguments:
    MARGIN: Integer margin around genes (required)
    SYNORTHO_DIR: Path to synorthogroups directory (e.g., synorthogroups/no_new_edges)

Options:
    --parsed-gff: Path to parsed_faa_gff file (default: ../parsed_faa_gff)
    --alignments-file: Path to pairwise alignments table (default: ../pairwise_alignments_table)
    --output-dir: Output directory for images (default: images_synorthogroups)
    --sogs: Comma-separated list of SOG IDs to draw (default: all)
    --min-species: Minimum species count for SOG to be drawn (default: 2)

Reads from:
    - SynOrthogroups.tsv: orthogroup assignments
    - parsed_faa_gff: coordinates of all genes
    - pairwise_alignments_table: alignment data

Writes to:
    - images_synorthogroups/{margin}/SOG*_*.svg: one figure per SOG
"""

from sys import argv
import argparse
import os
from matplotlib.pyplot import cm
from pygenomeviz import GenomeViz
from random import randint
import numpy as np
from subprocess import run
import pickle
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from os.path import isfile
from collections import defaultdict

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

def parse_synorthogroups_tsv(tsv_path):
    """
    Parse SynOrthogroups.tsv file.

    Returns:
        sogs: dict {sog_id: {species: [gene_ids]}}
        species_list: list of species IDs
    """
    sogs = {}
    species_list = []

    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        species_list = header[1:]  # Skip 'SynOrthogroup' column

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue

            sog_id = fields[0]
            sogs[sog_id] = {}

            for i, gene_str in enumerate(fields[1:]):
                if gene_str.strip():
                    species = species_list[i]
                    # Genes are comma-separated
                    genes = [g.strip() for g in gene_str.split(',')]
                    sogs[sog_id][species] = genes

    return sogs, species_list

def load_gene_coordinates_from_tsv(withcoords_tsv_path):
    """
    Load gene coordinates from SynOrthogroups_WithCoords.tsv.

    Returns:
        coords: dict {species: {chromosome: [(start, end, strand, gene_id)]}}
        gene_to_coord: dict {gene_id: (species, chromosome, start, end, strand)}
    """
    coords = {}
    gene_to_coord = {}

    with open(withcoords_tsv_path) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            sog_id, gene_id, species, chrom, start, end, strand, ref_gene = parts

            start, end = int(start), int(end)
            strand_int = 1 if strand == '+' else -1

            gene_to_coord[gene_id] = (species, chrom, start, end, strand_int)

            if species not in coords:
                coords[species] = {}
            if chrom not in coords[species]:
                coords[species][chrom] = []
            coords[species][chrom].append((start, end, strand_int, gene_id))

    return coords, gene_to_coord

parser = argparse.ArgumentParser(
    description='Draw SynOrthogroup synteny visualizations with pygenomeviz',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog="""
Examples:
  %(prog)s 50000 synorthogroups/no_new_edges
  %(prog)s 50000 synorthogroups/no_new_edges --sogs SOG0000001,SOG0000003
  %(prog)s 50000 synorthogroups/with_new_edges --min-species 5
    """
)
parser.add_argument('margin', type=int, help='Margin around genes (bp)')
parser.add_argument('synortho_dir', help='Path to synorthogroups directory')
parser.add_argument('--parsed-gff', default='../parsed_faa_gff', help='Path to parsed_faa_gff file (default: ../parsed_faa_gff)')
parser.add_argument('--alignments-file', default='../../utils/pairwise_alignments_table', help='Path to pairwise alignments table (default: ../../utils/pairwise_alignments_table)')
parser.add_argument('--output-dir', default='images_synorthogroups', help='Output directory for images (default: images_synorthogroups)')
parser.add_argument('--sogs', default=None, help='Comma-separated list of SOG IDs to draw (default: all)')
parser.add_argument('--min-species', type=int, default=2, help='Minimum species count for SOG to be drawn (default: 2)')

args = parser.parse_args()

margin = args.margin
synortho_tsv = f'{args.synortho_dir}/SynOrthogroups.tsv'

graph_version = os.path.basename(args.synortho_dir.rstrip('/'))
if graph_version not in ['no_new_edges', 'with_new_edges']:
    graph_version = 'unknown'

print(f'Loading SynOrthogroups from {synortho_tsv}...')
sogs, species_list = parse_synorthogroups_tsv(synortho_tsv)
print(f'  Found {len(sogs)} SynOrthogroups across {len(species_list)} species')

withcoords_tsv = f'{args.synortho_dir}/SynOrthogroups_WithCoords.tsv'
print(f'Loading gene coordinates from {withcoords_tsv}...')
coords, gene_to_coord = load_gene_coordinates_from_tsv(withcoords_tsv)
print(f'  Loaded coordinates for {len(gene_to_coord)} genes')

# Filter SOGs if specific list provided
if args.sogs:
    selected_sogs = set(args.sogs.split(','))
    sogs = {k: v for k, v in sogs.items() if k in selected_sogs}
    print(f'  Filtering to {len(sogs)} requested SOGs')

# Filter by minimum species count
sogs = {k: v for k, v in sogs.items() if len(v) >= args.min_species}
print(f'  After filtering (min {args.min_species} species): {len(sogs)} SOGs to draw')

# Load alignments
ori = {}
orientation = {}
alignments = {}
bib = {}

# Check if alignment file exists
has_alignments = os.path.exists(args.alignments_file)

if not has_alignments:
    print(f'Warning: alignments file {args.alignments_file} not found')
    print('  Figures will show genes only, without alignment links')

# Count alignments
no_colors = 0
if has_alignments:
    with open(args.alignments_file) as f:
        for line in f:
            if line[0] == '#':
                continue
            no_colors += 1

colors = cm.rainbow(np.linspace(0, 1, max(1, no_colors)))
used = {}

# Parse alignment file
if has_alignments:
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

            # Store alignments
            if (org1,chromo1,start1,end1) not in alignments:
                alignments[(org1,chromo1,start1,end1)] = []
            alignments[(org1,chromo1,start1,end1)].append([org2,chromo2,start2,end2,ori_alignment])

            if (org2,chromo2,start2,end2) not in alignments:
                alignments[(org2,chromo2,start2,end2)] = []
            alignments[(org2,chromo2,start2,end2)].append([org1,chromo1,start1,end1,ori_alignment])

            # Store orientation information
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

if has_alignments:
    print(f'Loaded {len(alignments)} alignment fragments')
else:
    print('Skipping alignments (file not found)')

# Process each SOG
output_path = f'{args.output_dir}/{graph_version}/{margin}'
run(f'mkdir -p {output_path}',shell=True)

for sog_idx, (sog_id, sog_genes) in enumerate(sogs.items()):
    print(f'\n=== Drawing {sog_id} ({sog_idx+1}/{len(sogs)}) ===')
    print(f'  Species: {len(sog_genes)}, Genes: {sum(len(genes) for genes in sog_genes.values())}')

    # Build coords structure for genes in this SOG
    sog_coords = {}
    sog_gene_set = set()

    for species, gene_list in sog_genes.items():
        if species not in coords:
            continue
        sog_coords[species] = {}

        for gene_id in gene_list:
            if gene_id not in gene_to_coord:
                continue

            sp, chrom, start, end, strand = gene_to_coord[gene_id]

            if chrom not in sog_coords[species]:
                sog_coords[species][chrom] = []

            sog_coords[species][chrom].append((start, end, strand, gene_id))
            sog_gene_set.add(gene_id)

    if len(sog_coords) == 0:
        print(f'  Skipping {sog_id}: no coordinates found')
        continue

    # Filter bib to only include alignments near genes in this SOG
    import copy
    sog_bib = {}

    for org, chroms in bib.items():
        if org not in sog_coords:
            continue

        for chromo, alignment_list in chroms.items():
            if chromo not in sog_coords[org]:
                continue

            for align_start, align_end, color_set in alignment_list:
                # Check if alignment is near any gene in this SOG
                for gene_start, gene_end, strand, gene_id in sog_coords[org][chromo]:
                    if align_start >= gene_start - margin and align_end <= gene_end + margin:
                        if org not in sog_bib:
                            sog_bib[org] = {}
                        if chromo not in sog_bib[org]:
                            sog_bib[org][chromo] = []
                        sog_bib[org][chromo].append((align_start, align_end, color_set))
                        break

    # Determine chromosome orientations
    # Pick first chromosome as reference
    ref_org = None
    ref_chromo = None
    for org in sog_coords:
        for chromo in sog_coords[org]:
            ref_org = org
            ref_chromo = chromo
            break
        if ref_org:
            break

    if not ref_org:
        print(f'  Skipping {sog_id}: no reference found')
        continue

    sog_ori = {}
    sog_ori[(ref_org, ref_chromo)] = 'forward'

    # Orient other chromosomes relative to reference
    for org, chroms in orientation.items():
        for chromo, targets in chroms.items():
            if (org, chromo) in sog_ori:
                continue

            # Check if this chromosome is in the SOG
            if org not in sog_coords or chromo not in sog_coords[org]:
                continue

            if ref_org not in targets or ref_chromo not in targets.get(ref_org, {}):
                continue

            oris = targets[ref_org][ref_chromo]
            if oris.count('forward') > oris.count('reverse'):
                sog_ori[(org, chromo)] = 'forward'
            else:
                sog_ori[(org, chromo)] = 'reverse'

    # Add strands to bib
    for org, chroms in sog_bib.items():
        for chromo, l in chroms.items():
            new_l = []
            for start,end,color_set in l:
                strands = set()
                if (org,chromo,start,end) in alignments:
                    alis = alignments[(org,chromo,start,end)]
                    for ali in alis:
                        org2,chromo2,start2,end2,orient = ali
                        if orient == '*':
                            strands.add('reverse')
                            continue
                        if (org,chromo) not in sog_ori:
                            continue
                        if (org2,chromo2) not in sog_ori:
                            continue
                        ori1 = sog_ori[(org,chromo)]
                        ori2 = sog_ori[(org2,chromo2)]
                        if (ori1 == ori2 and orient == 'forward') or (ori1 != ori2 and orient == 'reverse'):
                            strands.add('forward')
                        else:
                            strands.add('reverse')

                new_l.append((start,end,color_set,strands))
            sog_bib[org][chromo] = new_l

    # Merge overlapping alignment fragments
    color_sets = {}
    members = {}
    counter = 1

    for org, chroms in sog_bib.items():
        for chromo, l in chroms.items():
            l = list(overlaps(l))
            for start,end,color_set,strands in l:
                colors_taken = set()
                for color in color_set:
                    if color in members:
                        colors_taken.add(members[color])

                if len(colors_taken) == 0:
                    color_sets[counter] = color_set
                    for color in color_set:
                        members[color] = counter
                    counter += 1
                elif len(colors_taken) == 1:
                    taken = list(colors_taken)[0]
                    color_sets[taken].update(color_set)
                    for color in color_set:
                        members[color] = taken
                else:
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
            sog_bib[org][chromo] = l

    new_colors = iter(cm.rainbow(np.linspace(0, 1, len(color_sets) if len(color_sets) > 0 else 1)))
    final_colors = {}

    for org, chroms in sog_bib.items():
        for chromo, l in chroms.items():
            new_l = []
            for start,end,color_set,strands in l:
                # Skip if overlaps with gene itself
                cont = 0
                for coord in sog_coords[org][chromo]:
                    if (start >= coord[0] and start <= coord[1]) or (end >= coord[0] and end <= coord[1]) or (start <= coord[0] and end >= coord[1]) or (start >= coord[0] and end <= coord[1]):
                        cont = 1
                        break
                if cont == 1:
                    continue

                if len(color_set) > 0:
                    color_no = list(set([members[color] for color in color_set]))[0]
                    if color_no not in final_colors:
                        final_colors[color_no] = next(new_colors)
                    color = final_colors[color_no]
                else:
                    color = 'gray'

                new_l.append((start,end,color,strands))
            sog_bib[org][chromo] = new_l

    # Create figure
    gv = GenomeViz()
    shifts = {}
    genome_sizes = {}

    # Determine track order (use species_list order)
    org_chromo_order = []
    for species in species_list:
        if species in sog_coords:
            for chromo in sorted(sog_coords[species].keys()):
                org_chromo_order.append((species, chromo))

    # Assign default orientations
    for org, chromo in org_chromo_order:
        if (org, chromo) not in sog_ori:
            sog_ori[(org, chromo)] = 'forward'

    # Create tracks
    for org, chromo in org_chromo_order:
        if org not in sog_coords or chromo not in sog_coords[org]:
            continue

        # Follow draw_verbose.py exactly from line 755 onwards
        start = min(int(g[0]) for g in sog_coords[org][chromo]) - margin
        end = max(int(g[1]) for g in sog_coords[org][chromo]) + margin
        if start < 0:
            start = 0

        shift = start
        genome_size = end - start
        name = org+' '+chromo + ' ' + sog_ori[(org,chromo)]
        name = name.replace('_',' ')
        shifts[name] = shift

        gene_drawings = []
        for start_l,end_l,strand,hit_name in sog_coords[org][chromo]:
            start_l,end_l = int(start_l),int(end_l)
            if end_l < start or start_l > end:
                continue
            start_l -= shift
            end_l -= shift
            if sog_ori[(org,chromo)] == 'reverse':
                start_l,end_l = turn(start_l,end_l,[(0,genome_size)])
                if strand == -1:
                    strand = 1
                else:
                    strand = -1
            if start_l > genome_size or end_l > genome_size or start_l < 0 or end_l < 0:
                continue
            gene_drawings.append((max(0,start_l-margin),min(end_l+margin,genome_size)))

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
            if sog_ori[(org,chromo)] == 'reverse':
                seg_labs.append((segment.start,segment.end))
            else:
                segment.add_sublabel(f'{segment.start+shift}-{segment.end+shift}')
            genome_sizes[name].append((int(segment.start),int(segment.end)))

        # Add genes (like draw_verbose lines 832-864)
        for start_l,end_l,strand,hit_name in sog_coords[org][chromo]:
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
            if sog_ori[(org,chromo)] == 'reverse':
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
                except:
                    pass

        # Add alignment fragments (like draw_verbose lines 866-897)
        if org in sog_bib and chromo in sog_bib[org]:
            for start_l,end_l,color,strands in sog_bib[org][chromo]:
                if start_l > end_l:
                    start_l,end_l = end_l,start_l
                cont = 0
                for coord in sog_coords[org][chromo]:
                    if (start_l >= coord[0] and start_l <= coord[1]) or (end_l >= coord[0] and end_l <= coord[1]) or (start_l <= coord[0] and end_l >= coord[1]) or (start_l >= coord[0] and end_l <= coord[1]):
                        cont = 1
                        break
                if cont == 1:
                    continue
                start_l = start_l - shift
                end_l = end_l - shift
                if sog_ori[(org,chromo)] == 'reverse':
                    start_l,end_l = turn(start_l,end_l,genome_sizes[name])
                if 'reverse' in strands:
                    strand_l = -1
                else:
                    strand_l = 1
                label = ''
                strand_l = 1
                for segment in track.segments:
                    try:
                        segment.add_feature(start_l,end_l,facecolor=color,strand = strand_l,plotstyle='arrow',label=label)
                    except:
                        pass

    # Add links between tracks
    done = {}
    for key, ali_list in alignments.items():
        org, chromo, start1, end1 = key
        if (org, chromo) not in sog_ori:
            continue

        # Skip if overlaps with gene
        cont = 0
        if org in sog_coords and chromo in sog_coords[org]:
            for coord in sog_coords[org][chromo]:
                if (start1 >= coord[0] and start1 <= coord[1]) or (end1 >= coord[0] and end1 <= coord[1]):
                    cont = 1
                    break
        if cont == 1:
            continue

        name1 = f'{org} {chromo} {sog_ori[(org,chromo)]}'.replace('_', ' ')
        if name1 not in shifts:
            continue

        shift1 = shifts[name1]
        plot_start1 = start1 - shift1
        plot_end1 = end1 - shift1

        if sog_ori[(org, chromo)] == 'reverse':
            plot_start1, plot_end1 = turn(plot_start1, plot_end1, genome_sizes[name1])

        for org2, chromo2, start2, end2, ori_alignment in ali_list:
            if (org2, chromo2) not in sog_ori:
                continue

            if ((org,chromo,start1,end1),(org2,chromo2,start2,end2)) in done:
                continue
            else:
                done[((org,chromo,start1,end1),(org2,chromo2,start2,end2))] = 1
                done[((org2,chromo2,start2,end2),(org,chromo,start1,end1))] = 1

            # Skip if overlaps with gene
            cont = 0
            if org2 in sog_coords and chromo2 in sog_coords[org2]:
                for coord in sog_coords[org2][chromo2]:
                    if (start2 >= coord[0] and start2 <= coord[1]) or (end2 >= coord[0] and end2 <= coord[1]):
                        cont = 1
                        break
            if cont == 1:
                continue

            name2 = f'{org2} {chromo2} {sog_ori[(org2,chromo2)]}'.replace('_', ' ')
            if name2 not in shifts:
                continue

            shift2 = shifts[name2]
            plot_start2 = start2 - shift2
            plot_end2 = end2 - shift2

            if ori_alignment == 'reverse' and sog_ori[(org,chromo)] == sog_ori[(org2,chromo2)] or ori_alignment == 'forward' and sog_ori[(org,chromo)] != sog_ori[(org2,chromo2)]:
                plot_start2, plot_end2 = plot_end2, plot_start2

            if sog_ori[(org2, chromo2)] == 'reverse':
                plot_start2, plot_end2 = turn(plot_start2, plot_end2, genome_sizes[name2])

            try:
                seg1 = None
                for track in gv.feature_tracks:
                    if track.name != name1:
                        continue
                    for segment in track.segments:
                        segstart = int(segment.start)
                        segend = int(segment.end)
                        if plot_start1 < segstart or plot_end1 > segend:
                            continue
                        else:
                            seg1 = segment.name

                seg2 = None
                for track in gv.feature_tracks:
                    if track.name != name2:
                        continue
                    for segment in track.segments:
                        segstart = int(segment.start)
                        segend = int(segment.end)
                        if plot_start2 < segstart or plot_end2 > segend:
                            continue
                        else:
                            seg2 = segment.name

                if seg1 and seg2:
                    gv.add_link(
                        (name1, seg1, plot_start1, plot_end1),
                        (name2, seg2, plot_start2, plot_end2),
                        color="skyblue",
                        inverted_color="lime",
                        curve=True
                    )
            except:
                pass

    # Save figure
    fig = gv.plotfig()
    output_filename = f'{output_path}/{sog_id}_{len(sog_genes)}species_{sum(len(g) for g in sog_genes.values())}genes.svg'
    fig.savefig(output_filename)
    print(f'  Saved to {output_filename}')
    plt.close(fig)

print(f'\n=== Done! Figures saved to {output_path}/ ===')
