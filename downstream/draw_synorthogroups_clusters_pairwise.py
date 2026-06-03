#!/usr/bin/env python3

"""
Draw pairwise SynOrthogroup cluster comparisons.

This script clusters genes from synorthogroups by spatial proximity within
species, then draws ALL pairwise cluster comparisons between species where
clusters share sufficient SOGs. Each pairwise comparison is drawn as a
separate image with exactly 2 tracks.

Usage:
    ./draw_synorthogroups_clusters_pairwise.py --sog_file SynOrthogroups.tsv \\
        --coords_file SynOrthogroups_WithCoords.tsv \\
        --alignments_file pairwise_alignments_table \\
        --cluster_margin 100000 \\
        --min_shared_sogs 5 \\
        --min_cluster_size 10 \\
        --cores 4 \\
        --output_dir pairwise_clusters

Arguments:
    --sog_file: Path to SynOrthogroups.tsv (SOG assignments)
    --coords_file: Path to SynOrthogroups_WithCoords.tsv (gene coordinates)
    --alignments_file: Path to pairwise_alignments_table (for drawing links)
    --cluster_margin: Max distance (bp) between genes in same cluster (default: 100000)
    --min_shared_sogs: Min shared SOGs to draw cluster pair (default: 1)
    --min_cluster_size: Min genes per cluster to include (default: 1)
    --gene_margin: Margin around genes for visualization (default: 10000)
    --gap_threshold: Min gap size to create segment break (default: 200000)
    --pairs_file: Optional file with species pairs (format: species1<tab>species2)
    --cores: Number of parallel workers (default: 1)
    --output_dir: Output directory for images (default: pairwise_clusters)
    --show_labels: Display gene IDs on arrows (default: off)
"""

import argparse
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Set
from pygenomeviz import GenomeViz
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np


@dataclass
class RejectionRecord:
    """Record why a feature was not drawn"""
    feature_id: str
    feature_type: str  # 'gene' or 'alignment'
    reason: str
    details: str
    org: str = ""
    chromosome: str = ""
    coordinates: Tuple = ()
    extra_info: Dict = None


class FeatureRejectionLog:
    """Track all rejected features with detailed reasons"""

    def __init__(self):
        self.rejections: Dict[str, List[RejectionRecord]] = defaultdict(list)
        self.total_genes_processed = 0
        self.genes_drawn = 0
        self.total_alignments_attempted = 0
        self.alignments_drawn = 0

    def add_gene(self, gene_id: str, reason: str, details: str,
                 org: str = "", chrom: str = "", coords: Tuple = (), **extra):
        """Log a rejected gene"""
        record = RejectionRecord(
            feature_id=gene_id,
            feature_type='gene',
            reason=reason,
            details=details,
            org=org,
            chromosome=chrom,
            coordinates=coords,
            extra_info=extra if extra else None
        )
        self.rejections[reason].append(record)

    def add_alignment(self, link_id: str, reason: str, details: str, **extra):
        """Log a rejected alignment link"""
        record = RejectionRecord(
            feature_id=link_id,
            feature_type='alignment',
            reason=reason,
            details=details,
            extra_info=extra if extra else None
        )
        self.rejections[reason].append(record)

    def summary(self) -> str:
        """Generate summary report"""
        lines = []
        lines.append("\n" + "="*80)
        lines.append("FEATURE REJECTION SUMMARY")
        lines.append("="*80)
        lines.append(f"Genes: {self.genes_drawn}/{self.total_genes_processed} drawn "
                    f"({100*self.genes_drawn/max(1,self.total_genes_processed):.1f}%)")
        lines.append(f"Alignments: {self.alignments_drawn}/{self.total_alignments_attempted} drawn "
                    f"({100*self.alignments_drawn/max(1,self.total_alignments_attempted):.1f}%)")
        lines.append("")

        # Group by reason
        gene_reasons = [r for r in self.rejections if any(rec.feature_type == 'gene'
                       for rec in self.rejections[r])]
        alignment_reasons = [r for r in self.rejections if any(rec.feature_type == 'alignment'
                            for rec in self.rejections[r])]

        if gene_reasons:
            lines.append("GENE REJECTIONS:")
            for reason in sorted(gene_reasons):
                records = [r for r in self.rejections[reason] if r.feature_type == 'gene']
                lines.append(f"  {reason}: {len(records)} genes")

        if alignment_reasons:
            lines.append("\nALIGNMENT REJECTIONS:")
            for reason in sorted(alignment_reasons):
                records = [r for r in self.rejections[reason] if r.feature_type == 'alignment']
                lines.append(f"  {reason}: {len(records)} links")

        lines.append("="*80 + "\n")
        return "\n".join(lines)

    def write_tsv(self, filepath: str):
        """Write detailed rejection log"""
        with open(filepath, 'w') as f:
            f.write("feature_id\tfeature_type\treason\tdetails\torg\tchromosome\tcoordinates\n")
            for reason, records in self.rejections.items():
                for rec in records:
                    coords_str = f"{rec.coordinates}" if rec.coordinates else "N/A"
                    f.write(f"{rec.feature_id}\t{rec.feature_type}\t{reason}\t"
                           f"{rec.details}\t{rec.org}\t{rec.chromosome}\t{coords_str}\n")


def turn(start, end, segments):
    """
    Flip coordinates for reverse orientation.

    Always flips within the FULL track range (min to max of all segments),
    regardless of segment boundaries. This ensures coordinates remain
    track-relative (as required by pygenomeviz).

    Args:
        start: Start coordinate (track-relative)
        end: End coordinate (track-relative)
        segments: List of (start, end) tuples defining track segments

    Returns:
        Tuple of (new_start, new_end) in track-relative coordinates
    """
    # Calculate full track range
    max_e = max(seg[1] for seg in segments)
    min_s = min(seg[0] for seg in segments)
    len_line = max_e - min_s

    # Flip within full range
    new_end = len_line - start
    new_start = new_end - (end - start)

    return new_start, new_end


def parse_synorthogroups_tsv(tsv_path):
    """Parse SynOrthogroups.tsv file."""
    sogs = {}
    gene_to_sog = {}

    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        species_list = header[1:]

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue

            sog_id = fields[0]
            sogs[sog_id] = {}

            for i, gene_str in enumerate(fields[1:]):
                if gene_str.strip():
                    species = species_list[i]
                    genes = [g.strip() for g in gene_str.split(',')]
                    sogs[sog_id][species] = genes

                    for gene in genes:
                        gene_to_sog[gene] = sog_id

    return sogs, gene_to_sog


def load_gene_coordinates(coords_tsv_path):
    """Load gene coordinates from SynOrthogroups_WithCoords.tsv."""
    coords = {}
    gene_to_coord = {}

    with open(coords_tsv_path) as f:
        header = f.readline()

        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

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


def load_alignments(alignments_file):
    """Load pairwise alignments from file."""
    alignments = defaultdict(list)
    orientation = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

    if not os.path.exists(alignments_file):
        print(f'Warning: alignments file {alignments_file} not found')
        return alignments, orientation

    with open(alignments_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split()
            if len(fields) < 10:
                continue

            org1, chr1 = fields[1], fields[2]
            start1, end1 = int(fields[3]), int(fields[4])
            org2, chr2 = fields[5], fields[6]
            start2, end2 = int(fields[7]), int(fields[8])
            ori = fields[9]

            if start1 > end1:
                start1, end1 = end1, start1
            if start2 > end2:
                start2, end2 = end2, start2

            alignments[(org1, chr1)].append((start1, end1, org2, chr2, start2, end2, ori))
            alignments[(org2, chr2)].append((start2, end2, org1, chr1, start1, end1, ori))

            orientation[org1][chr1][org2][chr2].append(ori)
            reverse_ori = 'reverse' if ori == 'forward' else 'forward'
            orientation[org2][chr2][org1][chr1].append(reverse_ori)

    return alignments, orientation


def cluster_genes_on_chromosome(genes, margin):
    """Cluster genes by proximity on a single chromosome."""
    if not genes:
        return []

    sorted_genes = sorted(genes, key=lambda g: g[0])
    clusters = []
    current_cluster = [sorted_genes[0]]

    for gene in sorted_genes[1:]:
        prev_end = current_cluster[-1][1]
        curr_start = gene[0]

        if curr_start - prev_end <= margin:
            current_cluster.append(gene)
        else:
            clusters.append(current_cluster)
            current_cluster = [gene]

    clusters.append(current_cluster)
    return clusters


def cluster_within_species(coords, margin):
    """Cluster genes by proximity within each species."""
    clusters = {}
    cluster_ids = []

    for species in coords:
        clusters[species] = {}
        for chrom in coords[species]:
            genes = coords[species][chrom]
            chrom_clusters = cluster_genes_on_chromosome(genes, margin)
            clusters[species][chrom] = chrom_clusters

            for idx in range(len(chrom_clusters)):
                cluster_id = f"{species}_{chrom}_{idx}"
                cluster_ids.append(cluster_id)

    return clusters, cluster_ids


def filter_clusters_by_size(clusters: Dict, min_cluster_size: int) -> Dict:
    """Filter clusters to keep only those with >= min_cluster_size genes."""
    if min_cluster_size <= 1:
        return clusters

    filtered = {}
    removed_count = 0
    kept_count = 0

    for species in clusters:
        filtered[species] = {}
        for chrom in clusters[species]:
            filtered[species][chrom] = []
            for cluster in clusters[species][chrom]:
                if len(cluster) >= min_cluster_size:
                    filtered[species][chrom].append(cluster)
                    kept_count += 1
                else:
                    removed_count += 1

    print(f"  Cluster size filter: kept {kept_count}, removed {removed_count}")
    return filtered


def load_species_pairs(pairs_file: str) -> List[Tuple[str, str]]:
    """Load species pairs from file."""
    pairs = []
    with open(pairs_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) != 2:
                print(f"Warning: Line {line_num} has {len(parts)} fields, expected 2. Skipping.")
                continue

            species1, species2 = parts[0].strip(), parts[1].strip()
            pairs.append((species1, species2))

    print(f"  Loaded {len(pairs)} species pairs from {pairs_file}")
    return pairs


def generate_all_pairs(species_list: List[str]) -> List[Tuple[str, str]]:
    """Generate all unique species pairs."""
    pairs = []
    for i, sp1 in enumerate(species_list):
        for sp2 in species_list[i+1:]:
            pairs.append((sp1, sp2))

    print(f"  Generated {len(pairs)} species pairs from {len(species_list)} species")
    return pairs


def find_cluster_pairs(species1: str, species2: str,
                      clusters: Dict, gene_to_sog: Dict,
                      min_shared_sogs: int) -> List[Dict]:
    """Find all cluster pairs between two species that share >= min_shared_sogs."""
    cluster_pairs = []

    # Build SOG sets for species1 clusters
    sp1_cluster_sogs = {}
    for chrom1 in clusters.get(species1, {}):
        for idx1, genes1 in enumerate(clusters[species1][chrom1]):
            cluster_id1 = f"{species1}_{chrom1}_{idx1}"
            sogs1 = set()
            for gene in genes1:
                gene_id = gene[3]
                if gene_id in gene_to_sog:
                    sogs1.add(gene_to_sog[gene_id])
            sp1_cluster_sogs[cluster_id1] = {
                'chrom': chrom1,
                'idx': idx1,
                'genes': genes1,
                'sogs': sogs1
            }

    # Compare with species2 clusters
    for chrom2 in clusters.get(species2, {}):
        for idx2, genes2 in enumerate(clusters[species2][chrom2]):
            sogs2 = set()
            for gene in genes2:
                gene_id = gene[3]
                if gene_id in gene_to_sog:
                    sogs2.add(gene_to_sog[gene_id])

            # Check all species1 clusters for shared SOGs
            for cluster_id1, cluster1_data in sp1_cluster_sogs.items():
                sogs1 = cluster1_data['sogs']
                shared = sogs1 & sogs2

                if len(shared) >= min_shared_sogs:
                    cluster_pairs.append({
                        'species1': species1,
                        'chrom1': cluster1_data['chrom'],
                        'cluster1_idx': cluster1_data['idx'],
                        'genes1': cluster1_data['genes'],
                        'sogs1': sogs1,
                        'species2': species2,
                        'chrom2': chrom2,
                        'cluster2_idx': idx2,
                        'genes2': genes2,
                        'sogs2': sogs2,
                        'shared_sogs': shared,
                        'n_shared': len(shared)
                    })

    return cluster_pairs


def get_colormap(n_colors):
    """Get appropriate colormap based on number of colors needed."""
    if n_colors <= 10:
        return plt.cm.tab10
    elif n_colors <= 20:
        return plt.cm.tab20
    elif n_colors <= 40:
        colors = []
        for cmap_name in ['tab20', 'tab20b', 'tab20c']:
            cmap = plt.cm.get_cmap(cmap_name)
            for i in range(20):
                colors.append(cmap(i / 20))
        return ListedColormap(colors[:n_colors])
    else:
        return plt.cm.hsv


def assign_colors_global(sogs: Dict, gene_to_sog: Dict) -> Dict:
    """Assign colors to ALL SOGs in the dataset for consistency across images."""
    all_sog_ids = set(sogs.keys())
    n_sogs = len(all_sog_ids)
    print(f"  Assigning colors to {n_sogs} SOGs globally...")

    cmap = get_colormap(n_sogs)

    sog_colors = {}
    for i, sog_id in enumerate(sorted(all_sog_ids)):
        sog_colors[sog_id] = cmap(i / max(1, n_sogs))

    return sog_colors


def filter_alignments_for_cluster(alignments, track_ranges):
    """Filter alignments to only those within cluster ranges."""
    cluster_alignments = []
    seen = set()

    for (org1, chr1), (min1, max1) in track_ranges.items():
        if (org1, chr1) not in alignments:
            continue

        for align in alignments[(org1, chr1)]:
            start1, end1, org2, chr2, start2, end2, orientation = align

            if not (min1 <= start1 <= max1 and min1 <= end1 <= max1):
                continue

            if (org2, chr2) not in track_ranges:
                continue

            min2, max2 = track_ranges[(org2, chr2)]
            if not (min2 <= start2 <= max2 and min2 <= end2 <= max2):
                continue

            link_key = tuple(sorted([
                (org1, chr1, start1, end1),
                (org2, chr2, start2, end2)
            ]))
            if link_key in seen:
                continue
            seen.add(link_key)

            cluster_alignments.append((
                org1, chr1, start1, end1,
                org2, chr2, start2, end2,
                orientation
            ))

    return cluster_alignments


def draw_pairwise_cluster(cluster_pair: Dict, gene_to_sog: Dict,
                         sog_colors: Dict, alignments: Dict,
                         orientation: Dict, show_labels: bool,
                         output_dir: str, gene_margin: int,
                         gap_threshold: int, rejection_log: FeatureRejectionLog) -> Dict:
    """Draw a single pairwise cluster comparison."""

    # Extract data from cluster_pair
    species1 = cluster_pair['species1']
    chrom1 = cluster_pair['chrom1']
    genes1 = cluster_pair['genes1']

    species2 = cluster_pair['species2']
    chrom2 = cluster_pair['chrom2']
    genes2 = cluster_pair['genes2']

    # Create output directory structure
    pair_dir = os.path.join(output_dir, species1, species2)
    os.makedirs(pair_dir, exist_ok=True)

    # Build genes_by_species_chrom structure
    genes_by_species_chrom = {
        (species1, chrom1): genes1,
        (species2, chrom2): genes2
    }

    # Determine genomic ranges for each track
    track_ranges = {}
    for (species, chrom), genes in genes_by_species_chrom.items():
        min_start = min(g[0] for g in genes)
        max_end = max(g[1] for g in genes)
        track_ranges[(species, chrom)] = (min_start, max_end)

    # Filter alignments to cluster ranges
    cluster_alignments = filter_alignments_for_cluster(alignments, track_ranges)

    # Determine track orientations (first track = reference = forward)
    track_orientation = {
        (species1, chrom1): 'forward'
    }

    # Orient species2 track relative to species1
    if (species1 in orientation[species2][chrom2] and
        chrom1 in orientation[species2][chrom2][species1]):
        oris = orientation[species2][chrom2][species1][chrom1]
        if oris.count('forward') > oris.count('reverse'):
            track_orientation[(species2, chrom2)] = 'forward'
        else:
            track_orientation[(species2, chrom2)] = 'reverse'
    else:
        track_orientation[(species2, chrom2)] = 'forward'

    # Create pygenomeviz visualization
    gv = GenomeViz()
    shifts = {}
    genome_sizes = {}
    track_objects = {}
    actually_drawn_genes = set()

    # Add tracks
    for (org, chromo) in sorted(genes_by_species_chrom.keys()):
        genes = genes_by_species_chrom[(org, chromo)]
        min_start, max_end = track_ranges[(org, chromo)]

        shift = min_start
        genome_size = max_end - min_start

        # Build gene_drawings with margins
        gene_drawings = []
        for gene in genes:
            start_l, end_l, strand, gene_id = gene
            start_l, end_l = int(start_l), int(end_l)

            plot_start = start_l - shift
            plot_end = end_l - shift

            if track_orientation[(org, chromo)] == 'reverse':
                plot_start, plot_end = turn(plot_start, plot_end, [(0, genome_size)])

            if plot_start < 0 or plot_end < 0:
                rejection_log.add_gene(
                    gene_id, "NEGATIVE_COORDINATES",
                    f"Coordinates negative after transformation: start={plot_start}, end={plot_end}",
                    org, chromo, (plot_start, plot_end)
                )
                continue

            gene_drawings.append((max(0, plot_start - gene_margin), min(plot_end + gene_margin, genome_size)))

        if not gene_drawings:
            continue

        # Sort and build segments based on gaps
        gene_drawings = sorted(gene_drawings)
        target_ranges = []
        i = 1
        s1 = gene_drawings[0][0]
        e1 = gene_drawings[0][1]

        while i < len(gene_drawings):
            s2 = gene_drawings[i][0]
            e2 = gene_drawings[i][1]
            if s2 - e1 > gap_threshold:
                target_ranges.append((s1, e1))
                s1 = s2
            e1 = max(e1, e2)
            i += 1

        # Handle final segment - s1/e1 track current segment start/max-end
        target_ranges.append((s1, e1))

        name = f"{org} {chromo} {track_orientation[(org, chromo)]}"
        name = name.replace('_', ' ')
        shifts[name] = shift

        track = gv.add_feature_track(name=name, segments=target_ranges)

        genome_sizes[name] = []
        for seg_idx, segment in enumerate(track.segments):
            plot_seg_start = int(segment.start)
            plot_seg_end = int(segment.end)

            if track_orientation[(org, chromo)] == 'reverse':
                orig_start = genome_size - plot_seg_end
                orig_end = genome_size - plot_seg_start
                genomic_start = orig_start + shift
                genomic_end = orig_end + shift
                segment.add_sublabel(f'{genomic_start}-{genomic_end}')
            else:
                genomic_start = plot_seg_start + shift
                genomic_end = plot_seg_end + shift
                segment.add_sublabel(f'{genomic_start}-{genomic_end}')

            genome_sizes[name].append((plot_seg_start, plot_seg_end))

        track_objects[(org, chromo)] = (track, name)

        # Add gene features
        for gene in genes:
            orig_genomic_start, orig_genomic_end, strand, gene_id = gene
            sog_id = gene_to_sog.get(gene_id)
            if not sog_id:
                continue

            color = sog_colors.get(sog_id, 'gray')

            plot_start = orig_genomic_start - shift
            plot_end = orig_genomic_end - shift
            plot_strand = strand

            if track_orientation[(org, chromo)] == 'reverse':
                plot_start, plot_end = turn(plot_start, plot_end, genome_sizes[name])
                plot_strand = -1 if strand == 1 else 1

            if plot_start < 0 or plot_end < 0:
                rejection_log.add_gene(
                    gene_id, "NEGATIVE_COORDINATES_AFTER_TURN",
                    f"Coordinates negative after turn: start={plot_start}, end={plot_end}",
                    org, chromo, (plot_start, plot_end)
                )
                continue

            label = gene_id if show_labels else ''

            added = False
            if plot_start > plot_end:
                rejection_log.add_gene(
                    gene_id, "REVERSED_COORDINATES",
                    f"start={plot_start} > end={plot_end} after transformation",
                    org, chromo, (plot_start, plot_end)
                )
            else:
                for seg_idx, segment in enumerate(track.segments):
                    if not (segment.start <= plot_start < plot_end <= segment.end):
                        continue

                    try:
                        segment.add_feature(
                            plot_start, plot_end,
                            facecolor=color,
                            strand=plot_strand,
                            plotstyle='arrow',
                            label=label
                        )
                        actually_drawn_genes.add(gene_id)
                        added = True
                        rejection_log.total_genes_processed += 1
                        break
                    except (ValueError, TypeError) as e:
                        rejection_log.add_gene(
                            gene_id, "PYGENOMEVIZ_VALIDATION",
                            f"Segment {seg_idx} rejected: {str(e)}",
                            org, chromo, (plot_start, plot_end)
                        )
                        continue

            if not added and plot_start <= plot_end:
                rejection_log.total_genes_processed += 1

    # Add alignment links
    done = set()
    for align in cluster_alignments:
        org1, chr1, start1, end1, org2, chr2, start2, end2, ori_alignment = align

        link_key = ((org1, chr1, start1, end1), (org2, chr2, start2, end2))
        if link_key in done or (link_key[1], link_key[0]) in done:
            continue
        done.add(link_key)

        if (org1, chr1) not in track_objects or (org2, chr2) not in track_objects:
            continue

        track1, name1 = track_objects[(org1, chr1)]
        track2, name2 = track_objects[(org2, chr2)]

        # Transform coordinates for track 1
        shift1 = shifts[name1]
        plot_start1 = start1 - shift1
        plot_end1 = end1 - shift1
        if track_orientation[(org1, chr1)] == 'reverse':
            plot_start1, plot_end1 = turn(plot_start1, plot_end1, genome_sizes[name1])

        # Transform coordinates for track 2
        shift2 = shifts[name2]
        plot_start2 = start2 - shift2
        plot_end2 = end2 - shift2

        # Swap coordinates based on alignment orientation and track orientations
        if ((ori_alignment == 'reverse' and
             track_orientation[(org1, chr1)] == track_orientation[(org2, chr2)]) or
            (ori_alignment == 'forward' and
             track_orientation[(org1, chr1)] != track_orientation[(org2, chr2)])):
            plot_start2, plot_end2 = plot_end2, plot_start2

        if track_orientation[(org2, chr2)] == 'reverse':
            plot_start2, plot_end2 = turn(plot_start2, plot_end2, genome_sizes[name2])

        # Find segments
        seg1 = None
        for segment in track1.segments:
            if segment.start <= plot_start1 <= segment.end and segment.start <= plot_end1 <= segment.end:
                seg1 = segment.name
                break

        seg2 = None
        for segment in track2.segments:
            if segment.start <= plot_start2 <= segment.end and segment.start <= plot_end2 <= segment.end:
                seg2 = segment.name
                break

        if seg1 and seg2:
            rejection_log.total_alignments_attempted += 1
            try:
                gv.add_link(
                    (name1, seg1, plot_start1, plot_end1),
                    (name2, seg2, plot_start2, plot_end2),
                    color="skyblue",
                    inverted_color="lime",
                    curve=True
                )
                rejection_log.alignments_drawn += 1
            except Exception:
                pass

    # Add legend
    unique_sogs = set()
    for gene_id in actually_drawn_genes:
        if gene_id in gene_to_sog:
            unique_sogs.add(gene_to_sog[gene_id])

    legend_handles = [
        mpatches.Patch(color=sog_colors[sog], label=sog)
        for sog in sorted(unique_sogs)
        if sog in sog_colors
    ]

    # Save with descriptive filename
    n_genes1 = len(genes1)
    n_genes2 = len(genes2)
    n_shared = cluster_pair['n_shared']
    cluster1_idx = cluster_pair['cluster1_idx']
    cluster2_idx = cluster_pair['cluster2_idx']

    filename = (f"cluster_{cluster1_idx}_vs_{cluster2_idx}_"
                f"{n_shared}sogs_{n_genes1}+{n_genes2}genes.svg")
    output_path = os.path.join(pair_dir, filename)

    fig = gv.plotfig()
    if legend_handles:
        plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    fig.savefig(output_path, bbox_inches='tight')
    plt.close(fig)

    return {
        'species1': species1,
        'species2': species2,
        'output_path': output_path,
        'n_genes1': n_genes1,
        'n_genes2': n_genes2,
        'n_shared_sogs': n_shared,
        'success': True
    }


# Global variables for multiprocessing workers (set in main)
_worker_gene_to_sog = None
_worker_sog_colors = None
_worker_alignments = None
_worker_orientation = None
_worker_show_labels = None
_worker_output_dir = None
_worker_gene_margin = None
_worker_gap_threshold = None


def _draw_pairwise_worker(cluster_pair):
    """Worker function for parallel drawing (module-level for pickling)."""
    rejection_log = FeatureRejectionLog()
    try:
        return draw_pairwise_cluster(
            cluster_pair, _worker_gene_to_sog, _worker_sog_colors,
            _worker_alignments, _worker_orientation, _worker_show_labels,
            _worker_output_dir, _worker_gene_margin, _worker_gap_threshold,
            rejection_log
        )
    except Exception as e:
        return {
            'species1': cluster_pair['species1'],
            'species2': cluster_pair['species2'],
            'output_path': None,
            'success': False,
            'error': str(e)
        }


def draw_pairwise_clusters_parallel(cluster_pairs: List[Dict],
                                   gene_to_sog: Dict, sog_colors: Dict,
                                   alignments: Dict, orientation: Dict,
                                   show_labels: bool, output_dir: str,
                                   gene_margin: int, gap_threshold: int,
                                   n_workers: int) -> List[Dict]:
    """Draw cluster pairs in parallel using multiprocessing."""
    from multiprocessing import Pool

    print(f"\n=== Drawing {len(cluster_pairs)} cluster pairs with {n_workers} workers ===")

    if n_workers == 1:
        # Sequential processing
        results = []
        for i, pair in enumerate(cluster_pairs, 1):
            if i % 10 == 0 or i == len(cluster_pairs):
                print(f"Progress: {i}/{len(cluster_pairs)} pairs ({100*i/len(cluster_pairs):.1f}%)")

            rejection_log = FeatureRejectionLog()
            try:
                result = draw_pairwise_cluster(
                    pair, gene_to_sog, sog_colors, alignments, orientation,
                    show_labels, output_dir, gene_margin, gap_threshold, rejection_log
                )
                results.append(result)
            except Exception as e:
                results.append({
                    'species1': pair['species1'],
                    'species2': pair['species2'],
                    'output_path': None,
                    'success': False,
                    'error': str(e)
                })

        return results

    else:
        # Parallel processing - set global variables for workers
        global _worker_gene_to_sog, _worker_sog_colors, _worker_alignments
        global _worker_orientation, _worker_show_labels, _worker_output_dir
        global _worker_gene_margin, _worker_gap_threshold

        _worker_gene_to_sog = gene_to_sog
        _worker_sog_colors = sog_colors
        _worker_alignments = alignments
        _worker_orientation = orientation
        _worker_show_labels = show_labels
        _worker_output_dir = output_dir
        _worker_gene_margin = gene_margin
        _worker_gap_threshold = gap_threshold

        with Pool(n_workers) as pool:
            results = []
            for i, result in enumerate(pool.imap(_draw_pairwise_worker, cluster_pairs), 1):
                if i % 10 == 0 or i == len(cluster_pairs):
                    print(f"Progress: {i}/{len(cluster_pairs)} pairs ({100*i/len(cluster_pairs):.1f}%)")
                results.append(result)

        return results


def main():
    parser = argparse.ArgumentParser(
        description='Draw pairwise SynOrthogroup cluster comparisons',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--sog_file', required=True,
                       help='Path to SynOrthogroups.tsv')
    parser.add_argument('--coords_file', required=True,
                       help='Path to SynOrthogroups_WithCoords.tsv')
    parser.add_argument('--alignments_file', required=True,
                       help='Path to pairwise_alignments_table')
    parser.add_argument('--cluster_margin', type=int, default=100000,
                       help='Max distance (bp) between genes in same cluster (default: 100000)')
    parser.add_argument('--min_shared_sogs', type=int, default=1,
                       help='Min shared SOGs to draw cluster pair (default: 1)')
    parser.add_argument('--min_cluster_size', type=int, default=1,
                       help='Min genes per cluster to include (default: 1)')
    parser.add_argument('--gene_margin', type=int, default=10000,
                       help='Margin around genes for visualization (default: 10000)')
    parser.add_argument('--gap_threshold', type=int, default=200000,
                       help='Min gap size to create segment break (default: 200000)')
    parser.add_argument('--pairs_file', type=str, default=None,
                       help='Optional file with species pairs (format: species1<tab>species2)')
    parser.add_argument('--cores', type=int, default=1,
                       help='Number of parallel workers (default: 1)')
    parser.add_argument('--output_dir', default='pairwise_clusters',
                       help='Output directory (default: pairwise_clusters)')
    parser.add_argument('--show_labels', action='store_true',
                       help='Display gene IDs on arrows')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("\n=== Phase 1: Loading data ===")
    print(f"Loading SOGs from {args.sog_file}...")
    sogs, gene_to_sog = parse_synorthogroups_tsv(args.sog_file)
    print(f"  Found {len(sogs)} SOGs with {len(gene_to_sog)} genes")

    print(f"Loading coordinates from {args.coords_file}...")
    coords, gene_to_coord = load_gene_coordinates(args.coords_file)
    print(f"  Loaded coordinates for {len(gene_to_coord)} genes")

    print(f"Loading alignments from {args.alignments_file}...")
    alignments, orientation = load_alignments(args.alignments_file)
    print(f"  Loaded {sum(len(v) for v in alignments.values())} alignment entries")

    print("\n=== Phase 2: Clustering within species ===")
    print(f"Clustering genes with margin={args.cluster_margin}bp...")
    clusters, cluster_ids = cluster_within_species(coords, args.cluster_margin)
    print(f"  Created {len(cluster_ids)} clusters across all species")

    print("\n=== Phase 3: Filtering clusters by size ===")
    if args.min_cluster_size > 1:
        clusters = filter_clusters_by_size(clusters, args.min_cluster_size)

    print("\n=== Phase 4: Loading species pairs ===")
    if args.pairs_file:
        species_pairs = load_species_pairs(args.pairs_file)
    else:
        species_pairs = generate_all_pairs(list(coords.keys()))

    print("\n=== Phase 5: Assigning global SOG colors ===")
    sog_colors = assign_colors_global(sogs, gene_to_sog)

    print("\n=== Phase 6: Finding cluster pairs ===")
    all_cluster_pairs = []
    for species1, species2 in species_pairs:
        pairs = find_cluster_pairs(
            species1, species2, clusters, gene_to_sog, args.min_shared_sogs
        )
        all_cluster_pairs.extend(pairs)
        print(f"  {species1} <-> {species2}: {len(pairs)} cluster pairs")

    print(f"\nTotal cluster pairs to draw: {len(all_cluster_pairs)}")

    if len(all_cluster_pairs) == 0:
        print("No cluster pairs found meeting criteria. Exiting.")
        return

    print("\n=== Phase 7: Drawing pairwise clusters ===")
    results = draw_pairwise_clusters_parallel(
        all_cluster_pairs, gene_to_sog, sog_colors,
        alignments, orientation, args.show_labels, args.output_dir,
        args.gene_margin, args.gap_threshold, args.cores
    )

    print("\n=== Phase 8: Summary ===")
    successful = sum(1 for r in results if r['success'])
    failed = len(results) - successful
    print(f"Successfully drew {successful}/{len(results)} cluster pairs")
    if failed > 0:
        print(f"Failed: {failed} pairs")

    # Save summary file
    summary_file = os.path.join(args.output_dir, "pairwise_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Pairwise Cluster Drawing Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total cluster pairs drawn: {len(results)}\n")
        f.write(f"Successful: {successful}\n")
        f.write(f"Failed: {failed}\n\n")

        # Group by species pairs
        by_species = defaultdict(list)
        for r in results:
            key = (r['species1'], r['species2'])
            by_species[key].append(r)

        f.write("By Species Pair:\n")
        for (sp1, sp2), pair_results in sorted(by_species.items()):
            f.write(f"\n{sp1} <-> {sp2}: {len(pair_results)} clusters\n")

    print(f"Summary saved to {summary_file}")
    print(f"\n=== Done! Output saved to {args.output_dir}/ ===")


if __name__ == '__main__':
    main()
