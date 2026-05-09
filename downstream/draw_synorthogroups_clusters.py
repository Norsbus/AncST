#!/usr/bin/env python3

"""
Draw SynOrthogroup clusters by genomic proximity.

This script clusters genes from synorthogroups by spatial proximity within
species, then groups clusters across species based on shared SOG membership.
Each super-cluster (multi-species group) is drawn as a single image with
genes colored by SOG membership.

Usage:
    ./draw_synorthogroups_clusters.py --sog_file SynOrthogroups.tsv \\
        --coords_file SynOrthogroups_WithCoords.tsv \\
        --alignments_file pairwise_alignments_table \\
        --cluster_margin 100000 \\
        --min_shared_sogs 1 \\
        --gene_margin 10000 \\
        --gap_threshold 200000 \\
        --output_dir clusters_output \\
        --show_labels \\
        --interactive

Arguments:
    --sog_file: Path to SynOrthogroups.tsv (SOG assignments)
    --coords_file: Path to SynOrthogroups_WithCoords.tsv (gene coordinates)
    --alignments_file: Path to pairwise_alignments_table (for drawing links)
    --cluster_margin: Max distance (bp) between genes in same cluster (default: 100000)
    --min_shared_sogs: Min shared SOGs to connect clusters (default: 1)
    --gene_margin: Margin around genes for visualization (default: 10000)
    --gap_threshold: Min gap size to create segment break (default: 200000)
    --output_dir: Output directory for images (default: clusters_output)
    --show_labels: Display gene IDs on arrows (default: off)
    --interactive: Generate interactive HTML visualizations using Plotly
"""

import argparse
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
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
    """
    Parse SynOrthogroups.tsv file.

    Returns:
        sogs: dict {sog_id: {species: [gene_ids]}}
        gene_to_sog: dict {gene_id: sog_id}
    """
    sogs = {}
    gene_to_sog = {}

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
                    genes = [g.strip() for g in gene_str.split(',')]
                    sogs[sog_id][species] = genes

                    # Build reverse mapping
                    for gene in genes:
                        gene_to_sog[gene] = sog_id

    return sogs, gene_to_sog


def load_gene_coordinates(coords_tsv_path):
    """
    Load gene coordinates from SynOrthogroups_WithCoords.tsv.

    Returns:
        coords: dict {species: {chromosome: [(start, end, strand, gene_id)]}}
        gene_to_coord: dict {gene_id: (species, chromosome, start, end, strand)}
    """
    coords = {}
    gene_to_coord = {}

    with open(coords_tsv_path) as f:
        header = f.readline()  # Skip header

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
    """
    Load pairwise alignments from file.

    Format: no org1 chromo1 start1 end1 org2 chromo2 start2 end2 ori_alignment score snd1 snd2

    Returns:
        alignments: Dict[(org1, chr1): List[(start1, end1, org2, chr2, start2, end2, ori)]]
        orientation: Dict[org1][chr1][org2][chr2] = [ori_alignment, ...]
    """
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

            # Normalize coordinates
            if start1 > end1:
                start1, end1 = end1, start1
            if start2 > end2:
                start2, end2 = end2, start2

            # Store bidirectionally
            alignments[(org1, chr1)].append((start1, end1, org2, chr2, start2, end2, ori))
            alignments[(org2, chr2)].append((start2, end2, org1, chr1, start1, end1, ori))

            # Store orientation info
            orientation[org1][chr1][org2][chr2].append(ori)
            # Reverse direction for bidirectional lookup
            reverse_ori = 'reverse' if ori == 'forward' else 'forward'
            orientation[org2][chr2][org1][chr1].append(reverse_ori)

    return alignments, orientation


def cluster_genes_on_chromosome(genes, margin):
    """
    Cluster genes by proximity on a single chromosome.

    Args:
        genes: List of (start, end, strand, gene_id)
        margin: Maximum gap between genes in same cluster

    Returns:
        List of clusters, each cluster is list of genes
    """
    if not genes:
        return []

    sorted_genes = sorted(genes, key=lambda g: g[0])  # sort by start
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
    """
    Cluster genes by proximity within each species.

    Args:
        coords: dict {species: {chromosome: [(start, end, strand, gene_id)]}}
        margin: Maximum gap between genes in same cluster

    Returns:
        clusters: dict {species: {chromosome: [cluster0, cluster1, ...]}}
        cluster_ids: list of cluster IDs as "species_chromosome_index"
    """
    clusters = {}
    cluster_ids = []

    for species in coords:
        clusters[species] = {}
        for chrom in coords[species]:
            genes = coords[species][chrom]
            chrom_clusters = cluster_genes_on_chromosome(genes, margin)
            clusters[species][chrom] = chrom_clusters

            # Generate cluster IDs
            for idx in range(len(chrom_clusters)):
                cluster_id = f"{species}_{chrom}_{idx}"
                cluster_ids.append(cluster_id)

    return clusters, cluster_ids


class UnionFind:
    """Union-Find data structure for clustering."""

    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px != py:
            self.parent[px] = py

    def get_groups(self):
        groups = defaultdict(set)
        for e in self.parent:
            groups[self.find(e)].add(e)
        return list(groups.values())


def build_super_clusters(clusters, gene_to_sog, min_shared_sogs=1):
    """
    Group species-clusters into super-clusters based on shared SOGs.

    Args:
        clusters: Dict[species][chromosome][cluster_idx] = [genes]
        gene_to_sog: Dict[gene_id] = sog_id
        min_shared_sogs: Minimum number of shared SOGs to connect clusters

    Returns:
        List of super-clusters, each a set of cluster_ids
    """
    # Build cluster_id -> set of SOGs mapping
    cluster_sogs = {}
    all_cluster_ids = []

    for species in clusters:
        for chrom in clusters[species]:
            for idx, gene_list in enumerate(clusters[species][chrom]):
                cluster_id = f"{species}_{chrom}_{idx}"
                all_cluster_ids.append(cluster_id)

                # Extract SOGs from this cluster
                sogs_in_cluster = set()
                for gene in gene_list:
                    gene_id = gene[3]  # gene = (start, end, strand, gene_id)
                    if gene_id in gene_to_sog:
                        sogs_in_cluster.add(gene_to_sog[gene_id])
                cluster_sogs[cluster_id] = sogs_in_cluster

    # Initialize Union-Find
    uf = UnionFind(all_cluster_ids)

    # Connect clusters with sufficient shared SOGs
    for i, cluster1 in enumerate(all_cluster_ids):
        for cluster2 in all_cluster_ids[i+1:]:
            sogs1 = cluster_sogs[cluster1]
            sogs2 = cluster_sogs[cluster2]
            shared = sogs1 & sogs2  # Set intersection

            if len(shared) >= min_shared_sogs:
                uf.union(cluster1, cluster2)

    return uf.get_groups()


def filter_multi_species_clusters(super_clusters):
    """
    Filter to keep only super-clusters with genes from multiple species.

    Args:
        super_clusters: List of sets of cluster_ids

    Returns:
        drawable: List of super-clusters with multiple species
        single_species: List of cluster_ids from single-species clusters
    """
    drawable = []
    single_species_ids = []

    for super_cluster in super_clusters:
        species_set = set()
        for cluster_id in super_cluster:
            species = cluster_id.split('_')[0]
            species_set.add(species)

        if len(species_set) > 1:
            drawable.append(super_cluster)
        else:
            single_species_ids.extend(super_cluster)

    return drawable, single_species_ids


def get_colormap(n_colors):
    """
    Get appropriate colormap based on number of colors needed.

    Args:
        n_colors: Number of distinct colors needed

    Returns:
        matplotlib colormap
    """
    if n_colors <= 10:
        return plt.cm.tab10
    elif n_colors <= 20:
        return plt.cm.tab20
    elif n_colors <= 40:
        # Combine tab20 + tab20b + tab20c
        colors = []
        for cmap_name in ['tab20', 'tab20b', 'tab20c']:
            cmap = plt.cm.get_cmap(cmap_name)
            for i in range(20):
                colors.append(cmap(i / 20))
        return ListedColormap(colors[:n_colors])
    else:
        # Use HSV for large counts
        return plt.cm.hsv


def assign_colors(drawable_super_clusters, gene_to_sog, clusters):
    """
    Assign colors to SOGs based on frequency.

    Args:
        drawable_super_clusters: List of super-clusters to draw
        gene_to_sog: Dict mapping gene_id to sog_id
        clusters: Dict of all clusters

    Returns:
        Dict {sog_id: (R, G, B, A)}
    """
    # Collect all SOG IDs in drawable clusters
    all_sogs = set()

    for super_cluster in drawable_super_clusters:
        for cluster_id in super_cluster:
            species, chrom, idx = cluster_id.rsplit('_', 2)
            idx = int(idx)

            if species in clusters and chrom in clusters[species]:
                if idx < len(clusters[species][chrom]):
                    for gene in clusters[species][chrom][idx]:
                        gene_id = gene[3]
                        if gene_id in gene_to_sog:
                            all_sogs.add(gene_to_sog[gene_id])

    # Get appropriate colormap
    n_sogs = len(all_sogs)
    if n_sogs == 0:
        return {}

    cmap = get_colormap(n_sogs)

    # Assign colors
    sog_colors = {}
    for i, sog_id in enumerate(sorted(all_sogs)):
        sog_colors[sog_id] = cmap(i / max(1, n_sogs))

    return sog_colors


def filter_alignments_for_cluster(alignments, track_ranges):
    """
    Filter alignments to only those within cluster ranges.

    Args:
        alignments: Dict[(org, chr)] = [(start1, end1, org2, chr2, start2, end2, ori), ...]
        track_ranges: Dict[(org, chr): (min_start, max_end)]

    Returns:
        List of alignments to draw as links
    """
    cluster_alignments = []
    seen = set()

    for (org1, chr1), (min1, max1) in track_ranges.items():
        if (org1, chr1) not in alignments:
            continue

        for align in alignments[(org1, chr1)]:
            start1, end1, org2, chr2, start2, end2, orientation = align

            # Check if alignment is within cluster range for org1
            if not (min1 <= start1 <= max1 and min1 <= end1 <= max1):
                continue

            # Check if matching species/chromosome is in this cluster
            if (org2, chr2) not in track_ranges:
                continue

            # Check if alignment is within cluster range for org2
            min2, max2 = track_ranges[(org2, chr2)]
            if not (min2 <= start2 <= max2 and min2 <= end2 <= max2):
                continue

            # Avoid duplicates
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


def verify_plot_to_genomic(plot_start, plot_end, shift, genome_size, track_orientation):
    """
    INDEPENDENT reverse transformation: plot coordinates -> genomic coordinates.

    Based on LOGIC, not inverse of forward transformation:
    - Forward track: Plot is a window starting at genomic 'shift'
      -> genomic = plot + shift
    - Reverse track: Plot is flipped. Rightmost plot (genome_size) = leftmost genomic (shift)
      -> genomic = shift + (genome_size - plot_position)

    Args:
        plot_start, plot_end: Coordinates in plot space [0, genome_size]
        shift: Genomic coordinate where track starts
        genome_size: Length of track in plot space
        track_orientation: 'forward' or 'reverse'

    Returns:
        (genomic_start, genomic_end): Original genomic coordinates
    """
    if track_orientation == 'forward':
        # Plot 0 = genomic shift, plot increases with genomic position
        genomic_start = plot_start + shift
        genomic_end = plot_end + shift
    else:  # reverse
        # Plot 0 = genomic (shift + genome_size), plot increases as genomic decreases
        # Plot genome_size = genomic shift
        genomic_start = shift + (genome_size - plot_end)
        genomic_end = shift + (genome_size - plot_start)

    return int(genomic_start), int(genomic_end)


def verify_strand_transformation(orig_strand, plot_strand, track_orientation):
    """
    INDEPENDENT strand verification based on BIOLOGICAL/VISUAL LOGIC.

    Strands represent direction relative to genomic coordinates:
    - + strand: points in direction of INCREASING genomic coordinates
    - - strand: points in direction of DECREASING genomic coordinates

    When track is reversed:
    - We flip the view (genomic coordinates now increase right-to-left)
    - A + strand gene still points toward increasing coordinates
    - But visually it now appears to point the opposite direction
    - Therefore strand must be flipped for correct visual representation

    LOGIC (independent of implementation):
    - Forward track: Visual direction matches genomic direction
      -> plot_strand = orig_strand
    - Reverse track: Visual direction is opposite of genomic direction
      -> plot_strand = -orig_strand

    Args:
        orig_strand: Original strand from data (1 or -1)
        plot_strand: Strand used in plot (1 or -1)
        track_orientation: 'forward' or 'reverse'

    Returns:
        expected_plot_strand: What the plot strand SHOULD be based on logic
    """
    if track_orientation == 'forward':
        # Forward track: genomic and visual directions align
        expected_strand = orig_strand
    else:  # reverse
        # Reverse track: view is flipped, so strand must flip too
        # + becomes -, - becomes +
        expected_strand = -orig_strand

    return expected_strand


def draw_super_cluster(super_cluster_id, super_cluster, clusters, gene_to_sog,
                       sog_colors, alignments, orientation, show_labels, output_dir,
                       gene_margin, gap_threshold, rejection_log, track_order=None):
    """
    Draw a single super-cluster visualization.

    Args:
        super_cluster_id: Index of this super-cluster
        super_cluster: Set of cluster_ids in this super-cluster
        clusters: Dict of all clusters
        gene_to_sog: Dict mapping gene_id to sog_id
        sog_colors: Dict mapping sog_id to color
        alignments: Dict of alignment data
        orientation: Dict of orientation data
        show_labels: Whether to display gene IDs
        output_dir: Output directory for images
        gene_margin: Margin around genes for visualization
        gap_threshold: Min gap size to create segment break
        rejection_log: FeatureRejectionLog instance for tracking rejections
        track_order: Optional list of (species, chromosome) tuples specifying track order.
                    If None, uses alphabetical ordering.

    Returns:
        Tuple of (species_count, gene_count, actually_drawn_genes, output_filename)
    """
    # Initialize debug log
    debug_log = []
    debug_log.append("=" * 80)
    debug_log.append(f"CLUSTER {super_cluster_id} DEBUG LOG")
    debug_log.append("=" * 80)
    debug_log.append("")
    # Extract all genes by (species, chrom)
    genes_by_species_chrom = defaultdict(list)

    for cluster_id in super_cluster:
        # Parse cluster_id more carefully to handle species names with underscores
        # Format: species_chromosome_idx
        # Need to find the correct split point
        parts = cluster_id.split('_')

        # Try to find the split by checking which combinations exist in clusters
        found = False
        for split_idx in range(len(parts) - 2, 0, -1):
            test_species = '_'.join(parts[:split_idx])
            test_chrom = '_'.join(parts[split_idx:-1])
            test_idx_str = parts[-1]

            if test_species in clusters and test_chrom in clusters[test_species]:
                try:
                    idx = int(test_idx_str)
                    if idx < len(clusters[test_species][test_chrom]):
                        species = test_species
                        chrom = test_chrom
                        gene_list = clusters[species][chrom][idx]
                        genes_by_species_chrom[(species, chrom)].extend(gene_list)
                        found = True
                        break
                except ValueError:
                    continue

        if not found:
            debug_log.append(f"WARNING: Could not parse cluster_id: {cluster_id}")
            continue

    if not genes_by_species_chrom:
        return None

    # Determine genomic ranges for each track
    track_ranges = {}
    for (species, chrom), genes in genes_by_species_chrom.items():
        min_start = min(g[0] for g in genes)
        max_end = max(g[1] for g in genes)
        track_ranges[(species, chrom)] = (min_start, max_end)

    # Filter alignments to cluster ranges
    cluster_alignments = filter_alignments_for_cluster(alignments, track_ranges)

    # Determine track order (manual or alphabetical)
    if track_order is not None:
        # Use manual track order (filter to only those in this super-cluster)
        ordered_tracks = [t for t in track_order if t in genes_by_species_chrom]
    else:
        # Use alphabetical ordering
        ordered_tracks = sorted(genes_by_species_chrom.keys())

    # Determine track orientations (reference + majority vote)
    # Pick first track as reference (always 'forward')
    ref_org, ref_chromo = ordered_tracks[0]

    track_orientation = {}
    track_orientation[(ref_org, ref_chromo)] = 'forward'

    # Orient other tracks relative to reference by majority vote
    for (org, chromo) in genes_by_species_chrom.keys():
        if (org, chromo) in track_orientation:
            continue  # Already set (reference)

        # Check alignments to reference
        if (ref_org in orientation[org][chromo] and
            ref_chromo in orientation[org][chromo][ref_org]):
            oris = orientation[org][chromo][ref_org][ref_chromo]
            if oris.count('forward') > oris.count('reverse'):
                track_orientation[(org, chromo)] = 'forward'
            else:
                track_orientation[(org, chromo)] = 'reverse'
        else:
            # No alignments to reference, default to forward
            track_orientation[(org, chromo)] = 'forward'

    # Create pygenomeviz visualization
    gv = GenomeViz()
    shifts = {}
    genome_sizes = {}
    track_objects = {}
    actually_drawn_genes = set()  # Track genes that are ACTUALLY drawn

    # Add tracks with segment splitting logic (from draw_synorthogroups.py lines 511-558)
    for (org, chromo) in ordered_tracks:
        genes = genes_by_species_chrom[(org, chromo)]
        min_start, max_end = track_ranges[(org, chromo)]

        shift = min_start
        genome_size = max_end - min_start

        # Build gene_drawings with margins (like draw_synorthogroups.py)
        gene_drawings = []
        for gene in genes:
            start_l, end_l, strand, gene_id = gene
            start_l, end_l = int(start_l), int(end_l)

            # Apply shift
            plot_start = start_l - shift
            plot_end = end_l - shift

            # Apply orientation transformation
            if track_orientation[(org, chromo)] == 'reverse':
                plot_start, plot_end = turn(plot_start, plot_end, [(0, genome_size)])

            # Only check for negative coordinates (will be rejected by segment bounds later if needed)
            if plot_start < 0 or plot_end < 0:
                rejection_log.add_gene(
                    gene_id, "NEGATIVE_COORDINATES",
                    f"Coordinates negative after transformation: start={plot_start}, end={plot_end}",
                    org, chromo, (plot_start, plot_end)
                )
                continue

            gene_drawings.append((max(0, plot_start - gene_margin), min(plot_end + gene_margin, genome_size)))

        if not gene_drawings:
            # ALL genes for this track were rejected - log this
            print(f"  WARNING: Track {org} {chromo} skipped - all {len(genes)} genes rejected in pre-processing")
            for gene in genes:
                gene_id = gene[3]
                rejection_log.add_gene(
                    gene_id, "TRACK_SKIPPED_ALL_GENES_REJECTED",
                    f"Track {org} {chromo} skipped because all genes were rejected during coordinate transformation",
                    org, chromo, ()
                )
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

        # Track name includes orientation
        name = f"{org} {chromo} {track_orientation[(org, chromo)]}"
        name = name.replace('_', ' ')
        shifts[name] = shift

        # Add track with segments
        track = gv.add_feature_track(name=name, segments=target_ranges)

        # Store segment boundaries and add labels
        genome_sizes[name] = []
        debug_log.append(f"\n--- TRACK: {org} {chromo} ---")
        debug_log.append(f"  Orientation: {track_orientation[(org, chromo)]}")
        debug_log.append(f"  Shift (genomic start): {shift}")
        debug_log.append(f"  Genome size (plot space): {genome_size}")
        debug_log.append(f"  Number of segments: {len(track.segments)}")

        for seg_idx, segment in enumerate(track.segments):
            plot_seg_start = int(segment.start)
            plot_seg_end = int(segment.end)

            # Calculate genomic coordinates for label
            if track_orientation[(org, chromo)] == 'reverse':
                # For reverse tracks, need to un-turn coordinates before adding shift
                # turn() did: new_end = len_line - start, new_start = new_end - (end - start)
                # Inverse: orig_start = len_line - new_end, orig_end = len_line - new_start
                orig_start = genome_size - plot_seg_end
                orig_end = genome_size - plot_seg_start
                genomic_start = orig_start + shift
                genomic_end = orig_end + shift
                segment.add_sublabel(f'{genomic_start}-{genomic_end}')
            else:
                # For forward tracks, simply add shift
                genomic_start = plot_seg_start + shift
                genomic_end = plot_seg_end + shift
                segment.add_sublabel(f'{genomic_start}-{genomic_end}')

            genome_sizes[name].append((plot_seg_start, plot_seg_end))

            # Verify using independent function
            verify_start, verify_end = verify_plot_to_genomic(
                plot_seg_start, plot_seg_end, shift, genome_size, track_orientation[(org, chromo)]
            )

            match = (verify_start == genomic_start and verify_end == genomic_end)
            debug_log.append(f"\n  SEGMENT {seg_idx}:")
            debug_log.append(f"    Plot range: ({plot_seg_start}, {plot_seg_end})")
            debug_log.append(f"    Genomic (label): {genomic_start}-{genomic_end}")
            debug_log.append(f"    Genomic (verify): {verify_start}-{verify_end}")
            debug_log.append(f"    VERIFICATION: {'PASS' if match else 'FAIL'}")

        track_objects[(org, chromo)] = (track, name)

        # Add gene features with proper segment handling
        genes_on_track = []
        for gene in genes:
            orig_genomic_start, orig_genomic_end, strand, gene_id = gene
            sog_id = gene_to_sog.get(gene_id)
            if not sog_id:
                continue

            color = sog_colors.get(sog_id, 'gray')

            # Apply shift
            plot_start = orig_genomic_start - shift
            plot_end = orig_genomic_end - shift
            plot_strand = strand

            # Apply orientation transformation
            if track_orientation[(org, chromo)] == 'reverse':
                plot_start, plot_end = turn(plot_start, plot_end, genome_sizes[name])
                plot_strand = -1 if strand == 1 else 1

            # Check for negative coordinates
            if plot_start < 0 or plot_end < 0:
                rejection_log.add_gene(
                    gene_id, "NEGATIVE_COORDINATES_AFTER_TURN",
                    f"Coordinates negative after turn: start={plot_start}, end={plot_end}",
                    org, chromo, (plot_start, plot_end)
                )
                continue

            label = gene_id if show_labels else ''

            # Try to add to appropriate segment
            added = False
            drawn_segment_idx = None

            # Pre-validate: check reversed coordinates
            if plot_start > plot_end:
                rejection_log.add_gene(
                    gene_id, "REVERSED_COORDINATES",
                    f"start={plot_start} > end={plot_end} after transformation",
                    org, chromo, (plot_start, plot_end)
                )
                added = False
            else:
                for seg_idx, segment in enumerate(track.segments):
                    # Pre-validate: check segment bounds
                    if not (segment.start <= plot_start < plot_end <= segment.end):
                        # Gene doesn't fit this segment, try next
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
                        drawn_segment_idx = seg_idx
                        rejection_log.total_genes_processed += 1
                        break  # Successfully added
                    except ValueError as e:
                        # pygenomeviz coordinate validation failed
                        rejection_log.add_gene(
                            gene_id, "PYGENOMEVIZ_VALIDATION",
                            f"Segment {seg_idx} rejected: {str(e)}",
                            org, chromo, (plot_start, plot_end),
                            segment_range=(segment.start, segment.end)
                        )
                        continue  # Try next segment
                    except TypeError as e:
                        # Type mismatch (float vs int, wrong strand type, etc.)
                        rejection_log.add_gene(
                            gene_id, "TYPE_ERROR",
                            f"Segment {seg_idx} type error: {str(e)}",
                            org, chromo, (plot_start, plot_end),
                            types=f"start:{type(plot_start)}, end:{type(plot_end)}, strand:{type(plot_strand)}"
                        )
                        continue
                    except Exception as e:
                        # Unexpected error - should investigate
                        rejection_log.add_gene(
                            gene_id, "UNEXPECTED_ERROR",
                            f"Segment {seg_idx}: {type(e).__name__}: {str(e)}",
                            org, chromo, (plot_start, plot_end)
                        )
                        continue

            # After loop - if gene wasn't added
            if not added and plot_start <= plot_end:
                # Gene doesn't fit any segment
                seg_ranges = [(s.start, s.end) for s in track.segments]
                rejection_log.add_gene(
                    gene_id, "NO_CONTAINING_SEGMENT",
                    f"Gene at ({plot_start}, {plot_end}) doesn't fit any segment",
                    org, chromo, (plot_start, plot_end),
                    segment_ranges=seg_ranges
                )
                rejection_log.total_genes_processed += 1

            if added:
                # Verify coordinates using independent function
                verify_start, verify_end = verify_plot_to_genomic(
                    plot_start, plot_end, shift, genome_size, track_orientation[(org, chromo)]
                )

                # Verify strand using independent logic
                expected_strand = verify_strand_transformation(
                    strand, plot_strand, track_orientation[(org, chromo)]
                )

                coord_match = (verify_start == orig_genomic_start and verify_end == orig_genomic_end)
                strand_match = (plot_strand == expected_strand)

                genes_on_track.append({
                    'gene_id': gene_id,
                    'sog_id': sog_id,
                    'orig_genomic': (orig_genomic_start, orig_genomic_end),
                    'orig_strand': strand,
                    'plot': (plot_start, plot_end),
                    'plot_strand': plot_strand,
                    'verify_genomic': (verify_start, verify_end),
                    'expected_strand': expected_strand,
                    'segment_idx': drawn_segment_idx,
                    'coord_match': coord_match,
                    'strand_match': strand_match
                })

        # Log all genes for this track
        debug_log.append(f"\n  GENES on this track: {len(genes_on_track)}")
        for gene_info in genes_on_track:
            debug_log.append(f"\n    GENE: {gene_info['gene_id']} (SOG: {gene_info['sog_id']})")
            debug_log.append(f"      Original genomic: ({gene_info['orig_genomic'][0]}, {gene_info['orig_genomic'][1]}), strand={'+' if gene_info['orig_strand']==1 else '-'}")
            debug_log.append(f"      Plot coords: ({gene_info['plot'][0]}, {gene_info['plot'][1]}), strand={'+' if gene_info['plot_strand']==1 else '-'}")
            debug_log.append(f"      Segment: {gene_info['segment_idx']}")
            debug_log.append(f"      Verify genomic: ({gene_info['verify_genomic'][0]}, {gene_info['verify_genomic'][1]})")
            debug_log.append(f"      Verify strand: expected={'+' if gene_info['expected_strand']==1 else '-'}, got={'+' if gene_info['plot_strand']==1 else '-'}")
            debug_log.append(f"      VERIFICATION:")
            debug_log.append(f"        Coordinates: {'PASS' if gene_info['coord_match'] else 'FAIL'}")
            debug_log.append(f"        Strand: {'PASS' if gene_info['strand_match'] else 'FAIL'}")

    # Add alignment links between tracks with orientation handling
    debug_log.append("\n" + "=" * 80)
    debug_log.append("ALIGNMENT LINKS")
    debug_log.append("=" * 80)

    done = set()
    drawn_alignments = []
    for align in cluster_alignments:
        org1, chr1, start1, end1, org2, chr2, start2, end2, ori_alignment = align

        # Skip duplicates
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

        # Find which segments contain these coordinates
        seg1 = None
        for segment in track1.segments:
            segstart = int(segment.start)
            segend = int(segment.end)
            if segstart <= plot_start1 <= segend and segstart <= plot_end1 <= segend:
                seg1 = segment.name
                break

        seg2 = None
        for segment in track2.segments:
            segstart = int(segment.start)
            segend = int(segment.end)
            if segstart <= plot_start2 <= segend and segstart <= plot_end2 <= segend:
                seg2 = segment.name
                break

        # Add link if both segments found
        if seg1 and seg2:
            link_id = f"{org1}_{chr1}_{start1}-{end1}_to_{org2}_{chr2}_{start2}-{end2}"
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

                # Verify both endpoints
                verify1_start, verify1_end = verify_plot_to_genomic(
                    plot_start1, plot_end1, shifts[name1],
                    max(gs[1] for gs in genome_sizes[name1]),
                    track_orientation[(org1, chr1)]
                )
                verify2_start, verify2_end = verify_plot_to_genomic(
                    plot_start2, plot_end2, shifts[name2],
                    max(gs[1] for gs in genome_sizes[name2]),
                    track_orientation[(org2, chr2)]
                )

                # For alignments, coordinates might be swapped for drawing
                # Check if they match in either order
                match1 = ((verify1_start == start1 and verify1_end == end1) or
                         (verify1_start == end1 and verify1_end == start1))
                match2 = ((verify2_start == start2 and verify2_end == end2) or
                         (verify2_start == end2 and verify2_end == start2))

                drawn_alignments.append({
                    'org1': org1,
                    'chr1': chr1,
                    'orig1': (start1, end1),
                    'plot1': (plot_start1, plot_end1),
                    'verify1': (verify1_start, verify1_end),
                    'seg1': seg1,
                    'match1': match1,
                    'org2': org2,
                    'chr2': chr2,
                    'orig2': (start2, end2),
                    'plot2': (plot_start2, plot_end2),
                    'verify2': (verify2_start, verify2_end),
                    'seg2': seg2,
                    'match2': match2,
                    'orientation': ori_alignment
                })
            except ValueError as e:
                rejection_log.add_alignment(
                    link_id, "LINK_VALIDATION_FAILED",
                    f"pygenomeviz rejected link: {str(e)}",
                    org1=org1, org2=org2,
                    coords1=(plot_start1, plot_end1),
                    coords2=(plot_start2, plot_end2)
                )
            except TypeError as e:
                rejection_log.add_alignment(
                    link_id, "LINK_TYPE_ERROR",
                    f"Type mismatch: {str(e)}",
                    org1=org1, org2=org2
                )
            except Exception as e:
                rejection_log.add_alignment(
                    link_id, "LINK_UNEXPECTED_ERROR",
                    f"{type(e).__name__}: {str(e)}",
                    org1=org1, org2=org2
                )

    # Log all drawn alignments
    debug_log.append(f"\nTotal alignments drawn: {len(drawn_alignments)}")
    for idx, align_info in enumerate(drawn_alignments):
        debug_log.append(f"\n  ALIGNMENT {idx}:")
        debug_log.append(f"    Orientation: {align_info['orientation']}")
        debug_log.append(f"    Org1: {align_info['org1']} {align_info['chr1']}")
        debug_log.append(f"      Original genomic: ({align_info['orig1'][0]}, {align_info['orig1'][1]})")
        debug_log.append(f"      Plot coords: ({align_info['plot1'][0]}, {align_info['plot1'][1]})")
        debug_log.append(f"      Segment: {align_info['seg1']}")
        debug_log.append(f"      Verify genomic: ({align_info['verify1'][0]}, {align_info['verify1'][1]})")
        debug_log.append(f"      VERIFICATION: {'PASS' if align_info['match1'] else 'FAIL'}")
        debug_log.append(f"    Org2: {align_info['org2']} {align_info['chr2']}")
        debug_log.append(f"      Original genomic: ({align_info['orig2'][0]}, {align_info['orig2'][1]})")
        debug_log.append(f"      Plot coords: ({align_info['plot2'][0]}, {align_info['plot2'][1]})")
        debug_log.append(f"      Segment: {align_info['seg2']}")
        debug_log.append(f"      Verify genomic: ({align_info['verify2'][0]}, {align_info['verify2'][1]})")
        debug_log.append(f"      VERIFICATION: {'PASS' if align_info['match2'] else 'FAIL'}")

    # Add summary statistics
    debug_log.append("\n" + "=" * 80)
    debug_log.append("SUMMARY")
    debug_log.append("=" * 80)
    debug_log.append(f"Total genes actually drawn: {len(actually_drawn_genes)}")
    debug_log.append(f"Total alignments drawn: {len(drawn_alignments)}")

    # Save debug log to file
    debug_dir = f"{output_dir}/debug"
    os.makedirs(debug_dir, exist_ok=True)
    debug_file = f"{debug_dir}/cluster_{super_cluster_id}_debug.txt"
    with open(debug_file, 'w') as f:
        f.write('\n'.join(debug_log))

    # Add legend
    unique_sogs = set()
    for gene_id in actually_drawn_genes:
        if gene_id in gene_to_sog:
            unique_sogs.add(gene_to_sog[gene_id])

    legend_handles = [
        mpatches.Patch(color=sog_colors[sog], label=sog)
        for sog in sorted(unique_sogs)
        if sog in sog_colors  # Only include SOGs that have colors assigned
    ]

    # Save figure
    species_count = len(set(org for org, _ in genes_by_species_chrom.keys()))
    gene_count = sum(len(genes) for genes in genes_by_species_chrom.values())

    output_filename = f"{output_dir}/cluster_{super_cluster_id}_{species_count}species_{gene_count}genes.svg"

    fig = gv.plotfig()
    if legend_handles:
        plt.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1), loc='upper left')
    fig.savefig(output_filename, bbox_inches='tight')
    plt.close(fig)

    return (species_count, gene_count, actually_drawn_genes, output_filename)


def generate_not_drawn_report(single_species_cluster_ids, clusters, gene_to_sog,
                              gene_to_coord, all_drawn_genes, output_dir):
    """
    Generate report of genes not drawn.

    Args:
        single_species_cluster_ids: List of cluster IDs that were not drawn
        clusters: Dict of all clusters
        gene_to_sog: Dict mapping gene_id to sog_id
        gene_to_coord: Dict mapping gene_id to coordinates
        all_drawn_genes: Set of gene IDs that were actually drawn
        output_dir: Output directory
    """
    not_drawn_by_species = defaultdict(list)
    not_drawn_genes_set = set()  # Track to avoid double-counting

    # Genes from single-species clusters
    for cluster_id in single_species_cluster_ids:
        species, chrom, idx = cluster_id.rsplit('_', 2)
        idx = int(idx)

        if species not in clusters or chrom not in clusters[species]:
            continue
        if idx >= len(clusters[species][chrom]):
            continue

        for gene in clusters[species][chrom][idx]:
            gene_id = gene[3]
            if gene_id not in all_drawn_genes:
                not_drawn_by_species[species].append((gene_id, 'single-species cluster'))
                not_drawn_genes_set.add(gene_id)  # Track this gene

    # Genes that were in multi-species clusters but not actually drawn
    for gene_id in gene_to_coord:
        if gene_id not in all_drawn_genes and gene_id in gene_to_sog and gene_id not in not_drawn_genes_set:
            species = gene_to_coord[gene_id][0]
            not_drawn_by_species[species].append((gene_id, 'not drawn in visualization'))

    # Write report
    report_file = f"{output_dir}/elements_not_drawn.txt"
    with open(report_file, 'w') as f:
        f.write("# Elements Not Drawn\n")
        f.write("# Format: gene_id  chromosome  start  end  SOG_id  reason\n\n")

        total_not_drawn = 0
        for species in sorted(not_drawn_by_species.keys()):
            f.write(f"\n## {species}\n")
            for gene_id, reason in sorted(not_drawn_by_species[species]):
                if gene_id not in gene_to_coord:
                    continue
                coord = gene_to_coord[gene_id]
                sog_id = gene_to_sog.get(gene_id, 'unknown')
                f.write(f"{gene_id}\t{coord[1]}\t{coord[2]}\t{coord[3]}\t{sog_id}\t{reason}\n")
                total_not_drawn += 1

        f.write(f"\n# Total not drawn: {total_not_drawn}\n")

    print(f"  Not-drawn report: {report_file}")
    return total_not_drawn


def generate_stats_report(sogs, gene_to_sog, all_drawn_genes, total_not_drawn,
                         drawable_super_clusters, output_dir):
    """
    Generate statistics report.

    Args:
        sogs: Dict of all SOGs
        gene_to_sog: Dict mapping gene_id to sog_id
        all_drawn_genes: Set of gene IDs that were actually drawn
        total_not_drawn: Count of genes not drawn
        drawable_super_clusters: List of super-clusters that were drawn
        output_dir: Output directory
    """
    total_genes = len(gene_to_sog)
    drawn_count = len(all_drawn_genes)

    # Write stats
    stats_file = f"{output_dir}/clustering_stats.txt"
    with open(stats_file, 'w') as f:
        f.write("Clustering Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total synorthogroups: {len(sogs)}\n")
        f.write(f"Total genes: {total_genes}\n")
        f.write(f"  - Drawn in clusters: {drawn_count} ({100*drawn_count/max(1,total_genes):.1f}%)\n")
        f.write(f"  - Not drawn: {total_not_drawn} ({100*total_not_drawn/max(1,total_genes):.1f}%)\n")
        f.write(f"Images created: {len(drawable_super_clusters)}\n")

        if drawable_super_clusters and drawn_count > 0:
            avg_genes = drawn_count / len(drawable_super_clusters)
            f.write(f"Average genes per image: {avg_genes:.1f}\n")

    print(f"  Statistics report: {stats_file}")
    print(f"\nSummary:")
    print(f"  Total genes: {total_genes}")
    print(f"  Drawn: {drawn_count} ({100*drawn_count/max(1,total_genes):.1f}%)")
    print(f"  Not drawn: {total_not_drawn} ({100*total_not_drawn/max(1,total_genes):.1f}%)")
    print(f"  Images: {len(drawable_super_clusters)}")

    # Sanity check
    if drawn_count + total_not_drawn != total_genes:
        print(f"\n  WARNING: Accounting mismatch! {drawn_count} + {total_not_drawn} != {total_genes}")


def load_manual_selection(filepath):
    """
    Load manual track selection file.

    Format: genome<tab>chromosome per line

    Args:
        filepath: Path to manual selection file

    Returns:
        List of (genome, chromosome) tuples in order

    Raises:
        ValueError: If file format is incorrect
    """
    tracks = []

    if not os.path.exists(filepath):
        raise ValueError(f"ERROR: Manual selection file not found: {filepath}")

    with open(filepath) as f:
        for line_num, line in enumerate(f, 1):
            # Strip whitespace
            original_line = line
            line = line.strip()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            # Split by tab
            parts = line.split('\t')

            # Validate field count
            if len(parts) != 2:
                raise ValueError(
                    f"ERROR: Manual selection file format error at line {line_num}\n"
                    f"  Expected: genome<tab>chromosome\n"
                    f"  Found: {len(parts)} fields\n"
                    f"  Line: '{original_line.rstrip()}'"
                )

            genome = parts[0].strip()
            chrom = parts[1].strip()

            # Validate non-empty
            if not genome:
                raise ValueError(
                    f"ERROR: Manual selection file format error at line {line_num}\n"
                    f"  Genome field is empty\n"
                    f"  Line: '{original_line.rstrip()}'"
                )

            if not chrom:
                raise ValueError(
                    f"ERROR: Manual selection file format error at line {line_num}\n"
                    f"  Chromosome field is empty\n"
                    f"  Line: '{original_line.rstrip()}'"
                )

            tracks.append((genome, chrom))

    # Validate file is not empty
    if not tracks:
        raise ValueError(
            f"ERROR: Manual selection file contains no valid entries\n"
            f"  File: {filepath}\n"
            f"  (empty lines and comments starting with # are ignored)"
        )

    return tracks


def validate_and_build_manual_cluster(manual_tracks, clusters, coords):
    """
    Validate manual tracks and build a super-cluster from them.

    Args:
        manual_tracks: List of (genome, chromosome) tuples in order
        clusters: Dict of all clusters
        coords: Dict of gene coordinates (to check if species/chr exist)

    Returns:
        - manual_super_cluster: Set of cluster_ids from specified tracks
        - valid_tracks: List of (species, chr) that passed validation (in order)

    Raises:
        ValueError: If no valid tracks found
    """
    manual_cluster_ids = set()
    valid_tracks = []

    print("\n=== Validating Manual Track Selection ===")

    for species, chrom in manual_tracks:
        # Check if species exists in dataset
        if species not in coords:
            print(f"  WARNING: Species '{species}' not found in dataset. Skipping.")
            print(f"           Available species: {', '.join(sorted(coords.keys()))}")
            continue

        # Check if chromosome exists for this species
        if chrom not in coords[species]:
            print(f"  WARNING: Chromosome '{chrom}' not found for species '{species}'. Skipping.")
            available_chroms = ', '.join(sorted(coords[species].keys()))
            print(f"           Available chromosomes for {species}: {available_chroms}")
            continue

        # Check if there are any clusters for this (species, chrom)
        if species not in clusters or chrom not in clusters[species]:
            print(f"  WARNING: No clusters exist for {species} {chrom}. Skipping.")
            continue

        if len(clusters[species][chrom]) == 0:
            print(f"  WARNING: No clusters exist for {species} {chrom} (empty list). Skipping.")
            continue

        # Add all cluster_ids from this (species, chrom)
        n_clusters = len(clusters[species][chrom])
        for idx in range(n_clusters):
            cluster_id = f"{species}_{chrom}_{idx}"
            manual_cluster_ids.add(cluster_id)

        valid_tracks.append((species, chrom))
        print(f"  {species} {chrom}: {n_clusters} cluster(s)")

    # Validate at least one valid track
    if not valid_tracks:
        raise ValueError(
            "ERROR: No valid tracks found in manual selection file\n"
            "  All specified genome/chromosome combinations were either:\n"
            "  - Not found in the dataset\n"
            "  - Had no genes/clusters after filtering"
        )

    print(f"\nManual selection: {len(valid_tracks)} valid tracks, "
          f"{len(manual_cluster_ids)} total cluster(s)")

    return manual_cluster_ids, valid_tracks


def main():
    parser = argparse.ArgumentParser(
        description='Draw SynOrthogroup clusters by genomic proximity',
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
                       help='Min shared SOGs to connect clusters (default: 1)')
    parser.add_argument('--gene_margin', type=int, default=10000,
                       help='Margin around genes for visualization (default: 10000)')
    parser.add_argument('--gap_threshold', type=int, default=200000,
                       help='Min gap size to create segment break (default: 200000)')
    parser.add_argument('--output_dir', default='clusters_output',
                       help='Output directory (default: clusters_output)')
    parser.add_argument('--show_labels', action='store_true',
                       help='Display gene IDs on arrows')
    parser.add_argument('--interactive', action='store_true',
                       help='Generate interactive HTML visualizations (not yet implemented)')
    parser.add_argument('--manual_selection', type=str, default=None,
                       help='Manual track selection file (format: genome<tab>chromosome per line). '
                            'Draws only specified tracks in specified order. Bypasses super-clustering.')

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("=== Phase 1: Loading data ===")
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

    # Check if manual selection is provided
    if args.manual_selection:
        print("\n=== Phase 3: Manual Track Selection ===")
        # Load manual selection file (strict validation)
        manual_tracks = load_manual_selection(args.manual_selection)
        print(f"  Loaded {len(manual_tracks)} tracks from {args.manual_selection}")

        # Validate and build manual super-cluster
        manual_super_cluster, valid_tracks = validate_and_build_manual_cluster(
            manual_tracks, clusters, coords
        )

        # Create drawable list with single manual super-cluster
        drawable = [manual_super_cluster]
        track_order = valid_tracks  # Preserve order for drawing
        single_species_ids = []  # No single-species filtering in manual mode

        print("\n=== Phase 4: Assigning colors ===")
        sog_colors = assign_colors(drawable, gene_to_sog, clusters)
        print(f"  Assigned colors to {len(sog_colors)} SOGs")

    else:
        print("\n=== Phase 3: Building super-clusters ===")
        print(f"Grouping clusters with min_shared_sogs={args.min_shared_sogs}...")
        super_clusters = build_super_clusters(clusters, gene_to_sog, args.min_shared_sogs)
        print(f"  Created {len(super_clusters)} super-clusters")

        print("\n=== Phase 4: Filtering multi-species clusters ===")
        drawable, single_species_ids = filter_multi_species_clusters(super_clusters)
        print(f"  Drawable (multi-species): {len(drawable)}")
        print(f"  Single-species clusters: {len(single_species_ids)}")

        print("\n=== Phase 5: Assigning colors ===")
        sog_colors = assign_colors(drawable, gene_to_sog, clusters)
        print(f"  Assigned colors to {len(sog_colors)} SOGs")

        track_order = None  # Use alphabetical ordering

    print("\n=== Phase 6: Drawing super-clusters ===")
    # Initialize global rejection log
    rejection_log = FeatureRejectionLog()

    all_drawn_genes = set()
    for idx, super_cluster in enumerate(drawable):
        print(f"Drawing cluster {idx+1}/{len(drawable)}...", end=' ')
        result = draw_super_cluster(
            idx, super_cluster, clusters, gene_to_sog, sog_colors,
            alignments, orientation, args.show_labels, args.output_dir,
            args.gene_margin, args.gap_threshold, rejection_log,
            track_order=track_order if args.manual_selection else None
        )
        if result:
            species_count, gene_count, drawn_genes, filename = result
            all_drawn_genes.update(drawn_genes)
            print(f"OK ({species_count} species, {gene_count} total genes, {len(drawn_genes)} drawn)")
        else:
            print("SKIPPED (no genes found)")

    print("\n=== Phase 7: Generating reports ===")
    total_not_drawn = generate_not_drawn_report(single_species_ids, clusters, gene_to_sog,
                                                 gene_to_coord, all_drawn_genes, args.output_dir)

    print("\n=== Phase 8: Generating statistics ===")
    generate_stats_report(sogs, gene_to_sog, all_drawn_genes, total_not_drawn,
                         drawable, args.output_dir)

    print("\n=== Phase 9: Generating rejection reports ===")
    # Update rejection log counts
    rejection_log.genes_drawn = len(all_drawn_genes)

    # Print summary to stdout
    print(rejection_log.summary())

    # Write detailed TSV
    log_file = os.path.join(args.output_dir, "feature_rejection_log.tsv")
    rejection_log.write_tsv(log_file)
    print(f"Detailed rejection log written to: {log_file}")

    # Save summary to file too
    summary_file = os.path.join(args.output_dir, "feature_rejection_summary.txt")
    with open(summary_file, 'w') as f:
        f.write(rejection_log.summary())
    print(f"Summary saved to: {summary_file}")

    if args.interactive:
        print("\nNote: Interactive visualization not yet implemented")

    print(f"\n=== Done! Output saved to {args.output_dir}/ ===")


if __name__ == '__main__':
    main()
