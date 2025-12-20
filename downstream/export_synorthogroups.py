#!/usr/bin/env python3

"""
Export SynOrthogroups from Synthology Pipeline Output

Generates OrthoFinder-compatible outputs from synthology union graphs.
Works for both protein and nucleotide synthology runs.

Ensures complete gene accounting - all input genes are represented either in
SynOrthogroups or as unassigned genes (filtered or single-species singletons).

Usage:
    export_synorthogroups.py <synthology_output_dir> [options]

The output_dir argument should be your synthology OUTPUT directory
(e.g., synthology_out_prot), which contains final_union_graph_*.pickle files.
The script will create a synorthogroups/ subdirectory automatically.

Examples:
    cd synthology_out_prot
    export_synorthogroups.py .

    export_synorthogroups.py synthology_out_prot --graph-version no_new_edges
"""

import argparse
import os
import pickle
import sys
from collections import defaultdict
import networkx as nx


# =============================================================================
# Data Loading
# =============================================================================

def load_synthology_data(output_dir, graph_version):
    """
    Load synthology output data.

    Args:
        output_dir: Path to synthology output directory
        graph_version: 'no_new_edges' or 'with_new_edges'

    Returns:
        tuple: (graph, gene_species, parsed_gff, metadata)
    """
    graph_file = f"{output_dir}/final_union_graph_{graph_version}.pickle"

    if not os.path.exists(graph_file):
        raise FileNotFoundError(f"Graph file not found: {graph_file}")

    print(f"Loading {graph_file}...")
    with open(graph_file, 'rb') as f:
        graph, gene_species, metadata = pickle.load(f)

    # Load parsed GFF (contains all input genes with coordinates)
    # Try both prot and nucl naming conventions
    parsed_gff_file = None
    for candidate in ['parsed_faa_gff', 'parsed_fna_gff']:
        candidate_path = f"{output_dir}/{candidate}"
        if os.path.exists(candidate_path):
            parsed_gff_file = candidate_path
            break

    if not parsed_gff_file:
        raise FileNotFoundError(
            f"Parsed GFF file not found. Looked for:\n"
            f"  {output_dir}/parsed_faa_gff\n"
            f"  {output_dir}/parsed_fna_gff"
        )

    print(f"Loading {parsed_gff_file}...")
    with open(parsed_gff_file, 'rb') as f:
        parsed_gff = pickle.load(f)

    return graph, gene_species, parsed_gff, metadata


def get_all_input_genes(parsed_gff):
    """
    Extract all gene IDs from input GFF annotations.

    Args:
        parsed_gff: Dictionary {species: {chromosome: [(start, end, strand, gene_id)]}}

    Returns:
        dict: {gene_id: species} for all input genes
    """
    all_genes = {}
    for species, chroms in parsed_gff.items():
        for chrom, genes in chroms.items():
            for start, end, strand, gene_id in genes:
                all_genes[gene_id] = species
    return all_genes


# =============================================================================
# Gene Classification
# =============================================================================

def classify_genes(all_input_genes, graph, gene_species):
    """
    Classify all input genes into SynOrthogroups, singletons, or filtered.

    Args:
        all_input_genes: dict {gene_id: species} from input GFFs
        graph: NetworkX graph (final union graph)
        gene_species: dict {gene_id: species} from pipeline

    Returns:
        tuple: (multi_species_sogs, single_species_singletons, filtered_genes)
            - multi_species_sogs: list of sets (each set is a SOG with 2+ species)
            - single_species_singletons: list of gene IDs
            - filtered_genes: list of gene IDs removed during pipeline
    """
    genes_in_graph = set(graph.nodes())
    all_input_set = set(all_input_genes.keys())

    # Genes that were filtered out during pipeline
    filtered_genes = list(all_input_set - genes_in_graph)

    # Extract connected components from graph
    components = list(nx.connected_components(graph))

    multi_species_sogs = []
    single_species_singletons = []

    for comp in components:
        # Count species in this component
        species_in_comp = set()
        for gene in comp:
            if gene in gene_species:
                species_in_comp.add(gene_species[gene])

        if len(species_in_comp) >= 2:
            # Multi-species SynOrthogroup
            multi_species_sogs.append(comp)
        else:
            # Single-species component - treat as singletons
            single_species_singletons.extend(list(comp))

    return multi_species_sogs, single_species_singletons, filtered_genes


def is_single_copy_sog(sog, gene_species):
    """
    Check if a SynOrthogroup is single-copy (exactly 1 gene per species).

    Args:
        sog: set of gene IDs
        gene_species: dict {gene_id: species}

    Returns:
        bool: True if single-copy orthogroup
    """
    species_counts = defaultdict(int)
    for gene in sog:
        if gene in gene_species:
            species_counts[gene_species[gene]] += 1

    # Single-copy if all species have exactly 1 gene
    return all(count == 1 for count in species_counts.values())


# =============================================================================
# Output Writers
# =============================================================================

def write_synorthogroups_tsv(sogs, gene_species, species_list, output_path):
    """
    Write main SynOrthogroups.tsv file (OrthoFinder format).

    Format:
        - Tab-separated columns: SynOrthogroup | Species1 | Species2 | ...
        - Genes within same species: comma-separated
        - Empty cells for species not in SOG
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    with open(output_path, 'w') as f:
        # Header
        f.write("SynOrthogroup\t" + "\t".join(species_list) + "\n")

        # Sort SOGs by size (descending) for consistent ordering
        sorted_sogs = sorted(sogs, key=len, reverse=True)

        for sog_idx, sog in enumerate(sorted_sogs):
            sog_id = f"SOG{sog_idx:07d}"

            # Organize genes by species
            species_genes = defaultdict(list)
            for gene in sorted(sog):
                if gene in gene_species:
                    species_genes[gene_species[gene]].append(gene)

            # Build row
            row_parts = [sog_id]
            for species in species_list:
                if species in species_genes:
                    # Comma-separated genes for this species
                    row_parts.append(", ".join(sorted(species_genes[species])))
                else:
                    # Empty cell
                    row_parts.append("")

            f.write("\t".join(row_parts) + "\n")


def write_gene_count_tsv(sogs, gene_species, species_list, output_path):
    """
    Write SynOrthogroups.GeneCount.tsv (count matrix).

    Format:
        SynOrthogroup | Species1_count | Species2_count | ... | Total
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    with open(output_path, 'w') as f:
        # Header
        f.write("SynOrthogroup\t" + "\t".join(species_list) + "\tTotal\n")

        sorted_sogs = sorted(sogs, key=len, reverse=True)

        for sog_idx, sog in enumerate(sorted_sogs):
            sog_id = f"SOG{sog_idx:07d}"

            # Count genes per species
            species_counts = defaultdict(int)
            for gene in sog:
                if gene in gene_species:
                    species_counts[gene_species[gene]] += 1

            # Build row
            row_parts = [sog_id]
            for species in species_list:
                row_parts.append(str(species_counts.get(species, 0)))
            row_parts.append(str(len(sog)))

            f.write("\t".join(row_parts) + "\n")


def write_single_copy_orthologues(sogs, gene_species, output_path):
    """
    Write SynOrthogroups_SingleCopyOrthologues.txt.

    Lists SOG IDs that are single-copy (exactly 1 gene per species).
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    single_copy_sogs = []
    sorted_sogs = sorted(sogs, key=len, reverse=True)

    for sog_idx, sog in enumerate(sorted_sogs):
        if is_single_copy_sog(sog, gene_species):
            sog_id = f"SOG{sog_idx:07d}"
            single_copy_sogs.append(sog_id)

    with open(output_path, 'w') as f:
        for sog_id in single_copy_sogs:
            f.write(sog_id + "\n")


def write_unassigned_genes(all_input_genes, singletons, filtered, gene_species, output_path):
    """
    Write SynOrthogroups_UnassignedGenes.tsv.

    Includes:
        - Filtered genes (removed during pipeline)
        - Single-species singletons (in graph but only connected within species)
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    with open(output_path, 'w') as f:
        f.write("Gene_ID\tSpecies\tReason\n")

        # Filtered genes (not in final graph)
        for gene in sorted(filtered):
            species = all_input_genes.get(gene, "unknown")
            f.write(f"{gene}\t{species}\tfiltered_during_pipeline\n")

        # Single-species singletons (in graph but isolated)
        for gene in sorted(singletons):
            species = gene_species.get(gene, "unknown")
            f.write(f"{gene}\t{species}\tsingle_species_singleton\n")


def write_statistics_overall(all_input_genes, sogs, singletons, filtered, gene_species, output_path):
    """
    Write Statistics_Overall.tsv (summary statistics).
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    total_input = len(all_input_genes)
    genes_in_sogs = sum(len(sog) for sog in sogs)
    genes_unassigned_filtered = len(filtered)
    genes_unassigned_singletons = len(singletons)
    percentage_in_sogs = (genes_in_sogs / total_input * 100) if total_input > 0 else 0

    num_single_copy = sum(1 for sog in sogs if is_single_copy_sog(sog, gene_species))

    sog_sizes = [len(sog) for sog in sogs]
    mean_size = sum(sog_sizes) / len(sog_sizes) if sog_sizes else 0
    median_size = sorted(sog_sizes)[len(sog_sizes)//2] if sog_sizes else 0

    # Get unique species
    species_set = set(all_input_genes.values())
    num_species = len(species_set)

    with open(output_path, 'w') as f:
        f.write("Metric\tValue\n")
        f.write(f"Number of species\t{num_species}\n")
        f.write(f"Total input genes\t{total_input}\n")
        f.write(f"Genes in SynOrthogroups\t{genes_in_sogs}\n")
        f.write(f"Genes unassigned (filtered)\t{genes_unassigned_filtered}\n")
        f.write(f"Genes unassigned (singletons)\t{genes_unassigned_singletons}\n")
        f.write(f"Percentage in SynOrthogroups\t{percentage_in_sogs:.1f}\n")
        f.write(f"Number of SynOrthogroups\t{len(sogs)}\n")
        f.write(f"Number of single-copy SynOrthogroups\t{num_single_copy}\n")
        f.write(f"Mean SynOrthogroup size\t{mean_size:.1f}\n")
        f.write(f"Median SynOrthogroup size\t{median_size}\n")


def write_statistics_per_species(all_input_genes, sogs, singletons, filtered,
                                 gene_species, species_list, output_path):
    """
    Write Statistics_PerSpecies.tsv (per-species statistics).
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    # Count input genes per species
    species_input_counts = defaultdict(int)
    for gene, species in all_input_genes.items():
        species_input_counts[species] += 1

    # Count genes in SOGs per species
    species_in_sogs = defaultdict(int)
    for sog in sogs:
        for gene in sog:
            if gene in gene_species:
                species_in_sogs[gene_species[gene]] += 1

    # Count unassigned per species
    species_unassigned = defaultdict(int)
    for gene in filtered:
        species = all_input_genes.get(gene)
        if species:
            species_unassigned[species] += 1
    for gene in singletons:
        species = gene_species.get(gene)
        if species:
            species_unassigned[species] += 1

    # Count single-copy SOGs per species
    species_single_copy = defaultdict(int)
    for sog in sogs:
        if is_single_copy_sog(sog, gene_species):
            for gene in sog:
                if gene in gene_species:
                    species_single_copy[gene_species[gene]] += 1

    with open(output_path, 'w') as f:
        f.write("Species\tInput_Genes\tIn_SynOrthogroups\tUnassigned\tPercentage\tSingleCopy_SOGs\n")

        for species in species_list:
            input_count = species_input_counts[species]
            in_sogs = species_in_sogs[species]
            unassigned = species_unassigned[species]
            percentage = (in_sogs / input_count * 100) if input_count > 0 else 0
            single_copy = species_single_copy[species]

            f.write(f"{species}\t{input_count}\t{in_sogs}\t{unassigned}\t"
                   f"{percentage:.1f}\t{single_copy}\n")


def write_synorthogroups_with_coords(sogs, gene_species, parsed_gff, output_path):
    """
    Write SynOrthogroups_WithCoords.tsv (detailed table with coordinates).
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    gene_coords = {}
    for species, chroms in parsed_gff.items():
        for chrom, genes in chroms.items():
            for start, end, strand, gene_id in genes:
                gene_coords[gene_id] = (species, chrom, start, end, strand)

    sorted_sogs = sorted(sogs, key=len, reverse=True)

    with open(output_path, 'w') as f:
        f.write("SynOrthogroup\tGene_ID\tSpecies\tChromosome\tStart\tEnd\tStrand\tRef_Gene\n")

        for sog_idx, sog in enumerate(sorted_sogs):
            sog_id = f"SOG{sog_idx:07d}"

            sog_genes = []
            for gene in sog:
                if gene in gene_coords:
                    species, chrom, start, end, strand = gene_coords[gene]
                    sog_genes.append((species, chrom, start, gene, end, strand))

            sog_genes.sort()

            for species, chrom, start, gene, end, strand in sog_genes:
                ref_gene = gene.split('_', 3)[-1] if '_' in gene else gene
                strand_str = '+' if strand == 1 else '-'
                f.write(f"{sog_id}\t{gene}\t{species}\t{chrom}\t{start}\t{end}\t{strand_str}\t{ref_gene}\n")


def write_pickle_format(sogs, all_input_genes, singletons, filtered,
                        gene_species, parsed_gff, graph_version, output_path):
    """
    Write synorthogroups.pickle (Python-friendly format).
    """
    print(f"  Writing {os.path.basename(output_path)}...")

    sorted_sogs = sorted(sogs, key=len, reverse=True)

    synorthogroups = {}
    gene_to_sog = {}

    for sog_idx, sog in enumerate(sorted_sogs):
        sog_id = f"SOG{sog_idx:07d}"
        synorthogroups[sog_id] = sorted(list(sog))
        for gene in sog:
            gene_to_sog[gene] = sog_id

    unassigned = {}
    for gene in filtered:
        unassigned[gene] = 'filtered_during_pipeline'
    for gene in singletons:
        unassigned[gene] = 'single_species_singleton'

    coordinates = {}
    for species, chroms in parsed_gff.items():
        for chrom, genes in chroms.items():
            for start, end, strand, gene_id in genes:
                strand_str = '+' if strand == 1 else '-'
                coordinates[gene_id] = (species, chrom, start, end, strand_str)

    metadata = {
        'graph_version': graph_version,
        'total_input_genes': len(all_input_genes),
        'total_sogs': len(sogs),
        'genes_in_sogs': sum(len(sog) for sog in sogs),
        'genes_unassigned': len(unassigned)
    }

    data = {
        'synorthogroups': synorthogroups,
        'gene_to_sog': gene_to_sog,
        'unassigned': unassigned,
        'coordinates': coordinates,
        'metadata': metadata
    }

    with open(output_path, 'wb') as f:
        pickle.dump(data, f)


# =============================================================================
# Main Export Function
# =============================================================================

def export_synorthogroups_version(output_dir, graph_version):
    """
    Export SynOrthogroups for a specific graph version.

    Args:
        output_dir: Path to synthology output directory
        graph_version: 'no_new_edges' or 'with_new_edges'
    """
    print(f"\n{'='*80}")
    print(f"EXPORTING SYNORTHOGROUPS ({graph_version})")
    print(f"{'='*80}\n")

    # Load data
    graph, gene_species, parsed_gff, metadata = load_synthology_data(output_dir, graph_version)

    print(f"  Graph nodes: {graph.number_of_nodes()}")
    print(f"  Graph edges: {graph.number_of_edges()}")

    # Get all input genes
    all_input_genes = get_all_input_genes(parsed_gff)
    print(f"  Total input genes: {len(all_input_genes)}")

    # Classify genes
    sogs, singletons, filtered = classify_genes(all_input_genes, graph, gene_species)

    print(f"  Multi-species SynOrthogroups: {len(sogs)}")
    print(f"  Single-species singletons: {len(singletons)}")
    print(f"  Filtered genes: {len(filtered)}")

    # Verify complete accounting
    total_accounted = sum(len(sog) for sog in sogs) + len(singletons) + len(filtered)
    if total_accounted != len(all_input_genes):
        print(f"  WARNING: Gene accounting mismatch!")
        print(f"    Input: {len(all_input_genes)}, Accounted: {total_accounted}")
    else:
        print(f"  ✓ Complete gene accounting verified")

    # Create output directory
    synortho_dir = f"{output_dir}/synorthogroups/{graph_version}"
    os.makedirs(synortho_dir, exist_ok=True)

    # Get sorted species list
    species_list = sorted(set(all_input_genes.values()))

    # Write all output files
    print(f"\nWriting output files to {synortho_dir}/")

    write_synorthogroups_tsv(
        sogs, gene_species, species_list,
        f"{synortho_dir}/SynOrthogroups.tsv"
    )

    write_gene_count_tsv(
        sogs, gene_species, species_list,
        f"{synortho_dir}/SynOrthogroups.GeneCount.tsv"
    )

    write_single_copy_orthologues(
        sogs, gene_species,
        f"{synortho_dir}/SynOrthogroups_SingleCopyOrthologues.txt"
    )

    write_unassigned_genes(
        all_input_genes, singletons, filtered, gene_species,
        f"{synortho_dir}/SynOrthogroups_UnassignedGenes.tsv"
    )

    write_statistics_overall(
        all_input_genes, sogs, singletons, filtered, gene_species,
        f"{synortho_dir}/Statistics_Overall.tsv"
    )

    write_statistics_per_species(
        all_input_genes, sogs, singletons, filtered, gene_species, species_list,
        f"{synortho_dir}/Statistics_PerSpecies.tsv"
    )

    write_synorthogroups_with_coords(
        sogs, gene_species, parsed_gff,
        f"{synortho_dir}/SynOrthogroups_WithCoords.tsv"
    )

    write_pickle_format(
        sogs, all_input_genes, singletons, filtered, gene_species, parsed_gff, graph_version,
        f"{synortho_dir}/synorthogroups.pickle"
    )

    print(f"\n✓ Export complete: {synortho_dir}/")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Export SynOrthogroups from Synthology pipeline output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From within synthology output directory (recommended)
  cd synthology_out_prot
  export_synorthogroups.py .

  # From parent directory
  export_synorthogroups.py synthology_out_prot

  # Export only one graph version
  export_synorthogroups.py . --graph-version no_new_edges

Note:
  - Argument is the synthology OUTPUT directory (containing final_union_graph_*.pickle)
  - Creates synorthogroups/ subdirectory automatically
  - Output saved to: <output_dir>/synorthogroups/<graph_version>/
        """
    )

    parser.add_argument(
        'output_dir',
        help='Synthology output directory (must contain final_union_graph_*.pickle)'
    )

    parser.add_argument(
        '--graph-version',
        choices=['no_new_edges', 'with_new_edges', 'both'],
        default='both',
        help='Which graph version(s) to export (default: both)'
    )

    args = parser.parse_args()

    # Validate output directory
    if not os.path.isdir(args.output_dir):
        print(f"ERROR: Directory not found: {args.output_dir}")
        print(f"\nUsage: Provide the synthology OUTPUT directory (e.g., synthology_out_prot)")
        print(f"       NOT the synorthogroups directory (which will be created)")
        print(f"\nExample: cd synthology_out_prot && export_synorthogroups.py .")
        sys.exit(1)

    # Check for required files
    required_files = ['final_union_graph_no_new_edges.pickle']
    for req_file in required_files:
        if not os.path.exists(f"{args.output_dir}/{req_file}"):
            print(f"ERROR: Required file not found: {req_file}")
            print(f"\nThe directory '{args.output_dir}' does not appear to be a synthology output directory.")
            print(f"Expected files: final_union_graph_*.pickle, parsed_faa_gff (or parsed_fna_gff)")
            print(f"\nMake sure you're providing the synthology OUTPUT directory, not synorthogroups/")
            sys.exit(1)

    # Export requested version(s)
    if args.graph_version in ['no_new_edges', 'both']:
        export_synorthogroups_version(args.output_dir, 'no_new_edges')

    if args.graph_version in ['with_new_edges', 'both']:
        export_synorthogroups_version(args.output_dir, 'with_new_edges')

    print(f"\n{'='*80}")
    print("ALL EXPORTS COMPLETE")
    print(f"{'='*80}\n")


if __name__ == '__main__':
    main()
