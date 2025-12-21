#! /usr/bin/env python3

import pickle
import sys
import os
import networkx as nx
from collections import defaultdict, Counter

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 compute_element_stats.py <output_dir>")
        sys.exit(1)

    outd = sys.argv[1]

    # Load the final union graph (no new edges version)
    union_file = f'{outd}/final_union_graph_no_new_edges.pickle'
    if not os.path.exists(union_file):
        print(f"ERROR: {union_file} not found")
        sys.exit(1)

    print("Loading final union graph...")
    with open(union_file, 'rb') as f:
        data = pickle.load(f)
        if len(data) == 3:
            graph, gene_species, metadata = data
            print(f"Threshold: {metadata.get('threshold', 'N/A')}")
        else:
            graph, gene_species = data

    print(f"Loaded graph with {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges")
    print(f"Gene-species mapping: {len(gene_species)} genes")

    # Identify genes not in graph (filtered during pipeline)
    genes_in_graph = set(graph.nodes())
    genes_not_in_graph = set(gene_species.keys()) - genes_in_graph
    print(f"Genes filtered out during pipeline: {len(genes_not_in_graph)}")

    # Load filtering report if available
    filter_report_file = f'{outd}/genetic_elements_stats/filtering_report.pickle'
    filtering_reasons = {}
    if os.path.exists(filter_report_file):
        print(f"Loading filtering analysis from filtering_report.pickle...")
        with open(filter_report_file, 'rb') as f:
            filter_data = pickle.load(f)
            gene_status_all = filter_data.get('gene_status', {})

            # Categorize filtered genes by reason
            for gene in genes_not_in_graph:
                status = gene_status_all.get(gene, {})
                if not status.get('in_chains', False):
                    filtering_reasons[gene] = 'not_in_chains'
                elif not status.get('homology', False):
                    filtering_reasons[gene] = 'no_homology'
                elif not status.get('in_pairwise', False):
                    filtering_reasons[gene] = 'filtered_at_pairwise'
                elif not status.get('in_aligned', False):
                    filtering_reasons[gene] = 'filtered_at_alignment'
                elif not status.get('in_union_new', False):
                    filtering_reasons[gene] = 'filtered_at_union'
                else:
                    filtering_reasons[gene] = 'filtered_at_final_cograph'

        # Count by reason
        reason_counts = Counter(filtering_reasons.values())
        print(f"  Filtering reasons:")
        for reason, count in reason_counts.most_common():
            print(f"    - {reason}: {count} genes")
    else:
        print(f"  Run trace_filtered_genes.py first for detailed filtering analysis")

    # Get all species
    species_list = sorted(set(gene_species.values()))
    print(f"Species: {len(species_list)}")

    # Find connected components
    components = list(nx.connected_components(graph))
    print(f"Connected components: {len(components)}")

    # Classify elements
    element_stats = {
        'singletons': [],
        'single_copy_groups': [],
        'multi_copy_groups': [],
        'within_species_duplicates': [],
        'filtered_genes': list(genes_not_in_graph)
    }

    per_species_stats = {sp: {
        'total_elements': 0,  # Only genes in final graph
        'filtered': 0,  # Genes removed during pipeline
        'filtered_reasons': Counter(),  # Breakdown of filtering reasons
        'singletons': 0,
        'in_single_copy_groups': 0,
        'in_multi_copy_groups': 0,
        'self_duplicated': 0,
        'species_connectivity': Counter(),  # Count edges to other species
        'num_species_connected': Counter()  # Count how many elements connect to N species
    } for sp in species_list}

    # Count total elements per species (only those in graph)
    for gene, sp in gene_species.items():
        if gene in genes_in_graph:
            per_species_stats[sp]['total_elements'] += 1
        else:
            per_species_stats[sp]['filtered'] += 1
            # Track filtering reason if available
            if gene in filtering_reasons:
                per_species_stats[sp]['filtered_reasons'][filtering_reasons[gene]] += 1

    # Analyze each component
    for comp in components:
        comp_genes = list(comp)
        comp_size = len(comp_genes)

        # Get species composition
        comp_species_counts = Counter([gene_species[g] for g in comp_genes])
        num_species = len(comp_species_counts)

        # Classify component type
        is_singleton = comp_size == 1
        is_single_copy = all(count == 1 for count in comp_species_counts.values())
        is_multi_copy = any(count > 1 for count in comp_species_counts.values())

        if is_singleton:
            gene = comp_genes[0]
            sp = gene_species[gene]
            element_stats['singletons'].append(gene)
            per_species_stats[sp]['singletons'] += 1
        elif is_single_copy:
            element_stats['single_copy_groups'].append(comp_genes)
            for gene in comp_genes:
                sp = gene_species[gene]
                per_species_stats[sp]['in_single_copy_groups'] += 1
        else:  # multi-copy
            element_stats['multi_copy_groups'].append(comp_genes)
            for gene in comp_genes:
                sp = gene_species[gene]
                per_species_stats[sp]['in_multi_copy_groups'] += 1

                # Check if this gene itself is duplicated within species
                same_species_in_comp = [g for g in comp_genes if gene_species[g] == sp]
                if len(same_species_in_comp) > 1:
                    per_species_stats[sp]['self_duplicated'] += 1
                    if gene not in element_stats['within_species_duplicates']:
                        element_stats['within_species_duplicates'].append(gene)

        # Analyze connectivity for each gene in component
        for gene in comp_genes:
            sp = gene_species[gene]

            # Count connections to other species
            species_connected = set()
            for neighbor in graph.neighbors(gene):
                neighbor_sp = gene_species[neighbor]
                if neighbor_sp != sp:
                    per_species_stats[sp]['species_connectivity'][neighbor_sp] += 1
                    species_connected.add(neighbor_sp)

            # Track how many species this element connects to
            num_connected_species = len(species_connected)
            per_species_stats[sp]['num_species_connected'][num_connected_species] += 1

    # Print global statistics
    total_in_graph = graph.number_of_nodes()
    print("\n" + "="*60)
    print("GLOBAL ELEMENT STATISTICS")
    print("="*60)
    print(f"Total elements in input: {len(gene_species)}")
    print(f"Filtered during pipeline: {len(genes_not_in_graph)} ({100*len(genes_not_in_graph)/len(gene_species):.1f}%)")
    print(f"In final graph: {total_in_graph} ({100*total_in_graph/len(gene_species):.1f}%)")
    print(f"\nClassification of {total_in_graph} elements in final graph:")
    print(f"  Singletons: {len(element_stats['singletons'])} ({100*len(element_stats['singletons'])/total_in_graph:.1f}%)")
    print(f"  In single-copy groups: {sum(len(g) for g in element_stats['single_copy_groups'])} ({100*sum(len(g) for g in element_stats['single_copy_groups'])/total_in_graph:.1f}%)")
    print(f"  In multi-copy groups: {sum(len(g) for g in element_stats['multi_copy_groups'])} ({100*sum(len(g) for g in element_stats['multi_copy_groups'])/total_in_graph:.1f}%)")
    print(f"  Within-species duplicates*: {len(element_stats['within_species_duplicates'])}")

    print(f"\nOrthology groups:")
    print(f"  Single-copy groups: {len(element_stats['single_copy_groups'])}")
    print(f"  Multi-copy groups: {len(element_stats['multi_copy_groups'])}")

    # Print per-species statistics
    print("\n" + "="*60)
    print("PER-SPECIES ELEMENT STATISTICS")
    print("="*60)

    for sp in species_list:
        stats = per_species_stats[sp]
        total_with_filtered = stats['total_elements'] + stats['filtered']
        print(f"\n{sp}:")
        print(f"  Total in input: {total_with_filtered}")
        if stats['filtered'] > 0:
            print(f"  Filtered: {stats['filtered']}")
            if stats['filtered_reasons']:
                for reason, count in stats['filtered_reasons'].most_common():
                    print(f"    - {reason}: {count}")
        print(f"  In final graph: {stats['total_elements']}")
        if stats['total_elements'] > 0:
            print(f"    Singletons: {stats['singletons']} ({100*stats['singletons']/stats['total_elements']:.1f}%)")
            print(f"    In single-copy groups: {stats['in_single_copy_groups']} ({100*stats['in_single_copy_groups']/stats['total_elements']:.1f}%)")
            print(f"    In multi-copy groups: {stats['in_multi_copy_groups']} ({100*stats['in_multi_copy_groups']/stats['total_elements']:.1f}%)")
            if stats['self_duplicated'] > 0:
                print(f"      - Self-duplicated*: {stats['self_duplicated']}")

        # Top species connections
        if stats['species_connectivity']:
            top_partners = stats['species_connectivity'].most_common(3)
            print(f"  Top ortholog partners:")
            for partner_sp, count in top_partners:
                print(f"    - {partner_sp}: {count} edges")

    # Save results for plotting/tables
    results = {
        'element_stats': element_stats,
        'per_species_stats': per_species_stats,
        'species_list': species_list,
        'graph': graph,
        'gene_species': gene_species
    }

    stats_dir = f'{outd}/genetic_elements_stats'
    os.makedirs(stats_dir, exist_ok=True)
    output_file = f'{stats_dir}/element_statistics.pickle'
    with open(output_file, 'wb') as f:
        pickle.dump(results, f)
    print(f"\nSaved element statistics to {output_file}")

if __name__ == '__main__':
    main()
