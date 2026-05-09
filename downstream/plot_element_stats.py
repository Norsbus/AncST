#! /usr/bin/env python3

import pickle
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_element_stats.py <output_dir>")
        sys.exit(1)

    outd = sys.argv[1]
    stats_file = f'{outd}/genetic_elements_stats/element_statistics.pickle'

    if not os.path.exists(stats_file):
        print(f"ERROR: {stats_file} not found")
        print("Run compute_element_stats.py first")
        sys.exit(1)

    print("Loading element statistics...")
    with open(stats_file, 'rb') as f:
        data = pickle.load(f)

    element_stats = data['element_stats']
    per_species_stats = data['per_species_stats']
    species_list = data['species_list']

    # Create output directory
    plot_dir = f'{outd}/genetic_elements_stats/plots'
    os.makedirs(plot_dir, exist_ok=True)

    # Choose color palette based on number of species
    nspe = len(species_list)
    colormap = 'tab20' if nspe > 10 else 'tab10'
    colors = plt.colormaps[colormap](np.linspace(0, 1, nspe))

    # Plot 1: Global element classification
    print("\nGenerating global element classification plot...")
    total_elements = sum(per_species_stats[sp]['total_elements'] for sp in species_list)
    num_singletons = len(element_stats['singletons'])
    num_single_copy = sum(len(g) for g in element_stats['single_copy_groups'])
    num_multi_copy = sum(len(g) for g in element_stats['multi_copy_groups'])

    fig, ax = plt.subplots(figsize=(10, 6))
    categories = ['Singletons', 'Single-copy\ngroups', 'Multi-copy\ngroups']
    counts = [num_singletons, num_single_copy, num_multi_copy]
    bars = ax.bar(categories, counts, alpha=0.7, edgecolor='black')
    bars[0].set_color('lightgray')
    bars[1].set_color('skyblue')
    bars[2].set_color('salmon')

    ax.set_ylabel('Number of elements', fontsize=12)
    ax.set_title('Global element classification', fontsize=14)
    ax.grid(axis='y', alpha=0.3)

    # Add counts on top of bars
    for i, v in enumerate(counts):
        ax.text(i, v, str(v), ha='center', va='bottom', fontsize=11)

    plt.tight_layout()
    plt.savefig(f'{plot_dir}/global_element_classification.svg', format='svg', bbox_inches='tight')
    plt.savefig(f'{plot_dir}/global_element_classification.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved global classification plot")

    # Plot 2: Per-species element classification (stacked bar)
    print("\nGenerating per-species stacked bar plot...")
    fig, ax = plt.subplots(figsize=(14, 8))

    singletons = [per_species_stats[sp]['singletons'] for sp in species_list]
    single_copy = [per_species_stats[sp]['in_single_copy_groups'] for sp in species_list]
    multi_copy = [per_species_stats[sp]['in_multi_copy_groups'] for sp in species_list]

    x = np.arange(len(species_list))
    width = 0.6

    ax.bar(x, singletons, width, label='Singletons', alpha=0.8, color='lightgray', edgecolor='black')
    ax.bar(x, single_copy, width, bottom=singletons, label='Single-copy groups', alpha=0.8, color='skyblue', edgecolor='black')
    bottom = np.array(singletons) + np.array(single_copy)
    ax.bar(x, multi_copy, width, bottom=bottom, label='Multi-copy groups', alpha=0.8, color='salmon', edgecolor='black')

    ax.set_ylabel('Number of elements', fontsize=12)
    ax.set_xlabel('Species', fontsize=12)
    ax.set_title('Per-species element classification', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels([sp.replace('GCF_', '') for sp in species_list], rotation=45, ha='right', fontsize=9)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{plot_dir}/per_species_classification_stacked.svg', format='svg', bbox_inches='tight')
    plt.savefig(f'{plot_dir}/per_species_classification_stacked.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved per-species stacked bar plot")

    # Plot 3: Per-species subplots showing category breakdown
    print("\nGenerating per-species subplot grid...")
    ncols = 5
    nrows = (nspe + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, nrows * 4))
    axes = axes.flatten()

    for i, sp in enumerate(species_list):
        stats = per_species_stats[sp]
        categories = ['Singleton', 'Single-copy', 'Multi-copy']
        counts = [stats['singletons'], stats['in_single_copy_groups'], stats['in_multi_copy_groups']]

        bars = axes[i].bar(categories, counts, alpha=0.8, edgecolor='black')
        bars[0].set_color('lightgray')
        bars[1].set_color('skyblue')
        bars[2].set_color('salmon')

        axes[i].set_title(sp.replace('GCF_', ''), fontsize=10)
        axes[i].set_ylabel('Count', fontsize=9)
        axes[i].tick_params(axis='x', labelsize=8, rotation=45)
        axes[i].grid(axis='y', alpha=0.3)

    # Hide unused subplots
    for i in range(len(species_list), len(axes)):
        axes[i].axis('off')

    plt.suptitle('Per-species element classification breakdown', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/per_species_classification_grid.svg', format='svg', bbox_inches='tight')
    plt.savefig(f'{plot_dir}/per_species_classification_grid.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved per-species grid plot")

    # Plot 4: Species connectivity histogram (how many species does each element connect to)
    print("\nGenerating species connectivity histograms...")
    ncols = 5
    nrows = (nspe + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, nrows * 4))
    axes = axes.flatten()

    for i, sp in enumerate(species_list):
        stats = per_species_stats[sp]
        connectivity = stats['num_species_connected']

        if connectivity:
            # Get x-axis values (number of species connected to)
            max_connected = max(connectivity.keys())
            x_vals = list(range(max_connected + 1))
            y_vals = [connectivity.get(x, 0) for x in x_vals]

            axes[i].bar(x_vals, y_vals, alpha=0.7, edgecolor='black', color=colors[i])
            axes[i].set_title(sp.replace('GCF_', ''), fontsize=10)
            axes[i].set_xlabel('# other species', fontsize=9)
            axes[i].set_ylabel('# elements', fontsize=9)
            axes[i].grid(axis='y', alpha=0.3)
        else:
            axes[i].text(0.5, 0.5, 'No data', ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_title(sp.replace('GCF_', ''), fontsize=10)

    # Hide unused subplots
    for i in range(len(species_list), len(axes)):
        axes[i].axis('off')

    plt.suptitle('Per-species connectivity: How many other species does each element match?', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/per_species_connectivity_histogram.svg', format='svg', bbox_inches='tight')
    plt.savefig(f'{plot_dir}/per_species_connectivity_histogram.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved species connectivity histogram")

    # Plot 5: Top ortholog partners (bar chart instead of pie chart)
    print("\nGenerating top ortholog partners plot...")
    ncols = 5
    nrows = (nspe + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, nrows * 5))
    axes = axes.flatten()

    for i, sp in enumerate(species_list):
        stats = per_species_stats[sp]
        connectivity = stats['species_connectivity']

        if connectivity:
            # Get top 5 partners
            top_partners = connectivity.most_common(5)
            partner_names = [p.replace('GCF_', '') for p, _ in top_partners]
            partner_counts = [c for _, c in top_partners]

            y_pos = np.arange(len(partner_names))
            axes[i].barh(y_pos, partner_counts, alpha=0.7, edgecolor='black', color=colors[i])
            axes[i].set_yticks(y_pos)
            axes[i].set_yticklabels(partner_names, fontsize=8)
            axes[i].set_xlabel('# edges', fontsize=9)
            axes[i].set_title(sp.replace('GCF_', ''), fontsize=10)
            axes[i].grid(axis='x', alpha=0.3)
            axes[i].invert_yaxis()
        else:
            axes[i].text(0.5, 0.5, 'No partners', ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_title(sp.replace('GCF_', ''), fontsize=10)

    # Hide unused subplots
    for i in range(len(species_list), len(axes)):
        axes[i].axis('off')

    plt.suptitle('Per-species top ortholog partners (by edge count)', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/per_species_top_partners.svg', format='svg', bbox_inches='tight')
    plt.savefig(f'{plot_dir}/per_species_top_partners.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved top partners plot")

    print(f"\nAll plots saved to {plot_dir}/")

if __name__ == '__main__':
    main()
