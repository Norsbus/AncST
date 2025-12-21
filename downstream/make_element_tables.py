#! /usr/bin/env python3

import pickle
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 make_element_tables.py <output_dir>")
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
    table_dir = f'{outd}/genetic_elements_stats/tables'
    os.makedirs(table_dir, exist_ok=True)

    # Table 1: Global summary
    global_table = f'{table_dir}/global_element_summary.tex'
    total_elements = sum(per_species_stats[sp]['total_elements'] for sp in species_list)
    num_singletons = len(element_stats['singletons'])
    num_single_copy = sum(len(g) for g in element_stats['single_copy_groups'])
    num_multi_copy = sum(len(g) for g in element_stats['multi_copy_groups'])
    num_self_dup = len(element_stats['within_species_duplicates'])

    with open(global_table, 'w') as f:
        f.write("\\documentclass[11pt,a4paper]{article}\n")
        f.write("\\usepackage[utf8]{inputenc}\n")
        f.write("\\usepackage{booktabs}\n")
        f.write("\\usepackage[margin=1in]{geometry}\n")
        f.write("\\begin{document}\n")
        f.write("\\begin{table}[h]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{lr}\n")
        f.write("\\toprule\n")
        f.write("Category & Count \\\\\n")
        f.write("\\midrule\n")
        f.write(f"Total elements & {total_elements} \\\\\n")
        f.write(f"Singletons & {num_singletons} \\\\\n")
        f.write(f"In single-copy groups & {num_single_copy} \\\\\n")
        f.write(f"In multi-copy groups & {num_multi_copy} \\\\\n")
        f.write(f"Within-species duplicates* & {num_self_dup} \\\\\n")
        f.write("\\midrule\n")
        f.write(f"Single-copy ortholog groups & {len(element_stats['single_copy_groups'])} \\\\\n")
        f.write(f"Multi-copy ortholog groups & {len(element_stats['multi_copy_groups'])} \\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{Global element classification summary}\n")
        f.write("\\end{table}\n")
        f.write("\\end{document}\n")

    print(f"Saved global summary table: {global_table}")

    # Table 2: Per-species summary
    per_species_table = f'{table_dir}/per_species_element_summary.tex'
    with open(per_species_table, 'w') as f:
        f.write("\\documentclass[11pt,a4paper]{article}\n")
        f.write("\\usepackage[utf8]{inputenc}\n")
        f.write("\\usepackage{booktabs}\n")
        f.write("\\usepackage[margin=1in]{geometry}\n")
        f.write("\\begin{document}\n")
        f.write("\\begin{table}[h]\n")
        f.write("\\centering\n")
        f.write("\\small\n")
        f.write("\\begin{tabular}{lrrrrrr}\n")
        f.write("\\toprule\n")
        f.write("Species & Input & Filtered & Final & Singletons & Single-copy & Multi-copy \\\\\n")
        f.write("\\midrule\n")

        for sp in species_list:
            stats = per_species_stats[sp]
            sp_short = sp.replace('GCF_', '')
            total_input = stats['total_elements'] + stats.get('filtered', 0)
            filtered = stats.get('filtered', 0)
            f.write(f"{sp_short} & {total_input} & {filtered} & {stats['total_elements']} & {stats['singletons']} & ")
            f.write(f"{stats['in_single_copy_groups']} & {stats['in_multi_copy_groups']} \\\\\n")

        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{Per-species element classification. Input = elements in input data, Filtered = removed during pipeline, Final = in final orthology graph}\n")
        f.write("\\end{table}\n")
        f.write("\\end{document}\n")

    print(f"Saved per-species summary table: {per_species_table}")

    # Table 3: Top ortholog partners per species
    partner_table = f'{table_dir}/per_species_top_partners.tex'
    with open(partner_table, 'w') as f:
        f.write("\\documentclass[11pt,a4paper]{article}\n")
        f.write("\\usepackage[utf8]{inputenc}\n")
        f.write("\\usepackage{booktabs}\n")
        f.write("\\usepackage[margin=1in]{geometry}\n")
        f.write("\\begin{document}\n")
        f.write("\\begin{table}[h]\n")
        f.write("\\centering\n")
        f.write("\\small\n")
        f.write("\\begin{tabular}{llr}\n")
        f.write("\\toprule\n")
        f.write("Species & Top Partner & Edges \\\\\n")
        f.write("\\midrule\n")

        for sp in species_list:
            stats = per_species_stats[sp]
            sp_short = sp.replace('GCF_', '')
            if stats['species_connectivity']:
                top_partner, top_count = stats['species_connectivity'].most_common(1)[0]
                partner_short = top_partner.replace('GCF_', '')
                f.write(f"{sp_short} & {partner_short} & {top_count} \\\\\n")
            else:
                f.write(f"{sp_short} & -- & 0 \\\\\n")

        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{Top ortholog partner per species (by edge count)}\n")
        f.write("\\end{table}\n")
        f.write("\\end{document}\n")

    print(f"Saved top partners table: {partner_table}")

    print(f"\nAll tables saved to {table_dir}/")

if __name__ == '__main__':
    main()
