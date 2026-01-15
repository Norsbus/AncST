#! /usr/bin/env python3

from subprocess import run
from sys import argv

ad = argv[1]

run('rm -rf out && mkdir -p out out/experimental out/stable out/stable/synthology out/utils out/experimental/custom out/stable/MCScanX out/experimental/dialign out/experimental/Smore out/experimental/i-ADHoRe out/scripts_which_produced_these_results && cp README requirements.txt out',shell=True)
run('cp -r ../scaffolder out/',shell=True)
run('cp -r get_syn_regions.py run_synthology_nucl.py run_synthology_prot.py compute_ortho_stats.py make_ortho_tables.py plot_ortho_stats.py run_blast_on_regions.py blast_results_to_synthology.py extract_syntenic_fastas.py draw_verbose.py draw_synorthogroups.py generate_plotting_order.py export_synorthogroups.py refine_synorthogroups_with_phylogeny.py SYNORTHOGROUPS_README.md README_synthology compute_element_stats.py make_element_tables.py plot_element_stats.py out/stable/synthology',shell=True)
run('cp -r README_Smore make_smore_genes_and_prep_folder.py smore_anchors out/experimental/Smore',shell=True)
run('cp -r README_MCScanX orgs MCScanX.gff MCScanX.homology put_coords_into_MCScanX_gff_and_homology.py extract_MCScanX_orthologies_clusters.py extract_MCScanX_orthologies_pairwise.py phylogeny_from_mcscanx.py get_orders.py get_mc_blocks.py sc_mcdraw.py mcdraw.py convert_ncbi_idents.py root_NJ_tree_at_midpoint.py out/stable/MCScanX',shell=True)
run(f'cp -r README_utils pairwise_alignments_table replace_names.py mapping get_mapping.py clasp.x orgs compressed_maps_multis_to_one {ad} ../utils/small_meta NJTree.nwk UPGMATree.nwk convert_ncbi_idents.py root_NJ_tree_at_midpoint.py out/utils',shell=True)
run('cp -r README_dialign dialign.homology make_anc.py out/experimental/dialign',shell=True)
run('cp -r README_custom draw_simple.py draw_verbose.py find_pot_homologs.py generate_plotting_order.py out/experimental/custom',shell=True)
run('cp -r README_iadhore adhore_config adhore_pairwise.table adhore_gene_lists get_pairwise_orthologies_of_custom_elements.py put_coords_into_gene_lists.py out/experimental/i-ADHoRe',shell=True)
run('cp -r README_scripts make_gff3.py eval_dups_only_syn.py add_dups_to_anchors.py make_dup_free_ams.py make_chr_and_org_mapping_for_gff.py compress_maps_and_ignore_multis_with_dups.py get_gff_pairwise_succinct.py make_pw_table_from_aligned.py make_clusters.py phylogeny_from_anchors.py out/scripts_which_produced_these_results',shell=True)
run(f'cp -r gff out && cp README_gff out/gff/',shell=True)
run('chmod -R 755 out',shell=True)
