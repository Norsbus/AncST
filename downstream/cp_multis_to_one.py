#! /usr/bin/env python3

from subprocess import run

run('rm -rf out && mkdir -p out out/scaffolder out/utils out/custom out/MCScanX out/dialign out/Smore out/i-ADHoRe out/scripts_which_produced_these_results && cp README requirements.txt out',shell=True)
run('cp -r ../scaffolder out/scaffolder',shell=True)
run('cp -r README_Smore make_smore_genes_and_prep_folder.py smore_anchors out/Smore',shell=True)
run('cp -r README_MCScanX orgs MCScanX.gff MCScanX.homology put_coords_into_MCScanX_gff_and_homology.py extract_MCScanX_orthologies_clusters.py extract_MCScanX_orthologies_pairwise.py phylogeny_from_mcscanx.py get_orders.py get_mc_blocks.py mcdraw.py out/MCScanX',shell=True)
run('cp -r README_utils replace_names.py mapping get_mapping.py clasp.x find_pot_homologs.py orgs compressed_maps_multis_to_one/ ../anchors ../utils/small_meta NJTree.nwk UPGMATree.nwk out/utils',shell=True)
run('cp -r README_dialign dialign.homology make_anc.py out/dialign',shell=True)
run('cp -r README_custom draw_simple.py draw_verbose.py pairwise_alignments_table get_syn_regions.py out/custom',shell=True)
run('cp -r README_iadhore adhore_config adhore_pairwise.table adhore_gene_lists get_pairwise_orthologies_of_custom_elements.py put_coords_into_gene_lists.py out/i-ADHoRe',shell=True)
run('cp -r README_scripts make_gff3.py eval_dups_only_syn.py add_dups_to_anchors.py make_dup_free_ams.py make_chr_and_org_mapping_for_gff.py compress_maps_and_ignore_multis_with_dups.py get_gff_pairwise_succinct.py make_pw_table_from_aligned.py make_clusters.py phylogeny_from_anchors.py out/scripts_which_produced_these_results',shell=True)
run(f'cp -r README_gff gff out',shell=True)
run('chmod -R 755 out',shell=True)
