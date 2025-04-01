#! /usr/bin/env python3

from subprocess import run

run('rm -rf out && mkdir -p out out/utils out/custom out/MCScanX out/dialign out/Smore out/i-ADHoRe && cp README requirements.txt out',shell=True)
run('cp -r README_Smore make_smore_genes_and_prep_folder.py smore_anchors out/Smore',shell=True)
run('cp -r README_MCScanX MCScanX.gff MCScanX.homology put_coords_into_MCScanX_gff_and_homology.py extract_MCScanX_orthologies_clusters.py extract_MCScanX_orthologies_pairwise.py out/MCScanX',shell=True)
run('cp -r README_utils replace_names.py mapping get_mapping.py clasp.x find_pot_homologs.py orgs ../anchors ../utils/small_meta out/utils',shell=True)
run('cp -r README_dialign dialign.homology make_anc.py out/dialign',shell=True)
run('cp -r README_custom draw_simple.py draw_verbose.py pairwise_alignments_table get_syn_regions.py out/custom',shell=True)
run('cp -r README_iadhore adhore_config adhore_pairwise.table adhore_gene_lists get_pairwise_orthologies_of_custom_elements.py put_coords_into_gene_lists.py out/i-ADHoRe',shell=True)
run('chmod -R 755 out',shell=True)
