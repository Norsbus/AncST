#! /usr/bin/env python3

from subprocess import run

#run('rm -rf out && mkdir -p out out/utils out/custom out/MCScanX out/dialign out/Smore out/i-ADHoRe && cp README requirements.txt out',shell=True)
run(f'mkdir -p out_chr12_test4 out_chr12_test4/utils out_chr12_test4/custom out_chr12_test4/MCScanX out_chr12_test4/dialign out_chr12_test4/Smore out_chr12_test4/i-ADHoRe',shell=True)
run('cp -r README_Smore make_smore_genes_and_prep_folder.py smore_anchors out_chr12_test4/Smore',shell=True)
run('cp -r README_MCScanX MCScanX.gff MCScanX.homology put_coords_into_MCScanX_gff_and_homology.py extract_MCScanX_orthologies.py out_chr12_test4/MCScanX',shell=True)
run('cp -r README_utils replace_names.py mapping get_mapping.py clasp.x find_pot_homologs.py out_chr12_test4/utils',shell=True)
run('cp -r README_dialign dialign.homology make_anc.py out_chr12_test4/dialign',shell=True)
run('cp -r README_custom draw_synteny.py out_chr12_test4/custom',shell=True)
run('cp -r README_iadhore adhore_config adhore_families.table adhore_pairwise.table adhore_gene_lists out_chr12_test4/i-ADHoRe',shell=True)
run('chmod -R 755 out_chr12_test4',shell=True)
