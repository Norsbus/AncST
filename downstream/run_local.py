#!/usr/bin/env python3

from subprocess import run
from sys import argv

anchor_dir = argv[1]
run('chmod +x *py',shell=True)
run('./make_directories.py',shell=True)
run(f'./get_all_orgs.py {anchor_dir}',shell=True)
run(f'./make_chr_and_org_mapping_for_gff.py {anchor_dir}',shell=True)
run(f'./run_per_org_local.py {anchor_dir}',shell=True)
run(f'./get_gff_pairwise_succinct.py',shell=True)
run(f'./make_pw_table_from_aligned.py {anchor_dir}',shell=True)
run(f'./make_gff3.py {anchor_dir}',shell=True)
run(f'./phylogeny_from_anchors.py',shell=True)
run(f'./replace_names.py UPGMATree.nwk UPGMATree.nwk && replace_names.py NJTree.nwk NJTree.nwk',shell=True)
run('./cp_multis_to_one.py',shell=True)
