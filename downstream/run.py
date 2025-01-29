#!/usr/bin/env python3

from subprocess import run
from sys import argv

anchor_dir = argv[1]
run('./make_directories.py',shell=True)
run(f'./get_all_orgs.py {anchor_dir}',shell=True)
run(f'./make_chr_and_org_mapping_for_gff.py {anchor_dir}',shell=True)
run(f'./run_per_org.py {anchor_dir}',shell=True)
run('./make_clusters.py',shell=True)
run(f'./get_gff_and_homology_with_adhore.py',shell=True)
run(f'./write_clusters.py',shell=True)
run(f'./make_gap_free_clusters.py',shell=True)
run('./cp.py',shell=True)
