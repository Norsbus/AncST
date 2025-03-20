#! /usr/bin/env python3

from subprocess import run
from sys import argv

run(f'mkdir -p singles_out multi_out',shell=True)
refs = []
with open('refs.txt','r') as f:
    for line in f:
        refs.append(line.strip())
refs_list_str = str(refs[0])
for ref in refs[1:]:
    refs_list_str += f',{ref}'

target = argv[1]

for ref in refs:
    run(f'./single_scaff.py {ref} {target} 1>singles_out/out_single_scaff_ref_{ref}.txt',shell = True)

run(f'./merge_singles_for_blossom.py --refs {refs_list_str} --target {target}',shell = True)

lines = []
no_e = 0
with open(f'singles_out/blossom_v_infile_merged.txt') as f:
    first = 1
    for line in f:
        if first == 1:
            no_v = int(line.strip().split()[0])
            first = 0
            continue
        no_e += 1
        lines.append(line)

with open(f'singles_out/blossom5_in.txt','w') as f:
    f.write(f'{no_v} {no_e}\n')
    f.writelines(lines)

run(f'blossom5-v2.05.src/blossom5 -e singles_out/blossom5_in.txt -w blossom5-v2.05.src/out',shell=True)

run(f'./del_min_cycles_and_make_multi_ref_scaffolds.py {target}',shell = True)
