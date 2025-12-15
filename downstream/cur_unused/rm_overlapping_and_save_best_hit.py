#! /usr/bin/env python3

from sys import argv
from pprint import pprint

ranges = {}
with open(argv[1]) as f:
    for line in f:
        line = line.strip().split('\t')
        if len(line) == 8:
            org,chromo,start,end,orientation,hit,best_hits,score = line
            hit = f'{org}{chromo}{start}{end}'
        elif len(line) == 7:
            org,chromo,start,end,orientation,hit,score = line
        else:
            input('ERROR1')
        start = int(start)
        end = int(end)
        if start > end:
            start,end = end,start
        if org not in ranges:
            ranges[org] = {}
        if chromo not in ranges[org]:
            ranges[org][chromo] = []
        best_hit = hit.split('_')[-1]
        ranges[org][chromo].append((int(start),int(end),chromo,orientation,best_hit,hit,score))

new_ranges = {}
hits = {}
for org,bib in ranges.items():
    new_ranges[org] = {}
    for chromo,rang in bib.items():
        rang = sorted(list(set(rang)))
        to_del = []
        for i,r1 in enumerate(rang[:-1]):
            hits[r1] = f'{r1[4]}:{r1[6]}'
            for r2 in rang[i+1:]:
                if r2[0] > r1[1]:
                    hits[r2] = r2[4]
                    break
                else:
                    hits[r1] = hits[r1] + f'${r2[4]}:{r2[6]}'
                    hits[r2] = hits[r1]
                    if r1[1] - r1[0] > r2[1] -r2[0]:
                        to_del.append(r2)
                    else:
                        to_del.append(r1)

        if len(rang) == 1:
            hits[rang[0]] = f'{rang[0][4]}:{rang[0][6]}'

        new_ranges[org][chromo] = rang.copy()
        for r in set(to_del):
            new_ranges[org][chromo].remove(r)

ranges = new_ranges
with open(argv[2],'w+') as f:
    for org,bib in ranges.items():
        for chromo,rang in bib.items():
            for ele in rang:
                start,end,chromo,orientation,best_hit,hit,score = ele
                best_hit = hits[ele]
                f.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\t{hit}\t{best_hit}\t{score}\n')
