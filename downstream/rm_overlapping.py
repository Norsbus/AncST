#! /usr/bin/env python3

from sys import argv

ranges = {}
with open(argv[1]) as f:
    for line in f:
        line = line.split()
        if len(line) == 5:
            org,chromo,start,end,orientation = line
            hit = f'{org}{chromo}{start}{end}'
        elif len(line) == 6:
            org,chromo,start,end,orientation,hit = line
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
        ranges[org][chromo].append((int(start),int(end),chromo,orientation,hit))

new_ranges = {}
for org,bib in ranges.items():
    new_ranges[org] = {}
    for chromo,rang in bib.items():
        rang = sorted(list(set(rang)))
        to_del = set()
        for i,r1 in enumerate(rang[:-1]):
            for r2 in rang[i+1:]:
                if r2[0] > r1[1]:
                    break
                else:
                    if r1[1] - r1[0] > r2[1] -r2[0]:
                        to_del.add(r2)
                    else:
                        to_del.add(r1)

        new_ranges[org][chromo] = rang.copy()
        for r in to_del:
            new_ranges[org][chromo].remove(r)
    
ranges = new_ranges
with open(argv[2],'w+') as f:
    for org,bib in ranges.items():
        for chromo,rang in bib.items():
            for ele in rang:
                start,end,chromo,orientation,hit = ele
                f.write(f'{org}\t{chromo}\t{start}\t{end}\t{orientation}\t{hit}\n')
