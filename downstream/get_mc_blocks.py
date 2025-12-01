#! /usr/bin/env python3

import pickle
from sys import argv

borders = []
first = 1
with open(f'{argv[1]}','r') as f:
    for line in f:
        if line[0] == '#' and 'Alignment' not in line:
            continue
        elif 'Alignment' in line:
            if first == 1:
                first = 0
                matches = []
                maxi1 = [0,'']
                mini1 = [1e16,'']
                maxi2 = [0,'']
                mini2 = [1e16,'']
                continue
            else:
                if (int(matches[0][0].split('to')[1]) > int(matches[-1][0].split('to')[1])) or (int(matches[0][1].split('to')[1]) > int(matches[-1][1].split('to')[1])):
                    rev = 1
                else:
                    rev = 0
                first1 = mini1[1]
                last1 = maxi1[1]
                first2 = mini2[1]
                last2 = maxi2[1]
                borders.append(((first1,last1),(first2,last2),rev))
                matches = []
                maxi1 = [0,'']
                mini1 = [1e16,'']
                maxi2 = [0,'']
                mini2 = [1e16,'']
        else:
            line = line.split(':')[1].split()
            first = line[0].strip()
            second = line[1].strip()
            matches.append((first,second))
            start1 = int(first.split('ele')[1].split('to')[0])
            end1 = int(first.split('to')[1])
            start2 = int(second.split('ele')[1].split('to')[0])
            end2 = int(second.split('to')[1])
            if start1 < mini1[0]:
                mini1 = [start1,first]
            if end1 > maxi1[0]:
                maxi1 = [end1,first]
            if start2 < mini2[0]:
                mini2 = [start2,second]
            if end2 > maxi2[0]:
                maxi2 = [end2,second]

with open('naive_blocks.pickle','wb') as f:
    pickle.dump(borders,f)
