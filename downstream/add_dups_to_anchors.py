#! /usr/bin/env python3

import pickle
import os
from pprint import pprint
from bisect import bisect_left
from copy import deepcopy

def overlaps(x,y):
    if x[0] >= y[0] and x[0] <= y[1]:
        return True
    if x[1] >= y[0] and x[1] <= y[1]:
        return True
    if y[0] >= x[0] and y[0] <= x[1]:
        return True
    if y[1] >= x[0] and y[1] <= x[1]:
        return True
    return False

if __name__ == "__main__":

    orgs = []
    with open('orgs') as f:
        for line in f:
            orgs.append(line.strip())

    for org in orgs:
        print(f'======= {org} =======')

        with open('../anchors/candidates' + f'/{org}','rb') as f:
            bib = pickle.load(f)

        print(len(bib))

        ses = sorted([(i,i+x['end']-x['start']) for i,x in bib.items()])
        starts = [x[0] for x in ses]
        ends = [x[1] for x in ses]

        # Check if dups file exists for this organism
        dups_file = f'dups/{org}'
        if os.path.exists(dups_file):
            with open(dups_file,'rb') as f:
                dups = pickle.load(f)

            to_del = []
            for i,d in dups.items():
                end = i + d['end']-d['start']

                end_bigger_start_of_dup = bisect_left(ends,i)
                start_smaller_end_of_dup = bisect_left(starts,end)

                new_entry = {}
                new_entry['i'] = i
                new_entry['bib'] = deepcopy(d)

                cur_start = i
                cur_end = end
                for idx in range(max(0,end_bigger_start_of_dup-3),min(len(ends),start_smaller_end_of_dup+4)):
                    if overlaps((starts[idx],ends[idx]),(cur_start,cur_end)):
                        if cur_end >= starts[idx] and cur_end - starts[idx] < (ends[idx] - starts[idx])/10:
                            new_entry['bib']['end'] -= (cur_end - starts[idx] + 100)
                            cur_end = cur_end - (cur_end-starts[idx] + 100)
                        elif cur_start <= ends[idx] and ends[idx] - cur_start < (ends[idx] - starts[idx])/10:
                            new_entry['i'] = cur_start + (ends[idx] - cur_start + 100)
                            new_entry['bib']['start'] += (ends[idx] - cur_start + 100)
                            cur_start = cur_start + (ends[idx] - cur_start + 100)
                        to_del.append(starts[idx])
                # some cases at chromosome starts/ends
                if cur_end - cur_start < 300 or cur_end < cur_start or cur_end < 0 or new_entry['bib']['end'] < 0 or new_entry['bib']['start'] < 0 or new_entry['bib']['end'] - new_entry['bib']['start'] < 0:
                    continue
                bib[new_entry['i']] = new_entry['bib']
            for x in set(to_del):
                del bib[x]

            # check again not overlapping dups - so if everything has worked
            ses = sorted([(i,i+x['end']-x['start']) for i,x in bib.items() if 'dup' not in x])
            starts = [x[0] for x in ses]
            ends = [x[1] for x in ses]

            for i,d in dups.items():
                end = i + d['end']-d['start']
                end_bigger_start_of_dup = bisect_left(ends,i)
                start_smaller_end_of_dup = bisect_left(starts,end)
                for idx in range(end_bigger_start_of_dup-10,min(len(ends),start_smaller_end_of_dup+10)):
                    if starts[idx] == i:
                        continue
                    if overlaps((starts[idx],ends[idx]),(i,end)):
                        print(starts[idx],ends[idx],i,end)
                        pprint(bib[starts[idx]])
                        pprint(d)
                        print('error check non overlapping dups - candidates')
                        exit(1)
        else:
            print(f'No dups file found for {org}, skipping dups processing')

        # check again general overlap...
        ses = sorted([(i,i+x['end']-x['start']) for i,x in bib.items()])
        for enu,se in enumerate(ses[:-1]):
            if se[1] >= ses[enu+1][0]:
                print('error overlapping candidates after adding dups',se,ses[enu+1])
                exit(1)

        print(len(bib))
        with open('../anchors/candidates' + f'/{org}','wb') as f:
            pickle.dump(bib,f)
