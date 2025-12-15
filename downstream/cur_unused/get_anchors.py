#! /usr/bin/env python3

from sys import argv
from pprint import pprint
import pickle

if __name__ == "__main__":
    org,chromo,start,end,margin = argv[1:6]
    margin = int(margin)
    start = int(start)
    end = int(end)
    with open('/scr/k80san/karl//anchors_to_save/462/aligned_corrected/'+org,'rb') as f:
        anchors = pickle.load(f)

    for i,b in anchors.items():
        if b['chromosome'] == chromo:# and b['start'] >= start-margin and b['end'] <= end+margin:
            print(i)
            pprint(b)
            input()
