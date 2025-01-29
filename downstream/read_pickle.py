#! /usr/bin/env python3

import pickle
from sys import argv
from pprint import pprint
from statistics import mean,median,stdev

def read_pickle(path):
    with open(path,'rb') as f:
        x = pickle.load(f)
    for k,v in x.items():
        if v['chromosome'] == 'NT_033779.5' and v['start'] >= 5520000 and v['end'] <= 5524000 and 'GCF_025231255.1' in v['matches']:
            print('---')
            print(k)
            pprint(v['matches']['GCF_025231255.1'])
    return(0)

if __name__ == "__main__":
    read_pickle(argv[1])
