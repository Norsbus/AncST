#! /usr/bin/env python3

import pickle
from sys import argv
from pprint import pprint
from statistics import mean,median,stdev

def read_pickle(path):
    with open(path,'rb') as f:
        x = pickle.load(f)
    pprint(x[int(argv[2])])
    return(0)

if __name__ == "__main__":
    read_pickle(argv[1])
