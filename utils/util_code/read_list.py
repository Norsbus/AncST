#! /usr/bin/env python3

import pickle
from sys import argv

def read_pickle(path):
    with open(path,'rb') as f:
        x = pickle.load(f)

    if int(argv[2]) in x:
        print('in there')
    else:
        print('not in there')
    return(0)

if __name__ == "__main__":
    read_pickle(argv[1])
