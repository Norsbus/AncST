#! /usr/bin/env python3

import pickle
from sys import argv
from pprint import pprint

def read_pickle(path):
    with open(path,'rb') as f:
        x = pickle.load(f)
    if len(argv) == 3:
        try:
            pprint(x[int(argv[2])])
        except:
            input('why not')
            print(x)
    else:
        #for k,v in x.items():
        #    print(f'{k} : {v}')
        print(x)
    return(0)

if __name__ == "__main__":
    read_pickle(argv[1])
