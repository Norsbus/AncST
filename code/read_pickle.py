#! /usr/bin/env python3

import pickle
from sys import argv

def read_pickle(path):
    with open(path,'rb') as f:
        x = pickle.load(f)
    if len(argv) > 2:
        if argv[2] == 'multies':
            for i,d in x.items():
                if 1 in [x[0] for x in d['matches'].values()]:
                    print(i,d)
        elif argv[2] == 'chromos':
            no = 0
            mini = 0
            maxi = 1e15
            for i,d in x.items():
                if d['chromosome'] == argv[3]:
                    no += 1
                    if d['start'] > mini:
                        mini = d['start']
                    if d['end'] > maxi:
                        maxi = d['end']
            print(mini,maxi,no)

        else:
            try:
                print(x[int(argv[2])])
            except:
                #print('not there')
                #x = [x[0] for x in x]
                if int(argv[2]) in x:
                    print('in there')
                else:
                    print('not in there')
    else:
        print(len(list(x.keys())))
    return(0)

if __name__ == "__main__":
    read_pickle(argv[1])
