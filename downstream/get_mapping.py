#! /usr/bin/env python3

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

def get_mapping():
    org_mapping = {}
    chr_mapping = {}
    if os.path.exists(dir_path+'/../utils/mapping'):
        ff = dir_path+'/../utils/mapping'
    elif os.path.exists('mapping'):
        ff = dir_path+'/mapping'
    with open(ff) as f:
        org = False
        for line in f:
            if line[0] == '#':
                if 'species' in line:
                    org = True
                continue
            if org:
                org_mapping[line.strip().split()[0]] = line.strip().split()[1]
                org_mapping[line.strip().split()[1]] = line.strip().split()[0]
                org = False
                which_org = line.strip().split()[0]
                chr_mapping[which_org] = {}
                chr_mapping[org_mapping[which_org]] = {}
                continue
            else:
                chr_mapping[which_org][line.strip().split()[0]] = line.strip().split()[1]
                chr_mapping[which_org][line.strip().split()[1]] = line.strip().split()[0]
                chr_mapping[org_mapping[which_org]][line.strip().split()[0]] = line.strip().split()[1]
                chr_mapping[org_mapping[which_org]][line.strip().split()[1]] = line.strip().split()[0]
    return(org_mapping,chr_mapping)
