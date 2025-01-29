#! /usr/bin/python3

from os import listdir
import sys
sys.path.append('./utils/')
from get_mapping import get_mapping
from pprint import pprint
import pickle

org_mapping,chr_mapping = get_mapping()

with open('MCScanX/alloc_pairwise','rb') as f:
    alloc_m = pickle.load(f)

with open('i-ADHoRe/alloc_pairwise','rb') as f:
    alloc_a = pickle.load(f)

m_check = {}
for org,bib in alloc_m.items():
    m_check[org] = {}
    for element in bib:
        m_check[org][element] = 1
a_check = {}
for org,bib in alloc_a.items():
    a_check[org] = {}
    for element in bib:
        a_check[org][element] = 1
c_check = {}
custom_clusters = listdir('custom/clusters/50000')
for c in custom_clusters:
    with open(f'custom/clusters/50000/{c}') as f:
        for line in f:
            if 'species' in line:
                org = line.strip().split()[1]
                org_orig = org
                org = org_mapping[org]
                if org not in c_check:
                    c_check[org] = {}
            # element on chromosome seq453 start 2740 end 2792 strand 1
            if 'element on chromosome' in line:
                chromo = line.strip().split('element on chromosome')[1].split()[0]
                chromo = chr_mapping[org_orig][chromo]
                strand = line.strip().split('element on chromosome')[1].split()[6]
                start = line.strip().split('element on chromosome')[1].split()[2]
                end = line.strip().split('element on chromosome')[1].split()[4]
                if strand == '1':
                    strand = '+'
                else:
                    strand = '-'
                element = f'{chromo}CUSTOMele{start}to{end}{strand}'
                c_check[org][element] = 1
coords_check = {}
with open('coords') as f:
    for line in f:
        line = line.strip().split()
        if len(line) < 5:
            continue
        org,chromo,start,end,ori = line
        chromo = chr_mapping[org][chromo]
        org = org_mapping[org]
        if org not in coords_check:
            coords_check[org] = {}
        if ori == 'forward':
            strand = '+'
        elif ori == 'reverse':
            strand = '-'
        element = f'{chromo}CUSTOMele{start}to{end}{strand}'
        coords_check[org][element] = 1

sum_tot,sum_m,sum_a,sum_c = 0,0,0,0
out = open('latex_table','w')
for org,bib in coords_check.items():
    print(org,len(bib),len(m_check[org]),len(a_check[org]),len(c_check[org]))
    out.write(f'{org} & {len(bib)} & \hspace*{{0.42cm}} {len(m_check[org])} & {len(a_check[org])} & {len(c_check[org])}\\\\\n')
    #if org == '48org':
    #    for element,bib2 in bib.items():
    #        if element not in m_check[org]:
    #            print(f'element {element} not in m_check')
    #        else:
    #            pass
    #        if element not in a_check[org]:
    #            print(f'element {element} not in a_check')
    #        if element not in c_check[org]:
    #            print(f'element {element} not in c_check')
    #            input()
    #        else:
    #            pass
    sum_tot += len(bib)
    sum_m += len(m_check[org])
    sum_a += len(a_check[org])
    sum_c += len(c_check[org])
out.write('\hhline{|=|=|=|}')
out.write(f'sum & {sum_tot} & \hspace*{{0.42cm}} {sum_m} & {sum_a} & {sum_c}\\\\')
out.close()
