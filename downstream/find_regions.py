#! /usr/bin/env python3

import pickle,numpy
from subprocess import run
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pprint import pprint
from bisect import bisect_left,bisect_right
from pprint import pprint

def check_forward_ambiguity(match_bib):
    if match_bib['meta']['multiple matches out of tolerance range'] == 1 or match_bib['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
        return(True)
    return(False)
def check_reciprocal_ambiguity(org,org2,j):
    if am[org2][j]['matches'][org]['meta']['multiple matches out of tolerance range'] == 1 or am[org2][j]['matches'][org]['meta']['matches have ambiguous matches (tolerance/chromosome out of range)'] == 1:
        return(True)
    return(False)

def get_anchors(org,c,s,e,margin=50000):
    anchors = []
    for i,b in aligned.items():
        if b['chromosome'] == c and b['start'] >= s - margin and b['end'] <= e + margin:
            anchors.append(i)
    return anchors

def find_clusters(coords,pivot):
    result = {}
    oris = {}
    for org in orgs:
        oris[org] = {}
        result[org] = {}
    for c,s,e in coords:
        anchors = get_anchors(pivot,c,s,e,margin=50000)
        if len(anchors) == 0:
            print(f'no anchors for {pivot} and coords {c} {s} {e}')
        for a in anchors:
            for org2,bib2 in aligned[a]['matches'].items():
                if org2 not in orgs:
                    continue
                if check_forward_ambiguity(bib2):
                    continue
                for j,bib3 in bib2['matches'].items():
                    #if check_reciprocal_ambiguity(pivot,org2,j):
                    #    continue
                    if am[org2][j]['chromosome'] not in oris[org2]:
                        oris[org2][am[org2][j]['chromosome']] = []
                    if bib3['match is on other strand in other genome']:
                        oris[org2][am[org2][j]['chromosome']].append('reverse')
                    else:
                        oris[org2][am[org2][j]['chromosome']].append('forward')
                    if am[org2][j]['chromosome'] in result[org2]:
                        result[org2][am[org2][j]['chromosome']].add(int(am[org2][j]['start']+(am[org2][j]['end']-am[org2][j]['start'])/2))
                    else:
                        result[org2][am[org2][j]['chromosome']] = set([int(am[org2][j]['start']+(am[org2][j]['end']-am[org2][j]['start'])/2)])
    return result,oris


orgs = []
with open('orgs') as f:
    for line in f:
        orgs.append(line.strip())
name_mapping = {}
for org in orgs:
    name_mapping[org] = org

tables = {}
for org in orgs:
    tables[org] = open(f'tables_regions/{name_mapping[org]}.tsv','w')
    tables[org].write(f'species\tchromosome\tstart syntenic region\tend syntenic region\tnumber of anchor alignments to the reference species\n')

margin = 50000 #anchors
margin_seq_write = 100000
margin_regions_making = 500000
margin_regions = 50000

am = {}
for org in orgs:
    with open('HOMEDIR/anchors_to_save/462/candidates/' + org,'rb') as f:
        am[org] = pickle.load(f)

coords = {}
with open('../coords') as f:
    for line in f:
        if line[0] == '#':
            continue
        org,chromo,start,end,ori,hit_no = line.split()
        if org not in coords:
            coords[org] = []
        start = int(start)
        end = int(end)
        coords[org].append((chromo,start,end))

bib = {}
oris_g = {}
for org in orgs:
    bib[org] = {}
    oris_g[org] = {}

for pivot in coords:
    with open('HOMEDIR/anchors_to_save/462/aligned/' + pivot,'rb') as f:
        aligned = pickle.load(f)
    result,oris = find_clusters(coords[pivot],pivot)
    for org,r in result.items():
        for chromo,r2 in r.items():
            if chromo not in bib[org]:
                bib[org][chromo] = set()
                oris_g[org][chromo] = []
            bib[org][chromo].update(r2)
            oris_g[org][chromo] += oris[org][chromo]

oris = oris_g

regions = {}
for org in orgs:
    regions[org] = {}

for org,r in bib.items():
    for c,hits in r.items():
        if c not in regions[org]:
            regions[org][c] = [[]]
        l = sorted(hits)
        for i,x in enumerate(l[:-1]):
            regions[org][c][-1].append(x)
            if l[i+1] - x > margin_regions_making:
                regions[org][c].append([])
        if len(regions[org][c][-1]) > 0:
            if l[-1] - regions[org][c][-1][-1] > margin_regions_making:
                regions[org][c].append([l[-1]])
        else:
            regions[org][c][-1].append(l[-1])


never_true = 0
for org,bib2 in regions.items():
    print(org)
    seqs = []
    for c,l in bib2.items():
        print(c)
        region_no = 1
        for hits in l:
            if len(hits) < 5:
                continue
            print(hits)
            for quants in [(0.05,0.95)]:
                genome = SeqIO.parse('HOMEDIR/insects/stable_synteny/utils/genomes/' + org + '.fasta','fasta')
                for seq in genome:
                    if seq.id == c:
                        write_from = max(int(numpy.quantile(list(hits),quants[0]))-margin_seq_write,0)
                        write_to = min(int(numpy.quantile(list(hits),quants[1]))+margin_seq_write,len(seq))
                        print(write_from,write_to)
                        tables[org].write(f'{name_mapping[org]}\t{c}\t{write_from}\t{write_to}\t{len(hits)}\n')
                        run(f'mkdir -p fastas',shell=True)
                        run(f'mkdir -p fastas/{org}',shell=True)
                        run(f'mkdir -p fastas/{org}/{c}',shell=True)
                        seqs.append(SeqRecord(Seq(seq.seq[write_from:write_to]),id=f'{org}_{c}_{write_from}_{write_to}_forward',name=f"region_no_{region_no}_{org}",description=f"region_no_{region_no}_{org}"))
                        SeqIO.write(SeqRecord(Seq(seq.seq[write_from:write_to]),id=f'{org}_{c}_{write_from}_{write_to}_forward',name=f"region_no_{region_no}_{org}",description=f"region_no_{region_no}_{org}"),f'fastas/{org}/{c}/{region_no}.fasta','fasta')
                        region_no += 1
                        break
    if len(seqs) > 0:
        SeqIO.write(seqs,f'fastas/{org}/all.fasta','fasta')
for org in orgs:
    tables[org].close()
pprint(regions)
