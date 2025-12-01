#! /usr/bin/env python3


import sys
sys.path.append('../utils/')
from get_mapping import get_mapping
import pickle
from pygenomeviz import GenomeViz
from subprocess import run
from pprint import pprint
import argparse
from os.path import isfile

def turn(start,end,segments):
    max_e = 0
    min_s = 1e15
    for segment in segments:
        max_e = max(segment[1],max_e)
        min_s = min(segment[0],min_s)
    len_line = max_e - min_s
    new_end = len_line - start
    new_start = new_end - (end - start)
    return new_start, new_end

def get_gff():

    gff = {}
    with open('MCScanX.gff') as f:
        for line in f:
            chromo,name,start,end = line.strip().split()
            gff[name] = (chromo,start,end)

    return gff

if __name__ == "__main__":

    parser= argparse.ArgumentParser(description='Draw MCScanX chains. The order maximise pariwise aligned length greedily.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--mcscanx_file', required=True, type=str,nargs=1, default='MCScanX.collinearity',help="MCScanX collinearity file")
    parser.add_argument('--threshold_chr', required=True, type=int,nargs=1, default=1000000,help="Threshold size for Chromosome/Scaffold/Contig to include in drawing")
    args = parser.parse_args()
    mcscanx_file = args.mcscanx_file[0]
    threshold_chr = args.threshold_chr[0]
    
    run(f'mkdir -p singles_out',shell=True)
    run(f'mkdir -p simple_maps',shell=True)
    run(f'mkdir -p orders',shell=True)

    run(f'python3 get_mc_blocks.py {mcscanx_file}',shell=True)
    run(f'python3 get_orders.py {threshold_chr}',shell=True)
   
    tr = []
    with open('orders_orgs') as f:
        for line in f:
            target,ref = line.strip().split('\t')
            tr.append((target,ref))


    org_mapping,chr_mapping = get_mapping()

    gff = get_gff()
    # save if using multiple times
    #with open('parsed_gff','wb') as f:
    #    pickle.dump(gff,f)
    #with open('parsed_gff','rb') as f:
    #    gff = pickle.load(f)

    drawing_order = {}
    for target,ref in tr:
        drawing_order[target] = []
        if not isfile(f'orders/{target}'):
            continue
        with open(f'orders/{target}') as f:
            for line in f:
                seq, ori = line.strip().split('\t')
                drawing_order[target].append((seq,ori))

    small_meta = {}
    for org in [x[0] for x in tr]:
        with open(f'../utils/small_meta/{org}','rb') as f:
            small_meta[org] = pickle.load(f)
    with open(f'../utils/small_meta/{tr[0][1]}','rb') as f:
        small_meta[tr[0][1]] = pickle.load(f)

    gv = GenomeViz()
    new_lengths = {}
    cur = 0
    target_ranges = []
    draw_ch = {}
    draw_ch[tr[0][1]] = {}
    new_lengths[tr[0][1]] = {}
    orig_lengths = []
    tmporder = []
    for enu,ch in enumerate(small_meta[tr[0][1]][0]):
        if enu == 0:
            if small_meta[tr[0][1]][1][enu] < int(threshold_chr):
                continue
            else:
                draw_ch[tr[0][1]][ch] = 1
                target_ranges.append((cur,cur+small_meta[tr[0][1]][1][enu]))
                new_lengths[tr[0][1]][ch] = (cur,cur+small_meta[tr[0][1]][1][enu])
                orig_lengths.append((0,small_meta[tr[0][1]][1][enu]))
                cur = cur+small_meta[tr[0][1]][1][enu]
                tmporder.append(ch)
        elif small_meta[tr[0][1]][1][enu] - small_meta[tr[0][1]][1][enu-1] < int(threshold_chr):
            continue
        else:
            draw_ch[tr[0][1]][ch] = 1
            target_ranges.append((cur,cur+small_meta[tr[0][1]][1][enu] - small_meta[tr[0][1]][1][enu-1]))
            new_lengths[tr[0][1]][ch] = (cur,cur+small_meta[tr[0][1]][1][enu] - small_meta[tr[0][1]][1][enu-1])
            orig_lengths.append((small_meta[tr[0][1]][1][enu-1] + 1,small_meta[tr[0][1]][1][enu]))
            cur = cur+small_meta[tr[0][1]][1][enu] - small_meta[tr[0][1]][1][enu-1] + 1
            tmporder.append(ch)
    track = gv.add_feature_track(tr[0][1],segments=target_ranges)

    for enu,segment in enumerate(track.segments):
            segment.add_sublabel('')#f'{tmporder[enu]}')
    track.set_segment_sep(symbol="//")

    for org in [x[0] for x in tr]:
        cur = 0
        target_ranges = []
        draw_ch[org] = {}
        new_lengths[org] = {}
        orig_lengths = []
        tmporder = []
        for seq,ori in drawing_order[org]:
            enu = small_meta[org][0].index(seq)
            if enu == 0:
                if small_meta[org][1][enu] < int(threshold_chr):
                    continue
                else:
                    draw_ch[org][seq] = 1
                    target_ranges.append((cur,cur+small_meta[org][1][enu]))
                    new_lengths[org][seq] = (cur,cur+small_meta[org][1][enu])
                    orig_lengths.append((0,small_meta[org][1][enu]))
                    cur = cur+small_meta[org][1][enu]
                    tmporder.append(seq)
            elif small_meta[org][1][enu] - small_meta[org][1][enu-1] < int(threshold_chr):
                continue
            else:
                draw_ch[org][seq] = 1
                target_ranges.append((cur,cur+small_meta[org][1][enu] - small_meta[org][1][enu-1]))
                new_lengths[org][seq] = (cur,cur+small_meta[org][1][enu] - small_meta[org][1][enu-1])
                orig_lengths.append((max(0,small_meta[org][1][enu-1] + 1),small_meta[org][1][enu]))
                cur = cur+small_meta[org][1][enu] - small_meta[org][1][enu-1] + 1
                tmporder.append(seq)
        track = gv.add_feature_track(org,segments=target_ranges)

        for enu,segment in enumerate(track.segments):
            segment.add_sublabel('')#f'{tmporder[enu]}')
        track.set_segment_sep(symbol="//")
   

    draw = 0


    with open('naive_blocks.pickle','rb') as f:
        blocks = pickle.load(f)

    homology = {}
    for b in blocks:
        t1,t2,rev = b
        first1,last1 = t1
        first2,last2 = t2
        c11,s11,e11 = gff[first1]
        c12,s12,e12 = gff[last1]
        c21,s21,e21 = gff[first2]
        c22,s22,e22 = gff[last2]
        s11 = int(s11)
        e11 = int(e11)
        s12 = int(s12)
        e12 = int(e12)
        s21 = int(s21)
        e21 = int(e21)
        s22 = int(s22)
        e22 = int(e22)
        org1 = org_mapping[first1.split('ele')[0].split('chr')[0]]
        org2 = org_mapping[first2.split('ele')[0].split('chr')[0]]
        chromo1 = chr_mapping[org1][c11]
        chromo2 = chr_mapping[org2][c21]
        start1 = s11
        end1 = e12
        start2 = s21
        end2 = e22
        if org1 not in homology:
            homology[org1] = {}
        if chromo1 not in homology[org1]:
            homology[org1][chromo1] = {}
        if start1 not in homology[org1][chromo1]:
            homology[org1][chromo1][start1] = {'matches': {}}
        homology[org1][chromo1][start1]['matches'][org2] = (chromo2,start2,rev,end1)

        if org2 not in homology:
            homology[org2] = {}
        if chromo2 not in homology[org2]:
            homology[org2][chromo2] = {}
        if start2 not in homology[org2][chromo2]:
            homology[org2][chromo2][start2] = {'matches': {}}
        homology[org2][chromo2][start2]['matches'][org1] = (chromo1,start1,rev,end2)
    
    for org1,org2 in tr:
        with open(f'simple_maps/{org1}','rb') as f:
            simple_map = pickle.load(f)
        for seq,ori in drawing_order[org1]:
            if seq not in homology[org1]:
                continue
            enu = small_meta[org1][0].index(seq)
            for start1,bib in homology[org1][seq].items():
                if org2 in bib['matches']:
                    chromo2,start2,rev,end1 = bib['matches'][org2]
                    if chromo2 not in draw_ch[org2]:
                        continue
                    end2 = homology[org2][chromo2][start2]['matches'][org1][3]
                    nstart1 = new_lengths[org1][seq][0] + start1
                    nend1 = new_lengths[org1][seq][0] + end1
                    nstart2 = new_lengths[org2][chromo2][0] + start2
                    nend2 = new_lengths[org2][chromo2][0] + end2
                    if rev == 1:
                        nstart2,nend2 = nend2,nstart2
                    if chromo2 == simple_map[seq][0]:
                        col = "blue"
                        inv_col = "lime"
                        al = 0.1
                    else:
                        col = "magenta"
                        inv_col = "orange"
                        al = 0.9
                    for track in gv.feature_tracks:
                        if track.name != org1:
                            continue
                        for segment in track.segments:
                            segstart = int(segment.start)
                            segend = int(segment.end)
                            if nstart1 < segstart or nend1 > segend:
                                continue
                            else:
                                seg1 = segment.name
                    for track in gv.feature_tracks:
                        if track.name != org2:
                            continue
                        for segment in track.segments:
                            segstart = int(segment.start)
                            segend = int(segment.end)
                            if nstart2 < segstart or nend2 > segend:
                                continue
                            else:
                                seg2 = segment.name

                    gv.add_link((org1,seg1,nstart1,nend1),(org2,seg2,nstart2,nend2), color=col, inverted_color=inv_col, curve=True, alpha = al)
                    draw += 1

fig = gv.plotfig()
fig.savefig(f'mcscx_blocks.svg', dpi=300)
