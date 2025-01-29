#! /usr/bin/env python3

from sys import argv
from matplotlib.pyplot import cm
from pygenomeviz import GenomeViz
from random import randint
import numpy as np
from subprocess import run
from pprint import pprint
import pickle
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from os.path import isfile

### IMPORTANT NOTES $$$

# 1. since the chromsome/contig lengths is not explicitely known here, the drawing can be longer than a chromosome/contig

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

def overlaps(ranges):
    ranges = sorted(ranges)
    it = iter(ranges)
    try:
        curr_start, curr_stop, curr_color, curr_strands = next(it)
    except StopIteration:
        return
    for start, stop, color, strands in it:
        if curr_start <= start <= curr_stop:
            curr_stop = max(curr_stop, stop)
            curr_color.update(color)
            curr_strands.update(strands)
        else:
            yield curr_start, curr_stop, curr_color, curr_strands
            curr_start, curr_stop, curr_color, curr_strands = start, stop, color, strands
    yield curr_start, curr_stop, curr_color, curr_strands

if len(argv) > 2:
    ref = (argv[2],argv[3],argv[4])
else:
    ref = 0

margin = int(argv[1])

orgs = []
coords = {}
to_ori = set()
with open('../utils/coords_both_genomic_and_regional_nono') as f:
    for line in f:
        line = line.strip().split()
        if len(line) == 5:
            org,chromo,start,end,ori_gene = line
        elif len(line) == 6:
            org,chromo,start,end,ori_gene,hit_name = line
        if (chromo,org,start) == ref:
            ori[ref] = ori_gene
        if org not in coords:
            coords[org] = {}
            orgs.append(org)
        if chromo not in coords[org]:
            coords[org][chromo] = []
            to_ori.add((org,chromo))
        if ori_gene == 'forward':
            strand = 1
        else:
            strand = -1
        if ref == 0:
            ref = (org,chromo,start)
        coords[org][chromo].append((int(start),int(end),strand,hit_name))

run(f'mkdir -p tables_drawing',shell=True)
tables = {}
for org in orgs:
    tables[org] = open(f'tables_drawing/{org}.tsv','w')
    tables[org].write(f'species\tchromosome\tstart\tend\tstrand\tlabel\toriginal protein name (in ncbi annotation/proteome)\n')

tot_lens = {}
for org in orgs:
    tot_lens[org] = 0

ori = {}
bib = {}

no_colors = 0
with open('pairwise_alignments_table_only_relevant') as f:
    for line in f:
        if line[0] == '#':
            continue
        no_colors += 1

colors = cm.rainbow(np.linspace(0, 1, no_colors))
used = {}
orientation = {}
alignments = {}

pairs = {}

with open('pairwise_alignments_table_only_relevant') as f:
    for line in f:
        if line[0] == '#':
            continue
        no,org1,chromo1,start1,end1,org2,chromo2,start2,end2,ori_alignment,score,snd1,snd2 = line.split()
        if org1 not in coords or org2 not in coords or chromo1 not in coords[org1] or chromo2 not in coords[org2]:
            continue
        start1 = int(start1)
        end1 = int(end1)
        start2 = int(start2)
        end2 = int(end2)
        if start1 > end1:
            start1,end1 = end1,start1
        if start2 > end2:
            start2,end2 = end2,start2
        relevant1 = 0
        rel1 = []
        relevant2 = 0
        rel2 = []
        for t in coords[org1][chromo1]:
            start,end,strand,hit_name = t
            if start1 > start - margin and end1 < end + margin:
                relevant1 = 1
                rel1.append((org1,chromo1,t[0]))
        for t in coords[org2][chromo2]:
            start,end,strand,hit_name = t
            if start2 > start - margin and end2 < end + margin:
                relevant2 = 1
                rel2.append((org2,chromo2,t[0]))
        if relevant1 == 0 or relevant2 == 0:
            continue
        for ele1 in rel1:
            if ele1 not in pairs:
                pairs[ele1] = set()
            for ele2 in rel2:
                if ele2 not in pairs:
                    pairs[ele2] = set()
                pairs[ele1].add(ele2)
                pairs[ele2].add(ele1)

        if (org1,chromo1,start1,end1) not in alignments:
            alignments[(org1,chromo1,start1,end1)] = []
        alignments[(org1,chromo1,start1,end1)].append([org2,chromo2,start2,end2,ori_alignment])
        if (org2,chromo2,start2,end2) not in alignments:
            alignments[(org2,chromo2,start2,end2)] = []
        alignments[(org2,chromo2,start2,end2)].append([org1,chromo1,start1,end1,ori_alignment])
        if org1 not in bib:
            bib[org1] = {}
            orientation[org1] = {}
        if chromo1 not in bib[org1]:
            bib[org1][chromo1] = []
            orientation[org1][chromo1] = {}
        if org2 not in bib:
            bib[org2] = {}
            orientation[org2] = {}
        if chromo2 not in bib[org2]:
            bib[org2][chromo2] = []
            orientation[org2][chromo2] = {}
        if org1 not in orientation[org2][chromo2]:
            orientation[org2][chromo2][org1] = {}
        if chromo1 not in orientation[org2][chromo2][org1]:
            orientation[org2][chromo2][org1][chromo1] = []
        if org2 not in orientation[org1][chromo1]:
            orientation[org1][chromo1][org2] = {}
        if chromo2 not in orientation[org1][chromo1][org2]:
            orientation[org1][chromo1][org2][chromo2] = []
        ri = randint(0,no_colors-1)
        while ri in used:
            ri = randint(0,no_colors-1)
        used[ri] = 1
        color = ri
        bib[org1][chromo1].append((start1,end1,set([color])))
        bib[org2][chromo2].append((start2,end2,set([color])))
        orientation[org1][chromo1][org2][chromo2].append(ori_alignment)
        orientation[org2][chromo2][org1][chromo1].append(ori_alignment)

if (ref[0],ref[1],int(ref[2])) not in pairs or len(pairs[(ref[0],ref[1],int(ref[2]))]) == 0:
    print('exiting because focal element does not have alignments')
    exit(0)
else:
    to_draw = set()
    to_draw.add((ref[0],ref[1],int(ref[2])))
    new = 1
    cur_ele = (ref[0],ref[1],int(ref[2]))
    new_eles = pairs[cur_ele]
    while new > 0:
        new = 0
        ne = set()
        for cur_ele in new_eles:
            for ele2 in pairs[cur_ele]:
                if ele2 in to_draw:
                    continue
                else:
                    to_draw.add(ele2)
                    ne.add(ele2)
        new = len(new_eles)
        new_eles = ne

new_bib = {}
tdbib = {}
for org,chromo,start in to_draw:
    if org not in tdbib:
        tdbib[org] = {}
    if chromo not in tdbib[org]:
        tdbib[org][chromo] = []
    for startele,end,strand,hit_name in coords[org][chromo]:
        if start == startele:
            tdbib[org][chromo].append((start,end,strand))

for org,x in bib.items():
    if org not in tdbib:
        continue
    for chromo,l in x.items():
        if chromo not in tdbib[org]:
            continue
        for ele in l:
            for gene in tdbib[org][chromo]:
                if ele[0] >= gene[0]-margin and ele[1] <= gene[1]+margin:
                    if org not in new_bib:
                        new_bib[org] = {}
                    if chromo not in new_bib[org]:
                        new_bib[org][chromo] = []
                    new_bib[org][chromo].append(ele)

bib = new_bib

for org,x in orientation.items():
    for chromo,y in x.items():
        if (org,chromo) in ori:
            continue
        if org not in orientation[ref[0]][ref[1]]:
            continue
        if chromo not in orientation[ref[0]][ref[1]][org]:
            continue
        oris = orientation[ref[0]][ref[1]][org][chromo]
        if oris.count('forward') > oris.count('reverse'):
            ori[(org,chromo)] = 'forward'
        else:
            ori[(org,chromo)] = 'reverse'


for org,x in orientation.items():
    for chromo,y in x.items():
        if (org,chromo) in ori:
            continue
        for org2,z in y.items():
            for chromo2,oris in z.items():
                if (org2,chromo2) not in ori:
                    continue
                if org2 not in orientation[org][chromo]:
                    continue
                if chromo2 not in orientation[org][chromo][org2]:
                    continue
                oris = orientation[org][chromo][org2][chromo2]
                if oris.count('forward') > oris.count('reverse'):
                    if ori[(org2,chromo2)] == 'forward':
                        ori[(org,chromo)] = 'forward'
                    else:
                        ori[(org,chromo)] = 'reverse'
                else:
                    if ori[(org2,chromo2)] == 'forward':
                        ori[(org,chromo)] = 'reverse'
                    else:
                        ori[(org,chromo)] = 'forward'

first = 1
if isfile('try_to_draw_all'):
    while len(ori) != len(to_ori):
        if first == 1:
            for org,chromo in to_ori:
                if (org,chromo) not in ori:
                    ori[(org,chromo)] = 'forward'
                    first = 0
                    break
        else:
            first = 1
            for org,x in orientation.items():
                for chromo,y in x.items():
                    if (org,chromo) in ori:
                        continue
                    for org2,z in y.items():
                        for chromo2,oris in z.items():
                            if (org2,chromo2) not in ori:
                                continue
                            if org2 not in orientation[org][chromo]:
                                continue
                            if chromo2 not in orientation[org][chromo][org2]:
                                continue
                            oris = orientation[org][chromo][org2][chromo2]
                            if oris.count('forward') > oris.count('reverse'):
                                if ori[(org2,chromo2)] == 'forward':
                                    ori[(org,chromo)] = 'forward'
                                else:
                                    ori[(org,chromo)] = 'reverse'
                            else:
                                if ori[(org2,chromo2)] == 'forward':
                                    ori[(org,chromo)] = 'reverse'
                                else:
                                    ori[(org,chromo)] = 'forward'

for org,x in bib.items():
    for chromo,l in x.items():
        new_l = []
        for start,end,color_set in l:
            strands = set()
            alis = alignments[(org,chromo,start,end)]
            for ali in alis:
                org2,chromo2,start2,end2,orient = ali
                if orient == '*':
                    strands.add('reverse')
                    continue
                if (org,chromo) not in ori:
                    continue
                if (org2,chromo2) not in ori:
                    continue
                ori1 = ori[(org,chromo)]
                ori2 = ori[(org2,chromo2)]
                if (ori1 == ori2 and ali[-1] == 'forward') or (ori1 != ori2 and ali[-1] == 'reverse'):
                    strands.add('forward')
                else:
                    if ali[-1] == 'forward':
                        if ori1 == 'forward':
                            strands.add('forward')
                        else:
                            strands.add('reverse')
                    else:
                        strands.add('forward')
                        for enu,t in enumerate(alignments[(org2,chromo2,start2,end2)]):
                            a,b,c,d,e = t
                            if a == org and b == chromo and c == start and d == end:
                                alignments[(org2,chromo2,start2,end2)][enu][-1] = '*'
                        
            new_l.append((start,end,color_set,strands))
        bib[org][chromo] = new_l

color_sets = {}
members = {}
counter = 1

for org,x in bib.items():
    for chromo,l in x.items():
        l = list(overlaps(l))
        for start,end,color_set,strands in l:
            colors_taken = set()
            for color in color_set:
                if color in members:
                    colors_taken.add(members[color])
            if len(colors_taken) == 0:
                color_sets[counter] = color_set
                for color in color_set:
                    members[color] = counter
                counter += 1
            elif len(colors_taken) == 1:
                taken = list(colors_taken)[0]
                color_sets[taken].update(color_set)
                for color in color_set:
                    members[color] = taken
            else:
                color_sets[counter] = set()
                for color in colors_taken:
                    color_sets[counter].update(color_sets[color])
                    for color2 in color_sets[color]:
                        members[color2] = counter
                    del color_sets[color]
                color_sets[counter].update(color_set)
                for color in color_set:
                    members[color] = counter
                counter += 1
        bib[org][chromo] = l

new_colors = iter(cm.rainbow(np.linspace(0, 1, len(color_sets))))
colors = {}

for org,x in bib.items():
    for chromo,l in x.items():
        new_l = []
        for start,end,color_set,strands in l:
            cont = 0
            for coord in coords[org][chromo]:
                if (start >= coord[0] and start <= coord[1]) or (end >= coord[0] and end <= coord[1]) or (start <= coord[0] and end >= coord[1]) or (start >= coord[0] and end <= coord[1]):
                    cont = 1
                    break
            if cont == 1:
                continue
            if len(set([members[color] for color in color_set])) != 1:
                input('ERROR collecting all overlapping fragments and their colors')
            else:
                color_no = list(set([members[color] for color in color_set]))[0]
            if color_no not in colors:
                colors[color_no] = next(new_colors)
            color = colors[color_no]
            new_l.append((start,end,color,strands))
        bib[org][chromo] = new_l

color_counter = 1
used = {}
gv = GenomeViz()
shifts = {}
genome_sizes = {}

org_chromo_order = []
if isfile('plotting_order'):
    with open('plotting_order') as f:
        for line in f:
            org,chromo = line.strip().split()
            org_chromo_order.append((org,chromo))
else:
    for org,d in bib.items():
        for chromo in d:
            org_chromo_order.append((org,chromo))

new_order = []
for org,chromo in org_chromo_order:
    if org in bib and chromo in bib[org]:
        new_order.append((org,chromo))
org_chromo_order = new_order

for org,chromo in org_chromo_order:
    l = bib[org][chromo]
    if (org,chromo) not in ori or len(l) == 0:
        continue
    # first find elements and their size to make figures and get shift
    start = min(y[0] for y in l)
    end = max(y[1] for y in l)
    # i take anchors to actually check which genes are in the cluster
    if org not in coords or chromo not in coords[org]:
        continue
    start_min_ele = 1e15
    end_max_ele = 0
    for start_l,end_l,strand,hit_name in coords[org][chromo]:
        start_l,end_l = int(start_l),int(end_l)
        if (end_l < start and start-end_l > margin) or (start_l > end and start_l - end > margin):
            continue
        if start_l < start_min_ele:
            start_min_ele = start_l
        if end_l > end_max_ele:
            end_max_ele = end_l
    if start_min_ele == 1e15 and end_max_ele == 0:
        continue
    start = start_min_ele - margin
    local_margin_up = margin
    local_margin_down = margin
    if start < 0:
        local_margin_up = start_min_ele
        start = 0
    end = end_max_ele + margin
    if ori[(org,chromo)] == 'reverse':
        local_margin_up,local_margin_down = local_margin_down,local_margin_up
    shift = start
    genome_size = end - start
    name = org+' '+chromo + ' ' + ori[(org,chromo)]
    name = name.replace('_',' ')
    shifts[name] = shift

    gene_drawings = []
    for start_l,end_l,strand,hit_name in coords[org][chromo]:
        start_l,end_l = int(start_l),int(end_l)
        if end_l < start or start_l > end:
            continue
        start_l -= shift
        end_l -= shift
        if ori[(org,chromo)] == 'reverse':
            start_l,end_l = turn(start_l,end_l,[(0,genome_size)])
            if strand == -1:
                strand = 1
            else:
                strand = -1
        if start_l > genome_size or end_l > genome_size or start_l < 0 or end_l < 0:
            continue
        gene_drawings.append((max(0,start_l-local_margin_up),min(end_l+local_margin_down,genome_size)))

    gene_drawings = sorted(gene_drawings)
    target_ranges = []
    i = 1
    s1 = gene_drawings[0][0]
    e1 = gene_drawings[0][1]
    while i < len(gene_drawings):
        s2 = gene_drawings[i][0]
        e2 = gene_drawings[i][1]
        if s2 - e1 > 500000:
            target_ranges.append((s1,e1))
            s1 = s2
        e1 = e2
        i+=1
    if len(target_ranges) == 0:
        target_ranges = (gene_drawings[0][0],gene_drawings[-1][1])
    else:
        if gene_drawings[-1][1] >  target_ranges[-1][1]:
            if target_ranges[-1][0] != s1:
                target_ranges.append((s1,gene_drawings[-1][1]))
            else:
                target_ranges.append((gene_drawings[-1][0],gene_drawings[-1][1]))

    track = gv.add_feature_track(name,segments=target_ranges)
    seg_labs = []
    genome_sizes[name] = []
    for segment in track.segments:
        if ori[(org,chromo)] == 'reverse':
            seg_labs.append((segment.start,segment.end))
        else:
            segment.add_sublabel(f'{segment.start+shift}-{segment.end+shift}')
        genome_sizes[name].append((int(segment.start),int(segment.end)))
        tot_lens[org] += int(segment.end) - int(segment.start)
    if len(seg_labs) > 0:
        seg_labels = []
        seg_labs = sorted(seg_labs)
        tot_len = max(x[1] for x in seg_labs) - min(x[0] for x in seg_labs)
        cur_size = tot_len
        for i,t in enumerate(seg_labs[:-1]):
            s,e = t
            len_seg = e-s
            seg_labels.append(f'{cur_size} - {cur_size-len_seg}')
            cur_size = cur_size - len_seg - (seg_labs[i+1][0] - e)
        s = seg_labs[-1][0]
        e = seg_labs[-1][1]
        len_seg = e-s
        seg_labels.append(f'{cur_size} - {cur_size-len_seg}')
        for i,s in enumerate(track.segments):
            s.add_sublabel(seg_labels[i])
    track.set_segment_sep(symbol="//")

    for start_l,end_l,strand,hit_name in coords[org][chromo]:
        orig_start = start_l
        orig_end = end_l
        start_l,end_l = int(start_l),int(end_l)
        if end_l < start or start_l > end:
            continue
        start_l -= shift
        end_l -= shift
        if strand == -1:
            vpos = 'bottom'
            ymargin = 2.5
        else:
            vpos = 'top'
            ymargin = 0.8
        if ori[(org,chromo)] == 'reverse':
            start_l,end_l = turn(start_l,end_l,genome_sizes[name])
            if strand == -1:
                strand = 1
                vpos = 'top'
                ymargin = 0.8
            else:
                strand = -1
                vpos = 'bottom'
                ymargin = 2.5
            
        label = hit_name
        facecolor = 'black'
        for segment in track.segments:
            try:
                segment.add_feature(start_l,end_l,facecolor='black',label = label, text_kws=dict(color=facecolor,vpos=vpos,ymargin=ymargin,hpos='center'),strand = strand,plotstyle='arrow')
                tables[org].write(f'{org}\t{chromo}\t{orig_start}\t{orig_end}\t{strand}\t{label}\t{label_for_dict}\n')
            except:
                pass
    
    for start_l,end_l,color,strands in l:
        if start_l > end_l:
            start_l,end_l = end_l,start_l
        cont = 0
        for coord in coords[org][chromo]:
            if (start_l >= coord[0] and start_l <= coord[1]) or (end_l >= coord[0] and end_l <= coord[1]) or (start_l <= coord[0] and end_l >= coord[1]) or (start_l >= coord[0] and end_l <= coord[1]):
                cont = 1
                break
        if cont == 1:
            continue
        start_l = start_l - shift
        end_l = end_l - shift
        if ori[(org,chromo)] == 'reverse':
            start_l,end_l = turn(start_l,end_l,genome_sizes[name])
        if 'reverse' in strands:
            strand_l = -1
        else:
            strand_l = 1
        if tuple(color) not in used:
            used[tuple(color)] = color_counter
            color_counter += 1
        label = ''
        strand_l = 1
        added = 0
        for segment in track.segments:
            try:
                segment.add_feature(start_l,end_l,facecolor=color,strand = strand_l,plotstyle='arrow',label=label)
                added = 1
                no_anchors[org] += 1
            except:
                pass
        if added == 0:
            print(f'CANNOT ADD anchor of {org} {chromo} {start_l} {end_l} {strands}...probably out of range: genome_size = {genome_size}')

done = {}
for t1,l in alignments.items():
    org,chromo,start1,end1 = t1
    if (org,chromo) not in ori:
        continue
    orig_start1 = start1
    orig_end1 = end1
    orig_start1 = int(orig_start1)
    orig_end1 = int(orig_end1)
    cont = 0
    for coord in coords[org][chromo]:
        if (start1 >= coord[0] and start1 <= coord[1]) or (end1 >= coord[0] and end1 <= coord[1]) or (start1 <= coord[0] and end1 >= coord[1]) or (start1 >= coord[0] and end1 <= coord[1]):
            cont = 1
            break
    if cont == 1:
        continue
    name1 = org+' '+chromo+' '+ori[(org,chromo)]
    name1 = name1.replace('_',' ')
    if name1 not in shifts:
        continue
    shift1 = shifts[name1]
    start1 -= shift1
    end1 -= shift1
    if ori[(org,chromo)] == 'reverse':
        start1,end1 = turn(start1,end1,genome_sizes[name1])
    for org2,chromo2,start2,end2,ori_alignment in l:
        if (org2,chromo2) not in ori:
            continue
        orig_start2 = start2
        orig_end2 = end2
        orig_start2 = int(orig_start2)
        orig_end2 = int(orig_end2)
        if ((org,chromo,orig_start1,orig_end1),(org2,chromo2,orig_start2,orig_end2)) in done or ((org2,chromo2,orig_start2,orig_end2),(org,chromo,orig_start1,orig_end1)) in done:
            continue
        else:
            done[((org,chromo,orig_start1,orig_end1),(org2,chromo2,orig_start2,orig_end2))] = 1
            done[((org2,chromo2,orig_start2,orig_end2),(org,chromo,orig_start1,orig_end1))] = 1
        cont = 0
        for coord in coords[org2][chromo2]:
            if (start2 >= coord[0] and start2 <= coord[1]) or (end2 >= coord[0] and end2 <= coord[1]) or (start2 <= coord[0] and end2 >= coord[1]) or (start2 >= coord[0] and end2 <= coord[1]):
                cont = 1 
                break
        if cont == 1:
            continue
        name2 = org2+' '+chromo2+ ' ' + ori[(org2,chromo2)]
        name2 = name2.replace('_',' ')
        if name2 not in shifts:
            continue
        shift2 = shifts[name2]
        start2 -= shift2
        end2 -= shift2
        if ori_alignment == 'reverse' and ori[(org,chromo)] == ori[(org2,chromo2)] or ori_alignment == 'forward' and ori[(org,chromo)] != ori[(org2,chromo2)]:
            start2,end2 = end2,start2
        if ori[(org2,chromo2)] == 'reverse':
            start2,end2 = turn(start2,end2,genome_sizes[name2])

        try:
            for track in gv.feature_tracks:
                if track.name != name1:
                    continue
                for segment in track.segments:
                    segstart = int(segment.start)
                    segend = int(segment.end)
                    if start1 < segstart or end1 > segend:
                        continue
                    else:
                        seg1 = segment.name
            for track in gv.feature_tracks:
                if track.name != name2:
                    continue
                for segment in track.segments:
                    segstart = int(segment.start)
                    segend = int(segment.end)
                    if start2 < segstart or end2 > segend:
                        continue
                    else:
                        seg2 = segment.name
            gv.add_link((name1,seg1,start1,end1),(name2,seg2,start2,end2), color="skyblue", inverted_color="lime", curve=True)
        except:
            print(f'link FROM {name1} {start1} {end1} TO {name2} {start2} {end2} not adjacent')

for org in orgs:
    tables[org].close()

run(f'mkdir -p images_draw_verbose && mkdir -p images_draw_verbose/{margin} && mkdir -p images_draw_verbose/{margin}/',shell=True)
fig = gv.plotfig()
fig.savefig(f'images_draw_verbose/{margin}/{ref}_from_table_{argv[1].split("/")[-1]}_verbose_with_alignments_with_prot_labels.svg')
