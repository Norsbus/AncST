#! /usr/bin/env python3

import pickle
from statistics import mean,stdev,median
from numpy import quantile
from subprocess import run
from bisect import bisect_left
from pprint import pprint
from statistics import mean

""" collecting clusters going through all matches but in the end only saving wanted ones"""
""" most likely in the future there will be an error that some org does not have compressed
    map...if there are transitive matches which introduce new orgs that are not found by 
    looking only once in all matches of all orgs in orgs list they will not have a compressed 
    map since in get_all_orgs.py it is impossible to recursively go through matches without 
    heuristics in make_clusters.py """

def reconcile_overlapping():
    clusters_involved_glob = set()
    for org in orgs:
        l = sorted(list(cluster_i_bib[org].keys()))
        seen = dict()
        seen[cluster_i_bib[org][l[0]]] = set([l[0]])
        to_merge = set()
        cases = 0
        clusters_involved = set()
        containing = []
        for i,x in enumerate(l[1:]):
            if cluster_i_bib[org][x] in seen and cluster_i_bib[org][x] != cluster_i_bib[org][l[i-1]]:
                y = sorted(seen[cluster_i_bib[org][x]])
                idx1 = l.index(y[0])
                idx2 = l.index(y[-1])
                no_eles = idx2-idx1+1
                containing.append(no_eles)
                #pprint(clusters[cluster_i_bib[org][x]])
                #pprint(clusters[cluster_i_bib[org][l[i-1]]])
                #to_merge.add(())
                print(cluster_i_bib[org][x],cluster_i_bib[org][l[i-1]])
                cases += 1
                clusters_involved.add(cluster_i_bib[org][x])
                clusters_involved_glob.add(cluster_i_bib[org][x])
            elif cluster_i_bib[org][x] == cluster_i_bib[org][l[i-1]]:
                seen[cluster_i_bib[org][x]].add(x)
            else:
                seen[cluster_i_bib[org][x]] = set([x])
        print(org,len(l),len(clusters_involved))
    print(len(clusters_involved_glob))
    print(mean(containing))

    return 0

def majority_ori(l):
    f = l.count('forward')
    r = l.count('reverse')
    if f >= r:
        return('forward')
    else:
        return('reverse')

def check_rolling_ori(r_ori,ori):
    if r_ori == 'forward' and ori == 'forward':
        return('forward')
    if r_ori == 'reverse' and ori == 'reverse':
        return('forward')
    if r_ori == 'reverse' and ori == 'forward':
        return('reverse')
    if r_ori == 'forward' and ori == 'reverse':
        return('reverse')

def get_score(cluster):
    scores = []
    lens = []
    for org,s in cluster.items():
        for i in s:
            if i not in am[org]:
                continue
            for org2,bib in am[org][i]['matches'].items():
                if org2 in cluster:
                    for j,score_length_ori in bib.items():
                        scores.append(int(score_length_ori[0]))
                        lens.append(int(score_length_ori[1][1])-int(score_length_ori[1][0]))
    if len(scores) == 0:
        return('not valid due to ambiguous metadata when making compressed maps')
    return(int(mean(scores)),int(mean(lens)))


def get_matches_with_scores(org,i):
    matches = {}
    if i not in am[org]:
        problematic[org][i] = True
        return(matches)
    if i in problematic[org] or i in cluster_i_bib[org]:
        return(matches)
    for org2,bib in am[org][i]['matches'].items():
        if org2 not in orgs:
            continue
        matches[org2] = set()
        for j,score_length_ori in bib.items():
            if j not in am[org2]:
                problematic[org2][j] = True
                continue
            if j in problematic[org2] or j in cluster_i_bib[org2]:
                continue
            matches[org2].add((j,score_length_ori[0]))
    return matches
def get_matches(org,i):
    matches = {}
    if i not in am[org]:
        problematic[org][i] = True
        return(matches)
    if i in problematic[org] or i in cluster_i_bib[org]:
        return(matches)
    for org2,bib in am[org][i]['matches'].items():
        if org2 not in orgs:
            continue
        matches[org2] = set()
        for j,score_length_ori in bib.items():
            if j not in am[org2]:
                problematic[org2][j] = True
                continue
            if j in problematic[org2] or j in cluster_i_bib[org2]:
                continue
            matches[org2].add((j,score_length_ori[2]))
    return matches

def rec_collect_matches(start_org,start_i):
    orientations = {}
    orientations[start_org] = 'forward'
    bib = {}
    old = set()
    orig_matches = get_matches(start_org,start_i)
    new = set()
    diff = set()
    ranges = {}
    out = set()
    for org2,s in orig_matches.items():
        if org2 not in ranges:
            ranges[org2] = [1e15,0]
        if org2 not in orientations:
            orientations[org2] = []
        break_loop = False
        for t in s:
            j,ori = t
            if break_loop:
                break
            if j < ranges[org2][0]:
                ranges[org2][0] = j
            if j > ranges[org2][1]:
                ranges[org2][1] = j
            if ranges[org2][1] - ranges[org2][0] > 100000:
                out.add(org2)
                if org2 in bib:
                    del bib[org2]
                break_loop = True
                for j2 in s:
                    problematic[org2][j2] = True
                continue
            if org2 in bib:
                bib[org2].add(j)
            else:
                bib[org2] = set([j])
        
            orientations[org2].append(ori)
            new.add((org2,j))
            if (org2,j) not in old:
                diff.add((org2,j))

    while len(new) > len(old):
        for org2,j in diff:
            if org2 in out:
                continue
            if org2 != start_org:
                rolling_orientation = majority_ori(orientations[org2])
            else:
                rolling_orientation = 'forward'
            matches = get_matches(org2,j)
            break_loop = False
            for org2,s in matches.items():
                if org2 not in orientations:
                    orientations[org2] = []
                if break_loop:
                    break
                if org2 in out:
                    continue
                if org2 not in ranges:
                    ranges[org2] = [1e15,0]
                for t in s:
                    j,ori = t
                    if j < ranges[org2][0]:
                        ranges[org2][0] = j
                    if j > ranges[org2][1]:
                        ranges[org2][1] = j
                    if ranges[org2][1] - ranges[org2][0] > 100000:
                        out.add(org2)
                        if org2 in bib:
                            del bib[org2]
                        break_loop = True
                        for j2 in s:
                            problematic[org2][j2] = True
                        continue
                    if org2 != start_org:
                        orientations[org2].append(check_rolling_ori(rolling_orientation,ori))
                    if org2 in bib:
                        bib[org2].add(j)
                    else:
                        bib[org2] = set([j])
        old = new.copy()
        new = set()
        diff = set()
        for org2,s in bib.items():
            if org2 in out:
                continue
            for j in s:
                new.add((org2,j))
                if (org2,j) not in old:
                    diff.add((org2,j))

    not_valid = 0
    to_look = set()
    for o,s in bib.items():
        l = sorted(list(s))
        for iti,ele in enumerate(l[:-1]):
            ele_end = ele + am[o][ele]['end'] - am[o][ele]['start']
            ele2 = l[iti+1]
            if ele2-ele_end > 10000 or am[o][ele]['chromosome'] != am[o][ele2]['chromosome']:
                to_look.add(o)
                not_valid +=1
                break
    for o in to_look:
        #for i in bib[o]:
        #    problematic[o][i] = True
        del bib[o]
    #new_bib = {}
    #for o in to_look:
    #    new_bib[o] = set()
    #    tot = {}
    #    for i in bib[o]:
    #        tot[i] = 0
    #        matches = get_matches_with_scores(o,i)
    #        for org2,s2 in matches.items():
    #            for t in s2:
    #                j,score = t
    #                if org2 in bib and j in bib[org2]:
    #                    tot[i] += score
    #    reihe = sorted([(score,i) for i,score in tot.items()],reverse=True)
    #    for score,i in reihe:
    #        if len(new_bib[o]) == 0:
    #            new_bib[o].add(i)
    #        else:
    #            eles = new_bib[o].copy()
    #            for ele in eles:
    #                dis = abs(ele-i)
    #                if dis < 10000:
    #                    new_bib[o].add(i)
    #for o in bib:
    #    if o in to_look:
    #        continue
    #    if o not in new_bib:
    #        new_bib[o] = bib[o]
    #bib = new_bib
                
    #if not_valid == 0:
    if len(bib) > 0:
        return(bib,orientations)
    else:
        return(None)


if __name__ == "__main__":
    orgs = []
    problematic = {}
    for line in open('orgs'):
        orgs.append(line.strip())
    orgs_all_matches = []
    for line in open('orgs_all_matches'):
        orgs_all_matches.append(line.strip())
    
    am = {}

    for org in orgs_all_matches:
        with open('compressed_maps/'+org,'rb') as f:
            am[org] = pickle.load(f)

    clusters = {}
    count = 1
    cluster_i_bib = {}
    for org in orgs_all_matches:
        cluster_i_bib[org] = {}
        problematic[org] = {}

    for org in orgs:
        print('=========')
        print('///',org,'///')
        print(len(am[org]))
        for i,bib in am[org].items():
            if i in cluster_i_bib[org] or i in problematic[org]:
                continue
            res = rec_collect_matches(org,i)
            if res == None:
                continue
            cluster,orientations = res
            if len(cluster) < 2:
                continue
            score_cluster = {}
            for o,s in cluster.items():
                if len(s) == 0:
                    continue
                score_cluster[o] = s
            if len(score_cluster) < 2:
                continue
            else:
                res = get_score(score_cluster)
                if res == 'not valid due to ambiguous metadata when making compressed maps':
                    continue
                else:
                    rep_score,rep_length = res
            new_cluster = {}
            for o,s in score_cluster.items():
                if o not in orgs:
                    continue
                else:
                    new_cluster[o] = {'matches':s,'orientation':majority_ori(orientations[o])}
            if len(new_cluster) < 2:
                continue
            for o,bib in new_cluster.items():
                s = bib['matches']
                #for j in s:
                #    closest_lower = bisect_left(list(cluster_i_bib[o].keys()),j) - 1
                #    if closest_lower >= 0:
                #        print(clusters[cluster_i_bib[o][list(cluster_i_bib[o].keys())[closest_lower]]])
                #        print(new_cluster)
                #    else:
                #        print('smallest element')
                #    closest_upper = bisect_left(list(cluster_i_bib[o].keys()),j)
                #    if closest_upper < len(list(cluster_i_bib[o].keys())):
                #        print(clusters[cluster_i_bib[o][list(cluster_i_bib[o].keys())[closest_upper]]])
                #        print(new_cluster)
                #    else:
                #        print('largest element')
                for j in s:
                    cluster_i_bib[o][j] = count
            new_cluster['focal_orientation'] = org
            new_cluster['representative_score'] = rep_score
            new_cluster['representative_length'] = rep_length
            clusters[count] = new_cluster
            count += 1
        print('=========')
        print(len(problematic[org]))
        print(len(cluster_i_bib[org]))
        print(len(clusters))


    reconcile_overlapping()

    counter = {}
    with open(f'clusters','wb') as f:
        pickle.dump(clusters,f)
    with open(f'i_bib','wb') as f:
        pickle.dump(cluster_i_bib,f)
    for count,cluster in clusters.items():
        del cluster['representative_score']
        del cluster['representative_length']
        del cluster['focal_orientation']
        if len(cluster) in counter:
            counter[len(cluster)] += 1
        else:
            counter[len(cluster)] = 1
    for l,c in sorted(counter.items()):
        print(f'there are {c} clusters with {l} members')
    print('==================')
