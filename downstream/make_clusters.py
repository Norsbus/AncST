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


def get_matches(org,i,problematic):
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
        if len(matches[org2]) == 0:
            del matches[org2]
    return matches

def which_iss(org,iss,bib,max_range):
    iss = sorted(iss)
    scores = {}
    out_l = []
    for i in iss:
        scores[i] = 0
        if i not in am[org][i]:
            continue
        factor = max(len([j for j in iss if abs(j-i) < 20000]),1)
        for org2,s in bib.items():
            if org2 not in am[org][i]['matches']:
                continue
            for j in s:
                if j in am[org][i]['matches'][org2]:
                    scores[i] += factor * am[org][i]['matches'][org2][j][0]

    for i,s in sorted(scores.items(), key=lambda item: item[1], reverse=True):
        iss.remove(i)
        out_l.append(i)
        iss = sorted(iss)
        out_of_range = 0
        for iti,ele in enumerate(iss[:-1]):
            if ele not in am[org]:
                out_of_range = 1
                break
            ele_end = ele + am[org][ele]['end'] - am[org][ele]['start']
            ele2 = iss[iti+1]
            if ele2 not in am[org]:
                out_of_range = 1
                break
            if ele2-ele_end > max_range or am[org][ele]['chromosome'] != am[org][ele2]['chromosome']:
                out_of_range = 1
                break
            
        if out_of_range == 0:
            return(iss,out_l)

    return(iss,out_l)



            


def rec_collect_matches(start_org,start_i):
    orientations = {}
    orientations[start_org] = 'forward'
    max_range = 20000
    bib = {}
    old = set()
    out = {}
    for org in orgs:
        out[org] = {}
    orig_matches = get_matches(start_org,start_i,out)
    new = set()
    diff = set()
    #print('-------------------------------')
    #print(f'first matches from {start_org} {start_i}')
    #pprint(orig_matches)
    out_of_range = 0
    for org2,s in orig_matches.items():
        iss = sorted([x[0] for x in s])
        for iti,ele in enumerate(iss[:-1]):
            ele_end = ele + am[org2][ele]['end'] - am[org2][ele]['start']
            ele2 = iss[iti+1]
            if ele2-ele_end > max_range or am[org2][ele]['chromosome'] != am[org2][ele2]['chromosome']:
                out_of_range = 1
                break
        if out_of_range == 1:
            out_of_range = 0
            print(f'{org2} out of range 1')
            print(iss)
            iss,out_l = which_iss(org2,iss,bib,max_range)
            print(iss,out_l)
            for outer in out_l: 
                out[org2][outer] = 1
            if org2 in bib:
                if len(iss) == 0:
                    del bib[org2]
                else:
                    bib[org2] = set(iss)
            continue
        if org2 not in orientations:
            orientations[org2] = []
        for t in s:
            j,ori = t
            if org2 in bib:
                bib[org2].add(j)
            else:
                bib[org2] = set([j])
        
            orientations[org2].append(ori)
            new.add((org2,j))
            if (org2,j) not in old:
                diff.add((org2,j))

        max_end = am[org2][max(iss)]['end']
        min_start = am[org2][min(iss)]['start']
        if max_end - min_start > max_range:
            max_range = min(100000,max_end - min_start)

    while len(new) > len(old):
        out_of_range = 0
        for org2,j in diff:
            if j in out[org2]:
                continue
            if org2 != start_org:
                rolling_orientation = majority_ori(orientations[org2])
            else:
                rolling_orientation = 'forward'
            matches = get_matches(org2,j,out)
            #print(f'next matches from {org2} {j}')
            #print(matches)
            for org22,s in matches.items():
                iss = [x[0] for x in s]
                if org22 in bib:
                    iss += list(bib[org22])
                iss = sorted(iss)
                for iti,ele in enumerate(iss[:-1]):
                    ele_end = ele + am[org22][ele]['end'] - am[org22][ele]['start']
                    ele2 = iss[iti+1]
                    if ele2-ele_end > max_range or am[org22][ele]['chromosome'] != am[org22][ele2]['chromosome']:
                        out_of_range = 1
                        break
                if out_of_range == 1:
                    out_of_range = 0
                    print(f'{org22} out of range 2')
                    print(iss)
                    iss,out_l = which_iss(org22,iss,bib,max_range)
                    print(iss,out_l)
                    for outer in out_l: 
                        out[org22][outer] = 1
                    if org22 in bib:
                        if len(iss) == 0:
                            del bib[org22]
                        else:
                            bib[org22] = set(iss)
                    continue
                if org22 not in orientations:
                    orientations[org22] = []
                if org22 != start_org:
                    orientations[org22].append(check_rolling_ori(rolling_orientation,ori))
                for t in s:
                    j,ori = t
                    if org22 in bib:
                        bib[org22].add(j)
                    else:
                        bib[org22] = set([j])
                max_end = am[org22][max(iss)]['end']
                min_start = am[org22][min(iss)]['start']
                if max_end - min_start > max_range:
                    max_range = min(100000,max_end - min_start)

        old = new.copy()
        new = set()
        diff = set()
        for org2,s in bib.items():
            for j in s:
                if j in out[org2]:
                    continue
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
            if ele2-ele_end > max_range or am[o][ele]['chromosome'] != am[o][ele2]['chromosome']:
                to_look.add(o)
                not_valid +=1
                break
    for o in to_look:
        print(f'{o} out of range 3')
        print(sorted(list(bib[o])))
        del bib[o]
    if len(bib) > 0:
        return(bib,orientations)
    else:
        return(None)


if __name__ == "__main__":
    orgs = []
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

    for org in orgs:
        print('=========')
        print('///',org,'///')
        print(len(am[org]))
        for i,bib in am[org].items():
            if i in cluster_i_bib[org]:
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
            #with open(f'single_clusters_max_100000/cluster_{count}','wb') as f:
            #    pickle.dump(new_cluster,f)
        print('=========')
        print(len(cluster_i_bib[org]))
        print(len(clusters))


    #reconcile_overlapping()

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
