#! /usr/bin/env python3

import argparse
import os
import pickle
import bisect
from collections import defaultdict
from Bio import SeqIO
import networkx as nx
from subprocess import run
import sys
import multiprocessing as mp
import random

# =============================================================================
# COGRAPH EDITOR (requires tralda package: pip install tralda)
# =============================================================================

from tralda.datastructures.Tree import Tree, TreeNode
from tralda.cograph.Cograph import to_cograph

def edit_to_cograph(G, run_number=10, forbidden_edge_checker=None):
    

    ce = CographEditor(G)
    best_cotree = ce.cograph_edit()

    result_graph = to_cograph(best_cotree)

    # Copy edge attributes from original graph to result graph
    for u, v in result_graph.edges():
        if G.has_edge(u, v):
            result_graph[u][v].update(G[u][v])

    # Remove forbidden edges if checker provided
    if forbidden_edge_checker:
        edges_to_remove = []
        for u, v in result_graph.edges():
            if forbidden_edge_checker(u, v):
                edges_to_remove.append((u, v))

        if edges_to_remove:
            result_graph.remove_edges_from(edges_to_remove)

    return result_graph

class CographEditor:
    

    def __init__(self, G):
        if not isinstance(G, nx.Graph):
            raise TypeError('not a NetworkX Graph')

        self.G = G
        self.V = [v for v in G.nodes()]

        self.cotrees = []
        self.costs = []

        self.best_T = None
        self.best_cost = float('inf')

    def cograph_edit(self, run_number=10):
        

        for i in range(run_number):

            T = Tree(None)
            self.aux_counter = {}
            already_in_T = set()
            leaf_map = {}
            total_cost = 0

            if i > 0:
                random.shuffle(self.V)

            self._start_tree(T, already_in_T, leaf_map)

            if len(self.V) > 2:
                for x in self.V[2:]:
                    cost, x_node = self._insert(x, T, already_in_T, leaf_map)

                    current = x_node.parent
                    while current:
                        self.aux_counter[current] += 1
                        current = current.parent

                    leaf_map[x] = x_node
                    already_in_T.add(x)
                    total_cost += cost

            self.cotrees.append(T)
            self.costs.append(total_cost)

            if total_cost < self.best_cost:
                self.best_T = T
                self.best_cost = total_cost

            if self.best_cost <= 0:
                break

        return self.best_T

    def _start_tree(self, T, already_in_T, leaf_map):

        if len(self.V) == 0:
            raise RuntimeError('empty graph in cograph editing')
            return

        elif len(self.V) == 1:
            T.root = TreeNode(label=self.V[0])
            return

        v1, v2 = self.V[0], self.V[1]
        already_in_T.update([v1, v2])

        R = TreeNode(label='series')
        self.aux_counter[R] = 2
        T.root = R

        if self.G.has_edge(v1, v2):
            v1_node = TreeNode(label=v1)
            self.aux_counter[v1_node] = 1
            v2_node = TreeNode(label=v2)
            self.aux_counter[v2_node] = 1
            R.add_child(v1_node)
            R.add_child(v2_node)
        else:
            N = TreeNode(label='parallel')
            self.aux_counter[N] = 2
            R.add_child(N)
            v1_node = TreeNode(label=v1)
            self.aux_counter[v1_node] = 1
            v2_node = TreeNode(label=v2)
            self.aux_counter[v2_node] = 1
            N.add_child(v1_node)
            N.add_child(v2_node)

        leaf_map[v1] = v1_node
        leaf_map[v2] = v2_node

    def _insert(self, x, T, already_in_T, leaf_map):

        C_nh = {}
        C_nh_number = {}

        Nx_number = {}
        completion_forced = {}
        full = {}
        deletion_forced = {}

        Nx_counter = 0

        for v in self.G.neighbors(x):
            if v not in already_in_T:
                continue

            Nx_counter += 1

            v = leaf_map[v]
            C_nh[v]                = []
            C_nh_number[v]         = 0
            Nx_number[v]           = 1
            completion_forced[v]   = True
            full[v]                = True
            deletion_forced[v]     = False

            current = v
            while current.parent:
                if current.parent in C_nh:
                    C_nh[current.parent].append(current)
                    break
                else:
                    C_nh[current.parent] = [current]
                current = current.parent

        if self.aux_counter[T.root] == Nx_counter:
            R = T.root
            x_node = TreeNode(label=x)
            self.aux_counter[x_node] = 1
            R.add_child(x_node)
            leaf_map[x] = x_node
            return 0, x_node
        elif Nx_counter == 0:
            if len(T.root.children) == 1:
                N = T.root.children[0]
                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                N.add_child(x_node)
            else:
                R_old = T.root
                R_new = TreeNode(label='series')
                self.aux_counter[R_new] = self.aux_counter[R_old]
                N = TreeNode(label='parallel')
                self.aux_counter[N] = self.aux_counter[R_old]
                R_new.add_child(N)
                N.add_child(R_old)
                T.root = R_new

                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                N.add_child(x_node)
            return 0, x_node

        stack = [T.root]
        while stack:
            u = stack.pop()
            if not u.children:
                continue
            elif u not in C_nh_number:
                stack.append(u)
                C_nh_number[u] = len(C_nh[u])
                for nh_child in C_nh[u]:
                    stack.append(nh_child)
            else:
                completion_forced_counter = 0
                deletion_forced_counter = 0
                full_counter = 0
                Nx_number[u] = 0
                for nh_child in C_nh[u]:
                    Nx_number[u] += Nx_number[nh_child]
                    if completion_forced[nh_child]:
                        completion_forced_counter += 1
                    if deletion_forced[nh_child]:
                        deletion_forced_counter += 1
                    if full[nh_child]:
                        full_counter += 1
                full[u] = True if Nx_number[u] == self.aux_counter[u] else False

                if (full[u] or
                    (u.label == 'parallel' and len(u.children) == C_nh_number[u]) or
                    (u.label == 'series' and len(u.children) == completion_forced_counter)):
                    completion_forced[u] = True
                else:
                    completion_forced[u] = False

                if ((u.label == 'series' and full_counter == 0) or
                    (u.label == 'parallel' and C_nh_number[u] == deletion_forced_counter)):
                    deletion_forced[u] = True
                else:
                    deletion_forced[u] = False

        cost_above = {}
        C_mixed = {}
        C_full = {}

        stack = [T.root]
        while stack:
            u = stack.pop()
            C_mixed[u] = []
            C_full[u] = []
            for nh_child in C_nh[u]:
                if not full[nh_child]:
                    stack.append(nh_child)
                    C_mixed[u].append(nh_child)
                else:
                    C_full[u].append(nh_child)
            v = u.parent
            if not v:
                cost_above[u] = 0
            elif u.parent.label == 'parallel':
                cost_above[u] = cost_above[v] + Nx_number[v] - Nx_number[u]
            else:
                cost_above[u] = cost_above[v] + ((self.aux_counter[v] - Nx_number[v]) -
                                                 (self.aux_counter[u] - Nx_number[u]))

        mincost = {}

        for u in cost_above.keys():

            if len(u.children) == 2 and u.label == 'parallel':
                for i in range(2):
                    v, v2 = u.children[i], u.children[-1-i]
                    if v in completion_forced and completion_forced[v]:
                        Nx_v2 = Nx_number[v2] if (v2 in Nx_number) else 0
                        cost = cost_above[u] + (self.aux_counter[v] - Nx_number[v]) + Nx_v2
                        if (u not in mincost) or (cost < mincost[u][0]):
                            mincost[u] = (cost, [v])
            elif len(u.children) == 2 and u.label == 'series':
                for i in range(2):
                    v, v2 = u.children[i], u.children[-1-i]
                    if v not in deletion_forced or deletion_forced[v]:
                        Nx_v = Nx_number[v] if (v in Nx_number) else 0
                        Nx_v2 = Nx_number[v2] if (v2 in Nx_number) else 0
                        cost = cost_above[u] + Nx_v + (self.aux_counter[v2] - Nx_v2)
                        if (u not in mincost) or (cost < mincost[u][0]):
                            mincost[u] = (cost, [v2])

            elif len(u.children) >= 3:
                red = set()
                blue = set()

                if u.label == 'parallel':
                    if C_nh_number[u] == 1:
                        v = C_nh[u][0]
                        if completion_forced[v]:
                            cost = cost_above[u] + (self.aux_counter[v] - Nx_number[v])
                            mincost[u] = (cost, [v])
                        continue

                    for v in C_mixed[u]:
                        if Nx_number[v] >= self.aux_counter[v] - Nx_number[v]:
                            red.add(v)
                        else:
                            blue.add(v)

                    if C_nh_number[u] == len(u.children) and not blue:
                        current_min, current_min_node = float('inf'), None
                        for v in red:
                            diff = 2 * Nx_number[v] - self.aux_counter[v]
                            if diff < current_min:
                                current_min, current_min_node = diff, v
                        red.remove(current_min_node)
                        blue.add(current_min_node)

                    nb_filled = len(C_full[u]) + len(red)
                    if nb_filled < 2:
                        for i in range(2 - nb_filled):
                            current_min, current_min_node = float('inf'), None
                            for v in blue:
                                diff = self.aux_counter[v] - 2 * Nx_number[v]
                                if diff < current_min:
                                    current_min, current_min_node = diff, v
                            blue.remove(current_min_node)
                            red.add(current_min_node)
                else:
                    if len(u.children) - len(C_full[u]) == 1:
                        v = None
                        for child in u.children:
                            if (child not in full) or (not full[child]):
                                v = child
                                break
                        if (v not in deletion_forced) or (deletion_forced[v]):
                            Nx_v = Nx_number[v] if v in Nx_number else 0
                            cost = cost_above[u] + Nx_v
                            mincost[u] = (cost, C_full[u])
                        continue

                    for v in C_mixed[u]:
                        if self.aux_counter[v] - Nx_number[v] >= Nx_number[v]:
                            blue.add(v)
                        else:
                            red.add(v)

                    if not C_full[u] and not red:
                        current_min, current_min_node = float('inf'), None
                        for v in blue:
                            diff = self.aux_counter[v] - 2 * Nx_number[v]
                            if diff < current_min:
                                current_min, current_min_node = diff, v
                        blue.remove(current_min_node)
                        red.add(current_min_node)

                    nb_hollowed = len(u.children) - C_nh_number[u] + len(blue)
                    if nb_hollowed < 2:
                        for i in range(2 - nb_hollowed):
                            current_min, current_min_node = float('inf' ), None
                            for v in red:
                                diff = 2 * Nx_number[v] - self.aux_counter[v]
                                if diff < current_min:
                                    current_min, current_min_node = diff, v
                            red.remove(current_min_node)
                            blue.add(current_min_node)

                cost = cost_above[u]
                for v in red:
                    cost += (self.aux_counter[v] - Nx_number[v])
                for v in blue:
                    cost += Nx_number[v]

                mincost[u] = (cost, C_full[u] + list(red))

        insertion_mincost, settling_node = float('inf'), None
        for u in mincost.keys():
            if mincost[u][0] < insertion_mincost:
                insertion_mincost, settling_node = mincost[u][0], u

        u = settling_node
        filled = mincost[u][1]

        if u.label == 'parallel' and len(filled) == 1:
            w = filled[0]
            if w.is_leaf():
                new_node = TreeNode(label='series')
                self.aux_counter[new_node] = 1
                u.remove_child(w)
                u.add_child(new_node)
                new_node.add_child(w)

                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                new_node.add_child(x_node)
            else:
                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                w.add_child(x_node)

        elif (u.label == 'series' and
              len(u.children) - len(filled) == 1):
            set_A = set(filled)
            w = None
            for child in u.children:
                if child not in set_A:
                    w = child
                    break
            if w.is_leaf():
                new_node = TreeNode(label='parallel')
                self.aux_counter[new_node] = 1
                u.remove_child(w)
                u.add_child(new_node)
                new_node.add_child(w)

                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                new_node.add_child(x_node)
            else:
                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                w.add_child(x_node)

        else:
            filled_aux_counter = 0
            y = TreeNode(label=u.label)
            for a in filled:
                u.remove_child(a)
                y.add_child(a)
                filled_aux_counter += self.aux_counter[a]
            self.aux_counter[y] = filled_aux_counter

            if u.label == 'parallel':
                new_node = TreeNode(label='series')
                self.aux_counter[new_node] = filled_aux_counter
                u.add_child(new_node)

                new_node.add_child(y)
                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                new_node.add_child(x_node)
            else:
                par = u.parent
                if par is not None:
                    par.remove_child(u)
                    par.add_child(y)
                else:
                    T.root = y

                new_node = TreeNode(label='parallel')
                self.aux_counter[new_node] = 0
                y.add_child(new_node)
                new_node.add_child(u)
                x_node = TreeNode(label=x)
                self.aux_counter[x_node] = 1
                new_node.add_child(x_node)

                self.aux_counter[y] = self.aux_counter[u]
                self.aux_counter[u] -= filled_aux_counter
                self.aux_counter[new_node] = self.aux_counter[u]

        return insertion_mincost, x_node

# =============================================================================
# GENOMIC ALIGNMENT MODULE (align_genomic_simple_tuples)
# =============================================================================

# Scoring parameters
MATCH_SCORE = 2.0
GAP_PENALTY = -1.0
TANDEM_DUP_PENALTY = 1.0
SEG_DUP_PENALTY = 1.0
INVERSION_PENALTY = 1.0
ORIENTATION_FLIP_PENALTY = 0.1

MAX_TANDEM_LEN = 20
MAX_SEG_DUP_LEN = 10
MAX_INVERSION_LEN = 10

# Operation codes
MATCH = 1
GAP_SEQ1 = 2
GAP_SEQ2 = 3
TANDEM_DUP_1_TO_Q = 4
TANDEM_DUP_P_TO_1 = 5
SEG_DUP_SEQ1 = 6
SEG_DUP_SEQ2 = 7
INVERSION = 8
SEG_DUP_SEQ1_INVERTED = 9
SEG_DUP_SEQ2_INVERTED = 10

def parse_element_align(element):
    
    if isinstance(element, tuple):
        original_id, element_string = element
    else:
        original_id = None
        element_string = element

    if len(element_string) < 2:
        raise ValueError(f"Element too short: {element_string}")

    orientation = element_string[-1]
    base_name = element_string[:-1]

    if orientation not in ['+', '-']:
        raise ValueError(f"Element must end with + or -: {element_string}")

    return (original_id, base_name, orientation)

def parse_sequence_align(seq):
    
    return [parse_element_align(s) for s in seq]

def score_match_align(elem1, elem2):
    
    orig_id1, base1, ori1 = elem1
    orig_id2, base2, ori2 = elem2

    if base1 == base2:
        return MATCH_SCORE
    else:
        return float('-inf')

def score_tandem_dup_1_to_q_align(elem_a, seq_b_slice):
    
    q = len(seq_b_slice)
    orig_id_a, base_a, ori_a = elem_a

    orientation_mismatches = 0
    for elem_b in seq_b_slice:
        orig_id_b, base_b, ori_b = elem_b
        if base_a != base_b:
            return float('-inf')
        if ori_a != ori_b:
            orientation_mismatches += 1

    total_match = MATCH_SCORE * q
    penalty = TANDEM_DUP_PENALTY * (q - 1)
    penalty += ORIENTATION_FLIP_PENALTY * orientation_mismatches
    return total_match - penalty

def score_tandem_dup_p_to_1_align(seq_a_slice, elem_b):
    
    p = len(seq_a_slice)
    orig_id_b, base_b, ori_b = elem_b

    orientation_mismatches = 0
    for elem_a in seq_a_slice:
        orig_id_a, base_a, ori_a = elem_a
        if base_a != base_b:
            return float('-inf')
        if ori_a != ori_b:
            orientation_mismatches += 1

    total_match = MATCH_SCORE * p
    penalty = TANDEM_DUP_PENALTY * (p - 1)
    penalty += ORIENTATION_FLIP_PENALTY * orientation_mismatches
    return total_match - penalty

def score_segmental_dup_seq1_align(seq1, i, block_len, seq2, j, num_copies):
    
    target_block = seq2[j - block_len : j]

    orientation_mismatches = 0
    for copy_idx in range(num_copies):
        copy_start = i - (num_copies - copy_idx) * block_len
        copy_end = i - (num_copies - copy_idx - 1) * block_len
        block = seq1[copy_start : copy_end]

        for k in range(block_len):
            orig_id_b, base_b, ori_b = block[k]
            orig_id_t, base_t, ori_t = target_block[k]

            if base_b != base_t:
                return float('-inf')
            if ori_b != ori_t:
                orientation_mismatches += 1

    total_match = MATCH_SCORE * (num_copies * block_len)
    penalty = SEG_DUP_PENALTY * block_len * (num_copies - 1)
    penalty += ORIENTATION_FLIP_PENALTY * orientation_mismatches
    return total_match - penalty

def score_segmental_dup_seq2_align(seq1, i, seq2, j, block_len, num_copies):
    
    target_block = seq1[i - block_len : i]

    orientation_mismatches = 0
    for copy_idx in range(num_copies):
        copy_start = j - (num_copies - copy_idx) * block_len
        copy_end = j - (num_copies - copy_idx - 1) * block_len
        block = seq2[copy_start : copy_end]

        for k in range(block_len):
            orig_id_b, base_b, ori_b = block[k]
            orig_id_t, base_t, ori_t = target_block[k]

            if base_b != base_t:
                return float('-inf')
            if ori_b != ori_t:
                orientation_mismatches += 1

    total_match = MATCH_SCORE * (num_copies * block_len)
    penalty = SEG_DUP_PENALTY * block_len * (num_copies - 1)
    penalty += ORIENTATION_FLIP_PENALTY * orientation_mismatches
    return total_match - penalty

def score_inversion_align(seq1, i, inv_len, seq2, j):
    
    seq1_block = seq1[i - inv_len : i]
    seq2_block = seq2[j - inv_len : j]
    seq2_reversed = seq2_block[::-1]

    for k in range(inv_len):
        orig_id1, base1, ori1 = seq1_block[k]
        orig_id2, base2, ori2 = seq2_reversed[k]

        if base1 != base2:
            return float('-inf')

        if ori1 == ori2:
            return float('-inf')

    total_match = MATCH_SCORE * inv_len
    penalty = INVERSION_PENALTY * inv_len
    return total_match - penalty

def score_segmental_dup_seq1_inverted_align(seq1, i, block_len, seq2, j, num_copies):
    
    target_block = seq2[j - block_len : j]

    for copy_idx in range(num_copies):
        copy_start = i - (num_copies - copy_idx) * block_len
        copy_end = i - (num_copies - copy_idx - 1) * block_len
        block = seq1[copy_start : copy_end]
        block_reversed = block[::-1]

        for k in range(block_len):
            orig_id_b, base_b, ori_b = block_reversed[k]
            orig_id_t, base_t, ori_t = target_block[k]

            if base_b != base_t:
                return float('-inf')

            if ori_b == ori_t:
                return float('-inf')

    total_match = MATCH_SCORE * (num_copies * block_len)
    penalty = SEG_DUP_PENALTY * block_len * (num_copies - 1)
    penalty += INVERSION_PENALTY * block_len * num_copies
    return total_match - penalty

def score_segmental_dup_seq2_inverted_align(seq1, i, seq2, j, block_len, num_copies):
    
    target_block = seq1[i - block_len : i]

    for copy_idx in range(num_copies):
        copy_start = j - (num_copies - copy_idx) * block_len
        copy_end = j - (num_copies - copy_idx - 1) * block_len
        block = seq2[copy_start : copy_end]
        block_reversed = block[::-1]

        for k in range(block_len):
            orig_id_b, base_b, ori_b = block_reversed[k]
            orig_id_t, base_t, ori_t = target_block[k]

            if base_b != base_t:
                return float('-inf')

            if ori_b == ori_t:
                return float('-inf')

    total_match = MATCH_SCORE * (num_copies * block_len)
    penalty = SEG_DUP_PENALTY * block_len * (num_copies - 1)
    penalty += INVERSION_PENALTY * block_len * num_copies
    return total_match - penalty

class AlignGenomicSimpleTuples:
    

    @staticmethod
    def align(seq1_strings, seq2_strings):
        
        seq1 = parse_sequence_align(seq1_strings)
        seq2 = parse_sequence_align(seq2_strings)

        m = len(seq1)
        n = len(seq2)

        score = [[0 for _ in range(n+1)] for _ in range(m+1)]
        pointer = [[None for _ in range(n+1)] for _ in range(m+1)]

        for i in range(m+1):
            score[i][0] = GAP_PENALTY * i
            pointer[i][0] = (GAP_SEQ1, 1)

        for j in range(n+1):
            score[0][j] = GAP_PENALTY * j
            pointer[0][j] = (GAP_SEQ2, 1)

        pointer[0][0] = (0, 0)

        for i in range(1, m+1):
            for j in range(1, n+1):
                candidates = []

                match_sc = score[i-1][j-1] + score_match_align(seq1[i-1], seq2[j-1])
                candidates.append((match_sc, (MATCH, 1, 1)))

                gap_sc = score[i-1][j] + GAP_PENALTY
                candidates.append((gap_sc, (GAP_SEQ1, 1)))

                gap_sc = score[i][j-1] + GAP_PENALTY
                candidates.append((gap_sc, (GAP_SEQ2, 1)))

                for p in range(2, min(i+1, MAX_TANDEM_LEN+1)):
                    dup_sc = score[i-p][j-1] + score_tandem_dup_p_to_1_align(seq1[i-p:i], seq2[j-1])
                    if dup_sc == float('-inf'):
                        break
                    candidates.append((dup_sc, (TANDEM_DUP_P_TO_1, p, 1)))

                for q in range(2, min(j+1, MAX_TANDEM_LEN+1)):
                    dup_sc = score[i-1][j-q] + score_tandem_dup_1_to_q_align(seq1[i-1], seq2[j-q:j])
                    if dup_sc == float('-inf'):
                        break
                    candidates.append((dup_sc, (TANDEM_DUP_1_TO_Q, 1, q)))

                for block_len in range(2, min(j+1, MAX_SEG_DUP_LEN+1)):
                    for num_copies in range(2, min(i // block_len + 1, MAX_TANDEM_LEN+1)):
                        if i >= num_copies * block_len and j >= block_len:
                            seg_sc = score[i - num_copies*block_len][j - block_len] + \
                                     score_segmental_dup_seq1_align(seq1, i, block_len, seq2, j, num_copies)
                            if seg_sc == float('-inf'):
                                break
                            candidates.append((seg_sc, (SEG_DUP_SEQ1, num_copies*block_len, block_len, block_len, num_copies)))

                for block_len in range(2, min(i+1, MAX_SEG_DUP_LEN+1)):
                    for num_copies in range(2, min(j // block_len + 1, MAX_TANDEM_LEN+1)):
                        if i >= block_len and j >= num_copies * block_len:
                            seg_sc = score[i - block_len][j - num_copies*block_len] + \
                                     score_segmental_dup_seq2_align(seq1, i, seq2, j, block_len, num_copies)
                            if seg_sc == float('-inf'):
                                break
                            candidates.append((seg_sc, (SEG_DUP_SEQ2, block_len, num_copies*block_len, block_len, num_copies)))

                for inv_len in range(2, min(i+1, j+1, MAX_INVERSION_LEN+1)):
                    inv_sc = score[i - inv_len][j - inv_len] + \
                             score_inversion_align(seq1, i, inv_len, seq2, j)
                    if inv_sc > float('-inf'):
                        candidates.append((inv_sc, (INVERSION, inv_len, inv_len, inv_len)))

                for block_len in range(2, min(j+1, MAX_SEG_DUP_LEN+1)):
                    for num_copies in range(2, min(i // block_len + 1, MAX_TANDEM_LEN+1)):
                        if i >= num_copies * block_len and j >= block_len:
                            seg_sc = score[i - num_copies*block_len][j - block_len] + \
                                     score_segmental_dup_seq1_inverted_align(seq1, i, block_len, seq2, j, num_copies)
                            if seg_sc == float('-inf'):
                                break
                            candidates.append((seg_sc, (SEG_DUP_SEQ1_INVERTED, num_copies*block_len, block_len, block_len, num_copies)))

                for block_len in range(2, min(i+1, MAX_SEG_DUP_LEN+1)):
                    for num_copies in range(2, min(j // block_len + 1, MAX_TANDEM_LEN+1)):
                        if i >= block_len and j >= num_copies * block_len:
                            seg_sc = score[i - block_len][j - num_copies*block_len] + \
                                     score_segmental_dup_seq2_inverted_align(seq1, i, seq2, j, block_len, num_copies)
                            if seg_sc == float('-inf'):
                                break
                            candidates.append((seg_sc, (SEG_DUP_SEQ2_INVERTED, block_len, num_copies*block_len, block_len, num_copies)))

                best_score, best_op = max(candidates, key=lambda x: x[0])
                score[i][j] = best_score
                pointer[i][j] = best_op

        aligned1, aligned2, operations = AlignGenomicSimpleTuples._traceback(pointer, seq1, seq2)

        final_score = score[m][n]
        return aligned1, aligned2, operations, final_score

    @staticmethod
    def _traceback(pointer, seq1, seq2):
        
        m = len(seq1)
        n = len(seq2)
        i, j = m, n

        aligned1 = []
        aligned2 = []
        operations = []

        while i > 0 or j > 0:
            if i == 0:
                aligned1.append(None)
                aligned2.append(seq2[j-1])
                operations.append(('gap_seq1', None, seq2[j-1]))
                j -= 1
                continue

            if j == 0:
                aligned1.append(seq1[i-1])
                aligned2.append(None)
                operations.append(('gap_seq2', seq1[i-1], None))
                i -= 1
                continue

            op = pointer[i][j]
            op_type = op[0]

            if op_type == MATCH:
                elem1 = seq1[i-1]
                elem2 = seq2[j-1]
                aligned1.append(elem1)
                aligned2.append(elem2)

                orig_id1, base1, ori1 = elem1
                orig_id2, base2, ori2 = elem2
                if ori1 == ori2:
                    operations.append(('match', elem1, elem2, 'same_orientation'))
                else:
                    operations.append(('match', elem1, elem2, 'orientation_flip'))

                i -= 1
                j -= 1

            elif op_type == GAP_SEQ1:
                aligned1.append(seq1[i-1])
                aligned2.append(None)
                operations.append(('gap_seq2', seq1[i-1], None))
                i -= 1

            elif op_type == GAP_SEQ2:
                aligned1.append(None)
                aligned2.append(seq2[j-1])
                operations.append(('gap_seq1', None, seq2[j-1]))
                j -= 1

            elif op_type == TANDEM_DUP_P_TO_1:
                p = op[1]
                seq1_dups = seq1[i-p:i]
                elem2 = seq2[j-1]

                for elem1 in seq1_dups:
                    aligned1.append(elem1)
                    aligned2.append(elem2 if elem1 == seq1_dups[-1] else "DUP")

                operations.append(('tandem_dup_seq1', seq1_dups, elem2, p))
                i -= p
                j -= 1

            elif op_type == TANDEM_DUP_1_TO_Q:
                q = op[2]
                elem1 = seq1[i-1]
                seq2_dups = seq2[j-q:j]

                for elem2 in seq2_dups:
                    aligned1.append(elem1 if elem2 == seq2_dups[-1] else "DUP")
                    aligned2.append(elem2)

                operations.append(('tandem_dup_seq2', elem1, seq2_dups, q))
                i -= 1
                j -= q

            elif op_type == SEG_DUP_SEQ1:
                block_len = op[3]
                num_copies = op[4]
                seq2_block = seq2[j - block_len : j]

                seq1_blocks = []
                for copy_idx in range(num_copies):
                    copy_start = i - (num_copies - copy_idx) * block_len
                    copy_end = i - (num_copies - copy_idx - 1) * block_len
                    seq1_block = seq1[copy_start : copy_end]
                    seq1_blocks.append(seq1_block)

                    for k in range(block_len):
                        aligned1.append(seq1_block[k])
                        aligned2.append(seq2_block[k] if copy_idx == num_copies - 1 else "DUP")

                operations.append(('seg_dup_seq1', seq1_blocks, seq2_block, block_len, num_copies))
                i -= num_copies * block_len
                j -= block_len

            elif op_type == SEG_DUP_SEQ2:
                block_len = op[3]
                num_copies = op[4]
                seq1_block = seq1[i - block_len : i]

                seq2_blocks = []
                for copy_idx in range(num_copies):
                    copy_start = j - (num_copies - copy_idx) * block_len
                    copy_end = j - (num_copies - copy_idx - 1) * block_len
                    seq2_block = seq2[copy_start : copy_end]
                    seq2_blocks.append(seq2_block)

                    for k in range(block_len):
                        aligned1.append(seq1_block[k] if copy_idx == num_copies - 1 else "DUP")
                        aligned2.append(seq2_block[k])

                operations.append(('seg_dup_seq2', seq1_block, seq2_blocks, block_len, num_copies))
                i -= block_len
                j -= num_copies * block_len

            elif op_type == INVERSION:
                inv_len = op[3]
                seq1_block = seq1[i - inv_len : i]
                seq2_block = seq2[j - inv_len : j]
                seq2_reversed = seq2_block[::-1]

                for k in range(inv_len):
                    aligned1.append(seq1_block[k])
                    aligned2.append(seq2_reversed[k])

                operations.append(('inversion', seq1_block, seq2_block, inv_len))
                i -= inv_len
                j -= inv_len

            elif op_type == SEG_DUP_SEQ1_INVERTED:
                block_len = op[3]
                num_copies = op[4]
                seq2_block = seq2[j - block_len : j]

                seq1_blocks = []
                for copy_idx in range(num_copies):
                    copy_start = i - (num_copies - copy_idx) * block_len
                    copy_end = i - (num_copies - copy_idx - 1) * block_len
                    seq1_block = seq1[copy_start : copy_end]
                    seq1_blocks.append(seq1_block)
                    seq1_reversed = seq1_block[::-1]

                    for k in range(block_len):
                        aligned1.append(seq1_reversed[k])
                        aligned2.append(seq2_block[k] if copy_idx == num_copies - 1 else "DUP")

                operations.append(('seg_dup_seq1_inverted', seq1_blocks, seq2_block, block_len, num_copies))
                i -= num_copies * block_len
                j -= block_len

            elif op_type == SEG_DUP_SEQ2_INVERTED:
                block_len = op[3]
                num_copies = op[4]
                seq1_block = seq1[i - block_len : i]

                seq2_blocks = []
                for copy_idx in range(num_copies):
                    copy_start = j - (num_copies - copy_idx) * block_len
                    copy_end = j - (num_copies - copy_idx - 1) * block_len
                    seq2_block = seq2[copy_start : copy_end]
                    seq2_blocks.append(seq2_block)
                    seq2_reversed = seq2_block[::-1]

                    for k in range(block_len):
                        aligned1.append(seq1_block[k] if copy_idx == num_copies - 1 else "DUP")
                        aligned2.append(seq2_reversed[k])

                operations.append(('seg_dup_seq2_inverted', seq1_block, seq2_blocks, block_len, num_copies))
                i -= block_len
                j -= num_copies * block_len

            else:
                raise ValueError(f"Unknown operation: {op_type}")

        aligned1.reverse()
        aligned2.reverse()
        operations.reverse()

        return aligned1, aligned2, operations

# Create module-level align function
align_genomic_simple_tuples = AlignGenomicSimpleTuples()

# =============================================================================
# MAPPING FUNCTION
# =============================================================================

def get_mapping(mapping_file=None):
    
    org_mapping = {}
    chr_mapping = {}
    if mapping_file:
        ff = mapping_file
    elif os.path.exists('../../utils/mapping'):
        ff = '../../utils/mapping'
    elif os.path.exists('mapping'):
        ff = 'mapping'
    else:
        raise FileNotFoundError("No mapping file found")
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
    return(org_mapping, chr_mapping)

parser = argparse.ArgumentParser(description='Simple Synthology Pipeline (BLASTN)')
parser.add_argument('--cores', type=int, default=1, help="number of cores to run with")
parser.add_argument('--col', type=str, default="MCScanX.collinearity", help="name of MCScanX collinearity file")
parser.add_argument('--homology_threshold', type=int, default=40, help="BLAST bit score threshold")
parser.add_argument('--annotation_dir', type=str, default="annotations_for_synthology", help='directory with gff/ and seqs/ subdirs')
parser.add_argument('--new_blast', type=str, default="n", help='y to run BLAST, n to use cache')
parser.add_argument('--mapping_file', type=str, default=None, help='path to mapping file (optional, uses get_mapping() if not specified)')
parser.add_argument('--pairwise_alignments', type=str, default="../../utils/pairwise_alignments_table", help='path to pairwise alignments table file')
parser.add_argument('--use-clasp', type=str, default="y", help='y to use CLASP chaining for BLAST hits (default), n for raw BLAST')
parser.add_argument('--blast-mode', type=str, choices=['chains', 'all-vs-all'], default='chains',
                    help='chains: BLAST only within synteny (fast); all-vs-all: BLAST all genes for evaluation (slow)')
parser.add_argument('--blastn-path', type=str, default="blastn", help='path to blastn binary (default: blastn from PATH)')
parser.add_argument('--clasp-path', type=str, default="../../utils/clasp.x", help='path to clasp.x binary (default: ../../utils/clasp.x)')
args = parser.parse_args()

OUTPUT_DIR = "synthology_out_nucl"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/to_comp", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/to_comp_out", exist_ok=True)
if args.use_clasp == 'y':
    os.makedirs(f"{OUTPUT_DIR}/clasp_in", exist_ok=True)
    os.makedirs(f"{OUTPUT_DIR}/clasp_out", exist_ok=True)

BLAST_THRESHOLD = args.homology_threshold
MCSCANX_FILE = args.col
GFF_DIR = os.path.join(args.annotation_dir, 'gff/')
SEQ_DIR = os.path.join(args.annotation_dir, 'seqs/')
NEW_BLAST = args.new_blast == 'y'
NUM_CORES = args.cores
USE_CLASP = args.use_clasp == 'y'
BLAST_MODE = args.blast_mode
BLASTN_PATH = args.blastn_path
CLASP_PATH = args.clasp_path
PAIRWISE_ALIGNMENTS_FILE = args.pairwise_alignments

def parse_blast_output(blast_file):
    
    hits = []
    with open(blast_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) >= 12:
                hits.append((cols[0], cols[1], float(cols[11])))
                            # qseqid  sseqid   bitscore
    return hits

def parse_clasp_output(clasp_file):
    
    hits = []
    with open(clasp_file) as f:
        for line in f:
            cols = line.strip().split()
            if cols[0] == '#':
                continue
            hits.append((cols[1], cols[2], float(cols[7])))
                        # qseqid  sseqid   chained_score
    return hits

def validate_unique_gene_ids(parsed):
    
    for species in parsed:
        seen_ids = set()
        duplicates = []

        for chromo in parsed[species]:
            for start, end, strand, pid in parsed[species][chromo]:
                if pid in seen_ids:
                    duplicates.append(pid)
                else:
                    seen_ids.add(pid)

        if duplicates:
            raise ValueError(
                f"ERROR: Duplicate gene IDs found in {species}:\n"
                f"  {', '.join(sorted(set(duplicates))[:10])}"
                f"{' ...' if len(set(duplicates)) > 10 else ''}\n"
                f"  Total duplicates: {len(set(duplicates))}\n"
                f"  All gene IDs must be unique within each species."
            )

def filter_overlapping_chains(chains):
    
    from collections import defaultdict

    # Step 1: Flatten chains to list (avoid bidirectional duplicates)
    chain_list = []
    seen_keys = set()

    for org1 in chains:
        for chr1 in chains[org1]:
            for (mini1, maxi1) in chains[org1][chr1]:
                for org2 in chains[org1][chr1][(mini1, maxi1)]:
                    chr2, mini2, maxi2, orientation = chains[org1][chr1][(mini1, maxi1)][org2]

                    # Canonical key to avoid duplicates
                    chain_key = tuple(sorted([
                        (org1, chr1, mini1, maxi1),
                        (org2, chr2, mini2, maxi2)
                    ]))

                    if chain_key not in seen_keys:
                        seen_keys.add(chain_key)
                        chain_list.append({
                            'org1': org1, 'chr1': chr1, 'mini1': mini1, 'maxi1': maxi1,
                            'org2': org2, 'chr2': chr2, 'mini2': mini2, 'maxi2': maxi2,
                            'orientation': orientation
                        })

    print(f"  Total unique chains before filtering: {len(chain_list)}")

    # Step 2: Group by species pair
    species_pair_chains = defaultdict(list)
    for idx, chain in enumerate(chain_list):
        species_pair = tuple(sorted([chain['org1'], chain['org2']]))
        species_pair_chains[species_pair].append((idx, chain))

    # Step 3: Find overlapping chains
    chains_to_remove = set()
    one_sided_overlaps = []
    both_sided_overlaps = []

    for species_pair, chains_in_pair in species_pair_chains.items():
        sp1, sp2 = species_pair

        for i in range(len(chains_in_pair)):
            for j in range(i+1, len(chains_in_pair)):
                idx_a, ca = chains_in_pair[i]
                idx_b, cb = chains_in_pair[j]

                # Skip if either already marked for removal
                if idx_a in chains_to_remove or idx_b in chains_to_remove:
                    continue

                # Get regions for each species
                regions_a = [(ca['org1'], ca['chr1'], ca['mini1'], ca['maxi1']),
                            (ca['org2'], ca['chr2'], ca['mini2'], ca['maxi2'])]
                regions_b = [(cb['org1'], cb['chr1'], cb['mini1'], cb['maxi1']),
                            (cb['org2'], cb['chr2'], cb['mini2'], cb['maxi2'])]

                # Find regions for species1 and species2
                reg_s1_a = [r for r in regions_a if r[0] == sp1][0]
                reg_s1_b = [r for r in regions_b if r[0] == sp1][0]
                reg_s2_a = [r for r in regions_a if r[0] == sp2][0]
                reg_s2_b = [r for r in regions_b if r[0] == sp2][0]

                # Check overlap on species1
                ov_s1 = (reg_s1_a[1] == reg_s1_b[1] and  # same chromosome
                        reg_s1_a[2] <= reg_s1_b[3] and  # mini1_a <= maxi1_b
                        reg_s1_b[2] <= reg_s1_a[3])     # mini1_b <= maxi1_a

                # Check overlap on species2
                ov_s2 = (reg_s2_a[1] == reg_s2_b[1] and  # same chromosome
                        reg_s2_a[2] <= reg_s2_b[3] and  # mini2_a <= maxi2_b
                        reg_s2_b[2] <= reg_s2_a[3])     # mini2_b <= maxi2_a

                if ov_s1 and ov_s2:
                    # Both-sided overlap: keep longest
                    len_a = (ca['maxi1'] - ca['mini1']) + (ca['maxi2'] - ca['mini2'])
                    len_b = (cb['maxi1'] - cb['mini1']) + (cb['maxi2'] - cb['mini2'])

                    if len_a >= len_b:
                        chains_to_remove.add(idx_b)
                        both_sided_overlaps.append((idx_a, idx_b, 'kept_a'))
                    else:
                        chains_to_remove.add(idx_a)
                        both_sided_overlaps.append((idx_a, idx_b, 'kept_b'))

                elif ov_s1 or ov_s2:
                    # One-sided overlap: remove BOTH
                    chains_to_remove.add(idx_a)
                    chains_to_remove.add(idx_b)
                    one_sided_overlaps.append((idx_a, idx_b, ov_s1, ov_s2))

    # Report statistics
    print(f"  Found {len(both_sided_overlaps)} both-sided overlaps (kept longest)")
    print(f"  Found {len(one_sided_overlaps)} one-sided overlaps (removed both)")

    if one_sided_overlaps:
        print(f"\n  WARNING: {len(one_sided_overlaps)} one-sided overlaps detected!")
        print(f"  These chains overlap on one chromosome but not the other - suspicious!")
        print(f"  Both chains in each pair will be discarded.")

        # Show a few examples
        for i, (idx_a, idx_b, ov_s1, ov_s2) in enumerate(one_sided_overlaps[:3]):
            ca = chain_list[idx_a]
            cb = chain_list[idx_b]
            print(f"\n    Example {i+1}:")
            print(f"      Chain A: {ca['org1']}:{ca['chr1']}[{ca['mini1']}-{ca['maxi1']}] <-> {ca['org2']}:{ca['chr2']}[{ca['mini2']}-{ca['maxi2']}]")
            print(f"      Chain B: {cb['org1']}:{cb['chr1']}[{cb['mini1']}-{cb['maxi1']}] <-> {cb['org2']}:{cb['chr2']}[{cb['mini2']}-{cb['maxi2']}]")
            print(f"      Overlaps on species1: {ov_s1}, Overlaps on species2: {ov_s2}")

        if len(one_sided_overlaps) > 3:
            print(f"    ... and {len(one_sided_overlaps) - 3} more")

    # Step 4: Build filtered chain list
    filtered_chain_list = [chain for idx, chain in enumerate(chain_list) if idx not in chains_to_remove]

    print(f"  Chains after filtering: {len(filtered_chain_list)}")
    print(f"  Removed {len(chains_to_remove)} chains total")

    # Step 5: Rebuild chains dict with bidirectional entries
    filtered_chains = {}

    for chain in filtered_chain_list:
        org1 = chain['org1']
        chr1 = chain['chr1']
        mini1 = chain['mini1']
        maxi1 = chain['maxi1']
        org2 = chain['org2']
        chr2 = chain['chr2']
        mini2 = chain['mini2']
        maxi2 = chain['maxi2']
        orientation = chain['orientation']

        # Add org1 -> org2 direction
        if org1 not in filtered_chains:
            filtered_chains[org1] = {}
        if chr1 not in filtered_chains[org1]:
            filtered_chains[org1][chr1] = {}
        if (mini1, maxi1) not in filtered_chains[org1][chr1]:
            filtered_chains[org1][chr1][(mini1, maxi1)] = {}
        filtered_chains[org1][chr1][(mini1, maxi1)][org2] = (chr2, mini2, maxi2, orientation)

        # Add org2 -> org1 direction
        if org2 not in filtered_chains:
            filtered_chains[org2] = {}
        if chr2 not in filtered_chains[org2]:
            filtered_chains[org2][chr2] = {}
        if (mini2, maxi2) not in filtered_chains[org2][chr2]:
            filtered_chains[org2][chr2][(mini2, maxi2)] = {}
        filtered_chains[org2][chr2][(mini2, maxi2)][org1] = (chr1, mini1, maxi1, orientation)

    return filtered_chains

def run_blast_for_chain(args_tuple):
    i, chain, output_dir, seq_dir, new_blast, use_clasp, blastn_path, clasp_path = args_tuple

    # Define file paths
    blast_out_1 = f"{output_dir}/to_comp_out/{i}_1"
    blast_out_2 = f"{output_dir}/to_comp_out/{i}_2"

    # Check cache based on mode
    if not new_blast:
        if use_clasp:
            # Check if CLASP outputs exist
            clasp_files = [f'{output_dir}/clasp_out/{i}_{suffix}' for suffix in ['for_1', 'rev_1', 'for_2', 'rev_2']]
            if all(os.path.exists(f) for f in clasp_files):
                hits = []
                for suffix in ['for_1', 'rev_1', 'for_2', 'rev_2']:
                    hits.extend(parse_clasp_output(f'{output_dir}/clasp_out/{i}_{suffix}'))
                return i, hits
        else:
            # Check if BLAST outputs exist
            if os.path.exists(blast_out_1) and os.path.exists(blast_out_2):
                hits = []
                hits.extend(parse_blast_output(blast_out_1))
                hits.extend(parse_blast_output(blast_out_2))
                return i, hits

    # Prepare FASTA files
    faa1 = f"{output_dir}/to_comp/{i}_1.faa"
    faa2 = f"{output_dir}/to_comp/{i}_2.faa"

    pids1 = [pid for s, e, strand, pid in chain['genes1']]
    pids2 = [pid for s, e, strand, pid in chain['genes2']]

    seq_path1 = None
    seq_path2 = None
    for ext in ['.fasta', '.fna']:
        if os.path.exists(f"{seq_dir}/{chain['org1']}{ext}"):
            seq_path1 = f"{seq_dir}/{chain['org1']}{ext}"
        if os.path.exists(f"{seq_dir}/{chain['org2']}{ext}"):
            seq_path2 = f"{seq_dir}/{chain['org2']}{ext}"

    if not seq_path1 or not seq_path2:
        return i, []

    seqs1 = []
    for record in SeqIO.parse(seq_path1, "fasta"):
        if record.id.split()[0] in pids1:
            seqs1.append(record)

    seqs2 = []
    for record in SeqIO.parse(seq_path2, "fasta"):
        if record.id.split()[0] in pids2:
            seqs2.append(record)

    SeqIO.write(seqs1, faa1, "fasta")
    SeqIO.write(seqs2, faa2, "fasta")

    # Run BLAST in both directions
    run([blastn_path, "-query", faa1, "-subject", faa2, "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-out", blast_out_1, "-num_threads", "1","-word_size","4"], capture_output=True)
    run([blastn_path, "-query", faa2, "-subject", faa1, "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-out", blast_out_2, "-num_threads", "1","-word_size","4"], capture_output=True)

    if use_clasp:
        # Process BLAST 1 results (org1 → org2): split by strand
        forward_1 = []
        reverse_1 = []
        with open(blast_out_1) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 12:
                    # Check if reverse strand: sstart (cols[8]) > send (cols[9])
                    if int(cols[8]) > int(cols[9]):
                        cols[8],cols[9] = cols[9],cols[8]
                        line = '\t'.join(cols)
                        reverse_1.append(line+'\n')
                    else:
                        forward_1.append(line+'\n')

        # Process BLAST 2 results (org2 → org1): split by strand
        forward_2 = []
        reverse_2 = []
        with open(blast_out_2) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 12:
                    # Check if reverse strand: sstart (cols[8]) > send (cols[9])
                    if int(cols[8]) > int(cols[9]):
                        cols[8],cols[9] = cols[9],cols[8]
                        line = '\t'.join(cols)
                        reverse_2.append(line+'\n')
                    else:
                        forward_2.append(line+'\n')

        # Run CLASP on all 4 sets
        with open(f"{output_dir}/clasp_in/{i}_for_1","w") as f:
            f.writelines(forward_1)
        run(f'{clasp_path} -m -i {output_dir}/clasp_in/{i}_for_1 -c 7 8 9 10 12 -C 1 2 -l 0 -e 0 -o {output_dir}/clasp_out/{i}_for_1',shell=True)

        with open(f"{output_dir}/clasp_in/{i}_rev_1","w") as f:
            f.writelines(reverse_1)
        run(f'{clasp_path} -m -i {output_dir}/clasp_in/{i}_rev_1 -c 7 8 9 10 12 -C 1 2 -l 0 -e 0 -o {output_dir}/clasp_out/{i}_rev_1',shell=True)

        with open(f"{output_dir}/clasp_in/{i}_for_2","w") as f:
            f.writelines(forward_2)
        run(f'{clasp_path} -m -i {output_dir}/clasp_in/{i}_for_2 -c 7 8 9 10 12 -C 1 2 -l 0 -e 0 -o {output_dir}/clasp_out/{i}_for_2',shell=True)

        with open(f"{output_dir}/clasp_in/{i}_rev_2","w") as f:
            f.writelines(reverse_2)
        run(f'{clasp_path} -m -i {output_dir}/clasp_in/{i}_rev_2 -c 7 8 9 10 12 -C 1 2 -l 0 -e 0 -o {output_dir}/clasp_out/{i}_rev_2',shell=True)

        # Parse CLASP outputs
        hits = []
        for suffix in ['for_1', 'rev_1', 'for_2', 'rev_2']:
            hits.extend(parse_clasp_output(f'{output_dir}/clasp_out/{i}_{suffix}'))
    else:
        # Parse raw BLAST outputs directly
        hits = []
        hits.extend(parse_blast_output(blast_out_1))
        hits.extend(parse_blast_output(blast_out_2))

    return i, hits

def run_all_vs_all_blast(args_tuple):

    sp1, sp2, fasta1, fasta2, output_dir, use_clasp, blastn_path, clasp_path = args_tuple

    pair_id = f"{sp1}___{sp2}"
    blast_out_1 = f"{output_dir}/all_vs_all/blast_{pair_id}_1"
    blast_out_2 = f"{output_dir}/all_vs_all/blast_{pair_id}_2"

    # Run BLAST both directions
    run([blastn_path, "-query", fasta1, "-subject", fasta2, "-outfmt",
         "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
         "-out", blast_out_1, "-num_threads", "1", "-word_size", "4"], capture_output=True)
    run([blastn_path, "-query", fasta2, "-subject", fasta1, "-outfmt",
         "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
         "-out", blast_out_2, "-num_threads", "1", "-word_size", "4"], capture_output=True)

    if use_clasp:
        # Process and run CLASP (same as run_blast_for_chain)
        forward_1, reverse_1 = [], []
        with open(blast_out_1) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 12:
                    if int(cols[8]) > int(cols[9]):
                        cols[8], cols[9] = cols[9], cols[8]
                        reverse_1.append('\t'.join(cols) + '\n')
                    else:
                        forward_1.append(line)

        forward_2, reverse_2 = [], []
        with open(blast_out_2) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 12:
                    if int(cols[8]) > int(cols[9]):
                        cols[8], cols[9] = cols[9], cols[8]
                        reverse_2.append('\t'.join(cols) + '\n')
                    else:
                        forward_2.append(line)

        os.makedirs(f"{output_dir}/all_vs_all/clasp_in", exist_ok=True)
        os.makedirs(f"{output_dir}/all_vs_all/clasp_out", exist_ok=True)

        for suffix, data in [('for_1', forward_1), ('rev_1', reverse_1),
                              ('for_2', forward_2), ('rev_2', reverse_2)]:
            with open(f"{output_dir}/all_vs_all/clasp_in/{pair_id}_{suffix}", "w") as f:
                f.writelines(data)
            run(f'{clasp_path} -m -i {output_dir}/all_vs_all/clasp_in/{pair_id}_{suffix} '
                f'-c 7 8 9 10 12 -C 1 2 -l 0 -e 0 -o {output_dir}/all_vs_all/clasp_out/{pair_id}_{suffix}',
                shell=True)

        # Parse CLASP outputs
        hits = []
        for suffix in ['for_1', 'rev_1', 'for_2', 'rev_2']:
            hits.extend(parse_clasp_output(f'{output_dir}/all_vs_all/clasp_out/{pair_id}_{suffix}'))
    else:
        # Parse raw BLAST
        hits = []
        hits.extend(parse_blast_output(blast_out_1))
        hits.extend(parse_blast_output(blast_out_2))

    # Convert to dict
    edges_dict = {}
    for qseqid, sseqid, score in hits:
        edge = tuple(sorted([qseqid, sseqid]))
        edges_dict[edge] = max(edges_dict.get(edge, 0), score)

    return edges_dict

def build_graph_from_blast(args_tuple):
    i, chain, blast_hits, blast_threshold = args_tuple

    G = nx.Graph()

    for s, e, strand, pid in chain['genes1']:
        G.add_node(pid)
    for s, e, strand, pid in chain['genes2']:
        G.add_node(pid)

    # Add BLAST edges (NetworkX handles duplicate edges automatically - just keeps one)
    for qseqid, sseqid, bitscore in blast_hits:
        if bitscore >= blast_threshold:
            if G.has_node(qseqid) and G.has_node(sseqid):
                G.add_edge(qseqid, sseqid)

    edges_before = G.number_of_edges()
    original_edges = set(G.edges())

    def is_forbidden(u, v):
        return (u, v) not in original_edges and (v, u) not in original_edges

    G = edit_to_cograph(G, forbidden_edge_checker=is_forbidden)

    edges_after = G.number_of_edges()

    return G, chain, edges_before, edges_after

# =============================================================================
# HELPER FUNCTIONS FOR ALIGNMENT NODE AUGMENTATION
# =============================================================================

def find_alignments_in_range(alignments, ref_species, ref_chr, tgt_species,
                               ref_min, ref_max, tgt_min, tgt_max,
                               existing_genes, tolerance=10000, min_score_after_cut=50):
    
    # Check if this ref_species/chr/tgt_species combination exists
    if (ref_species not in alignments or
        ref_chr not in alignments[ref_species] or
        tgt_species not in alignments[ref_species][ref_chr]):
        return []

    aln_list = alignments[ref_species][ref_chr][tgt_species]

    # Expand ranges by tolerance
    ref_search_min = ref_min - tolerance
    ref_search_max = ref_max + tolerance
    tgt_search_min = tgt_min - tolerance
    tgt_search_max = tgt_max + tolerance

    # Use bisect to find alignments that overlap with ref range
    # Find first alignment that could overlap (start < search_max)
    start_idx = bisect.bisect_left([a[0] for a in aln_list], ref_search_min)

    valid_alignments = []

    for i in range(start_idx, len(aln_list)):
        aln = aln_list[i]
        ref_start, ref_end, tgt_chr, tgt_start, tgt_end, score, orientation = aln

        # Stop if we're past the search range
        if ref_start > ref_search_max:
            break

        # Check if alignment overlaps with search range
        if ref_end < ref_search_min:
            continue

        # Check if target coordinates are in range
        if not (tgt_search_min <= tgt_start <= tgt_search_max or
                tgt_search_min <= tgt_end <= tgt_search_max or
                (tgt_start <= tgt_search_min and tgt_end >= tgt_search_max)):
            continue

        # Cut alignment to fit within range (without tolerance for the cut)
        cut_ref_start = max(ref_start, ref_min)
        cut_ref_end = min(ref_end, ref_max)

        # Calculate how much was cut
        original_length = ref_end - ref_start
        cut_length = cut_ref_end - cut_ref_start
        length_cut = original_length - cut_length

        # Check if score is still good enough after cutting
        adjusted_score = score - (length_cut * 2)
        if adjusted_score < min_score_after_cut:
            continue

        # Check if it overlaps with existing genes
        overlaps = False
        for gene_start, gene_end, *_ in existing_genes:
            if not (cut_ref_end < gene_start or cut_ref_start > gene_end):
                overlaps = True
                break

        if overlaps:
            continue

        # Add valid alignment
        valid_alignments.append((
            ref_species, ref_chr, cut_ref_start, cut_ref_end,
            tgt_species, tgt_chr, tgt_start, tgt_end,
            adjusted_score, orientation
        ))

    return valid_alignments

def augment_graph_with_alignments(graph, region_info, alignments, gene_mapping,
                                    alignments_mapping, tolerance=10000):
    
    if not alignments:
        return graph

    ident, (org1, elements1), (org2, elements2) = region_info

    # Find coordinate ranges for each species
    if not elements1 or not elements2:
        return graph

    # Elements are (start, end, ori, gene_id)
    org1_starts = [e[0] for e in elements1]
    org1_ends = [e[1] for e in elements1]
    org1_min = min(org1_starts)
    org1_max = max(org1_ends)

    org2_starts = [e[0] for e in elements2]
    org2_ends = [e[1] for e in elements2]
    org2_min = min(org2_starts)
    org2_max = max(org2_ends)

    # Get chromosomes from gene_mapping
    org1_chr = None
    org2_chr = None

    for e in elements1:
        gene_id = e[3]
        if gene_id in gene_mapping:
            org1_chr = gene_mapping[gene_id][1]
            break

    for e in elements2:
        gene_id = e[3]
        if gene_id in gene_mapping:
            org2_chr = gene_mapping[gene_id][1]
            break

    if not org1_chr or not org2_chr:
        return graph

    # Search for alignments in one direction only (data is bidirectional)
    # Try org1 as reference first
    alns = find_alignments_in_range(
        alignments, org1, org1_chr, org2,
        org1_min, org1_max, org2_min, org2_max,
        elements1, tolerance=tolerance
    )

    # If no alignments found, try the reverse direction
    if not alns:
        alns = find_alignments_in_range(
            alignments, org2, org2_chr, org1,
            org2_min, org2_max, org1_min, org1_max,
            elements2, tolerance=tolerance
        )

    # Add alignment nodes to graph and populate alignments_mapping
    # IMPORTANT: Create TWO nodes per alignment (one for each species side)
    # This allows proper handling in NW alignment step
    alignments_found = 0
    for aln in alns:
        # Unpack alignment info
        (ref_species, ref_chr, ref_start, ref_end,
         tgt_species, tgt_chr, tgt_start, tgt_end,
         adjusted_score, orientation) = aln

        # Create separate node IDs for each species side
        aln_node_ref = f"aln_{ref_species}_{ref_chr}:{ref_start}-{ref_end}"
        aln_node_tgt = f"aln_{tgt_species}_{tgt_chr}:{tgt_start}-{tgt_end}"

        # Add both nodes to graph
        graph.add_node(aln_node_ref)
        graph.add_node(aln_node_tgt)

        # Connect them with an edge (they represent two sides of same alignment)
        graph.add_edge(aln_node_ref, aln_node_tgt)

        alignments_found += 2  # Added two nodes

        # Populate alignments_mapping for BOTH nodes
        # Format similar to gene_mapping: (org, chromo, start, end, ori, dict())
        # For alignment nodes, ori is stored in dict along with partner info and orientation
        alignments_mapping[aln_node_ref] = (
            ref_species, ref_chr, ref_start, ref_end, None,
            {'partner': aln_node_tgt, 'score': adjusted_score, 'orientation': orientation}
        )
        alignments_mapping[aln_node_tgt] = (
            tgt_species, tgt_chr, tgt_start, tgt_end, None,
            {'partner': aln_node_ref, 'score': adjusted_score, 'orientation': orientation}
        )

    # No need to re-run cograph editing after adding alignment nodes
    # Alignment nodes are only connected to their partner (2-node components)
    # Adding isolated 2-node components cannot violate cograph property

    return graph

def align_graph_worker(args):
    
    idx, G, chain, alignments, gene_mapping = args

    if G.number_of_nodes() == 0:
        empty_alignment_info = {
            'species1_elements': [],
            'species2_elements': [],
            'edges': []
        }
        return (G, chain, empty_alignment_info)

    # Create alignments_mapping for this graph
    alignments_mapping = {}

    # Add alignment nodes to graph
    region_info = (
        idx,
        (chain['org1'], chain['genes1']),
        (chain['org2'], chain['genes2'])
    )
    G = augment_graph_with_alignments(G, region_info, alignments, gene_mapping, alignments_mapping)

    # Use connected components to assign shared homolog IDs
    components = list(nx.connected_components(G))
    node_to_homolog = {}
    for comp_idx, comp in enumerate(components):
        shared_name = f"h{comp_idx}"
        for node in comp:
            node_to_homolog[node] = shared_name

    # Collect and sort elements for org1
    org1_elements = []
    for s, e, strand, pid in chain['genes1']:
        if pid in G.nodes():
            ori = '+' if strand in ['+', 1, '+1'] else '-'
            org1_elements.append(('gene', s, e, ori, pid))

    for aln_node, aln_info in alignments_mapping.items():
        species, chromo, start, end, _, extra = aln_info
        if species == chain['org1'] and aln_node in G.nodes():
            org1_elements.append(('alignment', start, end, None, aln_node))

    org1_elements.sort(key=lambda x: x[1])  # Sort by start position

    # Handle overlaps for org1: cut alignment nodes that overlap with genes
    org1_filtered = []
    for elem in org1_elements:
        if elem[0] == 'gene':
            org1_filtered.append(elem)
        elif elem[0] == 'alignment':
            _, aln_start, aln_end, _, aln_node = elem
            overlaps_removed = [(aln_start, aln_end)]
            for gene in [e for e in org1_elements if e[0] == 'gene']:
                _, gene_start, gene_end, _, _ = gene
                new_ranges = []
                for rng_start, rng_end in overlaps_removed:
                    if gene_end <= rng_start or gene_start >= rng_end:
                        new_ranges.append((rng_start, rng_end))
                    elif gene_start <= rng_start and gene_end >= rng_end:
                        pass
                    elif gene_start > rng_start and gene_end < rng_end:
                        new_ranges.append((rng_start, gene_start))
                        new_ranges.append((gene_end, rng_end))
                    elif gene_start <= rng_start:
                        new_ranges.append((gene_end, rng_end))
                    else:
                        new_ranges.append((rng_start, gene_start))
                overlaps_removed = new_ranges

            if overlaps_removed:
                longest = max(overlaps_removed, key=lambda x: x[1] - x[0])
                if longest[1] - longest[0] >= 50:
                    org1_filtered.append(('alignment', longest[0], longest[1], None, aln_node))

    # Same for org2
    org2_elements = []
    for s, e, strand, pid in chain['genes2']:
        if pid in G.nodes():
            ori = '+' if strand in ['+', 1, '+1'] else '-'
            org2_elements.append(('gene', s, e, ori, pid))

    for aln_node, aln_info in alignments_mapping.items():
        species, chromo, start, end, _, extra = aln_info
        if species == chain['org2'] and aln_node in G.nodes():
            org2_elements.append(('alignment', start, end, None, aln_node))

    org2_elements.sort(key=lambda x: x[1])

    # Handle overlaps for org2
    org2_filtered = []
    for elem in org2_elements:
        if elem[0] == 'gene':
            org2_filtered.append(elem)
        elif elem[0] == 'alignment':
            _, aln_start, aln_end, _, aln_node = elem
            overlaps_removed = [(aln_start, aln_end)]
            for gene in [e for e in org2_elements if e[0] == 'gene']:
                _, gene_start, gene_end, _, _ = gene
                new_ranges = []
                for rng_start, rng_end in overlaps_removed:
                    if gene_end <= rng_start or gene_start >= rng_end:
                        new_ranges.append((rng_start, rng_end))
                    elif gene_start <= rng_start and gene_end >= rng_end:
                        pass
                    elif gene_start > rng_start and gene_end < rng_end:
                        new_ranges.append((rng_start, gene_start))
                        new_ranges.append((gene_end, rng_end))
                    elif gene_start <= rng_start:
                        new_ranges.append((gene_end, rng_end))
                    else:
                        new_ranges.append((rng_start, gene_start))
                overlaps_removed = new_ranges

            if overlaps_removed:
                longest = max(overlaps_removed, key=lambda x: x[1] - x[0])
                if longest[1] - longest[0] >= 50:
                    org2_filtered.append(('alignment', longest[0], longest[1], None, aln_node))

    # Count '-' orientations in genes for orientation decision
    org1_minus_count = sum(1 for e in org1_filtered if e[0] == 'gene' and e[3] == '-')
    org2_minus_count = sum(1 for e in org2_filtered if e[0] == 'gene' and e[3] == '-')

    # Build tuple sequences with proper orientation for alignment nodes
    seq1 = []
    for elem in org1_filtered:
        elem_type, start, end, ori, node_id = elem
        if node_id in node_to_homolog:
            shared_name = node_to_homolog[node_id]
            if elem_type == 'gene':
                seq1.append((node_id, f"{shared_name}{ori}"))
            else:  # alignment node
                aln_info = alignments_mapping[node_id]
                aln_orientation = aln_info[5]['orientation']
                if aln_orientation == 'forward':
                    ori = '+'
                else:  # reverse
                    ori = '-' if org1_minus_count >= org2_minus_count else '+'
                seq1.append((node_id, f"{shared_name}{ori}"))

    seq2 = []
    for elem in org2_filtered:
        elem_type, start, end, ori, node_id = elem
        if node_id in node_to_homolog:
            shared_name = node_to_homolog[node_id]
            if elem_type == 'gene':
                seq2.append((node_id, f"{shared_name}{ori}"))
            else:  # alignment node
                aln_info = alignments_mapping[node_id]
                aln_orientation = aln_info[5]['orientation']
                if aln_orientation == 'forward':
                    ori = '+'
                else:  # reverse
                    ori = '-' if org2_minus_count > org1_minus_count else '+'
                seq2.append((node_id, f"{shared_name}{ori}"))

    if chain['orientation'] == 'reverse':
        seq2 = seq2[::-1]

    if not seq1 or not seq2:
        empty_alignment_info = {
            'species1_elements': [],
            'species2_elements': [],
            'edges': []
        }
        return (G, chain, empty_alignment_info)

    # Run alignment with tuples
    aligned1, aligned2, operations, score = align_genomic_simple_tuples.align(seq1, seq2)

    # Build NEW graph from alignment operations
    new_G = nx.Graph()

    for op in operations:
        op_type = op[0]

        if op_type == 'match':
            _, elem1, elem2, match_type = op
            if elem1 is not None and elem2 is not None:
                orig_id1, base1, ori1 = elem1
                orig_id2, base2, ori2 = elem2
                new_G.add_edge(orig_id1, orig_id2, event='match', match_type=match_type)

        elif op_type == 'tandem_dup_seq1':
            _, seq1_dups, elem2, p = op
            if elem2 is not None:
                orig_id2, base2, ori2 = elem2
                for dup_elem in seq1_dups:
                    orig_id1, base1, ori1 = dup_elem
                    new_G.add_edge(orig_id1, orig_id2, event='tandem_dup', direction='seq1_to_seq2')

        elif op_type == 'tandem_dup_seq2':
            _, elem1, seq2_dups, q = op
            if elem1 is not None:
                orig_id1, base1, ori1 = elem1
                for dup_elem in seq2_dups:
                    orig_id2, base2, ori2 = dup_elem
                    new_G.add_edge(orig_id1, orig_id2, event='tandem_dup', direction='seq2_to_seq1')

        elif op_type == 'inversion':
            _, seq1_block, seq2_block, inv_len = op
            for elem1, elem2 in zip(seq1_block, seq2_block[::-1]):
                orig_id1, base1, ori1 = elem1
                orig_id2, base2, ori2 = elem2
                new_G.add_edge(orig_id1, orig_id2, event='inversion')

        elif op_type == 'seg_dup_seq1':
            _, seq1_blocks, seq2_block, block_len, num_copies = op
            for seq1_block in seq1_blocks:
                for elem1, elem2 in zip(seq1_block, seq2_block):
                    orig_id1, base1, ori1 = elem1
                    orig_id2, base2, ori2 = elem2
                    new_G.add_edge(orig_id1, orig_id2, event='segmental_dup', direction='seq1')

        elif op_type == 'seg_dup_seq2':
            _, seq1_block, seq2_blocks, block_len, num_copies = op
            for seq2_block in seq2_blocks:
                for elem1, elem2 in zip(seq1_block, seq2_block):
                    orig_id1, base1, ori1 = elem1
                    orig_id2, base2, ori2 = elem2
                    new_G.add_edge(orig_id1, orig_id2, event='segmental_dup', direction='seq2')

        elif op_type == 'seg_dup_seq1_inverted':
            _, seq1_blocks, seq2_block, block_len, num_copies = op
            for seq1_block in seq1_blocks:
                for elem1, elem2 in zip(seq1_block[::-1], seq2_block):
                    orig_id1, base1, ori1 = elem1
                    orig_id2, base2, ori2 = elem2
                    new_G.add_edge(orig_id1, orig_id2, event='inverted_segmental_dup', direction='seq1')

        elif op_type == 'seg_dup_seq2_inverted':
            _, seq1_block, seq2_blocks, block_len, num_copies = op
            for seq2_block in seq2_blocks:
                for elem1, elem2 in zip(seq1_block, seq2_block[::-1]):
                    orig_id1, base1, ori1 = elem1
                    orig_id2, base2, ori2 = elem2
                    new_G.add_edge(orig_id1, orig_id2, event='inverted_segmental_dup', direction='seq2')

        elif op_type in ['gap_seq1', 'gap_seq2']:
            # Gaps don't create edges - genes become singletons
            pass

    # Collect alignment information before removing alignment nodes
    species1_elements = [elem[4] for elem in org1_filtered]  # elem[4] is the node_id
    species2_elements = [elem[4] for elem in org2_filtered]
    alignment_edges = [(u, v, dict(attrs)) for u, v, attrs in new_G.edges(data=True)]

    alignment_info = {
        'species1_elements': species1_elements,
        'species2_elements': species2_elements,
        'edges': alignment_edges
    }

    # Remove alignment nodes from new graph (they start with "aln_")
    alignment_nodes = [n for n in new_G.nodes() if isinstance(n, str) and n.startswith('aln_')]
    for node in alignment_nodes:
        new_G.remove_node(node)

    return (new_G, chain, alignment_info)

if __name__ == '__main__':
    print("="*80)
    print("STEP 1: PARSE MCSCANX")
    print("="*80)

    org_mapping, chr_mapping = get_mapping(args.mapping_file)

    chains_path = f"{OUTPUT_DIR}/parsed_mcscanx_chains"
    if os.path.isfile(chains_path):
        print("Loading cached chains...")
        with open(chains_path, 'rb') as f:
            chains = pickle.load(f)
    else:
        print(f"Parsing {MCSCANX_FILE}...")
        chains = {}
        mini1 = 'NA'
        elements1, elements2 = [], []

        with open(MCSCANX_FILE, 'r') as f:
            for line in f:
                if 'Alignment' in line:
                    if mini1 != 'NA' and len(elements1) > 0 and len(elements2) > 0:
                        if (elements1[0] > elements1[-1] and elements2[0] < elements2[-1]) or (elements1[0] < elements1[-1] and elements2[0] > elements2[-1]):
                            sign = 'reverse'
                        else:
                            sign = 'forward'
                        if (mini1, maxi1) not in chains[org1][chr1]:
                            chains[org1][chr1][(mini1, maxi1)] = {}
                        if org2 not in chains[org1][chr1][(mini1, maxi1)]:
                            chains[org1][chr1][(mini1, maxi1)][org2] = (chr2, mini2, maxi2, sign)
                        if (mini2, maxi2) not in chains[org2][chr2]:
                            chains[org2][chr2][(mini2, maxi2)] = {}
                        if org1 not in chains[org2][chr2][(mini2, maxi2)]:
                            chains[org2][chr2][(mini2, maxi2)][org1] = (chr1, mini1, maxi1, sign)

                    elements1, elements2 = [], []
                    mini1, maxi1 = 1e15, 0
                    mini2, maxi2 = 1e15, 0
                    org1 = line.split()[-2].split('&')[0].split('chr')[0]
                    org2 = line.split()[-2].split('&')[1].split('chr')[0]
                    chr1 = line.split()[-2].split('&')[0]
                    chr2 = line.split()[-2].split('&')[1]
                    org1 = org_mapping[org1]
                    org2 = org_mapping[org2]
                    chr1 = chr_mapping[org1][chr1]
                    chr2 = chr_mapping[org2][chr2]
                    if org1 not in chains:
                        chains[org1] = {}
                    if chr1 not in chains[org1]:
                        chains[org1][chr1] = {}
                    if org2 not in chains:
                        chains[org2] = {}
                    if chr2 not in chains[org2]:
                        chains[org2][chr2] = {}
                    continue

                if '#' in line or len(line) == 0:
                    continue

                line = line.split()
                first = line[-3]
                second = line[-2]
                elements1.append(first)
                elements2.append(second)

                start1 = int(line[-3].split('ele')[1].split('to')[0])
                end1 = int(line[-3].split('ele')[1].split('to')[1])
                start2 = int(line[-2].split('ele')[1].split('to')[0])
                end2 = int(line[-2].split('ele')[1].split('to')[1])

                mini1 = min(mini1, start1)
                maxi1 = max(maxi1, end1)
                mini2 = min(mini2, start2)
                maxi2 = max(maxi2, end2)

        if mini1 != 'NA' and len(elements1) > 0 and len(elements2) > 0:
            if (elements1[0] > elements1[-1] and elements2[0] < elements2[-1]) or (elements1[0] < elements1[-1] and elements2[0] > elements2[-1]):
                sign = 'reverse'
            else:
                sign = 'forward'
            if (mini1, maxi1) not in chains[org1][chr1]:
                chains[org1][chr1][(mini1, maxi1)] = {}
            if org2 not in chains[org1][chr1][(mini1, maxi1)]:
                chains[org1][chr1][(mini1, maxi1)][org2] = (chr2, mini2, maxi2, sign)
            if (mini2, maxi2) not in chains[org2][chr2]:
                chains[org2][chr2][(mini2, maxi2)] = {}
            if org1 not in chains[org2][chr2][(mini2, maxi2)]:
                chains[org2][chr2][(mini2, maxi2)][org1] = (chr1, mini1, maxi1, sign)

        with open(chains_path, 'wb') as f:
            pickle.dump(chains, f)

    print(f"Loaded chains for {len(chains)} species")

    # Filter overlapping chains
    print("\nFiltering overlapping MCScanX chains...")
    chains = filter_overlapping_chains(chains)

    print("\n" + "="*80)
    print("STEP 2: PARSE GFF/FAA")
    print("="*80)

    parsed_path = f"{OUTPUT_DIR}/parsed_faa_gff"
    if os.path.isfile(parsed_path):
        print("Loading cached annotations...")
        with open(parsed_path, 'rb') as f:
            parsed = pickle.load(f)
    else:
        print("Parsing GFF/FAA files...")
        parsed = {}

        for gff_file in os.listdir(GFF_DIR):
            if not gff_file.endswith('.gff'):
                continue

            species = gff_file.replace('.gff', '')
            seq_file = None
            for ext in ['.fasta', '.fna']:
                test_path = os.path.join(SEQ_DIR, species + ext)
                if os.path.exists(test_path):
                    seq_file = test_path
                    break

            if not seq_file:
                raise FileNotFoundError(f"GFF file exists but no sequence file (.fasta/.fna) found for: {species}")

            print(f"  {species}...")

            gene_ids = set()
            for record in SeqIO.parse(seq_file, "fasta"):
                gene_ids.add(record.id.split()[0])

            gene_intervals = defaultdict(lambda: defaultdict(list))

            with open(os.path.join(GFF_DIR, gff_file)) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 9:
                        continue

                    seqid, source, ftype, start, end, score, strand, phase, attrs = cols

                    start, end = int(start), int(end)

                    attr_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in attrs.split(";") if "=" in kv}

                    candidates = set()
                    for key, value in attr_dict.items():
                        candidates.add(value)
                        if ":" in value:
                            candidates.add(value.split(":", 1)[1])

                    matched_id = None
                    for candidate in candidates:
                        prefixed_candidate = f"{species}_{candidate}"
                        if prefixed_candidate in gene_ids:
                            matched_id = prefixed_candidate
                            break

                    if not matched_id:
                        continue

                    gene_intervals[seqid][matched_id].append((start, end, strand))

            parsed[species] = {}
            for seqid in gene_intervals:
                for gene_id in gene_intervals[seqid]:
                    intervals = gene_intervals[seqid][gene_id]
                    # Always merge intervals with same ID (handles both single and multi-exon features)
                    min_start = min(s for s, e, strand in intervals)
                    max_end = max(e for s, e, strand in intervals)
                    strand = intervals[0][2]

                    if seqid not in parsed[species]:
                        parsed[species][seqid] = []
                    parsed[species][seqid].append((min_start, max_end, strand, gene_id))

            for seqid in parsed[species]:
                coords = parsed[species][seqid]
                coords.sort()

                kept = []
                for start, end, strand, pid in coords:
                    overlap = False
                    for i, (ks, ke, kstrand, kid) in enumerate(kept):
                        if not (end < ks or start > ke):
                            overlap = True
                            new_len = end - start
                            old_len = ke - ks
                            if new_len > old_len:
                                kept[i] = (start, end, strand, pid)
                            break

                    if not overlap:
                        kept.append((start, end, strand, pid))

                parsed[species][seqid] = kept

        with open(parsed_path, 'wb') as f:
            pickle.dump(parsed, f)

    print(f"Loaded annotations for {len(parsed)} species")

    print("\nValidating unique gene IDs...")
    validate_unique_gene_ids(parsed)
    print("  ✓ All gene IDs are unique within each species")

    print("\n" + "="*80)
    print("STEP 3: PARSE ALIGNMENTS")
    print("="*80)

    alignments_path = f"{OUTPUT_DIR}/parsed_alignments.pickle"
    if os.path.isfile(alignments_path):
        print("Loading cached alignments...")
        with open(alignments_path, 'rb') as f:
            alignments = pickle.load(f)
    else:
        print(f"Parsing {PAIRWISE_ALIGNMENTS_FILE}...")
        temp_alns = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        with open(PAIRWISE_ALIGNMENTS_FILE) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue

                try:
                    ref_species = fields[1]
                    ref_chr = fields[2]
                    ref_start = int(fields[3])
                    ref_end = int(fields[4])
                    tgt_species = fields[5]
                    tgt_chr = fields[6]
                    tgt_start = int(fields[7])
                    tgt_end = int(fields[8])
                    orientation = fields[9]
                    score = float(fields[10])
                except (ValueError, IndexError) as e:
                    print(f"WARNING: Malformed alignment table line, skipping: {line.strip()[:100]} - Error: {e}")
                    continue

                temp_alns[ref_species][ref_chr][tgt_species].append((ref_start, ref_end, tgt_chr, tgt_start, tgt_end, score, orientation))

        alignments = {species: {chr: dict(tgt_dict) for chr, tgt_dict in chr_dict.items()} for species, chr_dict in temp_alns.items()}

        with open(alignments_path, 'wb') as f:
            pickle.dump(alignments, f)

    # Count all unique species (both reference and target)
    all_species = set(alignments.keys())
    for ref_species in alignments.values():
        for chr_data in ref_species.values():
            all_species.update(chr_data.keys())

    print(f"Loaded alignments: {len(alignments)} reference species, {len(all_species)} total species involved")

    if BLAST_MODE == 'all-vs-all':
        print("\n" + "="*80)
        print("STEP 3.5: ALL-VS-ALL BLAST/CLASP")
        print("="*80)

        all_vs_all_path = f"{OUTPUT_DIR}/clasp_all_vs_all.pickle" if USE_CLASP else f"{OUTPUT_DIR}/blast_all_vs_all.pickle"

        if not NEW_BLAST and os.path.isfile(all_vs_all_path):
            print("Loading cached all-vs-all results...")
            with open(all_vs_all_path, 'rb') as f:
                all_vs_all_edges = pickle.load(f)
        else:
            print("Running all-vs-all BLAST/CLASP...")

            # Extract all tRNA sequences per species
            species_fastas = {}
            for species in parsed:
                seq_file = None
                for ext in ['.fasta', '.fna']:
                    test_path = os.path.join(SEQ_DIR, species + ext)
                    if os.path.exists(test_path):
                        seq_file = test_path
                        break

                if seq_file:
                    all_pids = set()
                    for chromo in parsed[species]:
                        for s, e, strand, pid in parsed[species][chromo]:
                            all_pids.add(pid)

                    seqs = []
                    for record in SeqIO.parse(seq_file, "fasta"):
                        if record.id.split()[0] in all_pids:
                            seqs.append(record)

                    fasta_path = f"{OUTPUT_DIR}/all_vs_all/{species}.fasta"
                    os.makedirs(f"{OUTPUT_DIR}/all_vs_all", exist_ok=True)
                    SeqIO.write(seqs, fasta_path, "fasta")
                    species_fastas[species] = fasta_path

            # Create species pairs for all-vs-all
            all_species_list = sorted(species_fastas.keys())
            species_pairs = []
            for i, sp1 in enumerate(all_species_list):
                for sp2 in all_species_list[i+1:]:
                    species_pairs.append((sp1, sp2))

            print(f"  Running BLAST for {len(species_pairs)} species pairs using {NUM_CORES} cores...")

            # Run BLAST for each pair using multiprocessing
            with mp.Pool(processes=NUM_CORES) as pool:
                results = pool.map(run_all_vs_all_blast,
                    [(sp1, sp2, species_fastas[sp1], species_fastas[sp2], OUTPUT_DIR, USE_CLASP, BLASTN_PATH, CLASP_PATH)
                     for sp1, sp2 in species_pairs])

            # Merge results into single dict
            all_vs_all_edges = {}
            for edges_dict in results:
                all_vs_all_edges.update(edges_dict)

            # Save
            with open(all_vs_all_path, 'wb') as f:
                pickle.dump(all_vs_all_edges, f)

            print(f"  Completed all-vs-all: {len(all_vs_all_edges)} edges")

        print(f"\nLoaded {len(all_vs_all_edges)} all-vs-all edges")

    print("\n" + "="*80)
    print("STEP 4: MERGE CHAINS WITH GENES")
    print("="*80)

    merged_path = f"{OUTPUT_DIR}/merged_chains_genes"
    if os.path.isfile(merged_path):
        print("Loading cached merged data...")
        with open(merged_path, 'rb') as f:
            merged = pickle.load(f)
    else:
        print("Merging...")
        merged = []
        seen_pairs = set()  # Track (org1, chr1, mini1, maxi1, org2, chr2, mini2, maxi2) in canonical order
        total_before_dedup = 0
        skipped_duplicates = 0

        for org1 in chains:
            if org1 not in parsed:
                print(f"WARNING: Species {org1} in chains but not in parsed annotations - skipping")
                continue
            for chr1 in chains[org1]:
                if chr1 not in parsed[org1]:
                    print(f"WARNING: Chromosome {chr1} for {org1} in chains but not in parsed annotations - skipping")
                    continue
                for coords1 in chains[org1][chr1]:
                    for org2 in chains[org1][chr1][coords1]:
                        if org2 not in parsed:
                            print(f"WARNING: Species {org2} in chains but not in parsed annotations - skipping")
                            continue
                        chr2, mini2, maxi2, sign = chains[org1][chr1][coords1][org2]
                        if chr2 not in parsed[org2]:
                            print(f"WARNING: Chromosome {chr2} for {org2} in chains but not in parsed annotations - skipping")
                            continue

                        mini1, maxi1 = coords1

                        genes1 = []
                        for start, end, strand, pid in parsed[org1][chr1]:
                            if not (end < mini1 or start > maxi1):
                                genes1.append((start, end, strand, pid))

                        genes2 = []
                        for start, end, strand, pid in parsed[org2][chr2]:
                            if not (end < mini2 or start > maxi2):
                                genes2.append((start, end, strand, pid))

                        if genes1 and genes2:
                            total_before_dedup += 1

                            # Create canonical pair representation (smaller tuple first)
                            pair_key = tuple(sorted([
                                (org1, chr1, mini1, maxi1),
                                (org2, chr2, mini2, maxi2)
                            ]))

                            # Skip if we've already seen this pair in opposite direction
                            if pair_key in seen_pairs:
                                skipped_duplicates += 1
                                continue
                            seen_pairs.add(pair_key)

                            merged.append({
                                'org1': org1,
                                'chr1': chr1,
                                'genes1': genes1,
                                'org2': org2,
                                'chr2': chr2,
                                'genes2': genes2,
                                'orientation': sign
                            })

        with open(merged_path, 'wb') as f:
            pickle.dump(merged, f)

        print(f"  Before deduplication: {total_before_dedup} chains")
        print(f"  Skipped duplicates: {skipped_duplicates} chains")
        print(f"  After deduplication: {len(merged)} chains")
        if total_before_dedup > 0:
            dedup_ratio = skipped_duplicates / total_before_dedup
            print(f"  Deduplication ratio: {dedup_ratio:.2%} (should be ~50%)")
            if abs(dedup_ratio - 0.5) > 0.05:
                print(f"  WARNING: Deduplication ratio is not close to 50% - deduplication may not be working correctly!")

    print(f"\nCreated {len(merged)} merged chains")

    # Check for elements appearing in multiple chains for same species pair
    print("\nChecking for overlapping chains per species pair...")

    # Map (canonical_species_pair) -> {element_id: [chain_indices]}
    species_pair_elements = defaultdict(lambda: defaultdict(list))

    for chain_idx, chain in enumerate(merged):
        # Create canonical species pair key
        species_pair = tuple(sorted([chain['org1'], chain['org2']]))

        # Track all elements in this chain
        for s, e, strand, pid in chain['genes1']:
            species_pair_elements[species_pair][pid].append(chain_idx)
        for s, e, strand, pid in chain['genes2']:
            species_pair_elements[species_pair][pid].append(chain_idx)

    # Find elements appearing in multiple chains for same species pair
    total_overlapping_pairs = 0
    for species_pair, element_chains in species_pair_elements.items():
        elements_in_multiple = {elem: chains for elem, chains in element_chains.items() if len(chains) > 1}
        if elements_in_multiple:
            total_overlapping_pairs += 1
            print(f"  Species pair {species_pair[0]} <-> {species_pair[1]}:")
            print(f"    {len(elements_in_multiple)} elements appear in multiple chains")
            # Show a few examples
            for i, (elem, chain_indices) in enumerate(list(elements_in_multiple.items())[:3]):
                print(f"      {elem}: chains {chain_indices}")
            if len(elements_in_multiple) > 3:
                print(f"      ... ({len(elements_in_multiple) - 3} more)")

    if total_overlapping_pairs == 0:
        print("  ✓ No overlapping chains found - each element appears in only one chain per species pair")
    else:
        print(f"\n  Found {total_overlapping_pairs} species pairs with overlapping chains")
        print(f"  This is expected when MCScanX finds multiple syntenic regions between same species")

    # Build gene_mapping for alignment node augmentation
    print("\nBuilding gene_mapping...")
    gene_mapping = {}
    for species in parsed:
        for chromo in parsed[species]:
            for start, end, strand, pid in parsed[species][chromo]:
                gene_mapping[pid] = (species, chromo, start, end, strand, {})
    print(f"  Mapped {len(gene_mapping)} genes")

    print("\n" + "="*80)
    print("STEP 5: RUN BLAST")
    print("="*80)

    blast_results_path = f"{OUTPUT_DIR}/blast_results.pickle"
    if not NEW_BLAST and os.path.isfile(blast_results_path):
        print("Loading cached BLAST results...")
        with open(blast_results_path, 'rb') as f:
            blast_results = pickle.load(f)
    else:
        print(f"Running BLAST for {len(merged)} chains using {NUM_CORES} cores...")
        with mp.Pool(processes=NUM_CORES) as pool:
            blast_results = pool.map(run_blast_for_chain, [(i, chain, OUTPUT_DIR, SEQ_DIR, NEW_BLAST, USE_CLASP, BLASTN_PATH, CLASP_PATH) for i, chain in enumerate(merged)])

        with open(blast_results_path, 'wb') as f:
            pickle.dump(blast_results, f)
        print(f"  Completed BLAST for {len(merged)} chains")

    print("\n" + "="*80)
    print("STEP 6: BUILD PAIRWISE GRAPHS + EDIT TO COGRAPH")
    print("="*80)

    graphs_cache_path = f"{OUTPUT_DIR}/graphs.pickle"
    if os.path.isfile(graphs_cache_path):
        print("Loading cached graphs...")
        with open(graphs_cache_path, 'rb') as f:
            graphs, total_edges_before_cograph, total_edges_after_cograph, total_nodes = pickle.load(f)
    else:
        print(f"Building graphs from {len(merged)} chains using {NUM_CORES} cores...")

        blast_dict = {i: hits for i, hits in blast_results}

        graph_args = [(i, chain, blast_dict[i], BLAST_THRESHOLD) for i, chain in enumerate(merged)]

        with mp.Pool(processes=NUM_CORES) as pool:
            graph_results = pool.map(build_graph_from_blast, graph_args)

        graphs = []
        total_edges_before_cograph = 0
        total_edges_after_cograph = 0
        total_nodes = 0

        for G, chain, edges_before, edges_after in graph_results:
            total_nodes += G.number_of_nodes()
            total_edges_before_cograph += edges_before
            total_edges_after_cograph += edges_after
            graphs.append((G, chain))

        with open(graphs_cache_path, 'wb') as f:
            pickle.dump((graphs, total_edges_before_cograph, total_edges_after_cograph, total_nodes), f)

    print(f"\n  Built {len(graphs)} pairwise graphs")
    print(f"  Total nodes: {total_nodes}")
    print(f"  Edges before cograph editing: {total_edges_before_cograph}")
    print(f"  Edges after cograph editing: {total_edges_after_cograph}")
    print(f"  Edges removed by cograph editing: {total_edges_before_cograph - total_edges_after_cograph}")

    print("\n" + "="*80)
    print("STEP 7: RUN ALIGNMENT ON GRAPHS")
    print("="*80)

    # Prepare worker args
    worker_args = [(idx, G, chain, alignments, gene_mapping)
                   for idx, (G, chain) in enumerate(graphs)]

    # Run alignment in parallel
    with mp.Pool(NUM_CORES) as pool:
        aligned_graphs = pool.map(align_graph_worker, worker_args)

    print(f"\n  Aligned {len(aligned_graphs)} graphs")

    # Replace graphs with aligned_graphs and collect alignment info
    graphs = [(G, chain) for G, chain, _ in aligned_graphs]

    # Organize alignment results by species pair
    alignment_results = {}
    for idx, (G, chain, alignment_info) in enumerate(aligned_graphs):
        species_pair = tuple(sorted([chain['org1'], chain['org2']]))
        if species_pair not in alignment_results:
            alignment_results[species_pair] = {}
        alignment_results[species_pair][idx] = alignment_info

    # Save aligned graphs for statistics
    print("\nSaving aligned graphs...")
    with open(f"{OUTPUT_DIR}/aligned_graphs.pickle", 'wb') as f:
        pickle.dump(aligned_graphs, f)

    print("\n" + "="*80)
    print("STEP 8: BUILD UNION GRAPH")
    print("="*80)

    union = nx.Graph()
    for G, chain in graphs:
        union = nx.compose(union, G)

    print(f"Union graph: {union.number_of_nodes()} nodes, {union.number_of_edges()} edges")

    print("\n" + "="*80)
    print("STEP 9: EDIT UNION TO COGRAPH")
    print("="*80)

    if union.number_of_nodes() == 0:
        print("Union graph is empty, skipping cograph editing")
        gene_species = {}
        union_v1 = union.copy()
        union_v2 = union.copy()
    else:
        gene_species = {}
        for species in parsed:
            for chr in parsed[species]:
                for s, e, strand, pid in parsed[species][chr]:
                    gene_species[pid] = species

        original_union_edges = set(union.edges())

        print("\nVersion 1: Cograph can add edges (only same-species forbidden)...")
        def is_forbidden_only_same_species(u, v):
            if u in gene_species and v in gene_species:
                if gene_species[u] == gene_species[v]:
                    return True
            return False

        union_v1 = union.copy()
        edges_before_v1 = union_v1.number_of_edges()
        union_v1 = edit_to_cograph(union_v1, forbidden_edge_checker=is_forbidden_only_same_species)
        edges_after_v1 = union_v1.number_of_edges()

        print(f"  Edges before: {edges_before_v1}")
        print(f"  Edges after: {edges_after_v1}")
        print(f"  Net change: {edges_after_v1 - edges_before_v1:+d}")

        print("\nVersion 2: Cograph cannot add edges (same-species + non-original forbidden)...")
        def is_forbidden_no_new_edges(u, v):
            if u in gene_species and v in gene_species:
                if gene_species[u] == gene_species[v]:
                    return True
            if (u, v) not in original_union_edges and (v, u) not in original_union_edges:
                return True
            return False

        union_v2 = union.copy()
        edges_before_v2 = union_v2.number_of_edges()
        union_v2 = edit_to_cograph(union_v2, forbidden_edge_checker=is_forbidden_no_new_edges)
        edges_after_v2 = union_v2.number_of_edges()

        print(f"  Edges before: {edges_before_v2}")
        print(f"  Edges after: {edges_after_v2}")
        print(f"  Removed: {edges_before_v2 - edges_after_v2}")

    print("\nSaving both final union graphs...")
    with open(f"{OUTPUT_DIR}/final_union_graph_with_new_edges.pickle", 'wb') as f:
        pickle.dump((union_v1, gene_species), f)
    with open(f"{OUTPUT_DIR}/final_union_graph_no_new_edges.pickle", 'wb') as f:
        pickle.dump((union_v2, gene_species), f)

    print("\nSaving alignment results...")
    alignments_dir = f"{OUTPUT_DIR}/alignments"
    os.makedirs(alignments_dir, exist_ok=True)

    # Save global alignment results
    with open(f"{alignments_dir}/global.pickle", 'wb') as f:
        pickle.dump(alignment_results, f)
    print(f"  Saved global alignment results: {len(alignment_results)} species pairs")

    # Save per-species-pair results
    for species_pair, pair_alignments in alignment_results.items():
        species_pair_str = f"{species_pair[0]}___{species_pair[1]}"

        # Save pickle
        with open(f"{alignments_dir}/{species_pair_str}.pickle", 'wb') as f:
            pickle.dump(pair_alignments, f)

        # Save human-readable text file
        with open(f"{alignments_dir}/{species_pair_str}.txt", 'w') as f:
            f.write("="*80 + "\n")
            f.write(f"ALIGNMENT RESULTS: {species_pair[0]} <-> {species_pair[1]}\n")
            f.write("="*80 + "\n\n")
            f.write(f"Total alignments: {len(pair_alignments)}\n\n")

            for chain_idx in sorted(pair_alignments.keys()):
                aln_info = pair_alignments[chain_idx]
                f.write("-"*80 + "\n")
                f.write(f"Chain {chain_idx}\n")
                f.write("-"*80 + "\n\n")

                f.write(f"{species_pair[0]} elements ({len(aln_info['species1_elements'])}):\n")
                for elem in aln_info['species1_elements']:
                    f.write(f"  {elem}\n")

                f.write(f"\n{species_pair[1]} elements ({len(aln_info['species2_elements'])}):\n")
                for elem in aln_info['species2_elements']:
                    f.write(f"  {elem}\n")

                f.write(f"\nEdges ({len(aln_info['edges'])}):\n")
                for u, v, attrs in aln_info['edges']:
                    event = attrs.get('event', 'unknown')
                    direction = attrs.get('direction', '')
                    match_type = attrs.get('match_type', '')
                    f.write(f"  {u} <-> {v} [event={event}")
                    if direction:
                        f.write(f", direction={direction}")
                    if match_type:
                        f.write(f", match_type={match_type}")
                    f.write("]\n")

                f.write("\n")

        print(f"  Saved {species_pair_str}: {len(pair_alignments)} alignments")

    print("\n" + "="*80)
    print("STEP 10: EXTRACT ORTHOLOGY GROUPS")
    print("="*80)

    union = union_v2
    print("="*80)

    components = list(nx.connected_components(union))
    multi_species = []
    for comp in components:
        species_in_comp = set()
        for node in comp:
            if node in gene_species:
                species_in_comp.add(gene_species[node])
        if len(species_in_comp) >= 2:
            multi_species.append(comp)

    print(f"Total components: {len(components)}")
    print(f"Multi-species components: {len(multi_species)}")

    print("\n" + "="*80)
    print("DONE")
    print("="*80)
