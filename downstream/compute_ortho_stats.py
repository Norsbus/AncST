#! /usr/bin/env python3

import pickle
import networkx as nx
import sys
import os
from collections import defaultdict

def calc_rbh(edges,gspe):
    bhits = {}
    for (ga, gb), sc in edges.items():
        if ga not in gspe or gb not in gspe:
            continue
        spa = gspe[ga]
        spb = gspe[gb]
        if spa == spb:
            continue
        if ga not in bhits:
            bhits[ga] = {}
        if spb not in bhits[ga] or sc > bhits[ga][spb][1]:
            bhits[ga][spb] = (gb, sc)
        if gb not in bhits:
            bhits[gb] = {}
        if spa not in bhits[gb] or sc > bhits[gb][spa][1]:
            bhits[gb][spa] = (ga, sc)
    rbh = set()
    for ga in bhits:
        spa = gspe[ga]
        for spb, (gb, _) in bhits[ga].items():
            if spa in bhits.get(gb, {}):
                if bhits[gb][spa][0] == ga:
                    rbh.add(tuple(sorted([ga, gb])))
    return rbh

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 compute_ortho_stats.py <output_dir>")
        sys.exit(1)
    outd = sys.argv[1]
    print("="*80)
    print("BUILDING MASTER ORTHOLOGY EDGE GRAPH")
    print("="*80)
    ava_clasp = os.path.exists(f'{outd}/clasp_all_vs_all.pickle')
    ava_blast = os.path.exists(f'{outd}/blast_all_vs_all.pickle')
    ava_mode = ava_clasp or ava_blast
    use_clasp = ava_clasp or os.path.exists(f'{outd}/clasp_out')
    mode = 'all-vs-all' if ava_mode else 'chains'
    score_type = 'CLASP' if use_clasp else 'BLAST'
    print(f"Mode: {mode}")
    print(f"Score type: {score_type}")
    union_file = f'{outd}/final_union_graph_no_new_edges.pickle'
    if not os.path.exists(union_file):
        print(f"ERROR: {union_file} not found")
        sys.exit(1)
    with open(union_file, 'rb') as f:
        _, gspe = pickle.load(f)
    print(f"Loaded gene_species: {len(gspe)} genes")
    thresh = 50
    mgraph = nx.Graph()
    if ava_mode:
        print("\n--- ALL-VS-ALL MODE ---")
        ava_file = f'{outd}/clasp_all_vs_all.pickle' if use_clasp else f'{outd}/blast_all_vs_all.pickle'
        print(f"\nLoading {ava_file}...")
        with open(ava_file, 'rb') as f:
            all_edges = pickle.load(f)
        for (u, v), sc in all_edges.items():
            if sc >= thresh:
                if u in gspe and v in gspe:
                    if gspe[u] != gspe[v]:
                        mgraph.add_edge(u, v, score=sc, is_rbh=False, in_chains=False, after_pairwise_cograph=False, after_alignment=False, after_union_with_new_edges=False, after_union_no_new_edges=False)
        print(f"  Added {mgraph.number_of_edges()} cross-species edges above threshold")
        print("\nCalculating RBH...")
        rbh_set = calc_rbh(all_edges, gspe)
        for edge in rbh_set:
            if mgraph.has_edge(*edge):
                mgraph.edges[edge]['is_rbh'] = True
        print(f"  Found {len(rbh_set)} RBH edges")
        print("\nMarking edges in chains...")
        with open(f'{outd}/merged_chains_genes', 'rb') as f:
            merged = pickle.load(f)
        g2c = defaultdict(set)
        for i, ch in enumerate(merged):
            for s, e, strand, pid in ch['genes1']:
                g2c[pid].add((i, 'org1'))
            for s, e, strand, pid in ch['genes2']:
                g2c[pid].add((i, 'org2'))
        for u, v in mgraph.edges():
            in_opp = False
            for ci_a, si_a in g2c.get(u, []):
                for ci_b, si_b in g2c.get(v, []):
                    if ci_a == ci_b and si_a != si_b:
                        in_opp = True
                        break
                if in_opp:
                    break
            if in_opp:
                mgraph.edges[u, v]['in_chains'] = True
        in_ch_cnt = sum(1 for u, v, d in mgraph.edges(data=True) if d['in_chains'])
        print(f"  {in_ch_cnt} edges have both genes in opposing chains")
    else:
        print("\n--- CHAINS-ONLY MODE ---")
        print("\nLoading BLAST/CLASP results...")
        with open(f'{outd}/blast_results.pickle', 'rb') as f:
            blast_results = pickle.load(f)
        ch_edges = {}
        for i, hits in blast_results:
            for qid, sid, sc in hits:
                if qid in gspe and sid in gspe:
                    if gspe[qid] != gspe[sid]:
                        edge = tuple(sorted([qid, sid]))
                        ch_edges[edge] = max(ch_edges.get(edge, 0), sc)
        for (u, v), sc in ch_edges.items():
            if sc >= thresh:
                mgraph.add_edge(u, v, score=sc, is_rbh=False, in_chains=True, after_pairwise_cograph=False, after_alignment=False, after_union_with_new_edges=False, after_union_no_new_edges=False)
        print(f"  Added {mgraph.number_of_edges()} cross-species edges above threshold")
        print("\nCalculating RBH...")
        rbh_set = calc_rbh(ch_edges, gspe)
        for edge in rbh_set:
            if mgraph.has_edge(*edge):
                mgraph.edges[edge]['is_rbh'] = True
        print(f"  Found {len(rbh_set)} RBH edges")
    print("\nMarking edges after pairwise cograph...")
    with open(f'{outd}/graphs.pickle', 'rb') as f:
        graphs, _, _, _ = pickle.load(f)
    for G, ch in graphs:
        for u, v in G.edges():
            if mgraph.has_edge(u, v):
                mgraph.edges[u, v]['after_pairwise_cograph'] = True
    after_pcg = sum(1 for u, v, d in mgraph.edges(data=True) if d['after_pairwise_cograph'])
    print(f"  {after_pcg} edges remain after pairwise cograph")
    print("\nMarking edges after alignment...")
    with open(f'{outd}/aligned_graphs.pickle', 'rb') as f:
        aligned_graphs = pickle.load(f)
    for G, ch, _ in aligned_graphs:
        for u, v in G.edges():
            if mgraph.has_edge(u, v):
                mgraph.edges[u, v]['after_alignment'] = True
    after_aln = sum(1 for u, v, d in mgraph.edges(data=True) if d['after_alignment'])
    print(f"  {after_aln} edges remain after alignment")
    print("\nMarking edges after union cograph...")
    with open(f'{outd}/final_union_graph_with_new_edges.pickle', 'rb') as f:
        union_v1, _ = pickle.load(f)
    for u, v in union_v1.edges():
        if mgraph.has_edge(u, v):
            mgraph.edges[u, v]['after_union_with_new_edges'] = True
    with open(f'{outd}/final_union_graph_no_new_edges.pickle', 'rb') as f:
        union_v2, _ = pickle.load(f)
    for u, v in union_v2.edges():
        if mgraph.has_edge(u, v):
            mgraph.edges[u, v]['after_union_no_new_edges'] = True
    after_u1 = sum(1 for u, v, d in mgraph.edges(data=True) if d['after_union_with_new_edges'])
    after_u2 = sum(1 for u, v, d in mgraph.edges(data=True) if d['after_union_no_new_edges'])
    print(f"  {after_u1} edges in union (with new edges)")
    print(f"  {after_u2} edges in union (no new edges)")
    outf = f'{outd}/master_orthology_edge_graph.pickle'
    print(f"\nSaving master graph to {outf}...")
    with open(outf, 'wb') as f:
        pickle.dump({'graph': mgraph, 'gene_species': gspe, 'mode': mode, 'score_type': score_type, 'threshold': thresh}, f)
    print("\n" + "="*80)
    print("STATISTICS")
    print("="*80)
    def cnt_edges(attr_filter=None, rbh_only=False):
        cnt = 0
        for u, v, d in mgraph.edges(data=True):
            if rbh_only and not d['is_rbh']:
                continue
            if attr_filter is None or d.get(attr_filter, False):
                cnt += 1
        return cnt
    if mode == 'all-vs-all':
        print("\n=== ALL-VS-ALL STRATIFICATION ===")
        print(f"\n{'Stage':<45} {'All':<12} {'RBH':<12}")
        print("-" * 70)
        tot = cnt_edges()
        tot_rbh = cnt_edges(rbh_only=True)
        print(f"{'Initial (above threshold)':<45} {tot:<12} {tot_rbh:<12}")
        in_ch = cnt_edges('in_chains')
        in_ch_rbh = cnt_edges('in_chains', rbh_only=True)
        rem = tot - in_ch
        rem_rbh = tot_rbh - in_ch_rbh
        print(f"  {'Removed (not in chains)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After chain filter':<45} {in_ch:<12} {in_ch_rbh:<12}")
        pcg = cnt_edges('after_pairwise_cograph')
        pcg_rbh = cnt_edges('after_pairwise_cograph', rbh_only=True)
        rem = in_ch - pcg
        rem_rbh = in_ch_rbh - pcg_rbh
        print(f"  {'Removed (pairwise cograph)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After pairwise cograph':<45} {pcg:<12} {pcg_rbh:<12}")
        aln = cnt_edges('after_alignment')
        aln_rbh = cnt_edges('after_alignment', rbh_only=True)
        rem = pcg - aln
        rem_rbh = pcg_rbh - aln_rbh
        print(f"  {'Removed (alignment)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After alignment':<45} {aln:<12} {aln_rbh:<12}")
        un2 = cnt_edges('after_union_no_new_edges')
        un2_rbh = cnt_edges('after_union_no_new_edges', rbh_only=True)
        rem = aln - un2
        rem_rbh = aln_rbh - un2_rbh
        print(f"  {'Removed (union cograph no new)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'Final (no new edges)':<45} {un2:<12} {un2_rbh:<12}")
        print("\n\n=== CHAINS-ONLY SUBSET (from all-vs-all data) ===")
        print(f"\n{'Stage':<45} {'All':<12} {'RBH':<12}")
        print("-" * 70)
        print(f"{'Initial (in chains, above threshold)':<45} {in_ch:<12} {in_ch_rbh:<12}")
        def cnt_ch_edges(attr_filter=None, rbh_only=False):
            cnt = 0
            for u, v, d in mgraph.edges(data=True):
                if not d.get('in_chains', False):
                    continue
                if rbh_only and not d['is_rbh']:
                    continue
                if attr_filter is None or d.get(attr_filter, False):
                    cnt += 1
            return cnt
        pcg_ch = cnt_ch_edges('after_pairwise_cograph')
        pcg_ch_rbh = cnt_ch_edges('after_pairwise_cograph', rbh_only=True)
        rem = in_ch - pcg_ch
        rem_rbh = in_ch_rbh - pcg_ch_rbh
        print(f"  {'Removed (pairwise cograph)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After pairwise cograph':<45} {pcg_ch:<12} {pcg_ch_rbh:<12}")
        aln_ch = cnt_ch_edges('after_alignment')
        aln_ch_rbh = cnt_ch_edges('after_alignment', rbh_only=True)
        rem = pcg_ch - aln_ch
        rem_rbh = pcg_ch_rbh - aln_ch_rbh
        print(f"  {'Removed (alignment)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After alignment':<45} {aln_ch:<12} {aln_ch_rbh:<12}")
        un2_ch = cnt_ch_edges('after_union_no_new_edges')
        un2_ch_rbh = cnt_ch_edges('after_union_no_new_edges', rbh_only=True)
        rem = aln_ch - un2_ch
        rem_rbh = aln_ch_rbh - un2_ch_rbh
        print(f"  {'Removed (union cograph no new)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'Final (no new edges)':<45} {un2_ch:<12} {un2_ch_rbh:<12}")
    else:
        print("\n=== CHAINS-ONLY MODE ===")
        print(f"\n{'Stage':<45} {'All':<12} {'RBH':<12}")
        print("-" * 70)
        tot = cnt_edges()
        tot_rbh = cnt_edges(rbh_only=True)
        print(f"{'Initial (in chains, above threshold)':<45} {tot:<12} {tot_rbh:<12}")
        pcg = cnt_edges('after_pairwise_cograph')
        pcg_rbh = cnt_edges('after_pairwise_cograph', rbh_only=True)
        rem = tot - pcg
        rem_rbh = tot_rbh - pcg_rbh
        print(f"  {'Removed (pairwise cograph)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After pairwise cograph':<45} {pcg:<12} {pcg_rbh:<12}")
        aln = cnt_edges('after_alignment')
        aln_rbh = cnt_edges('after_alignment', rbh_only=True)
        rem = pcg - aln
        rem_rbh = pcg_rbh - aln_rbh
        print(f"  {'Removed (alignment)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'After alignment':<45} {aln:<12} {aln_rbh:<12}")
        un2 = cnt_edges('after_union_no_new_edges')
        un2_rbh = cnt_edges('after_union_no_new_edges', rbh_only=True)
        rem = aln - un2
        rem_rbh = aln_rbh - un2_rbh
        print(f"  {'Removed (union cograph no new)':<43} {rem:<12} {rem_rbh:<12}")
        print(f"{'Final (no new edges)':<45} {un2:<12} {un2_rbh:<12}")
    print("\n" + "="*80)
    print("PER-SPECIES EDGE COUNTS")
    print("="*80)
    def cnt_per_sp(attr_filter=None, rbh_only=False, chains_only=False):
        cnts = defaultdict(int)
        for u, v, d in mgraph.edges(data=True):
            if rbh_only and not d.get('is_rbh', False):
                continue
            if chains_only and not d.get('in_chains', False):
                continue
            if attr_filter is not None and not d.get(attr_filter, False):
                continue
            if u in gspe:
                cnts[gspe[u]] += 1
            if v in gspe:
                cnts[gspe[v]] += 1
        return cnts
    all_sp = sorted(set(gspe.values()))
    if mode == 'chains':
        print("\n=== CHAINS-ONLY MODE ===\n")
        init = cnt_per_sp()
        apcg = cnt_per_sp('after_pairwise_cograph')
        aaln = cnt_per_sp('after_alignment')
        fin = cnt_per_sp('after_union_no_new_edges')
        init_rbh = cnt_per_sp(rbh_only=True)
        apcg_rbh = cnt_per_sp('after_pairwise_cograph', rbh_only=True)
        aaln_rbh = cnt_per_sp('after_alignment', rbh_only=True)
        fin_rbh = cnt_per_sp('after_union_no_new_edges', rbh_only=True)
        print(f"{'Species':<20} {'Initial':<10} {'After PCG':<10} {'After Aln':<10} {'Final':<10}")
        print("-" * 70)
        for sp in all_sp:
            print(f"{sp:<20} {init.get(sp, 0):<10} {apcg.get(sp, 0):<10} {aaln.get(sp, 0):<10} {fin.get(sp, 0):<10}")
        print("\n=== CHAINS-ONLY MODE (RBH only) ===\n")
        print(f"{'Species':<20} {'Initial':<10} {'After PCG':<10} {'After Aln':<10} {'Final':<10}")
        print("-" * 70)
        for sp in all_sp:
            print(f"{sp:<20} {init_rbh.get(sp, 0):<10} {apcg_rbh.get(sp, 0):<10} {aaln_rbh.get(sp, 0):<10} {fin_rbh.get(sp, 0):<10}")
        per_sp_stats = {'species_list': all_sp, 'initial': init, 'after_pairwise_cograph': apcg, 'after_alignment': aaln, 'after_union_no_new_edges': fin, 'rbh': {'initial': init_rbh, 'after_pairwise_cograph': apcg_rbh, 'after_alignment': aaln_rbh, 'after_union_no_new_edges': fin_rbh}}
        with open(outf, 'rb') as f:
            mdata = pickle.load(f)
        mdata['per_species_stats'] = per_sp_stats
        with open(outf, 'wb') as f:
            pickle.dump(mdata, f)
    else:
        print("\n=== ALL-VS-ALL MODE ===\n")
        init = cnt_per_sp()
        in_ch_sp = cnt_per_sp('in_chains')
        apcg = cnt_per_sp('after_pairwise_cograph')
        aaln = cnt_per_sp('after_alignment')
        fin = cnt_per_sp('after_union_no_new_edges')
        init_rbh = cnt_per_sp(rbh_only=True)
        in_ch_sp_rbh = cnt_per_sp('in_chains', rbh_only=True)
        apcg_rbh = cnt_per_sp('after_pairwise_cograph', rbh_only=True)
        aaln_rbh = cnt_per_sp('after_alignment', rbh_only=True)
        fin_rbh = cnt_per_sp('after_union_no_new_edges', rbh_only=True)
        print(f"{'Species':<20} {'Initial':<10} {'In Chains':<10} {'After PCG':<10} {'After Aln':<10} {'Final':<10}")
        print("-" * 90)
        for sp in all_sp:
            print(f"{sp:<20} {init.get(sp, 0):<10} {in_ch_sp.get(sp, 0):<10} {apcg.get(sp, 0):<10} {aaln.get(sp, 0):<10} {fin.get(sp, 0):<10}")
        print("\n=== ALL-VS-ALL MODE (RBH only) ===\n")
        print(f"{'Species':<20} {'Initial':<10} {'In Chains':<10} {'After PCG':<10} {'After Aln':<10} {'Final':<10}")
        print("-" * 90)
        for sp in all_sp:
            print(f"{sp:<20} {init_rbh.get(sp, 0):<10} {in_ch_sp_rbh.get(sp, 0):<10} {apcg_rbh.get(sp, 0):<10} {aaln_rbh.get(sp, 0):<10} {fin_rbh.get(sp, 0):<10}")
        print("\n\n=== CHAINS-ONLY SUBSET (from all-vs-all data) ===\n")
        in_ch_ch = cnt_per_sp('in_chains', chains_only=True)
        apcg_ch = cnt_per_sp('after_pairwise_cograph', chains_only=True)
        aaln_ch = cnt_per_sp('after_alignment', chains_only=True)
        fin_ch = cnt_per_sp('after_union_no_new_edges', chains_only=True)
        in_ch_ch_rbh = cnt_per_sp('in_chains', rbh_only=True, chains_only=True)
        apcg_ch_rbh = cnt_per_sp('after_pairwise_cograph', rbh_only=True, chains_only=True)
        aaln_ch_rbh = cnt_per_sp('after_alignment', rbh_only=True, chains_only=True)
        fin_ch_rbh = cnt_per_sp('after_union_no_new_edges', rbh_only=True, chains_only=True)
        print(f"{'Species':<20} {'Initial':<10} {'After PCG':<10} {'After Aln':<10} {'Final':<10}")
        print("-" * 70)
        for sp in all_sp:
            print(f"{sp:<20} {in_ch_ch.get(sp, 0):<10} {apcg_ch.get(sp, 0):<10} {aaln_ch.get(sp, 0):<10} {fin_ch.get(sp, 0):<10}")
        print("\n=== CHAINS-ONLY SUBSET (RBH only) ===\n")
        print(f"{'Species':<20} {'Initial':<10} {'After PCG':<10} {'After Aln':<10} {'Final':<10}")
        print("-" * 70)
        for sp in all_sp:
            print(f"{sp:<20} {in_ch_ch_rbh.get(sp, 0):<10} {apcg_ch_rbh.get(sp, 0):<10} {aaln_ch_rbh.get(sp, 0):<10} {fin_ch_rbh.get(sp, 0):<10}")
        per_sp_stats = {'species_list': all_sp, 'initial': init, 'in_chains': in_ch_sp, 'after_pairwise_cograph': apcg, 'after_alignment': aaln, 'after_union_no_new_edges': fin, 'rbh': {'initial': init_rbh, 'in_chains': in_ch_sp_rbh, 'after_pairwise_cograph': apcg_rbh, 'after_alignment': aaln_rbh, 'after_union_no_new_edges': fin_rbh}, 'chains_only': {'initial': in_ch_ch, 'after_pairwise_cograph': apcg_ch, 'after_alignment': aaln_ch, 'after_union_no_new_edges': fin_ch, 'rbh': {'initial': in_ch_ch_rbh, 'after_pairwise_cograph': apcg_ch_rbh, 'after_alignment': aaln_ch_rbh, 'after_union_no_new_edges': fin_ch_rbh}}}
        with open(outf, 'rb') as f:
            mdata = pickle.load(f)
        mdata['per_species_stats'] = per_sp_stats
        with open(outf, 'wb') as f:
            pickle.dump(mdata, f)
    print("\n" + "="*80)
    print("DONE")
    print("="*80)

if __name__ == '__main__':
    main()
