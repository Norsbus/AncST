#! /usr/bin/env python3

import pickle
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_ortho_stats.py <output_dir>")
        sys.exit(1)
    outd = sys.argv[1]
    mfile = f'{outd}/master_orthology_edge_graph.pickle'
    if not os.path.exists(mfile):
        print(f"ERROR: {mfile} not found")
        print("Run compute_ortho_stats.py first")
        sys.exit(1)
    print("Loading master orthology edge graph...")
    with open(mfile, 'rb') as f:
        data = pickle.load(f)
    mgraph = data['graph']
    mode = data['mode']
    score_type = data['score_type']
    thresh = data['threshold']
    print(f"Mode: {mode}, Score: {score_type}, Threshold: {thresh}")
    pdir = f'{outd}/orthology_plots'
    os.makedirs(pdir, exist_ok=True)
    def cnt_edges(attr_filter=None, rbh_only=False):
        cnt = 0
        for u, v, d in mgraph.edges(data=True):
            if rbh_only and not d['is_rbh']:
                continue
            if attr_filter is None or d.get(attr_filter, False):
                cnt += 1
        return cnt
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
    if mode == 'all-vs-all':
        stg = ['Initial', 'In chains', 'After\npairwise\ncograph', 'After\nalignment', 'Final']
        acnt = [cnt_edges(), cnt_edges('in_chains'), cnt_edges('after_pairwise_cograph'), cnt_edges('after_alignment'), cnt_edges('after_union_no_new_edges')]
        rcnt = [cnt_edges(rbh_only=True), cnt_edges('in_chains', rbh_only=True), cnt_edges('after_pairwise_cograph', rbh_only=True), cnt_edges('after_alignment', rbh_only=True), cnt_edges('after_union_no_new_edges', rbh_only=True)]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.bar(stg, acnt, alpha=0.7, edgecolor='black')
        ax1.set_ylabel('Number of edges')
        ax1.set_title(f'All-vs-all: All edges (>{thresh})')
        ax1.grid(axis='y', alpha=0.3)
        for i, v in enumerate(acnt):
            ax1.text(i, v, str(v), ha='center', va='bottom')
        ax2.bar(stg, rcnt, alpha=0.7, edgecolor='black', color='orange')
        ax2.set_ylabel('Number of edges')
        ax2.set_title('All-vs-all: RBH only')
        ax2.grid(axis='y', alpha=0.3)
        for i, v in enumerate(rcnt):
            ax2.text(i, v, str(v), ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(f'{pdir}/all_vs_all_pipeline.svg', format='svg', bbox_inches='tight')
        plt.savefig(f'{pdir}/all_vs_all_pipeline.pdf', format='pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved all-vs-all plot")
        chstg = ['Initial\n(in chains)', 'After\npairwise\ncograph', 'After\nalignment', 'Final']
        chacnt = [cnt_ch_edges(), cnt_ch_edges('after_pairwise_cograph'), cnt_ch_edges('after_alignment'), cnt_ch_edges('after_union_no_new_edges')]
        chrcnt = [cnt_ch_edges(rbh_only=True), cnt_ch_edges('after_pairwise_cograph', rbh_only=True), cnt_ch_edges('after_alignment', rbh_only=True), cnt_ch_edges('after_union_no_new_edges', rbh_only=True)]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.bar(chstg, chacnt, alpha=0.7, edgecolor='black')
        ax1.set_ylabel('Number of edges')
        ax1.set_title(f'Chains-only: All edges (>{thresh})')
        ax1.grid(axis='y', alpha=0.3)
        for i, v in enumerate(chacnt):
            ax1.text(i, v, str(v), ha='center', va='bottom')
        ax2.bar(chstg, chrcnt, alpha=0.7, edgecolor='black', color='orange')
        ax2.set_ylabel('Number of edges')
        ax2.set_title('Chains-only: RBH only')
        ax2.grid(axis='y', alpha=0.3)
        for i, v in enumerate(chrcnt):
            ax2.text(i, v, str(v), ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(f'{pdir}/chains_only_pipeline.svg', format='svg', bbox_inches='tight')
        plt.savefig(f'{pdir}/chains_only_pipeline.pdf', format='pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved chains-only plot")
    else:
        stg = ['Initial', 'After\npairwise\ncograph', 'After\nalignment', 'Final']
        acnt = [cnt_edges(), cnt_edges('after_pairwise_cograph'), cnt_edges('after_alignment'), cnt_edges('after_union_no_new_edges')]
        rcnt = [cnt_edges(rbh_only=True), cnt_edges('after_pairwise_cograph', rbh_only=True), cnt_edges('after_alignment', rbh_only=True), cnt_edges('after_union_no_new_edges', rbh_only=True)]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.bar(stg, acnt, alpha=0.7, edgecolor='black')
        ax1.set_ylabel('Number of edges')
        ax1.set_title(f'Chains-only: All edges (>{thresh})')
        ax1.grid(axis='y', alpha=0.3)
        for i, v in enumerate(acnt):
            ax1.text(i, v, str(v), ha='center', va='bottom')
        ax2.bar(stg, rcnt, alpha=0.7, edgecolor='black', color='orange')
        ax2.set_ylabel('Number of edges')
        ax2.set_title('Chains-only: RBH only')
        ax2.grid(axis='y', alpha=0.3)
        for i, v in enumerate(rcnt):
            ax2.text(i, v, str(v), ha='center', va='bottom')
        plt.tight_layout()
        plt.savefig(f'{pdir}/chains_only_pipeline.svg', format='svg', bbox_inches='tight')
        plt.savefig(f'{pdir}/chains_only_pipeline.pdf', format='pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved chains-only plot")
    if 'per_species_stats' in data:
        print("\nGenerating per-species plots...")
        pss = data['per_species_stats']
        splist = pss['species_list']
        if mode == 'chains':
            stg = ['Initial', 'After\nPCG', 'After\nAlignment', 'Final']
            skeys = ['initial', 'after_pairwise_cograph', 'after_alignment', 'after_union_no_new_edges']
            x = np.arange(len(stg))
            width = 0.08
            nspe = len(splist)
            fig, ax = plt.subplots(figsize=(14, 8))
            colors = plt.cm.tab10(np.linspace(0, 1, nspe))
            for i, sp in enumerate(splist):
                cnts = [pss[k].get(sp, 0) for k in skeys]
                offset = (i - nspe / 2) * width
                ax.bar(x + offset, cnts, width, label=sp, color=colors[i], alpha=0.8, edgecolor='black', linewidth=0.5)
            ax.set_ylabel('Number of edges', fontsize=12)
            ax.set_title('Edge counts per species through pipeline stages', fontsize=14)
            ax.set_xticks(x)
            ax.set_xticklabels(stg)
            ax.legend(loc='upper right', fontsize=9, ncol=2)
            ax.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(f'{pdir}/per_species_pipeline.svg', format='svg', bbox_inches='tight')
            plt.savefig(f'{pdir}/per_species_pipeline.pdf', format='pdf', bbox_inches='tight')
            plt.close()
            print(f"  Saved per-species plot")
            fig, ax = plt.subplots(figsize=(10, 6))
            for i, sp in enumerate(splist):
                cnts = [pss[k].get(sp, 0) for k in skeys]
                ax.plot(stg, cnts, marker='o', label=sp, linewidth=2, markersize=6)
            ax.set_ylabel('Number of edges', fontsize=12)
            ax.set_xlabel('Pipeline stage', fontsize=12)
            ax.set_title('Edge counts per species through pipeline stages', fontsize=14)
            ax.legend(loc='upper right', fontsize=9, ncol=2)
            ax.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(f'{pdir}/per_species_pipeline_lines.svg', format='svg', bbox_inches='tight')
            plt.savefig(f'{pdir}/per_species_pipeline_lines.pdf', format='pdf', bbox_inches='tight')
            plt.close()
            print(f"  Saved per-species line plot")
            if 'rbh' in pss:
                fig, axes = plt.subplots(2, 5, figsize=(20, 8))
                axes = axes.flatten()
                for i, sp in enumerate(splist):
                    acnt = [pss[k].get(sp, 0) for k in skeys]
                    rcnt = [pss['rbh'][k].get(sp, 0) for k in skeys]
                    xpos = np.arange(len(stg))
                    w = 0.35
                    axes[i].bar(xpos - w/2, acnt, w, label='All edges', alpha=0.8, edgecolor='black')
                    axes[i].bar(xpos + w/2, rcnt, w, label='RBH only', alpha=0.8, edgecolor='black', color='orange')
                    axes[i].set_title(sp.replace('GCF_', ''), fontsize=10)
                    axes[i].set_xticks(xpos)
                    axes[i].set_xticklabels(stg, fontsize=8, rotation=45, ha='right')
                    axes[i].legend(fontsize=8)
                    axes[i].grid(axis='y', alpha=0.3)
                plt.suptitle('Per-species edge counts: All vs RBH', fontsize=14)
                plt.tight_layout()
                plt.savefig(f'{pdir}/per_species_all_vs_rbh.svg', format='svg', bbox_inches='tight')
                plt.savefig(f'{pdir}/per_species_all_vs_rbh.pdf', format='pdf', bbox_inches='tight')
                plt.close()
                print(f"  Saved per-species All vs RBH comparison plot")
        else:
            stg = ['Initial', 'In\nChains', 'After\nPCG', 'After\nAlignment', 'Final']
            skeys = ['initial', 'in_chains', 'after_pairwise_cograph', 'after_alignment', 'after_union_no_new_edges']
            x = np.arange(len(stg))
            width = 0.08
            nspe = len(splist)
            fig, ax = plt.subplots(figsize=(16, 8))
            colors = plt.cm.tab10(np.linspace(0, 1, nspe))
            for i, sp in enumerate(splist):
                cnts = [pss[k].get(sp, 0) for k in skeys]
                offset = (i - nspe / 2) * width
                ax.bar(x + offset, cnts, width, label=sp, color=colors[i], alpha=0.8, edgecolor='black', linewidth=0.5)
            ax.set_ylabel('Number of edges', fontsize=12)
            ax.set_title('Edge counts per species (all-vs-all mode)', fontsize=14)
            ax.set_xticks(x)
            ax.set_xticklabels(stg)
            ax.legend(loc='upper right', fontsize=9, ncol=2)
            ax.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(f'{pdir}/per_species_all_vs_all.svg', format='svg', bbox_inches='tight')
            plt.savefig(f'{pdir}/per_species_all_vs_all.pdf', format='pdf', bbox_inches='tight')
            plt.close()
            print(f"  Saved per-species all-vs-all plot")
            if 'chains_only' in pss:
                stg = ['Initial\n(in chains)', 'After\nPCG', 'After\nAlignment', 'Final']
                skeys = ['initial', 'after_pairwise_cograph', 'after_alignment', 'after_union_no_new_edges']
                x = np.arange(len(stg))
                fig, ax = plt.subplots(figsize=(14, 8))
                for i, sp in enumerate(splist):
                    cnts = [pss['chains_only'][k].get(sp, 0) for k in skeys]
                    offset = (i - nspe / 2) * width
                    ax.bar(x + offset, cnts, width, label=sp, color=colors[i], alpha=0.8, edgecolor='black', linewidth=0.5)
                ax.set_ylabel('Number of edges', fontsize=12)
                ax.set_title('Edge counts per species (chains-only subset)', fontsize=14)
                ax.set_xticks(x)
                ax.set_xticklabels(stg)
                ax.legend(loc='upper right', fontsize=9, ncol=2)
                ax.grid(axis='y', alpha=0.3)
                plt.tight_layout()
                plt.savefig(f'{pdir}/per_species_chains_only.svg', format='svg', bbox_inches='tight')
                plt.savefig(f'{pdir}/per_species_chains_only.pdf', format='pdf', bbox_inches='tight')
                plt.close()
                print(f"  Saved per-species chains-only plot")
    print(f"\nPlots saved to {pdir}/")

if __name__ == '__main__':
    main()
