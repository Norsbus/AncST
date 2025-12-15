#! /usr/bin/env python3

import pickle
import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 make_ortho_tables.py <output_dir>")
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
    tdir = f'{outd}/orthology_tables'
    os.makedirs(tdir, exist_ok=True)
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
        fname = f'{tdir}/all_vs_all_pipeline.tex'
        with open(fname, 'w') as f:
            f.write("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=1in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\begin{tabular}{lrrrr}\n\\toprule\n")
            f.write(f"Stage & \\multicolumn{{2}}{{c}}{{All edges (>{thresh})}} & \\multicolumn{{2}}{{c}}{{RBH only}} \\\\\n\\cmidrule(lr){{2-3}} \\cmidrule(lr){{4-5}}\n & Count & Removed & Count & Removed \\\\\n\\midrule\n")
            sdata = []
            sdata.append(('Initial (above threshold)', cnt_edges(), cnt_edges(rbh_only=True)))
            sdata.append(('After chain filter', cnt_edges('in_chains'), cnt_edges('in_chains', rbh_only=True)))
            sdata.append(('After pairwise cograph', cnt_edges('after_pairwise_cograph'), cnt_edges('after_pairwise_cograph', rbh_only=True)))
            sdata.append(('After alignment', cnt_edges('after_alignment'), cnt_edges('after_alignment', rbh_only=True)))
            sdata.append(('Final (no new edges)', cnt_edges('after_union_no_new_edges'), cnt_edges('after_union_no_new_edges', rbh_only=True)))
            for i, (stage, acnt, rcnt) in enumerate(sdata):
                if i == 0:
                    f.write(f"{stage} & {acnt} & -- & {rcnt} & -- \\\\\n")
                else:
                    pall = sdata[i-1][1]
                    prbh = sdata[i-1][2]
                    rall = pall - acnt
                    rrbh = prbh - rcnt
                    f.write(f"{stage} & {acnt} & {rall} & {rcnt} & {rrbh} \\\\\n")
            f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{All-vs-all {score_type} edge filtering through pipeline}}\n\\end{{table}}\n\\end{{document}}\n")
        print(f"  Saved {fname}")
        fname = f'{tdir}/chains_only_pipeline.tex'
        with open(fname, 'w') as f:
            f.write("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=1in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\begin{tabular}{lrrrr}\n\\toprule\n")
            f.write(f"Stage & \\multicolumn{{2}}{{c}}{{All edges (>{thresh})}} & \\multicolumn{{2}}{{c}}{{RBH only}} \\\\\n\\cmidrule(lr){{2-3}} \\cmidrule(lr){{4-5}}\n & Count & Removed & Count & Removed \\\\\n\\midrule\n")
            sdata = []
            sdata.append(('Initial (in chains)', cnt_ch_edges(), cnt_ch_edges(rbh_only=True)))
            sdata.append(('After pairwise cograph', cnt_ch_edges('after_pairwise_cograph'), cnt_ch_edges('after_pairwise_cograph', rbh_only=True)))
            sdata.append(('After alignment', cnt_ch_edges('after_alignment'), cnt_ch_edges('after_alignment', rbh_only=True)))
            sdata.append(('Final (no new edges)', cnt_ch_edges('after_union_no_new_edges'), cnt_ch_edges('after_union_no_new_edges', rbh_only=True)))
            for i, (stage, acnt, rcnt) in enumerate(sdata):
                if i == 0:
                    f.write(f"{stage} & {acnt} & -- & {rcnt} & -- \\\\\n")
                else:
                    pall = sdata[i-1][1]
                    prbh = sdata[i-1][2]
                    rall = pall - acnt
                    rrbh = prbh - rcnt
                    f.write(f"{stage} & {acnt} & {rall} & {rcnt} & {rrbh} \\\\\n")
            f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{Chains-only {score_type} edge filtering through pipeline}}\n\\end{{table}}\n\\end{{document}}\n")
        print(f"  Saved {fname}")
    else:
        fname = f'{tdir}/chains_only_pipeline.tex'
        with open(fname, 'w') as f:
            f.write("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=1in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\begin{tabular}{lrrrr}\n\\toprule\n")
            f.write(f"Stage & \\multicolumn{{2}}{{c}}{{All edges (>{thresh})}} & \\multicolumn{{2}}{{c}}{{RBH only}} \\\\\n\\cmidrule(lr){{2-3}} \\cmidrule(lr){{4-5}}\n & Count & Removed & Count & Removed \\\\\n\\midrule\n")
            sdata = []
            sdata.append(('Initial (in chains)', cnt_edges(), cnt_edges(rbh_only=True)))
            sdata.append(('After pairwise cograph', cnt_edges('after_pairwise_cograph'), cnt_edges('after_pairwise_cograph', rbh_only=True)))
            sdata.append(('After alignment', cnt_edges('after_alignment'), cnt_edges('after_alignment', rbh_only=True)))
            sdata.append(('Final (no new edges)', cnt_edges('after_union_no_new_edges'), cnt_edges('after_union_no_new_edges', rbh_only=True)))
            for i, (stage, acnt, rcnt) in enumerate(sdata):
                if i == 0:
                    f.write(f"{stage} & {acnt} & -- & {rcnt} & -- \\\\\n")
                else:
                    pall = sdata[i-1][1]
                    prbh = sdata[i-1][2]
                    rall = pall - acnt
                    rrbh = prbh - rcnt
                    f.write(f"{stage} & {acnt} & {rall} & {rcnt} & {rrbh} \\\\\n")
            f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{{score_type} edge filtering through pipeline (chains-only mode)}}\n\\end{{table}}\n\\end{{document}}\n")
        print(f"  Saved {fname}")
    if 'per_species_stats' in data:
        print("\nGenerating per-species tables...")
        pss = data['per_species_stats']
        splist = pss['species_list']
        if mode == 'chains':
            fname = f'{tdir}/per_species_pipeline.tex'
            with open(fname, 'w') as f:
                f.write("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=1in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\small\n\\begin{tabular}{lrrrrrrrr}\n\\toprule\n")
                f.write("Species & \\multicolumn{2}{c}{Initial} & \\multicolumn{2}{c}{After PCG} & \\multicolumn{2}{c}{After Aln} & \\multicolumn{2}{c}{Final} \\\\\n\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\n & Count & Removed & Count & Removed & Count & Removed & Count & Removed \\\\\n\\midrule\n")
                skeys = ['initial', 'after_pairwise_cograph', 'after_alignment', 'after_union_no_new_edges']
                for sp in splist:
                    cnts = [pss[k].get(sp, 0) for k in skeys]
                    sname = sp.replace('GCF_', '')
                    f.write(f"{sname} & {cnts[0]} & --")
                    for i in range(1, len(cnts)):
                        rem = cnts[i-1] - cnts[i]
                        f.write(f" & {cnts[i]} & {rem}")
                    f.write(" \\\\\n")
                f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{Per-species edge counts through pipeline (chains-only mode, {score_type} >{thresh})}}\n\\end{{table}}\n\\end{{document}}\n")
            print(f"  Saved {fname}")
            if 'rbh' in pss:
                fname = f'{tdir}/per_species_pipeline_rbh.tex'
                with open(fname, 'w') as f:
                    f.write("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=1in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\small\n\\begin{tabular}{lrrrrrrrr}\n\\toprule\n")
                    f.write("Species & \\multicolumn{2}{c}{Initial} & \\multicolumn{2}{c}{After PCG} & \\multicolumn{2}{c}{After Aln} & \\multicolumn{2}{c}{Final} \\\\\n\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\n & Count & Removed & Count & Removed & Count & Removed & Count & Removed \\\\\n\\midrule\n")
                    for sp in splist:
                        cnts = [pss['rbh'][k].get(sp, 0) for k in skeys]
                        sname = sp.replace('GCF_', '')
                        f.write(f"{sname} & {cnts[0]} & --")
                        for i in range(1, len(cnts)):
                            rem = cnts[i-1] - cnts[i]
                            f.write(f" & {cnts[i]} & {rem}")
                        f.write(" \\\\\n")
                    f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{Per-species RBH edge counts through pipeline (chains-only mode, {score_type} >{thresh})}}\n\\end{{table}}\n\\end{{document}}\n")
                print(f"  Saved {fname}")
        else:
            fname = f'{tdir}/per_species_all_vs_all.tex'
            with open(fname, 'w') as f:
                f.write("\\documentclass[11pt,a4paper,landscape]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=0.5in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\tiny\n\\begin{tabular}{lrrrrrrrrrr}\n\\toprule\n")
                f.write("Species & \\multicolumn{2}{c}{Initial} & \\multicolumn{2}{c}{In Chains} & \\multicolumn{2}{c}{After PCG} & \\multicolumn{2}{c}{After Aln} & \\multicolumn{2}{c}{Final} \\\\\n\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11}\n & Count & Removed & Count & Removed & Count & Removed & Count & Removed & Count & Removed \\\\\n\\midrule\n")
                skeys = ['initial', 'in_chains', 'after_pairwise_cograph', 'after_alignment', 'after_union_no_new_edges']
                for sp in splist:
                    cnts = [pss[k].get(sp, 0) for k in skeys]
                    sname = sp.replace('GCF_', '')
                    f.write(f"{sname} & {cnts[0]} & --")
                    for i in range(1, len(cnts)):
                        rem = cnts[i-1] - cnts[i]
                        f.write(f" & {cnts[i]} & {rem}")
                    f.write(" \\\\\n")
                f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{Per-species edge counts (all-vs-all mode, {score_type} >{thresh})}}\n\\end{{table}}\n\\end{{document}}\n")
            print(f"  Saved {fname}")
            if 'chains_only' in pss:
                fname = f'{tdir}/per_species_chains_only.tex'
                with open(fname, 'w') as f:
                    f.write("\\documentclass[11pt,a4paper]{article}\n\\usepackage[utf8]{inputenc}\n\\usepackage{booktabs}\n\\usepackage[margin=1in]{geometry}\n\\begin{document}\n\\begin{table}[h]\n\\centering\n\\small\n\\begin{tabular}{lrrrrrrrr}\n\\toprule\n")
                    f.write("Species & \\multicolumn{2}{c}{Initial} & \\multicolumn{2}{c}{After PCG} & \\multicolumn{2}{c}{After Aln} & \\multicolumn{2}{c}{Final} \\\\\n\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}\n & Count & Removed & Count & Removed & Count & Removed & Count & Removed \\\\\n\\midrule\n")
                    skeys = ['initial', 'after_pairwise_cograph', 'after_alignment', 'after_union_no_new_edges']
                    for sp in splist:
                        cnts = [pss['chains_only'][k].get(sp, 0) for k in skeys]
                        sname = sp.replace('GCF_', '')
                        f.write(f"{sname} & {cnts[0]} & --")
                        for i in range(1, len(cnts)):
                            rem = cnts[i-1] - cnts[i]
                            f.write(f" & {cnts[i]} & {rem}")
                        f.write(" \\\\\n")
                    f.write(f"\\bottomrule\n\\end{{tabular}}\n\\caption{{Per-species edge counts (chains-only subset, {score_type} >{thresh})}}\n\\end{{table}}\n\\end{{document}}\n")
                print(f"  Saved {fname}")
    print(f"\nTables saved to {tdir}/")

if __name__ == '__main__':
    main()
