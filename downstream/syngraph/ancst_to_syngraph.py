#!/usr/bin/env python3
"""
Convert AncST / MCScanX data into syngraph-compatible per-species TSV files.

Two approaches:
  mcscanx  – MCScanX collinearity blocks (whole blocks or individual anchors)
  anchors  – native AncST anchor pickle files

Clustering pipeline: Union-Find -> optional MCL refinement -> species-count filter.
Use --no-trans to require fully connected clusters (cliques) instead of
transitively connected components.

Filtering:
  --include-species  restrict to specific species (file or comma-separated)
  --include-chrs     restrict to specific chromosomes (file or comma-separated)

Filters are applied before clustering so only relevant data influences marker
grouping.

Usage examples:
  # From AncST anchors
  ./ancst_to_syngraph.py --approach anchors \\
      --anchor-dir out/anchors --mapping out/utils/mapping \\
      --min-score 200 --min-species 13 --no-trans \\
      --output-dir syngraph_input

  # From MCScanX collinearity, filtered to sex chromosomes
  ./ancst_to_syngraph.py --approach mcscanx \\
      --collinearity MCScanX.collinearity \\
      --mcscanx-mode blocks --block-merge-strategy shared_anchors \\
      --include-chrs sex_chrs.txt --min-species 7 \\
      --output-dir syngraph_input_sex_chrs

  # Filter file format (one entry per line):
  #   include-species: 1org  (or accessions for anchors approach)
  #   include-chrs:    1orgchr12
"""

import argparse
import os
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


# ---------------------------------------------------------------------------
# Element name parsing (for MCScanX approach)
# ---------------------------------------------------------------------------

ELEMENT_RE = re.compile(r'^(\d+org)(chr\d+)ele(\d+)to(\d+)$')


def parse_element(name):
    """Parse '1orgchr1ele26113050to26114676' -> dict with org, chr, start, end."""
    m = ELEMENT_RE.match(name)
    if not m:
        raise ValueError(f"Cannot parse element name: {name}")
    return {
        'org': m.group(1),
        'chr': m.group(2),
        'start': int(m.group(3)),
        'end': int(m.group(4)),
    }


def _extract_org_from_orgchr(orgchr):
    """'1orgchr1' -> '1org', '10orgchr3' -> '10org'"""
    m = re.match(r'^(\d+org)', orgchr)
    return m.group(1) if m else orgchr


# ---------------------------------------------------------------------------
# .collinearity parser
# ---------------------------------------------------------------------------

ALIGNMENT_RE = re.compile(
    r'^## Alignment (\d+): score=([\d.]+) e_value=([\d.eE+-]+) N=(\d+) (\S+)&(\S+) (plus|minus)'
)


def parse_collinearity(filepath):
    """
    Parse an MCScanX .collinearity file.

    Returns a list of block dicts:
        {
            'block_id': int, 'score': float, 'e_value': float, 'n_pairs': int,
            'org1_chr': str, 'org2_chr': str, 'orientation': str,
            'pairs': [(element1_str, element2_str), ...]
        }
    """
    blocks = []
    current_block = None

    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')

            m = ALIGNMENT_RE.match(line)
            if not m and line.startswith('#'):
                continue
            if m:
                if current_block is not None:
                    blocks.append(current_block)
                current_block = {
                    'block_id': int(m.group(1)),
                    'score': float(m.group(2)),
                    'e_value': float(m.group(3)),
                    'n_pairs': int(m.group(4)),
                    'org1_chr': m.group(5),
                    'org2_chr': m.group(6),
                    'orientation': m.group(7),
                    'pairs': [],
                }
                continue

            if current_block is not None and ':' in line:
                parts = line.split(':', 1)[1].split()
                if len(parts) >= 2:
                    current_block['pairs'].append((parts[0], parts[1]))

    if current_block is not None:
        blocks.append(current_block)

    return blocks


# ---------------------------------------------------------------------------
# Mapping file parser
# ---------------------------------------------------------------------------

def parse_mapping(filepath):
    """
    Parse the AncST utils/mapping file.

    Returns:
        org_mapping:  dict  accession <-> alias  (bidirectional)
        chr_mapping:  dict  {species_key: {accession <-> alias}}  (bidirectional)
    """
    org_mapping = {}
    chr_mapping = {}

    with open(filepath) as f:
        reading_species = False
        current_accession = None

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                if 'species' in line:
                    reading_species = True
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            if reading_species:
                accession, alias = parts[0], parts[1]
                org_mapping[accession] = alias
                org_mapping[alias] = accession
                current_accession = accession
                chr_mapping[accession] = {}
                chr_mapping[alias] = {}
                reading_species = False
            else:
                if current_accession is None:
                    continue
                acc, alias = parts[0], parts[1]
                alias_org = org_mapping[current_accession]
                chr_mapping[current_accession][acc] = alias
                chr_mapping[current_accession][alias] = acc
                chr_mapping[alias_org][acc] = alias
                chr_mapping[alias_org][alias] = acc

    return org_mapping, chr_mapping


# ---------------------------------------------------------------------------
# Filter helpers
# ---------------------------------------------------------------------------

def parse_filter_arg(value):
    """
    Parse a filter argument that can be either a file path (one entry per line)
    or a comma-separated list.

    Returns a set of strings.
    """
    if value is None:
        return None
    if os.path.isfile(value):
        entries = set()
        with open(value) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    entries.add(line)
        return entries
    return set(v.strip() for v in value.split(',') if v.strip())


def filter_blocks_by_species(blocks, include_species):
    """Keep only blocks where both species are in include_species."""
    filtered = []
    for b in blocks:
        org1 = _extract_org_from_orgchr(b['org1_chr'])
        org2 = _extract_org_from_orgchr(b['org2_chr'])
        if org1 in include_species and org2 in include_species:
            filtered.append(b)
    return filtered


def filter_blocks_by_chrs(blocks, include_chrs):
    """Keep only blocks where both chromosomes are in include_chrs."""
    filtered = []
    for b in blocks:
        if b['org1_chr'] in include_chrs and b['org2_chr'] in include_chrs:
            filtered.append(b)
    return filtered


# ---------------------------------------------------------------------------
# Overlap filtering (adapted from AncST synthology pipeline)
# ---------------------------------------------------------------------------

def _compute_block_spans(blocks):
    """Compute genomic spans per (org, chr) for each block."""
    block_spans = {}
    for b in blocks:
        spans = {}
        for elem1, elem2 in b['pairs']:
            for elem in (elem1, elem2):
                p = parse_element(elem)
                key = (p['org'], p['chr'])
                if key in spans:
                    cur_start, cur_end = spans[key]
                    spans[key] = (min(cur_start, p['start']), max(cur_end, p['end']))
                else:
                    spans[key] = (p['start'], p['end'])
        block_spans[b['block_id']] = spans
    return block_spans


def filter_overlapping_blocks(blocks):
    """
    Filter overlapping collinearity blocks.

    Each block is pairwise (2 species). For two blocks sharing the same species pair:
    - Both-sided overlap (overlap on both species' chromosomes): keep the longer block
    - One-sided overlap (overlap on one species but not the other): remove BOTH
    """
    block_spans = _compute_block_spans(blocks)

    block_species_pair = {}
    block_regions = {}
    for b in blocks:
        org1 = _extract_org_from_orgchr(b['org1_chr'])
        org2 = _extract_org_from_orgchr(b['org2_chr'])
        sp_pair = tuple(sorted([org1, org2]))
        block_species_pair[b['block_id']] = sp_pair
        regions = {}
        for (org, chrom), (start, end) in block_spans[b['block_id']].items():
            regions[org] = (chrom, start, end)
        block_regions[b['block_id']] = regions

    species_pair_blocks = defaultdict(list)
    for b in blocks:
        species_pair_blocks[block_species_pair[b['block_id']]].append(b)

    blocks_to_remove = set()
    n_both_sided = 0
    n_one_sided = 0

    for sp_pair, sp_blocks in species_pair_blocks.items():
        sp1, sp2 = sp_pair
        for i in range(len(sp_blocks)):
            for j in range(i + 1, len(sp_blocks)):
                bid_a = sp_blocks[i]['block_id']
                bid_b = sp_blocks[j]['block_id']
                if bid_a in blocks_to_remove or bid_b in blocks_to_remove:
                    continue

                reg_a = block_regions[bid_a]
                reg_b = block_regions[bid_b]

                def overlaps(org):
                    if org not in reg_a or org not in reg_b:
                        return False
                    chr_a, s_a, e_a = reg_a[org]
                    chr_b, s_b, e_b = reg_b[org]
                    return chr_a == chr_b and s_a <= e_b and s_b <= e_a

                ov_sp1 = overlaps(sp1)
                ov_sp2 = overlaps(sp2)

                if ov_sp1 and ov_sp2:
                    len_a = sum(e - s for (s, e) in block_spans[bid_a].values())
                    len_b = sum(e - s for (s, e) in block_spans[bid_b].values())
                    if len_a >= len_b:
                        blocks_to_remove.add(bid_b)
                    else:
                        blocks_to_remove.add(bid_a)
                    n_both_sided += 1
                elif ov_sp1 or ov_sp2:
                    blocks_to_remove.add(bid_a)
                    blocks_to_remove.add(bid_b)
                    n_one_sided += 1

    filtered = [b for b in blocks if b['block_id'] not in blocks_to_remove]
    print(f"  Overlap filtering: {len(blocks)} blocks -> {len(filtered)} blocks")
    print(f"    Both-sided overlaps (kept longest): {n_both_sided}")
    print(f"    One-sided overlaps (removed both): {n_one_sided}")
    if n_one_sided > 0:
        print(f"    WARNING: {n_one_sided} one-sided overlaps detected")
    return filtered


# ---------------------------------------------------------------------------
# Union-Find
# ---------------------------------------------------------------------------

class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1

    def components(self):
        groups = defaultdict(set)
        for x in self.parent:
            groups[self.find(x)].add(x)
        return list(groups.values())


# ---------------------------------------------------------------------------
# Clique finding (for --no-trans mode)
# ---------------------------------------------------------------------------

def _find_maximal_cliques(adj, nodes):
    """Bron-Kerbosch with pivoting. Returns all maximal cliques (size >= 2)."""
    cliques = []
    nodes = set(nodes)

    def bron_kerbosch(R, P, X):
        if not P and not X:
            if len(R) >= 2:
                cliques.append(set(R))
            return
        pivot = max(P | X, key=lambda v: len(adj[v] & P))
        for v in list(P - adj[pivot]):
            neighbors_v = adj[v]
            bron_kerbosch(R | {v}, P & neighbors_v, X & neighbors_v)
            P.remove(v)
            X.add(v)

    bron_kerbosch(set(), nodes, set())
    return cliques


# ---------------------------------------------------------------------------
# MCL clustering
# ---------------------------------------------------------------------------

def run_mcl_external(edges, inflation):
    """Write ABC format, call mcl binary, parse output."""
    tmpdir = tempfile.mkdtemp(prefix='syngraph_mcl_')
    abc_path = os.path.join(tmpdir, 'graph.abc')
    mci_path = os.path.join(tmpdir, 'graph.mci')
    tab_path = os.path.join(tmpdir, 'graph.tab')
    out_path = os.path.join(tmpdir, 'graph.mcl')

    with open(abc_path, 'w') as f:
        for n1, n2, w in edges:
            f.write(f"{n1}\t{n2}\t{w}\n")

    try:
        r1 = subprocess.run(
            ['mcxload', '-abc', abc_path, '--stream-mirror',
             '-o', mci_path, '-write-tab', tab_path],
            capture_output=True, text=True,
        )
        if r1.returncode != 0:
            raise RuntimeError(
                f"mcxload failed (exit {r1.returncode}):\n{r1.stderr.strip()}")
        r2 = subprocess.run(
            ['mcl', mci_path, '-I', str(inflation),
             '-o', out_path, '-use-tab', tab_path],
            capture_output=True, text=True,
        )
        if r2.returncode != 0:
            raise RuntimeError(
                f"mcl failed (exit {r2.returncode}):\n{r2.stderr.strip()}")
        clusters = []
        with open(out_path) as f:
            for line in f:
                nodes = line.strip().split('\t')
                if nodes and nodes[0]:
                    clusters.append(set(nodes))
        return clusters
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def run_mcl_python(edges, inflation):
    """Fallback: use the markov-clustering Python package."""
    import markov_clustering as mcl_lib
    import scipy.sparse as sp

    nodes = set()
    for n1, n2, _ in edges:
        nodes.add(n1)
        nodes.add(n2)
    node_list = sorted(nodes, key=str)
    node_idx = {n: i for i, n in enumerate(node_list)}
    n = len(node_list)

    rows, cols, data = [], [], []
    for n1, n2, w in edges:
        i, j = node_idx[n1], node_idx[n2]
        rows.extend([i, j])
        cols.extend([j, i])
        data.extend([w, w])

    matrix = sp.csc_matrix((data, (rows, cols)), shape=(n, n))
    result = mcl_lib.run_mcl(matrix, inflation=inflation)
    mcl_clusters = mcl_lib.get_clusters(result)

    return [set(node_list[i] for i in cluster) for cluster in mcl_clusters]


def run_mcl(edges, inflation=2.0, mcl_impl='auto'):
    """
    Run MCL clustering.
    mcl_impl: 'auto' (external first, Python fallback), 'external', or 'python'
    """
    if not edges:
        return []

    if mcl_impl == 'external':
        if not (shutil.which('mcl') and shutil.which('mcxload')):
            raise RuntimeError("mcl/mcxload binaries not found in PATH")
        return run_mcl_external(edges, inflation)

    if mcl_impl == 'python':
        try:
            import markov_clustering  # noqa: F401
            return run_mcl_python(edges, inflation)
        except ImportError:
            raise RuntimeError(
                "markov-clustering package not installed: "
                "pip install markov-clustering scipy")

    # auto: try external first
    if shutil.which('mcl') and shutil.which('mcxload'):
        return run_mcl_external(edges, inflation)
    try:
        import markov_clustering  # noqa: F401
        return run_mcl_python(edges, inflation)
    except ImportError:
        raise RuntimeError(
            "Neither mcl binary nor markov-clustering package found.\n"
            "Install mcl:  conda install -c bioconda mcl  OR  brew install mcl\n"
            "Or install:   pip install markov-clustering scipy"
        )


# ---------------------------------------------------------------------------
# Shared clustering pipeline: UF -> MCL/cliques -> filter
# ---------------------------------------------------------------------------

def cluster_pipeline(edges, inflation=2.0, use_mcl=True, min_species=2,
                     node_to_org=None, mcl_impl='auto', no_trans=False):
    """
    edges:       [(node1, node2, weight), ...]
    node_to_org: callable  node -> org string  (for species counting)
    no_trans:    if True, find maximal cliques instead of transitive components

    Returns list of sets (each set = one cluster of node IDs).
    """
    if not edges:
        return []

    uf = UnionFind()
    for n1, n2, _ in edges:
        uf.union(n1, n2)
    components = uf.components()

    if no_trans:
        adj = defaultdict(set)
        for n1, n2, _ in edges:
            adj[n1].add(n2)
            adj[n2].add(n1)
        refined = []
        for comp in components:
            if len(comp) <= 2:
                refined.append(comp)
                continue
            cliques = _find_maximal_cliques(adj, comp)
            if cliques:
                refined.extend(cliques)
    elif use_mcl:
        comp_edges_map = defaultdict(list)
        for n1, n2, w in edges:
            comp_edges_map[uf.find(n1)].append((n1, n2, w))
        refined = []
        for comp in components:
            if len(comp) <= 2:
                refined.append(comp)
                continue
            root = uf.find(next(iter(comp)))
            comp_edges = comp_edges_map.get(root, [])
            if not comp_edges:
                refined.append(comp)
                continue
            mcl_clusters = run_mcl(comp_edges, inflation, mcl_impl)
            if mcl_clusters:
                refined.extend(mcl_clusters)
            else:
                refined.append(comp)
    else:
        refined = components

    if node_to_org is not None and min_species > 1:
        filtered = []
        for cluster in refined:
            species = {node_to_org(n) for n in cluster}
            if len(species) >= min_species:
                filtered.append(cluster)
        return filtered

    return refined


# ---------------------------------------------------------------------------
# Cluster -> marker conversion (shared by both approaches)
# ---------------------------------------------------------------------------

def clusters_to_markers(clusters, node_info_fn, include_chrs=None):
    """
    Convert clusters to marker dicts.

    node_info_fn: callable(node) -> {'org': str, 'chr': str, 'start': int, 'end': int}
                  or None if node has no info
    include_chrs: optional set of chromosomes (XorgchrY format) to keep.
                  Applied after clustering as a final filter.
    """
    markers = []
    for idx, cluster in enumerate(clusters, 1):
        marker_id = f"marker_{idx:06d}"
        elements_by_org = defaultdict(list)
        for node in cluster:
            info = node_info_fn(node)
            if info is not None:
                elements_by_org[info['org']].append(info)

        per_species = {}
        for org, elems in elements_by_org.items():
            chr_counts = defaultdict(int)
            for e in elems:
                chr_counts[e['chr']] += 1
            best_chr = max(chr_counts, key=chr_counts.get)
            full_chr = f"{org}{best_chr}"

            if include_chrs is not None and full_chr not in include_chrs:
                continue

            chr_elems = [e for e in elems if e['chr'] == best_chr]
            per_species[org] = {
                'chr': best_chr,
                'start': min(e['start'] for e in chr_elems),
                'end': max(e['end'] for e in chr_elems),
            }

        if not per_species:
            continue

        markers.append({
            'marker_id': marker_id,
            'elements': per_species,
            'n_elements': len(cluster),
            'n_species': len(per_species),
        })

    return markers


# ---------------------------------------------------------------------------
# Approach: MCScanX (individual anchors or whole blocks)
# ---------------------------------------------------------------------------

def approach_mcscanx_anchors(blocks, min_species=2, inflation=2.0, use_mcl=True,
                             mcl_impl='auto', no_trans=False, include_chrs=None):
    """Cluster individual MCScanX anchor elements from collinearity blocks."""
    edges = []
    for block in blocks:
        w = block['score']
        for elem1, elem2 in block['pairs']:
            edges.append((elem1, elem2, w))

    def node_to_org(elem):
        return parse_element(elem)['org']

    clusters = cluster_pipeline(edges, inflation, use_mcl, min_species,
                                node_to_org, mcl_impl, no_trans)

    def node_info_fn(node):
        try:
            p = parse_element(node)
            return p
        except ValueError:
            return None

    return clusters_to_markers(clusters, node_info_fn, include_chrs)


def approach_mcscanx_blocks(blocks, strategy='shared_anchors', min_species=2,
                            inflation=2.0, use_mcl=True, mcl_impl='auto',
                            no_trans=False, include_chrs=None):
    """Cluster whole MCScanX collinearity blocks, then produce markers."""
    block_by_id = {}
    for b in blocks:
        block_by_id[b['block_id']] = b
        block_by_id[str(b['block_id'])] = b

    if strategy == 'shared_anchors':
        edges = _block_edges_shared_anchors(blocks)
    elif strategy == 'overlapping_regions':
        edges = _block_edges_overlapping(blocks)
    else:
        raise ValueError(f"Unknown block merge strategy: {strategy}")

    print(f"  Built {len(edges)} block edges")

    clusters = cluster_pipeline(edges, inflation, use_mcl, min_species=0,
                                node_to_org=None, mcl_impl=mcl_impl,
                                no_trans=no_trans)
    print(f"  Clustering produced {len(clusters)} clusters")

    # Add singleton blocks
    connected_bids = set()
    for n1, n2, _ in edges:
        connected_bids.add(n1)
        connected_bids.add(n2)
    for b in blocks:
        if b['block_id'] not in connected_bids:
            clusters.append({b['block_id']})

    # Convert block clusters to markers
    markers = []
    marker_num = 0
    for cluster in clusters:
        elements_by_org = defaultdict(list)
        for bid in cluster:
            if bid not in block_by_id:
                continue
            b = block_by_id[bid]
            for elem1, elem2 in b['pairs']:
                for elem in (elem1, elem2):
                    p = parse_element(elem)
                    elements_by_org[p['org']].append(p)

        per_species = {}
        for org, elems in elements_by_org.items():
            chr_counts = defaultdict(int)
            for e in elems:
                chr_counts[e['chr']] += 1
            best_chr = max(chr_counts, key=chr_counts.get)
            full_chr = f"{org}{best_chr}"
            if include_chrs is not None and full_chr not in include_chrs:
                continue
            chr_elems = [e for e in elems if e['chr'] == best_chr]
            per_species[org] = {
                'chr': best_chr,
                'start': min(e['start'] for e in chr_elems),
                'end': max(e['end'] for e in chr_elems),
            }

        if len(per_species) < min_species:
            continue

        marker_num += 1
        markers.append({
            'marker_id': f"marker_{marker_num:06d}",
            'elements': per_species,
            'n_elements': sum(len(v) for v in elements_by_org.values()),
            'n_species': len(per_species),
        })

    return markers


def _block_edges_shared_anchors(blocks):
    """Connect MCScanX blocks that share at least one anchor element."""
    element_to_blocks = defaultdict(set)
    for b in blocks:
        for elem1, elem2 in b['pairs']:
            element_to_blocks[elem1].add(b['block_id'])
            element_to_blocks[elem2].add(b['block_id'])

    edge_set = set()
    for _elem, bids in element_to_blocks.items():
        bid_list = sorted(bids)
        for i in range(len(bid_list)):
            for j in range(i + 1, len(bid_list)):
                pair = (bid_list[i], bid_list[j])
                if pair not in edge_set:
                    edge_set.add(pair)

    return [(b1, b2, 1.0) for b1, b2 in edge_set]


def _block_edges_overlapping(blocks):
    """Connect MCScanX blocks whose genomic spans overlap on the same (org, chr)."""
    block_spans = {}
    for b in blocks:
        spans = defaultdict(lambda: (float('inf'), 0))
        for elem1, elem2 in b['pairs']:
            for elem in (elem1, elem2):
                p = parse_element(elem)
                key = (p['org'], p['chr'])
                cur = spans[key]
                spans[key] = (min(cur[0], p['start']), max(cur[1], p['end']))
        block_spans[b['block_id']] = dict(spans)

    orgchr_to_blocks = defaultdict(list)
    for bid, spans in block_spans.items():
        for key, (start, end) in spans.items():
            orgchr_to_blocks[key].append((bid, start, end))

    edge_set = set()
    for _key, blist in orgchr_to_blocks.items():
        blist.sort(key=lambda x: x[1])
        for i in range(len(blist)):
            bid_i, _start_i, end_i = blist[i]
            for j in range(i + 1, len(blist)):
                bid_j, start_j, _end_j = blist[j]
                if start_j > end_i:
                    break
                pair = (min(bid_i, bid_j), max(bid_i, bid_j))
                edge_set.add(pair)

    return [(b1, b2, 1.0) for b1, b2 in edge_set]


# ---------------------------------------------------------------------------
# Approach: AncST anchors (pickle files)
# ---------------------------------------------------------------------------

def discover_anchor_orgs(anchor_dir):
    """
    Auto-discover organism anchor files under anchor_dir.

    Supports layouts:
      1. anchor_dir/aligned_succinct/{accession}   (AncST default)
      2. anchor_dir/{accession}/{aligned_succinct}  (per-org subdirs)
      3. anchor_dir/{accession}                     (flat files)
    """
    orgs = []
    paths = {}

    aligned_dir = os.path.join(anchor_dir, 'aligned_succinct')
    if os.path.isdir(aligned_dir):
        for entry in sorted(os.listdir(aligned_dir)):
            full = os.path.join(aligned_dir, entry)
            if os.path.isfile(full):
                orgs.append(entry)
                paths[entry] = full
        if orgs:
            return orgs, paths

    for entry in sorted(os.listdir(anchor_dir)):
        candidate = os.path.join(anchor_dir, entry, 'aligned_succinct')
        if os.path.isfile(candidate):
            orgs.append(entry)
            paths[entry] = candidate

    if orgs:
        return orgs, paths

    for entry in sorted(os.listdir(anchor_dir)):
        full = os.path.join(anchor_dir, entry)
        if os.path.isfile(full) and not entry.startswith('.'):
            orgs.append(entry)
            paths[entry] = full

    return orgs, paths


def approach_anchors(anchor_dir, org_mapping, chr_mapping, min_species=2,
                     inflation=2.0, use_mcl=True, min_score=0, mcl_impl='auto',
                     no_trans=False, include_species=None, include_chrs=None):
    """
    Load AncST aligned_succinct anchor files, build graph, cluster.

    Anchor files use NCBI accession names for organisms and chromosomes.
    org_mapping and chr_mapping (from parse_mapping) convert these to
    Xorg/XorgchrY format for syngraph output.
    """
    orgs, pkl_paths = discover_anchor_orgs(anchor_dir)

    if not orgs:
        print("Error: no anchor files found.", file=sys.stderr)
        return []

    print(f"  Discovered {len(orgs)} organisms")

    # Build accession -> alias maps
    acc_to_alias = {}   # e.g. GCA_019343175.1 -> 1org
    chr_to_alias = {}   # e.g. (1org, CM033270.1) -> chr1
    for acc, alias in org_mapping.items():
        if re.match(r'^\d+org$', alias):
            acc_to_alias[acc] = alias
    for species_key, chrmap in chr_mapping.items():
        alias = org_mapping.get(species_key, species_key)
        if not re.match(r'^\d+org$', alias):
            alias = org_mapping.get(alias, alias)
        for chr_acc, chr_alias in chrmap.items():
            m = re.match(r'^(\d+org)?(chr\d+)$', chr_alias)
            if m:
                # Strip org prefix: mapping has '9orgchr1213' but we
                # need 'chr1213' since org is prepended later in output
                chr_to_alias[(alias, chr_acc)] = m.group(2)

    # Load anchor files, applying species filter
    am = {}
    for org_acc in orgs:
        org_alias = acc_to_alias.get(org_acc)
        if org_alias is None:
            print(f"  Warning: no mapping for {org_acc}, skipping", file=sys.stderr)
            continue
        if include_species is not None and org_alias not in include_species:
            continue
        with open(pkl_paths[org_acc], 'rb') as f:
            am[org_acc] = pickle.load(f)
        print(f"  Loaded {org_alias} ({org_acc}): {len(am[org_acc])} anchors")

    if not am:
        print("Error: no anchor files loaded.", file=sys.stderr)
        return []

    def convert_chr(org_acc, chr_acc):
        org_alias = acc_to_alias.get(org_acc, org_acc)
        return chr_to_alias.get((org_alias, chr_acc), chr_acc)

    # Pre-populate node_info, applying chromosome filter
    node_info = {}
    unconverted_chrs = set()
    for org_acc, anchors in am.items():
        org_alias = acc_to_alias[org_acc]
        for i, data in anchors.items():
            chr_alias = convert_chr(org_acc, data['chromosome'])
            if chr_alias == data['chromosome']:
                unconverted_chrs.add((org_alias, data['chromosome']))
                continue
            full_chr = f"{org_alias}{chr_alias}"
            if include_chrs is not None and full_chr not in include_chrs:
                continue
            node_id = f"{org_alias}_anchor_{i}"
            node_info[node_id] = {
                'org': org_alias,
                'chr': chr_alias,
                'start': data['start'],
                'end': data['end'],
            }
    if unconverted_chrs:
        print(f"Error: {len(unconverted_chrs)} chromosomes not found in mapping:",
              file=sys.stderr)
        for org, chrom in sorted(unconverted_chrs):
            print(f"  {org}: {chrom}", file=sys.stderr)
        sys.exit(1)

    # Build edges
    edge_best = {}
    for org_acc, anchors in am.items():
        org_alias = acc_to_alias[org_acc]
        for i, data in anchors.items():
            node1 = f"{org_alias}_anchor_{i}"
            if node1 not in node_info:
                continue
            for org2_acc, matches in data.get('matches', {}).items():
                if org2_acc not in am:
                    continue
                org2_alias = acc_to_alias.get(org2_acc)
                if org2_alias is None:
                    continue
                for j, match_data in matches.items():
                    score = int(match_data[0])
                    if score < min_score:
                        continue
                    node2 = f"{org2_alias}_anchor_{j}"
                    # Populate node_info for anchors only referenced via matches
                    if node2 not in node_info:
                        if j in am[org2_acc]:
                            chr2_alias = convert_chr(
                                org2_acc, am[org2_acc][j]['chromosome'])
                            full_chr2 = f"{org2_alias}{chr2_alias}"
                            if include_chrs is not None and full_chr2 not in include_chrs:
                                continue
                            node_info[node2] = {
                                'org': org2_alias,
                                'chr': chr2_alias,
                                'start': am[org2_acc][j]['start'],
                                'end': am[org2_acc][j]['end'],
                            }
                        else:
                            continue
                    key = (node1, node2) if node1 <= node2 else (node2, node1)
                    if key not in edge_best or score > edge_best[key]:
                        edge_best[key] = score

    edges = [(k[0], k[1], v) for k, v in edge_best.items()]
    del edge_best

    print(f"  Built {len(edges)} edges from {len(node_info)} anchor nodes")

    def node_to_org(node):
        return node_info[node]['org'] if node in node_info else node.split('_')[0]

    clusters = cluster_pipeline(edges, inflation, use_mcl, min_species,
                                node_to_org, mcl_impl, no_trans)

    def node_info_fn(node):
        return node_info.get(node)

    return clusters_to_markers(clusters, node_info_fn)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_syngraph_tsvs(markers, output_dir, include_species=None):
    """Write per-species syngraph TSV files."""
    os.makedirs(output_dir, exist_ok=True)

    per_species = defaultdict(list)
    for m in markers:
        for org, coords in m['elements'].items():
            if include_species is not None and org not in include_species:
                continue
            per_species[org].append(
                (m['marker_id'], f"{org}{coords['chr']}",
                 coords['start'], coords['end'])
            )

    for org, rows in per_species.items():
        rows.sort(key=lambda r: (r[1], r[2]))
        filepath = os.path.join(output_dir, f"{org}.tsv")
        with open(filepath, 'w') as f:
            for marker_id, chrom, start, end in rows:
                f.write(f"{marker_id}\t{chrom}\t{start}\t{end}\n")

    return set(per_species.keys())


def write_cluster_details(markers, output_dir):
    """Write marker_clusters.txt with full cluster details."""
    filepath = os.path.join(output_dir, 'marker_clusters.txt')
    with open(filepath, 'w') as f:
        f.write("marker_id\tn_species\tn_elements\tspecies\tcoordinates\n")
        for m in markers:
            species_list = sorted(m['elements'].keys())
            coords_list = []
            for org in species_list:
                c = m['elements'][org]
                coords_list.append(f"{org}{c['chr']}:{c['start']}-{c['end']}")
            f.write(f"{m['marker_id']}\t{m['n_species']}\t{m['n_elements']}\t"
                    f"{','.join(species_list)}\t{','.join(coords_list)}\n")


def write_cluster_stats(markers, output_dir, species_written):
    """Write cluster_stats.txt summary."""
    filepath = os.path.join(output_dir, 'cluster_stats.txt')
    with open(filepath, 'w') as f:
        f.write(f"Total markers: {len(markers)}\n")
        f.write(f"Species with output: {', '.join(sorted(species_written))}\n")
        if markers:
            species_counts = defaultdict(int)
            for m in markers:
                for org in m['elements']:
                    species_counts[org] += 1
            f.write("\nMarkers per species:\n")
            for org in sorted(species_counts):
                f.write(f"  {org}: {species_counts[org]}\n")

            n_species_dist = defaultdict(int)
            for m in markers:
                n_species_dist[m['n_species']] += 1
            f.write("\nCluster species-count distribution:\n")
            for n in sorted(n_species_dist):
                f.write(f"  {n} species: {n_species_dist[n]} markers\n")

            sizes = [m['n_elements'] for m in markers]
            f.write(f"\nCluster size: min={min(sizes)}, max={max(sizes)}, "
                    f"mean={sum(sizes)/len(sizes):.1f}\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert AncST/MCScanX data to syngraph input TSVs.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Input source
    parser.add_argument('--approach', required=True,
                        choices=['mcscanx', 'anchors'],
                        help='Data source: mcscanx (.collinearity) or '
                             'anchors (AncST pickle files)')
    parser.add_argument('--collinearity',
                        help='Path to MCScanX .collinearity file '
                             '(required for --approach mcscanx)')
    parser.add_argument('--anchor-dir',
                        help='Directory containing AncST anchor files '
                             '(required for --approach anchors)')
    parser.add_argument('--mapping',
                        help='Path to AncST utils/mapping file '
                             '(required for --approach anchors)')

    # MCScanX options
    parser.add_argument('--mcscanx-mode', default='blocks',
                        choices=['blocks', 'anchors'],
                        help='MCScanX clustering mode: whole blocks or '
                             'individual anchor elements (default: blocks)')
    parser.add_argument('--block-merge-strategy', default='shared_anchors',
                        choices=['shared_anchors', 'overlapping_regions'],
                        help='Block merging strategy (default: shared_anchors)')
    parser.add_argument('--filter-overlaps', action='store_true',
                        help='Filter overlapping MCScanX blocks before clustering')

    # Clustering options
    parser.add_argument('--min-species', type=int, default=2,
                        help='Minimum species per marker (default: 2)')
    parser.add_argument('--mcl-inflation', type=float, default=2.0,
                        help='MCL inflation parameter (default: 2.0)')
    parser.add_argument('--no-mcl', action='store_true',
                        help='Skip MCL, use Union-Find components only')
    parser.add_argument('--no-trans', action='store_true',
                        help='Require fully connected clusters (cliques) instead '
                             'of transitively connected components')
    parser.add_argument('--mcl-impl', default='auto',
                        choices=['auto', 'external', 'python'],
                        help='MCL implementation (default: auto)')
    parser.add_argument('--min-score', type=int, default=0,
                        help='Minimum BLAST score for anchor edges (default: 0)')

    # Filtering
    parser.add_argument('--include-species',
                        help='Restrict to these species. Either a file (one alias '
                             'per line, e.g. "1org") or comma-separated list')
    parser.add_argument('--include-chrs',
                        help='Restrict to these chromosomes. Either a file (one '
                             'per line, e.g. "1orgchr12") or comma-separated list. '
                             'For mcscanx: filters blocks before clustering. '
                             'For anchors: filters anchors before clustering.')

    # Output
    parser.add_argument('--output-dir', default='syngraph_input',
                        help='Output directory (default: syngraph_input)')

    args = parser.parse_args()

    # Validate inputs
    if args.approach == 'mcscanx' and not args.collinearity:
        parser.error("--collinearity is required for --approach mcscanx")
    if args.approach == 'anchors':
        if not args.anchor_dir:
            parser.error("--anchor-dir is required for --approach anchors")
        if not args.mapping:
            parser.error("--mapping is required for --approach anchors "
                         "(anchor files use NCBI accessions that need converting)")

    for path, name in [(args.collinearity, '--collinearity'),
                       (args.mapping, '--mapping'),
                       (args.anchor_dir, '--anchor-dir')]:
        if path is not None and not os.path.exists(path):
            parser.error(f"{name} path does not exist: {path}")

    include_species = parse_filter_arg(args.include_species)
    include_chrs = parse_filter_arg(args.include_chrs)

    if include_species:
        print(f"Species filter: {len(include_species)} species")
    if include_chrs:
        print(f"Chromosome filter: {len(include_chrs)} chromosomes")

    use_mcl = not args.no_mcl
    no_trans = args.no_trans

    # Run the chosen approach
    if args.approach == 'mcscanx':
        print(f"Parsing collinearity file: {args.collinearity}")
        blocks = parse_collinearity(args.collinearity)
        print(f"  Found {len(blocks)} synteny blocks, "
              f"{sum(len(b['pairs']) for b in blocks)} anchor pairs")

        if args.filter_overlaps:
            blocks = filter_overlapping_blocks(blocks)

        if include_species:
            before = len(blocks)
            blocks = filter_blocks_by_species(blocks, include_species)
            print(f"  Species filter: {before} -> {len(blocks)} blocks")

        if include_chrs:
            before = len(blocks)
            blocks = filter_blocks_by_chrs(blocks, include_chrs)
            print(f"  Chromosome filter: {before} -> {len(blocks)} blocks")

        if args.mcscanx_mode == 'anchors':
            print(f"Clustering individual MCScanX anchors "
                  f"(MCL={'on' if use_mcl else 'off'}, "
                  f"inflation={args.mcl_inflation}, "
                  f"min_species={args.min_species}, no_trans={no_trans})")
            markers = approach_mcscanx_anchors(
                blocks, args.min_species, args.mcl_inflation,
                use_mcl, args.mcl_impl, no_trans, include_chrs)
        else:
            print(f"Clustering MCScanX blocks "
                  f"(strategy={args.block_merge_strategy}, "
                  f"MCL={'on' if use_mcl else 'off'}, "
                  f"inflation={args.mcl_inflation}, "
                  f"min_species={args.min_species}, no_trans={no_trans})")
            markers = approach_mcscanx_blocks(
                blocks, args.block_merge_strategy, args.min_species,
                args.mcl_inflation, use_mcl, args.mcl_impl, no_trans,
                include_chrs)

    elif args.approach == 'anchors':
        print(f"Loading anchors from: {args.anchor_dir}")
        org_mapping, chr_mapping = parse_mapping(args.mapping)
        if not org_mapping:
            print(f"Error: mapping file contains no organism entries: "
                  f"{args.mapping}", file=sys.stderr)
            sys.exit(1)
        markers = approach_anchors(
            args.anchor_dir, org_mapping, chr_mapping,
            args.min_species, args.mcl_inflation, use_mcl,
            args.min_score, args.mcl_impl, no_trans,
            include_species, include_chrs)

    print(f"\nGenerated {len(markers)} markers")

    species_written = write_syngraph_tsvs(markers, args.output_dir,
                                          include_species)
    write_cluster_details(markers, args.output_dir)
    write_cluster_stats(markers, args.output_dir, species_written)

    print(f"Output written to: {args.output_dir}/")
    for org in sorted(species_written):
        print(f"  {org}.tsv")
    print(f"  marker_clusters.txt")
    print(f"  cluster_stats.txt")


if __name__ == '__main__':
    main()
