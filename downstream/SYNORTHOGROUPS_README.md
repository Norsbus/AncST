# SynOrthogroups Export Tool

## Overview

`export_synorthogroups.py` generates OrthoFinder-compatible output files from synthology pipeline results.

## Key Features

- ✅ **Works for both protein and nucleotide synthology runs**
- ✅ **Complete gene accounting** - All input genes represented (either in SOGs or as unassigned)
- ✅ **Dual graph version support** - Exports both `no_new_edges` and `with_new_edges` versions
- ✅ **OrthoFinder-compatible formats** - Can use downstream OrthoFinder tools
- ✅ **Standalone script** - No modifications to existing pipeline

## Usage

### Basic Usage

```bash
# From within synthology output directory
cd synthology_out_prot
export_synorthogroups.py .

# Or specify full path
export_synorthogroups.py /path/to/synthology_out_prot
```

### Options

```bash
# Export both graph versions (default)
export_synorthogroups.py . --graph-version both

# Export only one version
export_synorthogroups.py . --graph-version no_new_edges
export_synorthogroups.py . --graph-version with_new_edges
```

## Output Structure

```
synthology_out_prot/
└── synorthogroups/
    ├── no_new_edges/
    │   ├── SynOrthogroups.tsv
    │   ├── SynOrthogroups.GeneCount.tsv
    │   ├── SynOrthogroups_SingleCopyOrthologues.txt
    │   ├── SynOrthogroups_UnassignedGenes.tsv
    │   ├── SynOrthogroups_WithCoords.tsv
    │   ├── Statistics_Overall.tsv
    │   ├── Statistics_PerSpecies.tsv
    │   └── synorthogroups.pickle
    └── with_new_edges/
        └── [same files]
```

## Output Files

### 1. SynOrthogroups.tsv

Main orthogroups table (OrthoFinder format):
- Tab-separated columns: `SynOrthogroup | Species1 | Species2 | ...`
- Genes within species: comma-separated
- Empty cells for species not in SOG

Example:
```
SynOrthogroup	GCF_000001215.4	GCF_016746365.2	...
SOG0000001	hit-42_..., hit-11_...	hit-7_...	...
SOG0000002	hit-48_...		...
```

### 2. SynOrthogroups.GeneCount.tsv

Count matrix:
```
SynOrthogroup	GCF_000001215.4	GCF_016746365.2	...	Total
SOG0000001	2	1	...	8
SOG0000002	1	0	...	3
```

### 3. SynOrthogroups_SingleCopyOrthologues.txt

Lists SOG IDs with exactly 1 gene per species (1:1:1:1:... orthologs):
```
SOG0000003
SOG0000004
```

### 4. SynOrthogroups_UnassignedGenes.tsv

**Critical for completeness** - Lists genes not in multi-species SOGs:
```
Gene_ID	Species	Reason
hit-8_GCF_030788295.1_NP_001097507.1	GCF_030788295.1	filtered_during_pipeline
hit-999_GCF_018153835.1_...	GCF_018153835.1	single_species_singleton
```

Reasons:
- `filtered_during_pipeline`: Gene was in input GFFs but removed during cograph/alignment steps
- `single_species_singleton`: Gene is in final graph but only connected to same-species genes

### 5. SynOrthogroups_WithCoords.tsv

**Self-contained coordinate information** - Detailed table with genomic coordinates for all genes:
```
SynOrthogroup	Gene_ID	Species	Chromosome	Start	End	Strand	Ref_Gene
SOG0000000	hit-32_GCF_000001215.4_NP_647920.1	GCF_000001215.4	NT_037436.4	4634469	4634789	-	NP_647920.1
SOG0000000	hit-1_GCF_000001215.4_NP_001097507.1	GCF_000001215.4	NT_037436.4	4677848	4678114	-	NP_001097507.1
```

This file makes synorthogroups output self-contained - no need to reference original pipeline files for coordinates. Used by `draw_synorthogroups.py` for visualization.

### 6. Statistics_Overall.tsv

Summary statistics:
```
Metric	Value
Number of species	11
Total input genes	102
Genes in SynOrthogroups	101
Genes unassigned (filtered)	1
Genes unassigned (singletons)	0
Percentage in SynOrthogroups	99.0
Number of SynOrthogroups	5
Number of single-copy SynOrthogroups	2
Mean SynOrthogroup size	20.2
Median SynOrthogroup size	12
```

### 7. Statistics_PerSpecies.tsv

Per-species breakdown:
```
Species	Input_Genes	In_SynOrthogroups	Unassigned	Percentage	SingleCopy_SOGs
GCF_000001215.4	9	8	1	88.9	2
GCF_016746365.2	10	10	0	100.0	2
```

### 8. synorthogroups.pickle

Python-friendly format for programmatic access:
```python
import pickle

with open('synorthogroups.pickle', 'rb') as f:
    data = pickle.load(f)

# data['synorthogroups'] - dict {SOG_ID: [gene_ids]}
# data['gene_to_sog'] - dict {gene_id: SOG_ID}
# data['unassigned'] - dict {gene_id: reason}
# data['coordinates'] - dict {gene_id: (species, chrom, start, end, strand)}
# data['metadata'] - summary information
```

## Validation

### Tested On

1. **synthology_new** (protein, buffered, TBLASTN on syntenic regions)
   - Input: 102 genes
   - SynOrthogroups: 101 genes (5 SOGs)
   - Unassigned: 1 filtered gene
   - ✓ 100% accounting verified

2. **halos_synthology_run_dedupped** (protein, BLASTP whole genome)
   - Input: 88 genes
   - SynOrthogroups: 88 genes (4 SOGs)
   - Unassigned: 0 genes
   - ✓ 100% accounting verified

### Results Match Earlier Analysis

Component structures verified to match manual analysis:

**Buffered run (no_new_edges):**
- Component sizes: [50, 18, 12, 11, 10] ✓
- Single-copy SOGs: SOG0000003, SOG0000004 ✓

**BLASTP run (no_new_edges):**
- Component sizes: [58, 12, 11, 7] ✓
- Single-copy SOGs: SOG0000002, SOG0000003 ✓

## Important Notes

### Complete Gene Accounting

The tool ensures **all input genes are represented**:
- `Genes in SynOrthogroups` + `Genes unassigned (filtered)` + `Genes unassigned (singletons)` = `Total input genes`

This is critical for publication - no genes are "lost" without documentation.

### Graph Versions

The synthology pipeline produces two union graphs:

1. **no_new_edges**: More conservative, fewer components
2. **with_new_edges**: Adds edges during alignment, more fragmented

Example from buffered run:
- no_new_edges: 5 SOGs
- with_new_edges: 8 SOGs (more fragmented)

Always export both to compare!

### Nucleotide Synthology

Works identically for nucleotide runs - looks for `parsed_fna_gff` instead of `parsed_faa_gff`.

## Integration with OrthoFinder Tools

The output formats are OrthoFinder-compatible, so you can use:
- OrthoFinder visualization tools
- Comparative genomics analysis scripts
- Any tool expecting OrthoFinder output

## Visualizing SynOrthogroups

Use `draw_synorthogroups.py` to create synteny visualizations:

```bash
# Basic usage
cd synthology_out_prot
draw_synorthogroups.py 50000 synorthogroups/no_new_edges

# Draw specific SOGs
draw_synorthogroups.py 50000 synorthogroups/no_new_edges --sogs SOG0000001,SOG0000003

# Filter by minimum species count
draw_synorthogroups.py 50000 synorthogroups/with_new_edges --min-species 5
```

Features:
- Automatically loads coordinates from `SynOrthogroups_WithCoords.tsv` (self-contained)
- Creates one SVG per SOG showing genes across species
- Optional alignment links if `pairwise_alignments_table` available
- Margin parameter controls genomic context shown around genes

Output saved to `images_synorthogroups/<margin>/SOG*_<N>species_<M>genes.svg`

## Refining SynOrthogroups with Phylogenetic Analysis

SynOrthogroups are derived from connected components in the synthology union graph. This means they include all genes that are directly or transitively connected, which can include both orthologs and paralogs. For SOGs with multiple genes from the same species (paralogs), phylogenetic refinement can help distinguish true orthologs from paralogs.

### Overview

`refine_synorthogroups_with_phylogeny.py` uses state-of-the-art phylogenetic methods to:
1. Build gene trees for each SOG using multiple sequence alignment and maximum likelihood
2. Reconcile gene trees with the species tree to identify duplication vs. speciation events
3. Generate refined ortholog groups based on evolutionary relationships

### Dependencies

Install required tools:

```bash
# MAFFT (multiple sequence alignment)
conda install -c bioconda mafft

# IQ-TREE 2 (gene tree reconstruction)
conda install -c bioconda iqtree

# Optional: GeneRax (gene/species tree reconciliation)
# See: https://github.com/BenoitMorel/GeneRax

# Optional: Notung (alternative reconciliation tool)
# Download from: http://www.cs.cmu.edu/~durand/Notung/
```

### Basic Usage

```bash
# Refine SOGs with gene tree reconstruction only
refine_synorthogroups_with_phylogeny.py \
    synthology_out_prot/synorthogroups/no_new_edges \
    --faa-dir annotations/proteins \
    --species-tree ../../utils/NJTree.nwk

# With GeneRax reconciliation
refine_synorthogroups_with_phylogeny.py \
    synthology_out_prot/synorthogroups/no_new_edges \
    --faa-dir annotations/proteins \
    --species-tree ../../utils/NJTree.nwk \
    --reconciliation-tool generax

# Test on first 10 SOGs
refine_synorthogroups_with_phylogeny.py \
    synthology_out_prot/synorthogroups/no_new_edges \
    --faa-dir annotations/proteins \
    --species-tree ../../utils/NJTree.nwk \
    --max-sogs 10
```

### What Gets Refined?

The tool automatically identifies SOGs with paralogs (multiple genes from the same species):

```
Found 47 SynOrthogroups
SOGs with paralogs (multi-copy): 12
SOGs without paralogs: 35
```

Only SOGs with paralogs are processed (single-copy SOGs don't need refinement).

### Output Structure

```
synorthogroups/no_new_edges/
└── phylogeny_refinement/
    ├── refinement_summary.txt
    ├── refinement_results.pickle
    ├── SOG0000001/
    │   ├── SOG0000001.faa            # Extracted sequences
    │   ├── SOG0000001.aln            # MAFFT alignment
    │   ├── SOG0000001.treefile       # IQ-TREE gene tree
    │   ├── SOG0000001.iqtree         # IQ-TREE log
    │   └── generax/                  # Reconciliation results (if enabled)
    ├── SOG0000005/
    │   └── ...
    └── ...
```

### Methods

#### Gene Tree Reconstruction

Uses IQ-TREE 2 with:
- **ModelFinder Plus (MFP)**: Automatically selects best substitution model
- **Ultrafast Bootstrap (UFBoot)**: 1000 replicates for branch support
- **Maximum Likelihood**: State-of-the-art accuracy

#### Reconciliation

Two tools supported:

**GeneRax** (recommended):
- Joint maximum likelihood inference
- DTL model (Duplication-Transfer-Loss)
- Co-estimates gene tree and reconciliation

**Notung**:
- Fast parsimony-based reconciliation
- Good alternative when GeneRax is unavailable

### Interpreting Results

The `refinement_summary.txt` file reports:

```
SynOrthogroups Phylogenetic Refinement Summary
================================================================================

Gene trees built: 12
Trees reconciled: 12
Failed: 0
```

For each SOG:
- **Gene tree** (`.treefile`): Shows evolutionary relationships among genes
- **Reconciled tree** (if enabled): Maps duplication/speciation events
- **IQ-TREE log** (`.iqtree`): Model selection, bootstrap support, tree statistics

### Performance Notes

- **MAFFT**: Fast (~1-10 seconds per SOG)
- **IQ-TREE**: Moderate (10-60 seconds per SOG depending on size)
- **GeneRax**: Slow (1-10 minutes per SOG)

For large analyses (100+ SOGs with paralogs), consider:
1. Running with `--reconciliation-tool none` first to get gene trees
2. Using `--threads` to parallelize IQ-TREE
3. Running reconciliation separately on SOGs of interest

### State-of-the-Art Methods (Dec 2024 - Jan 2025)

This tool uses cutting-edge phylogenetic methods:

**Gene Trees**:
- IQ-TREE 2 (Minh et al. 2020): Best accuracy/speed balance
- Alternative: FastTree (faster but less accurate)

**Reconciliation**:
- GeneRax (Morel et al. 2020): Joint ML inference
- AleRax (NEW 2024): Probabilistic co-estimation
- Notung (Chen et al. 2000): Fast parsimony

**Comparison to OrthoFinder**:
- OrthoFinder 2.4+ uses hierarchical orthogroups (HOGs) approach
- Reported 12-20% more accurate than graph-based methods
- Synthology adds synteny constraint, complementary to HOGs

## Troubleshooting

### "File not found" errors

Ensure you're in a synthology output directory containing:
- `final_union_graph_no_new_edges.pickle`
- `final_union_graph_with_new_edges.pickle`
- `parsed_faa_gff` (or `parsed_fna_gff`)

### Gene accounting mismatch

If you see:
```
WARNING: Gene accounting mismatch!
  Input: 102, Accounted: 101
```

This indicates a bug - please report with the specific run directory.

## Future Enhancements

Potential additions (not currently implemented):
- Pairwise ortholog tables (Species1__v__Species2.tsv)
- ~~Gene tree inference~~ ✓ Implemented (see Phylogenetic Refinement section)
- ~~Duplication event detection~~ ✓ Implemented via reconciliation
- FASTA sequence files per orthogroup
- Automated creation of refined ortholog groups from reconciliation results

## Citation

If using SynOrthogroups output in publications, cite both:
1. The synthology pipeline (your publication)
2. OrthoFinder (since we use their output format)

## Contact

For issues or questions about this tool, please open an issue in the repository.
