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
