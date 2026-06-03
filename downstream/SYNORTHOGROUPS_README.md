# SynOrthogroups Export Tool

`export_synorthogroups.py` generates OrthoFinder-compatible output files from
synthology pipeline results. It works for both protein and nucleotide
synthology runs, accounts for all input genes (either in SOGs or as
unassigned), supports both `no_new_edges` and `with_new_edges` graph versions,
and runs as a standalone script with no modifications to the existing pipeline.

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
hit-790_Solanum_lycopersicum_Solyc07g042630.4.1	Solanum_lycopersicum	filtered_out_due_to_overlap_with_larger_seq
hit-8_GCF_030788295.1_NP_001097507.1	GCF_030788295.1	filtered_during_pipeline
hit-999_GCF_018153835.1_...	GCF_018153835.1	single_species_singleton
```

Reasons:
- `filtered_out_due_to_overlap_with_larger_seq`: Gene overlapped with a larger gene during GFF parsing and was removed (see `filter_overlapping_sequences.py` to pre-filter)
- `filtered_during_pipeline`: Gene was in parsed GFFs but not in any synteny chains
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
Total input genes (in parsed GFF)	102
Genes filtered during GFF parsing (overlaps)	38
Genes in SynOrthogroups	101
Genes unassigned (not in synteny chains)	1
Genes unassigned (single-species singletons)	0
Percentage in SynOrthogroups	99.0
Number of SynOrthogroups	5
Number of single-copy SynOrthogroups	2
Mean SynOrthogroup size	20.2
Median SynOrthogroup size	12
```

### 7. Statistics_PerSpecies.tsv

Per-species breakdown:
```
Species	Input_Genes	Overlap_Filtered	In_SynOrthogroups	Unassigned	Percentage	SingleCopy_SOGs
GCF_000001215.4	9	2	8	1	88.9	2
GCF_016746365.2	10	0	10	0	100.0	2
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
   - 100% accounting verified

2. **halos_synthology_run_dedupped** (protein, BLASTP whole genome)
   - Input: 88 genes
   - SynOrthogroups: 88 genes (4 SOGs)
   - Unassigned: 0 genes
   - 100% accounting verified

### Results Match Earlier Analysis

Component structures verified to match manual analysis:

**Buffered run (no_new_edges):**
- Component sizes: [50, 18, 12, 11, 10]
- Single-copy SOGs: SOG0000003, SOG0000004

**BLASTP run (no_new_edges):**
- Component sizes: [58, 12, 11, 7]
- Single-copy SOGs: SOG0000002, SOG0000003

## Important Notes

### Complete Gene Accounting

The tool ensures **all genes are accounted for**:
- Genes that made it into parsed GFF: `Genes in SynOrthogroups` + `Genes unassigned (filtered)` + `Genes unassigned (singletons)` = `Total input genes (in parsed GFF)`
- Genes filtered during GFF parsing: Reported separately as `Genes filtered during GFF parsing (overlaps)`

All filtered genes are documented in `SynOrthogroups_UnassignedGenes.tsv` with specific reasons.

This is critical for publication - no genes are "lost" without documentation.

**Note:** If many genes are filtered due to overlaps, consider using `filter_overlapping_sequences.py` (see Troubleshooting section below) to pre-filter your input files and keep only the longest isoforms.

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

### Many genes filtered due to overlaps

If `Statistics_Overall.tsv` shows many genes in "Genes filtered during GFF parsing (overlaps)", you likely have:
- Multiple isoforms per gene
- Alternative transcripts
- Overlapping gene annotations

**Solution:** Pre-filter your input files to keep only the longest isoform per locus:

```bash
# Create filtered input files
filter_overlapping_sequences.py input_annotations/ filtered_annotations/

# Then run synthology with filtered files
run_synthology_prot.py \
  --annotation-dir filtered_annotations/ \
  --output-dir synthology_out_prot \
  --mcscanx-dir mcscanx_results/ \
  --pairwise-alignments pairwise_alignments_table
```

This will reduce redundancy and give cleaner synorthogroups.

## Citation

If using SynOrthogroups output in publications, cite both:
1. The synthology pipeline (your publication)
2. OrthoFinder (since we use their output format)

## Contact

For issues or questions about this tool, please open an issue in the repository.
