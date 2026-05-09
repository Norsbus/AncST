# Project Name

> AncST is an alignment engine designed to produce reciprocal best hits based on their uniqueness within their own genome. It works efficiently using string statistics and can be universally applied to any genome(s). The downstream results contain useful data and scripts to determine macrosynteny with MCScanX, use the anchors for homology-based scaffolding, and provides pipelines to infer orthology relationships with the help of conserved genomic position.

---

## Requirements

- Operating system: Linux/MacOS
- Python: 3.10+
- External tools: GenMap, macle, BLAST, clasp.x, NCBI datasets
- Memory / disk: see GenMap and macle - those are the memory bottlenecks

## Installation

```bash
./install.sh
```

Run `./install.sh --help` to see all available options.

## Usage

```bash
./run_pipeline.sh --species examples/orgs.example --cores 8 --max-mem 16000
```

Run `./run_pipeline.sh --help` to see all available options.

## Input

Describe the input format(s):

- **Species file** — plain text, one genome identifier per line. They must correspond to {name} of genomes found in utils/genomes/ (create dir if it doesnt exist!) like utils/genomes/{name}.fasta. If they are NCBI accessions and not found in utils/genomes/ the script will automatically download them if NCBI datasets binary is available in path.

## Output

What the pipeline produces and where it lives:

- `anchors/<organism>/` — per-genome anchor candidates
- `anchors/aligned/<organism>` and `anchors/aligned_succinct/<organism>` — pairwise anchor alignments
- `downstream/...` — directory with useful downstream data and script, see READMEs there. An archive with the main results will be created there which is most useful. Use the installation script there to install additional requirements to fully use the downstream features.

## Configuration

Key configuration files:

- `template/pipeline_config.yaml` — pipeline parameters (window sizes, score thresholds).
- `code/config.yaml` — Snakemake configuration.
- `template/genmap_params.txt` / `template/macle_params.txt` / `template/dups_params.txt`— per-organism parameter sets (generated at run time if not provided, otherwise see examples/)

## Webserver/Tutorial
There is a webserver at [https://anchored.bioinf.uni-leipzig.de/](https://anchored.bioinf.uni-leipzig.de/) implementing this pipeline where users can run jobs on our server. On this website there also exist a tutorials for downstream application of AncST anchors and respective homology-based scaffolding. Moreover, precomputed anchor candidates and publicly available sets of computed anchors can be downloaded.

## Citation

If you use AncST, please cite:

1. Käther K, Lemke S, Stadler PF. (2023). Annotation-free Identification of Potential Synteny Anchors. In *International Work-Conference on Bioinformatics and Biomedical Engineering (IWWBIO'23)*, Lect. Notes Comp. Sci., **13919**, 217–230. <https://doi.org/10.1007/978-3-031-34953-9_17>
2. Käther KK, Remmel A, Lemke S, Stadler PF. (2025). Unbiased anchors for reliable genome-wide synteny detection. *Algorithms for Molecular Biology*, **20**, 5. <https://doi.org/10.1186/s13015-025-00275-9>
3. Käther KK, Gatter T, Lemke S, Stadler PF. (2025). Anchors for Homology-Based Scaffolding. *bioRxiv* preprint, 2025.04.28.650980. <https://doi.org/10.1101/2025.04.28.650980>

## License

This project is released under the [GNU General Public License v3.0](LICENSE).

## Contact

- karl.kaether@posteo.de
- karl@bioinf.uni-leipzig.de
