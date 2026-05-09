#!/usr/bin/env bash
#SBATCH --job-name=genmap_tick_GCF_016920785.2
#SBATCH --cpus-per-task=30
#SBATCH --mem=64000
#SBATCH --time=3-00:00:00
#SBATCH --error stderr_map_tick_GCF_016920785.2
#SBATCH --output stdout_map_tick_GCF_016920785.2
eval "$(conda shell.bash hook)"
conda activate /homes/biertank/karl//miniconda3/envs/snakemake
genmap map -fl -K 13 -E 0 -I ../genmap_indices/tick_GCF_016920785.2 -O ../genmap_out/tick_GCF_016920785.2/13_0 -r
