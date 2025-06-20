#!/usr/bin/env bash
#SBATCH --job-name=genmap_GCA_030067965.1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64000
#SBATCH --time=3-00:00:00
#SBATCH --error log/gm_map/stderr_map_GCA_030067965.1
#SBATCH --output log/gm_map/stdout_map_GCA_030067965.1
eval "$(conda shell.bash hook)"
conda activate /homes/biertank/karl//miniconda3/envs/snakemake
genmap map -fl -K 50 -E 4 -I /scr/k80san/karl/AncST_server/AncST/template/../utils/genmap_indices/GCA_030067965.1 -O /scr/k80san/karl/AncST_server/AncST/template/../utils/genmap_out/GCA_030067965.1/50_4 -r && touch touch/GCA_030067965.1_50_4_map_done
