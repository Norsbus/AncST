#!/usr/bin/env bash
#SBATCH --job-name=macle_Cm_GCF_945859705.1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --time=3-00:00:00
#SBATCH --error log/macle_Cm/stderr_map_GCF_945859705.1
#SBATCH --output log/macle_Cm/stdout_map_GCF_945859705.1
eval "$(conda shell.bash hook)"
conda activate /homes/biertank/karl//miniconda3/envs/snakemake
macle -i /scr/k80san/karl/new_AncST_server/AncST/template/../utils/macle_indices/GCF_945859705.1 -w 1000 -k 50 > /scr/k80san/karl/new_AncST_server/AncST/template/../utils/macle_out/GCF_945859705.1/1000_50.txt && touch touch/GCF_945859705.1_1000_50_macle_Cm_done
