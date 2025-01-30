#!/usr/bin/env bash
#SBATCH --job-name=genmap_GCF_026898425.1
#SBATCH --cpus-per-task=30
#SBATCH --mem=64000
#SBATCH --time=3-00:00:00
#SBATCH --error log/gm_index/stderr_index_GCF_026898425.1
#SBATCH --output log/gm_index/stdout_index_GCF_026898425.1
eval "$(conda shell.bash hook)"
conda activate CONDAHOMEDIR/miniconda3/envs/snakemake
genmap index -F /scr/k80san/karl/AncST_test/template/../utils/genomes/GCF_026898425.1.fasta -I /scr/k80san/karl/AncST_test/template/../utils/genmap_indices/GCF_026898425.1 && touch touch/GCF_026898425.1_index_done
