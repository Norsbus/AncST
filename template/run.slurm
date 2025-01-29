#!/usr/bin/env bash
#SBATCH --job-name=metadata_GCF_026898425.1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --time=3-00:00:00
#SBATCH --error log/metadata/stderr_metadata_GCF_026898425.1
#SBATCH --output log/metadata/stdout_metadata_GCF_026898425.1
    eval "$(conda shell.bash hook)"
    conda activate CONDAHOMEDIR/miniconda3/envs/snakemake
    ./make_blastdb_and_metadata.py GCF_026898425.1 && touch touch/GCF_026898425.1_metadata_done
    