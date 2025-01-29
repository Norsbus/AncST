#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=32000
#SBATCH --time=7-00:00:00
#SBATCH --output=log/out_main_%j
#SBATCH --error=log/err_main_%j
#SBATCH --partition=r0bornb

module load Miniconda3/4.12.0
eval "$(conda shell.bash hook)"
conda activate snakemake
./run_checks.py
