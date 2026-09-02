#!/bin/bash
#SBATCH --job-name=download_genes
#SBATCH --output=/no_backup/rg/ileahy/logs/download_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/download_%A_%a.err
#SBATCH --array=0-10
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Input files
TAXID_FILE="annotations.tsv"
GENES_FILE="genes.txt"
OUTDIR="/no_backup/rg/ileahy/mammals/readthrough_0926"

# Pass the array task ID
bash download_genes.sh \
    "$TAXID_FILE" \
    "$GENES_FILE" \
    "$OUTDIR" \
    "$SLURM_ARRAY_TASK_ID"