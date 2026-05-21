#!/bin/bash
#SBATCH --job-name=download_genes
#SBATCH --output=/no_backup/rg/ileahy/logs/download_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/download_%A_%a.err
#SBATCH --array=1-393%4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Input files
TAXID_FILE="taxids.txt"
GENES_FILE="genes.txt"
OUTDIR="/no_backup/rg/ileahy/mammals/readthrough_candidates"

# Get the taxid corresponding to this array task
TAXID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TAXID_FILE")

echo "Running taxid: $TAXID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"

# Run your script
bash download_genes.sh "$TAXID" "$GENES_FILE" "$OUTDIR"