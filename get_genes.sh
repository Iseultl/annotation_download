#!/bin/bash
#SBATCH --job-name=gff_async
#SBATCH --output=/no_backup/rg/ileahy/logs/gff_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/gff_%A_%a.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-11
#SBATCH --nice=10000

set -euo pipefail

TABLE="/no_backup/rg/ileahy/mammals/SEPHS2_locations_mapped.tsv"
CHUNK_SIZE=20

START=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
END=$((START + CHUNK_SIZE))

singularity exec ~/singularities/python.sif python get_genes.py \
    --tsv_file "${TABLE}" \
    --start "${START}" \
    --end "${END}"