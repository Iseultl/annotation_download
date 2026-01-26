#!/bin/bash
#SBATCH --job-name=download_genes
#SBATCH --output=/no_backup/rg/ileahy/logs/download_genes_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/download_genes_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-11%3
#SBATCH --nice=10000

set -euo pipefail

TABLE="/no_backup/rg/ileahy/mammals/SEPHS2_locations_mapped.tsv"
ASSEMBLY_URLS="/no_backup/rg/ileahy/mammals/assembly_urls.json"
GENOME_DIR="/no_backup/rg/ileahy/mammals"

CHUNK_SIZE=20

START=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
END=$((START + CHUNK_SIZE))

singularity exec ~/singularities/python.sif python get_genes.py \
    --tsv_file "${TABLE}" \
    --assembly_urls "${ASSEMBLY_URLS}" \
    --start "${START}" \
    --end "${END}" \
    --outdir "${GENOME_DIR}"