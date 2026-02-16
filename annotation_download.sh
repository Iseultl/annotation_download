#!/bin/bash
#SBATCH --job-name=gff_async
#SBATCH --output=/no_backup/rg/ileahy/logs/gff_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/gff_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-3
#SBATCH --nice=10000

set -euo pipefail

TABLE="/no_backup/rg/ileahy/fish/annotations_report.tsv"
CHUNK_SIZE=20

START=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
END=$((START + CHUNK_SIZE))



singularity exec ~/singularities/python.sif python async_gff_downloader.py "${TABLE}" "${START}" "${END}"
