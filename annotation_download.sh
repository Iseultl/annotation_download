#!/bin/bash
#SBATCH --job-name=gff_async
#SBATCH --output=/no_backup/rg/ileahy/logs/gff_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/gff_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1,2,3,4,13,14,15
#SBATCH --nice=10000

set -euo pipefail

TABLE="/no_backup/rg/ileahy/mammals/annotations_20260115_090352.tsv"
CHUNK_SIZE=20

START=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
END=$((START + CHUNK_SIZE))



singularity exec ~/singularities/python.sif python async_gff_downloader.py "${TABLE}" "${START}" "${END}"
