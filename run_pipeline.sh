#!/bin/bash
#SBATCH --job-name=download_genes
#SBATCH --output=/no_backup/rg/ileahy/logs/download_genes_%A_%a.out
#SBATCH --error=/no_backup/rg/ileahy/logs/download_genes_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-105
#SBATCH --nice=10000

set -euo pipefail

TABLE="/no_backup/rg/ileahy/fish/annotations_report.tsv"
CHUNK_SIZE=5

START=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE))
END=$((START + CHUNK_SIZE))


python annotation_download/pipeline.py \
  --annotations-report ${TABLE} \
  --gene 'SEPHS2|SPS2|SEPHS3' \
  --outdir /no_backup/rg/ileahy/fish/pipeline_out \
  --start ${START} \
  --end ${END} \
  --gffread-container quay.io/biocontainers/gffread:0.12.7--h077b44d_6 \
  --container-runtime singularity \
  --concurrency 2