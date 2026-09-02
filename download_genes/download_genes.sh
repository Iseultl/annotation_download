#!/bin/bash
set -euo pipefail

url_file=$1
genes=$2
output_dir=$3
task_id=$4

# Get the line corresponding to this SLURM array task
line=$(sed -n "${task_id}p" "${url_file}")

# Exit if the line does not exist
if [[ -z "${line}" ]]; then
    echo "No input found for array task ${task_id}"
    exit 0
fi

# Read the columns from this line
IFS=$'\t' read -r \
    taxid \
    organism_name \
    annotation_id \
    assembly_accession \
    annotation_url \
    assembly_url \
    busco_complete \
    busco_duplicated \
    busco_single_copy \
    busco_fragmented \
    busco_missing <<< "${line}"

echo "========================================"
echo "SLURM task: ${task_id}"
echo "Processing: ${organism_name}"
echo "TaxID: ${taxid}"
echo "Assembly: ${assembly_accession}"
echo "========================================"

# Replace spaces with underscores
species_name="${organism_name// /_}"

species_dir="${output_dir}/${species_name}_${taxid}"

mkdir -p "${species_dir}"

# Download files
(
    cd "${species_dir}" || exit 1

    singularity exec $HOME/singularities/python.sif \
        python $HOME/git/gitlab/annotation_download/download_genes/download_file.py \
        --taxid "${taxid}" \
        --annotation-url "${annotation_url}" \
        --fasta-url "${assembly_url}" \
        --retry-log "download_retry.tsv"
)

echo "species_dir=${species_dir}"
ls -lah "${species_dir}"

#Annocli alias match
annocli alias "${species_dir}/annotation.gff.gz" "${species_dir}/annotation.fasta.gz" --output "${species_dir}/annotation.aliasMatch.gff.gz"
gunzip -c "${species_dir}/annotation.aliasMatch.gff.gz" > "${species_dir}/annotation.aliasMatch.gff"

# Filter annotation
singularity exec $HOME/singularities/python.sif \
    python $HOME/git/gitlab/annotation_download/download_genes/filter_for_gene.py \
    --gff "${species_dir}/annotation.aliasMatch.gff" \
    --genes "${genes}" \
    --output "${species_dir}/filtered.gff"

# Create transcripts
if [[ -s "${species_dir}/filtered.gff" ]]; then

    gunzip -c "${species_dir}/annotation.fasta.gz" \
        > "${species_dir}/genome.fa"

    singularity exec $HOME/singularities/gffread.sif \
        gffread \
        -w "${species_dir}/transcripts.fa" \
        -g "${species_dir}/genome.fa" \
        "${species_dir}/filtered.gff" \
        --w-add 2000
fi

# Remove temporary/downloaded files
find "${species_dir}" -type f \
    ! -name "filtered.gff" \
    ! -name "transcripts.fa" \
    -delete