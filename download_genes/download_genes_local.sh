#!/bin/bash

set -euo pipefail

# ============================================================
# Usage:
#
# ./download_genes_local.sh <url_file> <genes> <output_dir>
#
# Example:
# ./download_genes_local.sh \
#     annotations.tsv \
#     genes.txt \
#     /Users/iseult/Desktop/Mammalian_selenoproteins
#
# ============================================================

url_file=$1
genes=$2
output_dir=$3

# Path to your local annotation_download repository
ANNOTATION_DOWNLOAD="$HOME/gitlab/annotations_download/annotation_download"

# Scripts
DOWNLOAD_SCRIPT="${ANNOTATION_DOWNLOAD}/download_genes/download_file.py"
FILTER_SCRIPT="${ANNOTATION_DOWNLOAD}/download_genes/filter_for_gene.py"


# ============================================================
# Check inputs
# ============================================================

if [[ ! -f "$url_file" ]]; then
    echo "ERROR: URL file not found:"
    echo "$url_file"
    exit 1
fi

if [[ ! -f "$genes" ]]; then
    echo "ERROR: Genes file not found:"
    echo "$genes"
    exit 1
fi

if [[ ! -f "$DOWNLOAD_SCRIPT" ]]; then
    echo "ERROR: Download script not found:"
    echo "$DOWNLOAD_SCRIPT"
    exit 1
fi

if [[ ! -f "$FILTER_SCRIPT" ]]; then
    echo "ERROR: Filter script not found:"
    echo "$FILTER_SCRIPT"
    exit 1
fi


mkdir -p "$output_dir"


# ============================================================
# Process each organism
# ============================================================

while IFS=$'\t' read -r \
    taxid \
    organism_name \
    annotation_id \
    annotation_url \
    assembly_accession \
    assembly_url \
    busco_complete \
    busco_duplicated \
    busco_single_copy \
    busco_fragmented \
    busco_missing
do

    # Skip empty lines
    [[ -z "$taxid" ]] && continue

    echo
    echo "========================================"
    echo "Processing: ${organism_name}"
    echo "TaxID: ${taxid}"
    echo "Assembly: ${assembly_accession}"
    echo "========================================"


    # --------------------------------------------------------
    # Create species directory
    # --------------------------------------------------------

    species_name="${organism_name// /_}"

    species_dir="${output_dir}/${species_name}_${taxid}"

    mkdir -p "$species_dir"

    # --------------------------------------------------------
    # Skip if already processed
    # --------------------------------------------------------

    if [[ -s "${species_dir}/transcripts.fa" ]]; then
	echo "========================================"
    	echo "Already processed: ${organism_name}"
    	echo "Found: ${species_dir}/transcripts.fa"
    	echo "Skipping..."
    	echo "========================================"
    	continue
    fi

    # --------------------------------------------------------
    # Download annotation + genome
    # --------------------------------------------------------

    echo "Downloading annotation and genome..."

    (
        cd "$species_dir"

        python "$DOWNLOAD_SCRIPT" \
            --taxid "$taxid" \
            --annotation-url "$annotation_url" \
            --fasta-url "$assembly_url" \
            --retry-log "download_retry.tsv"
    )


    echo "Download complete"


    # --------------------------------------------------------
    # Annocli alias match
    # --------------------------------------------------------

    echo "Running annocli alias match..."

    annocli alias \
        "${species_dir}/annotation.gff.gz" \
        "${species_dir}/annotation.fasta.gz" \
        --output "${species_dir}/annotation.aliasMatch.gff.gz"


    gunzip -c \
        "${species_dir}/annotation.aliasMatch.gff.gz" \
        > "${species_dir}/annotation.aliasMatch.gff"


    # --------------------------------------------------------
    # Filter annotation
    # --------------------------------------------------------

    echo "Filtering annotation for requested genes..."

    if ! python "$FILTER_SCRIPT" \
        --gff "${species_dir}/annotation.aliasMatch.gff" \
        --genes "$genes" \
        --output "${species_dir}/filtered.gff"
    then

        echo "No matching genes found for ${species_dir}"
        rm -rf "$species_dir"
        continue

    fi


    # --------------------------------------------------------
    # Create transcripts
    # --------------------------------------------------------

    if [[ -s "${species_dir}/filtered.gff" ]] && [[ -s "${species_dir}/annotation.aliasMatch.gff.gz.aliasMappings.tsv" ]]; then

        echo "Creating transcript sequences..."

        gunzip -c \
            "${species_dir}/annotation.fasta.gz" \
            > "${species_dir}/genome.fa"


        if ! command -v gffread &> /dev/null; then
            echo "ERROR: gffread is not installed or not in PATH"
            exit 1
        fi


        gffread \
            -w "${species_dir}/transcripts.fa" \
            -g "${species_dir}/genome.fa" \
            "${species_dir}/filtered.gff" \
            --w-add 2000

    else

        echo "No genes found in filtered annotation"
        rm -rf "$species_dir"
        continue

    fi


    # --------------------------------------------------------
    # Remove temporary/downloaded files
    # --------------------------------------------------------

    echo "Cleaning up..."

    find "$species_dir" -type f \
        ! -name "filtered.gff" \
        ! -name "transcripts.fa" \
        -delete


    echo "Finished: ${organism_name}"

done < "$url_file"


echo
echo "========================================"
echo "All organisms processed"
echo "========================================"
