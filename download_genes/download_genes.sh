#!/bin/bash
set -euo pipefail
taxid=$1
genes=$2
output_dir=$3

# Create the output directory
mkdir -p "${output_dir}"

#Functions to check if assemblies are present in directory 
has_fasta() {
    local dir=$1
    compgen -G "${dir}/*.fa*" > /dev/null || \
    compgen -G "${dir}/*.fna*" > /dev/null || \
    compgen -G "${dir}/*.fasta*" > /dev/null
}

has_gff() {
    local dir=$1
    compgen -G "${dir}/*.gff*" > /dev/null || \
    compgen -G "${dir}/*.gtf*" > /dev/null
}

# Create the descendants list to iterate through
annocli download --ref-only --taxids "${taxid}" --mode links | cut -d ' ' -f3 | cut -d '/' -f2,3 > "${output_dir}/descendants_list.txt"
dirs=$(cat "${output_dir}/descendants_list.txt")

for i in ${dirs}; do
    echo "${i}"
    taxon=$(echo "${i}" | cut -d '/' -f1 | awk -F'_' '{print $NF}')
    # Download the annotations for the taxid 
    if ! has_fasta "${output_dir}/${i}" || ! has_gff "${output_dir}/${i}"; then
        annocli download --ref-only --taxids "${taxon}" --add_asm --fix_alias --output "${output_dir}"
    fi
    rm -f "${output_dir}/${i}/"*.aliasMappings.tsv 2>/dev/null || true
    # Filter the annotations for the genes of interest
    singularity exec ~/singularities/python.sif python filter_for_gene.py --gff "$(ls "${output_dir}/${i}/"*.aliasMatch.*)" --genes "${genes}" --output "${output_dir}/${i}/filtered.gff"
    # Create the transcripts file with gffread
    if [ -s "${output_dir}/${i}/filtered.gff" ]; then
        genome=$(echo ${output_dir}/${i}/*.f*)
        echo "${genome}"
        gunzip -c "$genome" > "${output_dir}/${i}/genome.fa"
        singularity exec ~/singularities/gffread.sif gffread -w "${output_dir}/${i}/transcripts.fa" -g "${output_dir}/${i}/genome.fa" "${output_dir}/${i}/filtered.gff" --w-add 600
    fi
    # Remove the OG annotations and genome file
    find "${output_dir}/${i}" -type f \
        ! -name "filtered.gff" \
        ! -name "transcripts.fa" \
        -delete
done
