#!/bin/bash
set -euo pipefail
shopt -s nullglob

usage() {
    cat <<'EOF'
Usage: bash scrap.sh <genes_file> <output_dir> [--force]

Rebuild missing filtered.gff and transcripts.fa without re-downloading.
Only directories that already contain aliasMatch GFF/GTF and a genome FASTA
are processed. Existing outputs are skipped unless --force is provided.
EOF
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
    usage
    exit 1
fi

genes=$1
output_dir=$2
force="false"

if [[ $# -eq 3 ]]; then
    case "$3" in
        --force|--all)
            force="true"
            ;;
        *)
            usage
            exit 1
            ;;
    esac
fi

if [[ ! -f "$genes" ]]; then
    echo "Genes file not found: $genes" >&2
    exit 1
fi

if [[ ! -d "$output_dir" ]]; then
    echo "Output directory not found: $output_dir" >&2
    exit 1
fi

find_alias_gff() {
    local dir=$1
    local candidates=()
    local f

    candidates=(
        "$dir"/*.aliasMatch.gff*
        "$dir"/*.aliasMatch.gff3*
        "$dir"/*.aliasMatch.gtf*
        "$dir"/*.aliasMatch.*
    )

    for f in "${candidates[@]}"; do
        if [[ -f "$f" ]]; then
            if [[ "$f" == *".gff"* || "$f" == *".gtf"* ]]; then
                echo "$f"
                return 0
            fi
        fi
    done

    return 1
}

find_genome_fasta() {
    local dir=$1
    local candidates=()
    local f

    candidates=(
        "$dir"/*.fa*
        "$dir"/*.fna*
        "$dir"/*.fasta*
    )

    for f in "${candidates[@]}"; do
        if [[ -f "$f" && "$(basename "$f")" != "genome.fa" ]]; then
            echo "$f"
            return 0
        fi
    done

    return 1
}

process_dir() {
    local dir=$1
    local filtered="$dir/filtered.gff"
    local transcripts="$dir/transcripts.fa"

    if [[ "$force" != "true" && -s "$filtered" && -s "$transcripts" ]]; then
        echo "Skipping (outputs exist): $dir"
        return 0
    fi

    local gff
    if ! gff=$(find_alias_gff "$dir"); then
        echo "Skipping (missing aliasMatch GFF): $dir"
        return 0
    fi

    rm -f "$dir"/*.aliasMappings.tsv 2>/dev/null || true

    echo "Filtering: $dir"
    if ! singularity exec ~/singularities/python.sif \
        python filter_for_gene.py --gff "$gff" --genes "$genes" --output "$filtered"; then
        echo "Filter failed: $dir"
        return 0
    fi

    if [[ ! -s "$filtered" ]]; then
        echo "Filtered GFF is empty: $dir"
        return 0
    fi

    local genome
    if ! genome=$(find_genome_fasta "$dir"); then
        echo "Skipping (missing genome FASTA): $dir"
        return 0
    fi

    local genome_tmp="$dir/genome.fa"
    if [[ "$genome" == *.gz ]]; then
        if ! gunzip -c "$genome" > "$genome_tmp"; then
            echo "Failed to unzip genome: $dir"
            return 0
        fi
    else
        if ! cp "$genome" "$genome_tmp"; then
            echo "Failed to copy genome: $dir"
            return 0
        fi
    fi

    echo "Creating transcripts: $dir"
    if ! singularity exec ~/singularities/gffread.sif \
        gffread -w "$transcripts" -g "$genome_tmp" "$filtered" --w-add 600; then
        echo "gffread failed: $dir"
        rm -f "$genome_tmp"
        return 0
    fi

    rm -f "$genome_tmp"

    if [[ ! -s "$transcripts" ]]; then
        echo "Transcripts file is empty: $dir"
        return 0
    fi

    find "$dir" -type f \
        ! -name "filtered.gff" \
        ! -name "transcripts.fa" \
        -delete
}

gff_dirs=()
while IFS= read -r -d '' gff; do
    gff_dirs+=("$(dirname "$gff")")
done < <(find "$output_dir" -type f -name "*.aliasMatch.*" -print0)

if [[ ${#gff_dirs[@]} -eq 0 ]]; then
    echo "No aliasMatch files found under: $output_dir"
    exit 0
fi

unique_dirs=$(printf '%s\n' "${gff_dirs[@]}" | sort -u)
while IFS= read -r dir; do
    [[ -z "$dir" ]] && continue
    process_dir "$dir"
done <<< "$unique_dirs"