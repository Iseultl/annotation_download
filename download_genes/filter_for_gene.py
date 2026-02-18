#!/usr/bin/env python3

import sys
import gzip
from pathlib import Path
import argparse

# --------------------------------------------------
# Utilities
# --------------------------------------------------
def open_gff(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_attributes(attr_str):
    """Convert GFF3 attribute column into dict"""
    attrs = {}
    for item in attr_str.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.lower()] = v
    return attrs


def load_genes(arg):
    """Accept comma list OR file"""
    p = Path(arg)
    if p.exists():
        return {l.strip().lower() for l in p.read_text().splitlines() if l.strip()}
    return {g.strip().lower() for g in arg.split(",") if g.strip()}


# --------------------------------------------------
# Core filtering logic
# --------------------------------------------------
def filter_gff(gff_path, genes):

    orig_lines = []
    lower_lines = []
    with open_gff(gff_path) as f:
        for line in f:
            if not line.startswith("#"):
                stripped = line.rstrip("\n")
                orig_lines.append(stripped)
                lower_lines.append(stripped.lower())

    # pass 1 - find matching gene IDs
    gene_ids = set()

    for line in lower_lines:
        cols = line.split("\t")
        if len(cols) < 9:
            continue

        feature = cols[2]
        attrs = parse_attributes(cols[8])

        if feature == "gene":
            gene_name = attrs.get("gene") or attrs.get("name")
            if gene_name in genes:   # exact match
                if "id" in attrs:
                    gene_ids.add(attrs["id"])

    if not gene_ids:
        return False, []

    # pass 2 - collect all descendants
    keep_ids = set(gene_ids)
    changed = True

    while changed:
        changed = False
        for line in lower_lines:
            cols = line.split("\t")
            if len(cols) < 9:
                continue

            attrs = parse_attributes(cols[8])
            parents = attrs.get("parent", "").split(",")

            if any(p in keep_ids for p in parents):
                if "id" in attrs and attrs["id"] not in keep_ids:
                    keep_ids.add(attrs["id"])
                    changed = True

    # pass 3 - output original-case lines
    out_lines = []

    for orig, lower in zip(orig_lines, lower_lines):
        cols = lower.split("\t")
        if len(cols) < 9:
            continue

        attrs = parse_attributes(cols[8])
        parents = attrs.get("parent", "").split(",")

        if (
            attrs.get("id") in keep_ids
            or any(p in keep_ids for p in parents)
        ):
            out_lines.append(orig)

    return True, out_lines


# --------------------------------------------------
# Main
# --------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Filter a GFF file for gene(s) of interest")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--gff", type=str, help="Path to a (alias-matched) GFF/GFF3 file (.gz ok)")
    parser.add_argument("--genes", required=True, type=str, help="Comma-separated list of genes to filter for")
    parser.add_argument(
        "--output",
        default="filtered.gff",
        help="Output path (default: filtered.gff next to the input GFF)",
    )
    args = parser.parse_args()

    genes = load_genes(args.genes)

    gff_path = Path(args.gff)

    if not gff_path.exists():
        print(f"GFF file not found: {gff_path}")
        raise SystemExit(2)

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = gff_path.parent / out_path
    if out_path.exists():
        out_path.unlink()

    print("Genes:", ", ".join(sorted(genes)))
    print("Input:", gff_path)
    print("Output:", out_path)
    print()

    found, lines = filter_gff(gff_path, genes)
    if not found:
        print("No matching genes found in GFF file.")
        raise SystemExit(2)

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Wrote:", out_path)

if __name__ == "__main__":
    main()