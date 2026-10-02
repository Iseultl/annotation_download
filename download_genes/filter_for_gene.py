#!/usr/bin/env python3

import sys
import gzip
from pathlib import Path
import argparse
import re
from urllib.parse import unquote

# --------------------------------------------------
# Utilities
# --------------------------------------------------
def open_gff(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def parse_attributes(attr_str):
    if attr_str is None:
        return {}

    if not isinstance(attr_str, str):
        return {}

    attr_str = attr_str.strip()

    if not attr_str or attr_str == ".":
        return {}

    attrs = {}
    items = re.split(r';(?=(?:[^"]*"[^"]*")*[^"]*$)', attr_str)

    for item in items:
        item = item.strip()

        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            key = key.strip().lower()
            value = value.strip()

            # URL-decode GFF3 values
            value = unquote(value)

            # Remove surrounding quotes if present
            if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
                value = value[1:-1]
        else:
            match = re.match(
                r'^\s*([^\s]+)\s+(.*?)\s*$',
                item
            )

            if not match:
                key = item.strip().lower()
                value = True
            else:
                key, value = match.groups()
                key = key.strip().lower()
                value = value.strip()

                # Remove surrounding quotes
                if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
                    value = value[1:-1]
                value = unquote(value)
        if key in attrs:
            if isinstance(attrs[key], list):
                attrs[key].append(value)
            else:
                attrs[key] = [attrs[key], value]
        else:
            attrs[key] = value
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
def _as_list(value):
    """Return an attribute value as a list."""
    if value is None:
        return []

    if isinstance(value, list):
        return value

    return [value]


def filter_gff(gff_path, genes):

    # --------------------------------------------------
    # Read file
    # --------------------------------------------------

    orig_lines = []
    parsed_records = []

    with open_gff(gff_path) as f:
        for line in f:

            if line.startswith("#"):
                continue

            stripped = line.rstrip("\n")

            cols = stripped.split("\t")

            if len(cols) < 9:
                continue

            feature = cols[2].lower()
            attrs = parse_attributes(cols[8])

            orig_lines.append(stripped)

            parsed_records.append({
                "line": stripped,
                "feature": feature,
                "attrs": attrs,
            })

    target_gene_ids = set()
    target_gene_names = set()

    for record in parsed_records:

        if record["feature"] != "gene":
            continue

        attrs = record["attrs"]

        # Possible names for the biological gene name
        candidate_names = []

        for key in (
            "gene_name",
            "gene",
            "name",
            "gene_symbol",
        ):
            candidate_names.extend(_as_list(attrs.get(key)))

        # Possible gene identifiers
        candidate_ids = []

        for key in (
            "gene_id",
            "id",
        ):
            candidate_ids.extend(_as_list(attrs.get(key)))

        # Normalise
        candidate_names = {
            str(x).strip().lower()
            for x in candidate_names
            if x is not None
        }

        candidate_ids = {
            str(x).strip().lower()
            for x in candidate_ids
            if x is not None
        }

        # Check either gene name or gene ID
        matched = (
            bool(candidate_names & genes)
            or bool(candidate_ids & genes)
        )

        if matched:

            target_gene_names.update(candidate_names)
            target_gene_ids.update(candidate_ids)

    # --------------------------------------------------
    # If no genes were found, provide useful diagnostics
    # --------------------------------------------------

    if not target_gene_ids and not target_gene_names:

        print("No matching genes found.")

        # Print a few example gene annotations so that
        # we can see what the file actually contains.
        print("\nExample gene annotations from the input:")

        shown = 0

        for record in parsed_records:

            if record["feature"] != "gene":
                continue

            attrs = record["attrs"]

            print(
                "  gene_id=",
                attrs.get("gene_id"),
                "gene_name=",
                attrs.get("gene_name"),
                "id=",
                attrs.get("id"),
                "name=",
                attrs.get("name"),
            )

            shown += 1

            if shown >= 5:
                break

        return False, []

    print("Matched gene IDs:", sorted(target_gene_ids))
    print("Matched gene names:", sorted(target_gene_names))

    keep_ids = set(target_gene_ids)

    keep_transcript_ids = set()

    changed = True

    while changed:

        changed = False

        for record in parsed_records:

            attrs = record["attrs"]

            # ------------------------------------------
            # GTF
            # ------------------------------------------

            gene_ids = {
                str(x).strip().lower()
                for x in _as_list(attrs.get("gene_id"))
                if x is not None
            }

            transcript_ids = {
                str(x).strip().lower()
                for x in _as_list(attrs.get("transcript_id"))
                if x is not None
            }

            # If this annotation belongs directly to one
            # of our target genes, retain its transcript ID.
            if gene_ids & target_gene_ids:

                for tx_id in transcript_ids:
                    if tx_id not in keep_transcript_ids:
                        keep_transcript_ids.add(tx_id)
                        changed = True

            # ------------------------------------------
            # GFF3
            # ------------------------------------------

            parents = {
                str(x).strip().lower()
                for x in _as_list(attrs.get("parent"))
                if x is not None
            }

            ids = {
                str(x).strip().lower()
                for x in _as_list(attrs.get("id"))
                if x is not None
            }

            # A record whose parent is already retained
            # should also be retained.
            if parents & keep_ids:

                for obj_id in ids:

                    if obj_id not in keep_ids:
                        keep_ids.add(obj_id)
                        changed = True

                    # If this is a transcript, retain it
                    # explicitly as well.
                    if record["feature"] in (
                        "transcript",
                        "mrna",
                        "rna",
                    ):
                        keep_transcript_ids.add(obj_id)


    out_lines = []

    for record in parsed_records:

        attrs = record["attrs"]

        # ----------------------------------------------
        # GTF relationship
        # ----------------------------------------------

        gene_ids = {
            str(x).strip().lower()
            for x in _as_list(attrs.get("gene_id"))
            if x is not None
        }

        transcript_ids = {
            str(x).strip().lower()
            for x in _as_list(attrs.get("transcript_id"))
            if x is not None
        }

        gtf_match = (
            bool(gene_ids & target_gene_ids)
            or bool(transcript_ids & keep_transcript_ids)
        )

        # ----------------------------------------------
        # GFF3 relationship
        # ----------------------------------------------

        ids = {
            str(x).strip().lower()
            for x in _as_list(attrs.get("id"))
            if x is not None
        }

        parents = {
            str(x).strip().lower()
            for x in _as_list(attrs.get("parent"))
            if x is not None
        }

        gff_match = (
            bool(ids & keep_ids)
            or bool(parents & keep_ids)
        )

        if gtf_match or gff_match:
            out_lines.append(record["line"])

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
