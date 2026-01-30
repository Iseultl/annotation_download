#!/usr/bin/env python3
import pandas as pd

# Load the TSVs
ann = pd.read_csv("annotations_20260115_090352.tsv", sep="\t")
sephs2 = pd.read_csv("SEPHS2_locations_mapped.tsv", sep="\t")

# Get unique accessions present in SEPHS2 file
sephs2_accessions = set(sephs2["assembly_accession"].dropna())

# Filter annotations to keep only missing ones
ann_missing = ann[~ann["assembly_accession"].isin(sephs2_accessions)]

# Write output
ann_missing.to_csv(
    "annotations_missing_SEPHS2.tsv",
    sep="\t",
    index=False
)