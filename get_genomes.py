#!/usr/bin/env python3

import csv
import os
import requests

ANNOTRIEVE_API = "https://genome.crg.es/annotrieve/api/v0"
OUTDIR = "genomes"

os.makedirs(OUTDIR, exist_ok=True)

def get_genome_fasta_url(assembly_accession):
    """Query Annotrieve for genome FASTA URL"""
    url = f"{ANNOTRIEVE_API}/assemblies/{assembly_accession}"
    print(url)
    r = requests.get(url, headers={"Accept": "application/json"})
    print(r)
    if r.status_code != 200:
        raise RuntimeError(f"Annotrieve error for {assembly_accession}")

    data = r.json()

    return data["download_url"]

def download_file(url, outfile):
    if os.path.exists(outfile):
        print(f"[SKIP] {outfile}")
        return

    print(f"[DOWNLOAD] {outfile}")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(outfile, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

def main(tsv_file):
    assemblies = set()

    # Collect unique assemblies
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row["Assembly_accession"].strip()
            if acc:
                assemblies.add(acc)

    print(f"Found {len(assemblies)} unique assemblies")

    # Download genomes
    for acc in sorted(assemblies):
        try:
            fasta_url = get_genome_fasta_url(acc)
            outfile = os.path.join(OUTDIR, f"{acc}.genome.fa.gz")
            download_file(fasta_url, outfile)
        except Exception as e:
            print(f"[ERROR] {acc}: {e}")

if __name__ == "__main__":
    main("/Users/iseult/gitlab/annotations_download/data/test_locations_mapped.tsv")
