#!/usr/bin/env python3

import os
import csv
import gzip
import requests
import time
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import json

# --------------------------------------------------
# Download genome file
# --------------------------------------------------
def get_genome_fasta_url(assembly_accession):
    try:
        return ASSEMBLY_URLS[assembly_accession]
    except KeyError:
        raise RuntimeError(f"No URL for {assembly_accession}")

def download_file(url, outfile, retries=5, chunk_size=1024*1024):
    tmp = outfile + ".part"

    headers = {}
    pos = 0

    if os.path.exists(tmp):
        pos = os.path.getsize(tmp)
        headers["Range"] = f"bytes={pos}-"
        print(f"[RESUME] {outfile} from byte {pos}")

    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, headers=headers, timeout=60) as r:
                r.raise_for_status()
                mode = "ab" if pos else "wb"
                with open(tmp, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)

            os.rename(tmp, outfile)
            print(f"[DONE] {outfile}")
            return outfile

        except Exception as e:
            print(f"[RETRY {attempt}/{retries}] {e}")

    raise RuntimeError(f"Failed downloading after {retries} attempts: {url}")


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def extract_gene_name(attributes):
    for field in attributes.split(";"):
        if field.startswith("Name="):
            return field.replace("Name=", "")
    return "UNKNOWN"

def gene_output_path(row):
    species = row["Species_dir"].strip()
    seqid = row["Mapped_seqid"].strip()

    if seqid == "NOT_FOUND":
        seqid = row["seqid"].strip()

    start = int(row["start"])
    end = int(row["end"])
    attributes = row.get("attributes", "")

    start = max(1, start - 200)
    end = end + 200

    gene_name = extract_gene_name(attributes)

    species_dir = os.path.join(GENE_DIR, species)

    return os.path.join(
        species_dir,
        f"{gene_name}_{seqid}_{start}_{end}.fa"
    )
# --------------------------------------------------
# Core logic
# --------------------------------------------------
def extract_gene(genome_fasta_gz, row):
    species = row["Species_dir"].strip()
    seqid = row["Mapped_seqid"].strip()
    start = int(row["start"])
    end = int(row["end"])
    strand = row["strand"].strip()
    attributes = row.get("attributes", "")

    if seqid == "NOT_FOUND":
        seqid = row["seqid"].strip()

    assert end > start

    # Expand region
    start = max(1, start - 200)
    end = end + 200

    gene_name = extract_gene_name(attributes)

    # Load genome
    with gzip.open(genome_fasta_gz, "rt") as handle:
        genome = SeqIO.to_dict(
            SeqIO.parse(handle, "fasta"),
            key_function=lambda rec: rec.id.split()[0]
        )
        print(f"[DEBUG] Genome keys: {list(genome.keys())[:10]}")
    # Find chromosome
    if seqid not in genome:
        print(f"[WARN] {seqid} not found in genome: {genome_fasta_gz}")
        return

    seq = genome[seqid].seq[start - 1:end]

    if strand == "-":
        seq = seq.reverse_complement()

    species_dir = os.path.join(GENE_DIR, species)
    os.makedirs(species_dir, exist_ok=True)

    outfile = os.path.join(
        species_dir,
        f"{gene_name}_{seqid}_{start}_{end}.fa"
    )

    with open(outfile, "w") as out:
        out.write(f">{gene_name}|{seqid}:{start}-{end}({strand})\n")
        out.write(str(seq) + "\n")

    print(f"[WRITE] {outfile}")
    return outfile

def process_tsv(tsv_file, start=None, end=None):
    with open(tsv_file) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for i, row in enumerate(reader):
            if start is not None and i < start:
                continue
            if end is not None and i >= end:
                break
            
            output_path = gene_output_path(row)
            if os.path.exists(output_path):
                print(f"[SKIP] {output_path} already exists")
                continue
            
            assembly = row["Assembly_accession"].strip()

            genome_gz = os.path.join(
                GENOME_DIR,
                f"{assembly}.genome.fa.gz"
            )

            print(f"[ASSEMBLY] {assembly}")
            
            try:
                genome_url = get_genome_fasta_url(assembly)
                download_file(genome_url, genome_gz)
                extract_gene(genome_gz, row)

            finally:
                # ALWAYS clean up, even if extract_gene fails
                if os.path.exists(genome_gz):
                    os.remove(genome_gz)
                    print(f"[CLEANUP] removed {genome_gz}")

# --------------------------------------------------
# Main
# --------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract genes from genomes")
    parser.add_argument("--tsv_file", help="TSV file with gene annotation information")
    parser.add_argument("--assembly_urls", help="JSON file with assembly URLs")
    parser.add_argument("--start", type=int, help="Start index")
    parser.add_argument("--end", type=int, help="End index")
    parser.add_argument("--outdir", help="Output directory")
    
    args = parser.parse_args()
    
    ANNOTRIEVE_API = "https://genome.crg.es/annotrieve/api/v0"
    GENOME_DIR = args.outdir + "/genomes"
    GENE_DIR = args.outdir + "/genes"

    os.makedirs(GENOME_DIR, exist_ok=True)
    os.makedirs(GENE_DIR, exist_ok=True)
    
    with open(args.assembly_urls) as f:
        ASSEMBLY_URLS = json.load(f)

    process_tsv(args.tsv_file, args.start, args.end)

# Command for running script 
"""
python get_genes.py --tsv_file data/genes.tsv --assembly_urls data/assembly_urls.json --outdir data --start 0 --end 10
"""

