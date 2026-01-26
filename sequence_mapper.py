import requests
import csv
import re
import argparse

ANNOTRIEVE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
"""
def sequence_mapper(target_sequence, assembly_accession):
    """
#    This function helps to map the seqid id to all the possible sequences names by using the FTP server of the NCBI (assembly_report.txt)

#    Args:
#        target_sequence: the sequence id to map
#        assembly_accession: the assembly accession to use

#    Returns:
#        a list of possible sequence names


def sequence_mapper(target_sequence, assembly_accession):
    """
    Map TSV seqid (e.g. '15') to RefSeq/GenBank accessions
    using NCBI assembly_report.txt
    """

    base = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

    acc_path = "/".join([
        assembly_accession[0:3],
        assembly_accession[4:7],
        assembly_accession[7:10],
        assembly_accession[10:13]
    ])

    assembly_dir_url = f"{base}/{acc_path}/"

    # List directory contents
    r = requests.get(assembly_dir_url)
    if r.status_code != 200:
        raise RuntimeError(f"Cannot access FTP directory for {assembly_accession}")

    # Find the assembly directory
    match = re.search(rf'({assembly_accession}_[^"/]+)/', r.text)
    if not match:
        raise RuntimeError(f"Assembly directory not found for {assembly_accession}")

    assembly_dir = match.group(1)

    report_url = (
        f"{base}/{acc_path}/{assembly_dir}/"
        f"{assembly_dir}_assembly_report.txt"
    )

    report = requests.get(report_url)
    
    if report.status_code != 200:
        raise RuntimeError(f"Cannot download assembly report for {assembly_accession}")

    matches = []

    for line in report.text.splitlines():
        
        if line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 7:
            continue

        (
            sequence_name,
            sequence_role,
            assigned_molecule,
            molecule_type,
            genbank_acc,
            relationship,
            refseq_acc,
            *_
        ) = fields

        if (sequence_name == target_sequence) or (assigned_molecule == target_sequence) or (genbank_acc == target_sequence):
            matches.append({
                "sequence_name": sequence_name,
                "genbank": None if genbank_acc == "na" else genbank_acc,
                "refseq": None if refseq_acc == "na" else refseq_acc,
            })

    return matches
        
def run_mapper(tsv_file, output_file):
    with open(tsv_file) as fin, open(output_file, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames + ["Mapped_seqid"] + ['sequence_name'], delimiter="\t")
        writer.writeheader()

        for row in reader:
            # Skip empty rows
            if not row["seqid"].strip():
                continue

            mappings = sequence_mapper(
                target_sequence=row["seqid"],
                assembly_accession=row["Assembly_accession"]
            )
            if not mappings:
                row["Mapped_seqid"] = "NOT_FOUND"
                row["sequence_name"] = ""
            else:
                # Prefer RefSeq accession, fallback to GenBank
                best = mappings[0]
                row["Mapped_seqid"] = best["refseq"] or best["genbank"] 
                row["sequence_name"] = best["sequence_name"]

            writer.writerow(row)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input TSV file", required=True)
    parser.add_argument("--output", help="Output TSV file", required=True)
    args = parser.parse_args()
    
    run_mapper(args.input, args.output)
    
# Code for running the script 
"""
python annotation_download/sequence_mapper.py --input data/SEPHS2_locations.tsv --output data/SEPHS2_locations_mapped.tsv
"""

