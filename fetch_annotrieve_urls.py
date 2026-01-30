# fetch_annotrieve_urls.py
# The purpose of this script is to get the download URLs for assembly genomes
# listed in the SEPHS2_locations_mapped.tsv file.
import csv
import requests
import json

ANNOTRIEVE_API = "https://genome.crg.es/annotrieve/api/v0"

urls = {}

with open("/Users/iseult/gitlab/annotations_download/data/missing_SEPHS2_locations_mapped.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        asm = row["Assembly_accession"].strip()
        species = row["Species_dir"].strip()
        if asm in urls:
            continue

        r = requests.get(
            f"{ANNOTRIEVE_API}/assemblies/{asm}",
            headers={"Accept": "application/json"},
            timeout=20
        )
        r.raise_for_status()
        urls[asm] = r.json()["download_url"]

with open("/Users/iseult/gitlab/annotations_download/data/missing_assembly_urls.json", "w") as out:
    json.dump(urls, out, indent=2)
