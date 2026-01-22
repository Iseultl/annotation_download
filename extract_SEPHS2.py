import gzip
from pathlib import Path
import csv
import re
import zlib
import re

# =====================
# Paths
# =====================
BASE_DIR = Path("/no_backup/rg/ileahy/mammals")
OUTPUT_FILE = Path("/no_backup/rg/ileahy/mammals/SEPHS2_locations.tsv")

# =====================
# Output columns
# =====================
# GFF columns: seqid, source, type, start, end, score, strand, phase, attributes
OUTPUT_HEADER = [
    "seqid", "source", "type", "start", "end", "score",
    "strand", "phase", "attributes",
    "Annotation_ID", "Assembly_accession", "Database",
    "Species_dir"
]

# =====================
# Parse Attributes
# =====================
def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key.lower()] = val
    return attrs


# =====================
# Function to parse GFF
# =====================
def parse_gff_for_gene(gff_file, gene_aliases=None):
    """
    Returns the first row containing the gene_name in the attributes column.
    
    GENE_KEYS = ("Name", "gene", "gene_name", "Alias")

    if gene_aliases is None:
        gene_aliases = {
            "sephs2"
        }
    
    alias_pattern = "|".join(re.escape(a) for a in gene_aliases)

    pattern = re.compile(
        rf"(?:^|;)(?:{'|'.join(GENE_KEYS)})=({alias_pattern})(?:;|$)",
        re.IGNORECASE
    )
    """
    try:
        open_func = gzip.open if gff_file.suffix == ".gz" else open
        with open_func(gff_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
            
                if ("=SEPHS2" in fields[8]) or ("=sephs2" in fields[8]) or ("=SPS2" in fields[8]):
                    print(f"Found SEPHS2 in {gff_file}")
                    return fields  
        return None  # not found
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as e:
        print(f"[CORRUPT FILE] Skipping {gff_file}: {e}")
        return None

# =====================
# Iterate species directories
# =====================
rows_to_write = []

for species_dir in BASE_DIR.iterdir():
    if not species_dir.is_dir():
        continue

    # Extract metadata from README.txt
    readme_file = species_dir / "README.txt"
    annotation_id = assembly_acc = database = ""
    if readme_file.exists():
        with open(readme_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith("Annotation ID"):
                    annotation_id = line.split(":", 1)[1].strip()
                elif line.startswith("Assembly accession"):
                    assembly_acc = line.split(":", 1)[1].strip()
                elif line.startswith("Database"):
                    database = line.split(":", 1)[1].strip()

    # Search GFF files in directory
    gff_files = list(species_dir.glob("*.gff*"))
    sep_row = None
    for gff_file in gff_files:
        sep_row = parse_gff_for_gene(gff_file, "SEPHS2")
        if sep_row:
            break  # stop at first file with SEPHS2

    if sep_row:
        # Add metadata columns
        sep_row.extend([annotation_id, assembly_acc, database, species_dir.name])
        rows_to_write.append(sep_row)
    else:
        # No SEPHS2 found: create empty row with only final column
        empty_row = [""] * 12 + [species_dir.name]
        rows_to_write.append(empty_row)

# =====================
# Write output
# =====================
with open(OUTPUT_FILE, "w", newline="") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(OUTPUT_HEADER)
    writer.writerows(rows_to_write)

print(f"Done! SEPHS2 locations written to {OUTPUT_FILE}")
