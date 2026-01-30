Download mammalian SEPHS2 genes

This repository contains scripts to download mammalian genomes with SEPHS2 genes.
Step 1: Download annotations with async_gff_downloader.py <- run via annotation_download.sh
Step 2: Extract SEPHS2 gene annotations with extract_SEPHS2.py -> SEPHS2_locations.tsv
Step 3: Ensure the chromosomes are in the correct format with sequence_mapper.py
Step 4: Download genomes with fetch_annotrieve_urls.py -> assembly_urls.json
Step 5: Download genomes and extract gene of interest with get_genes.py -> genome asssemblies are automatically removed after run