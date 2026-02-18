Update the script so that it can download all the genes given a gene name and taxonomic group.

This repository contains scripts to download mammalian genomes with SEPHS2 genes.
Step 1: Download annotations with async_gff_downloader.py <- run via annotation_download.sh
Step 2: Extract SEPHS2 gene annotations with extract_SEPHS2.py -> SEPHS2_locations.tsv
Step 3: Ensure the chromosomes are in the correct format with sequence_mapper.py
Step 4: Get genome online locations with fetch_annotrieve_urls.py -> assembly_urls.json
Step 5: Download genomes and extract gene of interest with get_genes.py -> genome asssemblies are automatically removed after run - run with get_genes.sh script