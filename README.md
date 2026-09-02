This repository contains scripts to download mammalian genomes with SEPHS2 genes.

To download transcripts via download_genes.sh provide a list of taxids and gene names in the provided text files: taxids.txt and genes.txt

Download genes steps
1. The links for the gff and genome files are created by running fetch_annotrieve_urls.py
2. Next run the download_genes.sh script via SLURM with run_download_genes.slurm, edit the script to provide an output directory

Download genes scripts uses local locations for gffread and python containers. Annocli must be available on the command line. 