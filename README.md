This repository contains scripts to download the genes listed in the genes.txt file. Thus it requires that these genes are well annotated in the gff file. We extract the transcript sequence +2000 nucleotides upstream and downstream to ensure full extraction of the UTR regions. 

To download transcripts via download_genes.sh provide a list of taxids and gene names in the provided text files: taxids.txt and genes.txt

Download genes steps
1. The links for the gff and genome files are created by running fetch_annotrieve_urls.py. A single annotation and genome file is chosen for each taxid based on the best overall BUSCO scores.
2. Next run the download_genes.sh script via SLURM with run_download_genes.slurm, edit the script to provide an output directory. 

NOTE: download_genes.sh uses local locations for gffread and python containers. Annocli must be available on the command line. 