# 3-DownloadAllViruses

Here I download the CDS fastas of all of the available complete viral genomes through RefSeq. I then calculate the codon frequency and RSCU of each transcript and total genome. Then, I calculate the log2 RSCU ratio between viral transcripts and genome and the human/mosquito genomes. Lastly, I calculate the Spearman correlation between the ratios and the codon optimality of human/mosquito, producing a table that can be used for future analysis.

## TXT

+ **./metadata.txt**: data downloaded from [here](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download) that contains virus names, realm, host, etc.

**./host_metadata.txt**: cleaned version of metadata.txt that contains only the virus name, realm, and host

**./complete_assembly.txt**: data downloaded from [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt) that contains virus names and urls to download from RefSeq

**./urls.txt**: list of all of the urls needed to download the viral CDS fastas from RefSeq

## Snakemake

**./Snakefile**: file defining all of the rules of the analysis using Snakemake

## Scripts

**./scripts/download_metadata.sh**: Bash script to download the metadata.txt and complete_assembly.txt files

**./scripts/clean_assembly.R**: R script to clean the complete_assembly.txt file to use for analysis

**./scripts/create_host_metadata.R**: R script to create the host_metadata.txt file from the metadata.txt file

**./scripts/get_cds_fastas.sh**: Bash script to extract urls from the complete_assembly.txt file, save them to the urls.txt file, and use them to download the CDS fastas

**./scripts/change_file_names.sh**: Bash script to change the names of the fasta files to the names of the viruses

**./scripts/create_codon_count_tables.R**: R script to calculate codon frequency of each virus transcript and total genome

**./scripts/concatenate_codon_counts.R**: R script to concatenate the codon frequency tables of all viruses (both for individual transcripts and total genomes) into one file

**./scripts/create_rscu_tables.R**: R script to calculate the RSCU of each virus transcript and total genome

**./scripts/concatenate_rscu_tables.R**: R script to concatenate the RSCU tables of all viruses (both for individual transcripts and total genomes) into one file

**./scripts/calculate_rscu_ratios.R**: R script to calculate the log2 of the ratio between the RSCUs of the viral transcripts/genome and the RSCUs of the human/mosquito genomes

**./scripts/calculate_correlations.R**: R script to calculate the Spearman correlation between viral genome RSCU and human/mosquito genome RSCU

## Data

+ **./data/fastas/*.fna.gz**: all CDS fasta files of viral genomes from RefSeq

+ **./data/codon_tables/codon_counts/*.csv**: all codon frequency tables for viral transcripts

+ **./data/codon_tables/codon_counts_total/*.csv**: all codon frequency tables for viral genomes

+ **./data/codon_tables/all_codon_counts.csv**: table of all viral transcripts' codon frequency

+ **./data/codon_tables/all_codon_counts_total.csv**: table of all viral genomes' codon frequency

+ **./data/rscu_tables/rscu/*.csv**: all RSCU tables for viral transcripts

+ **./data/rscu_tables/rscu_total/*.csv**: all RSCU tables for viral genomes

+ **./data/rscu_tables/all_rscu.csv**: table of all viral transcripts' RSCU

+ **./data/rscu_tables/all_rscu_total.csv**: table of all viral genomes' RSCU

+ **./data/rscu_ratios/all_rscu_ratio_human.csv**: table of all viral transcripts' RSCU ratio relative to human genome

+ **./data/rscu_ratios/all_rscu_ratio_mosquito.csv**: table of all viral transcripts' RSCU ratio relative to mosquito genome

+ **./data/rscu_ratios/all_rscu_ratio_total_human.csv**: table of all viral genomes' RSCU ratio relative to human genome

+ **./data/rscu_ratios/all_rscu_ratio_total_mosquito.csv**: table of all viral genomes' RSCU ratio relative to mosquito genome

+ **./data/rscu_correlations/all_correlations_total_human.csv**: table of all viral genomes' Spearman correlation between RSCU ratio and human codon optimality

+ **./data/rscu_correlations/all_correlations_total_mosquito.csv**: table of all viral genomes' Spearman correlation between RSCU ratio and mosquito codon optimality 
