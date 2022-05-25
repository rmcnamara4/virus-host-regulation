#!/bin/bash

# script to run all of the analysis scripts

# run preprocessing steps
bash ./scripts/download_metadata.sh
Rscript ./scripts/create_virus_composition.R
Rscript ./scripts/create_host_metadata.R
Rscript ./scripts/clean_assembly.R
bash ./scripts/get_cds_fastas.sh
bash ./scripts/change_file_names.sh

# run snakemake analysis
snakemake --cores 25 concatenate_codon_counts
snakemake --cores 25 concatenate_rscu_tables
snakemake --cores 25 calculate_rscu_ratios
snakemake --cores 25 calculate_correlations
