#!/bin/bash

# script to run all of the analysis scripts

# run preprocessing steps
bash download_metadata.sh
Rscript create_host_metadata.R
Rscript clean_assembly.R
bash get_cds_fastas.sh
bash change_file_names.sh

# run snakemake analysis
snakemake --cores 25 all
