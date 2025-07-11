#!/bin/bash

# script to run all of the analysis
# run bash script to download data
bash download_data.sh

# run R script to extract passage 9 data
Rscript extract_passage_9.R

# run R script to create the RSCU fold change data.frame
Rscript create_rscu_fc_df.R

# run R script to create the condensed mutation data.frame
Rscript create_mutation_dfs.R

# create a folder for the figures
mkdir -p figures
