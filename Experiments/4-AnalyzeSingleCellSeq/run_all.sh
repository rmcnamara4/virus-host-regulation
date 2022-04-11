#!/bin/bash

# script to run all of the analysis
# run script to download data
bash download_data.sh

# run script to calculate CPMs
mkdir -p data/cpms
Rscript calculate_cpms.R
