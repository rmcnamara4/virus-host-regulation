#!/bin/bash

# script to run all data collection scripts and script to create supplemental table 1
# run bash script to collect data
bash ./collect_data.sh

# run R script to create supplemental table 1
Rscript ./create_supplemental_table_1.R
