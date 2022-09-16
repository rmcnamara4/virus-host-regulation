#!/bin/bash

# script to run all data collection scripts and script to create supplemental table 3
# run bash script to collect data
bash ./collect_data.sh

# run Rscript to create supplemental table 3
Rscript create_supplemental_table_3.R
