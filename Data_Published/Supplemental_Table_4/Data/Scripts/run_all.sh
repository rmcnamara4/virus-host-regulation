#!/bin/bash

# script to run all data collection scripts and create supplemental table 4
# run bash script to collect required data
bash ./collect_data.sh

# run Rscript to create supplemental table 4
Rscript create_supplemental_table_4.R
