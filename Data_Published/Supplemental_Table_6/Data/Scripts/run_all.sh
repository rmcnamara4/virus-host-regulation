#!/bin/bash

# script to run all data collection scripts and create supplemental table 6
# run bash script to collect required data
bash ./collect_data.sh

# run Rscript to create supplemental table 6
Rscript ./create_supplemental_table_6.R
