#!/bin/bash

# script to run data collection and create of supplemental table 2
# run bash script to collect data
bash ./collect_data.sh

# run Rscript to create supplemental table 2
Rscript ./create_supplemental_table_2.R
