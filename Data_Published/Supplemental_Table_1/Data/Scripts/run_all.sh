#!/bin/bash

# run bash script to collect data
bash ./collect_data.sh

# run R script to create supplemental table 1
Rscript ./create_supplemental_table_1.R
