#!/bin/bash

# script to run all of the data collection/creation scripts
# run data collection script
bash collect_data.sh

# run script to create viral amino acid fold change mosquito table
Rscript create_viral_aa_fc_mosquito.R
