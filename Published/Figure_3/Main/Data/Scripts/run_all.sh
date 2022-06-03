#!/bin/bash

# script to run all of the data creation scripts
# run script to collect optimality data
bash collect_opt.sh

# run script to create viral RSCU fold change mosquito table
Rscript create_viral_rscu_fc_mosquito.R
