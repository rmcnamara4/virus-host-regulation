#!/bin/bash 

# script to run all of the data creation/collection 
# run collection script 

bash collect_data.sh
Rscript create_viral_aa_fc_human.R
