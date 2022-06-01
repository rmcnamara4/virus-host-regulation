#!/bin/bash

# script to run all scripts for data collection 
# run script to gather optimality data 
bash gather_optimality_data.sh

# run script to create mosquito codon frequency ratio table 
Rscript create_mosquito_codon_frequency_ratio.R
