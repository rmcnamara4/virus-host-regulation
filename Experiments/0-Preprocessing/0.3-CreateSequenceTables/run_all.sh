#!/bin/bash

# script to run R scripts for creating sequence tables
# convert git lfs pointers to original data 
git lfs pull

# make dir to put the tables
mkdir -p data

# run script for making human and mosquito sequence tables
Rscript make_human_mosquito_seq_tables.R

# run script for making dengue isolates sequence tables
Rscript make_isolates_seq_tables.R
