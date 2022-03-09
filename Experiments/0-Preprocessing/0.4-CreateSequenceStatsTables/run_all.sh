#!/bin/bash

# script to run R scripts for creating sequence stats tables
# make dir to put the tables
mkdir -p data

# run script to create sequence stats tables
Rscript make_sequence_stats_tables.R
