#!/bin/bash

# script to concatenate all of the codon optimalities we know
# create folder to put data in 
mkdir -p data 

# run concatenate_optimalities.R 
Rscript concatenate_optimalities.R
