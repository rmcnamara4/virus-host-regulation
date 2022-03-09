#!/bin/bash

# script to run all of the collection and filtering of the human CSC data

# call get_human_csc.sh to download the CSC data 
bash get_human_csc.sh

# call gather_human_endo_csc.R script to filter for and average the endogenous CSC data 
Rscript gather_human_endo_csc.R 
