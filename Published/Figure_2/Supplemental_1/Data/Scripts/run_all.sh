#!/bin/bash

# script to run all scripts for data collection
# run script to gather optimality data
bash gather_optimality_data.sh

# run script to create mosquito codon frequency ratio table
Rscript create_mosquito_codon_frequency_ratio.R

# run script to create EGFP reporters codon frequency ratio table
Rscript create_egfp_reporters_codon_frequency_ratio.R

# run script to create 1nt out of frame reporters codon frequency ratio table
Rscript create_nt_out_of_frame_reporters_codon_frequency_ratio.R
