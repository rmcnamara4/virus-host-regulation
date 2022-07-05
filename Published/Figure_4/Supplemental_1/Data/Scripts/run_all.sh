#!/bin/bash

# script to run all of the data creation/collection scripts
# run collect_data.sh script
bash collect_data.sh

# run create_dengue_strain_16681_codon_comp.R
Rscript create_dengue_strain_16681_codon_comp.R
