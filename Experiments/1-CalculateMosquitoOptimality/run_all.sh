#!/bin/bash 

# script to run all of the bash scripts in order 
# run script to collect required files 
bash get_files.sh

# run script to split gff file 
bash split_gff.sh

# run script to map the introns 
bash count_introns.sh
