#!/bin/bash 

# script to run all of the data creation/collection for Figure 1 (main)

Rscript create_viral_rscu_fc_human.R
bash collect_human_opt.sh
Rscript create_dengue_isolates_h-relative_rscu_fc.R
Rscript create_dengue_2_iso_and_human_codon_freq.R
