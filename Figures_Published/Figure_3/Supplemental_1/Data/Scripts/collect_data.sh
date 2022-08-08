#!/bin/bash

# script to collect data from the Main data folder
# collect codon optimality data
cp ../../../Main/Data/codon_optimalities.csv ..

# collect viral RSCU fold change mosquito table
cp ../../../Main/Data/viral_rscu_fc_mosquito.csv ..

# collect viral RSCU fold change human table
cp ../../../Main/Data/viral_rscu_fc_human.csv ..

# collect Dengue isolates RSCU fc mosquito table
cp ../../../Main/Data/dengue_isolates_m-relative_rscu_fc.csv ..
