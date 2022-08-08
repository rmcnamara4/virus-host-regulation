#!/bin/bash

# script for collecting the required data
# collect RSCU fold change data for Dengue 2 strain 16681
cp ../../../../../Experiments/5-MutationAnalysis/data/rscu_fc.csv ..

# collect 9th passage point mutation fitness data
cp ../../../../../Experiments/5-MutationAnalysis/data/point_mutations_p9.csv ..

# collect human optimality
cp ../../../../../Experiments/2-ConcatenateOptimalities/data/codon_optimalities.csv ..

# collect codon groups
cp ../../../../../Experiments/5-MutationAnalysis/data/codon_groups.xlsx ..
