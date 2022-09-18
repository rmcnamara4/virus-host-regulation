#!/bin/bash

# script to collect required data
# collect synonymous mrf data by species
cp ../../../../Experiments/5-MutationAnalysis/data/syn_mut_mrf_per_species_p9.csv ..

# collect RSCU fold change data
cp ../../../../Experiments/5-MutationAnalysis/data/rscu_fc.csv ..

# collect optimality data
cp ../../../../Experiments/2-ConcatenateOptimalities/data/codon_optimalities.csv ..
