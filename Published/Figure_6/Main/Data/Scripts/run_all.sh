#!/bin/bash

# script to run all data collection/creation scripts
# run data collection script
bash collect_data.sh

# run script to create optimality and RSCU fold change table
Rscript create_opt_rscu_data.R

# run script to create condensed point mutations per species p9 table
Rscript create_syn_mutations_per_species_p9.R

# run script to create syn mutations MRF per rep p9 table
Rscript create_syn_mutations_mrf_per_rep_p9.R

# run script to create syn mutations MRF per species p9 table
Rscript create_syn_mutations_mrf_per_species_p9.R

# run script to create syn mutations fitness class percent per species p9 table
Rscript create_syn_mutations_fitness_class_percent_per_species_p9.R

# run script to create syn mutations MRF codon groups per species p9 table
Rscript create_syn_mutations_mrf_codon_groups_per_species_p9.R
