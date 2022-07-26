################################################################################
# script to create syn MRF codon groups per species for passage 9
################################################################################

# load libraries
library(tidyverse)
library(reshape2)

# source helper functions
source('./helpers.R')

# load data 
codon_groups = load_codon_groups()
syn_mutations_mrf_per_species = read.table(
  '../syn_mutations_mrf_per_species_p9.csv', 
  sep = ',', header = TRUE
)

# add host column to codon groups df 
codon_groups$Human$host = 'Human'
codon_groups$Mosquito$host = 'Mosquito'

# rbind to one df
codon_groups = Reduce('rbind', codon_groups)

# melted df
melted_codon_groups = melt(codon_groups, id = 'host')

# add a group column to link the pairs of codon groups
melted_codon_groups$group = ifelse(melted_codon_groups$variable %in% c('frequently used', 'infrequently used'), 1, 
                                   ifelse(melted_codon_groups$variable %in% c('preferentially used', 'non_preferentially_used'), 2, 3))

# rename columns
names(melted_codon_groups) = c('host', 'codon_class', 'mutcodon', 'group')

# merge df with the MRF values 
syn_mrf_codon_groups_per_species = merge(melted_codon_groups, syn_mutations_mrf_per_species, by = c('host', 'mutcodon'))

# save
write.table(
  syn_mrf_codon_groups_per_species, 
  '../syn_mutations_mrf_codon_groups_per_species_p9.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################