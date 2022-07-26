################################################################################
# script to create MRF of synonymous mutations to each codon per species for 
# passage 9
################################################################################

# load libraries
library(tidyverse)
library(hash)

# source helper functions
source('./helpers.R')

# load data 
opt_rscu_data = load_opt_rscu_data()
mutation_data = load_p9_mutation_data()

# prepare mutation data and merge with optimality and RSCU fold change data 
syn_mutations_mrf_per_species = mutation_data %>%
  filter(muttype == 'Syn') %>%
  select(host, mutcodon, mutRes, wrel, wrel.ciLower, wrel.ciUpper) %>%
  group_by(host, mutcodon, mutRes) %>%
  summarize(
    mean_wrel = mean(wrel), 
    mean_wrel_lower = mean(wrel.ciLower), 
    mean_wrel_upper = mean(wrel.ciUpper)
  ) %>%
  ungroup() %>%
  inner_join(opt_rscu_data, by = c('host', 'mutcodon')) %>%
  select(-wtcodon) %>%
  rename(
    mutcodon_csc = csc, 
    mutcodon_rscu_fc = rscu_fc
  ) %>%
  mutate(
    fitness_class = assign_fitness_class(.)
  ) %>%
  rename(
    mut_aa = mutRes
  )

# convert amino acid letter to three letter abbreviation 
syn_mutations_mrf_per_species$mut_aa = vapply(syn_mutations_mrf_per_species$mut_aa, function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)

# save
write.table(
  syn_mutations_mrf_per_species, 
  '../syn_mutations_mrf_per_species_p9.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################