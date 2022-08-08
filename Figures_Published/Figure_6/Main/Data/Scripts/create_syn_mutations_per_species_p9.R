################################################################################
# script to create syn point mutations per species data set for passage 9
################################################################################

# load libraries 
library(tidyverse)
library(hash)

# source helper functions 
source('./helpers.R')

# load data 
opt_rscu_data = load_opt_rscu_data()
mutation_data = load_p9_mutation_data()

# prepare the mutation data and merge with opt and rscu data 
syn_mutations_per_species = mutation_data %>%
  filter(muttype == 'Syn') %>%
  select(host, wtcodon, mutcodon, wtRes, mutRes, wrel, wrel.ciLower, wrel.ciUpper) %>%
  group_by(host, wtcodon, mutcodon, wtRes, mutRes) %>%
  summarize(
    mean_wrel = mean(wrel), 
    sd_wrel = sd(wrel), 
    mean_wrel_lower = mean(wrel.ciLower), 
    mean_wrel_upper = mean(wrel.ciUpper),
    n = n()
  ) %>%
  slice(1) %>%
  ungroup() %>%
  full_join(opt_rscu_data, by = c('host', 'wtcodon')) %>%
  select(-mutcodon.y) %>%
  rename(
    wtcodon_csc = csc, 
    wtcodon_rscu_fc = rscu_fc, 
    mutcodon = mutcodon.x
  ) %>%
  inner_join(opt_rscu_data, by = c('host', 'mutcodon')) %>%
  select(-wtcodon.y) %>%
  rename(
    mutcodon_csc = csc, 
    mutcodon_rscu_fc = rscu_fc, 
    wtcodon = wtcodon.x
  ) %>%
  mutate(
    delta_csc = .$mutcodon_csc - .$wtcodon_csc, 
    delta_rscu_fc = .$mutcodon_rscu_fc - .$wtcodon_rscu_fc, 
    fitness_class = assign_fitness_class(.)
  ) %>%
  rename(
    wt_aa = wtRes, 
    mut_aa = mutRes
  )

# use hash to convert single letter amino acid to three letter abbreviation
syn_mutations_per_species$wt_aa = vapply(syn_mutations_per_species$wt_aa, function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)
syn_mutations_per_species$mut_aa = vapply(syn_mutations_per_species$mut_aa, function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)

# save 
write.table(
  syn_mutations_per_species, 
  '../syn_mutations_per_species_p9.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################













