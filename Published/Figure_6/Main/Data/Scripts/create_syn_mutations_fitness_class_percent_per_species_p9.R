################################################################################
# script to create syn mutations fitness class percent per species data set for 
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
syn_mutations_fitness_class_percent_per_species = mutation_data %>%
  filter(muttype == 'Syn') %>%
  select(host, wtcodon, mutcodon, wtRes, mutRes, status) %>%
  group_by(host, wtcodon, mutcodon, wtRes, mutRes, status) %>%
  summarize(
    n = n() 
  ) %>%
  ungroup() %>%
  group_by(host, wtcodon, mutcodon, wtRes, mutRes) %>%
  summarize(
    status = status, 
    n = n, 
    tot = sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    perc = 100 * (n / tot) 
  ) %>%
  rename(
    wt_aa = wtRes,
    mut_aa = mutRes, 
    fitness_class = status
  ) %>%
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
    delta_rscu_fc = .$mutcodon_rscu_fc - .$wtcodon_rscu_fc
  )

# convert amino acid letter to three letter abbreviation 
syn_mutations_fitness_class_percent_per_species$wt_aa = vapply(syn_mutations_fitness_class_percent_per_species$wt_aa, 
                                                                function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)
syn_mutations_fitness_class_percent_per_species$mut_aa = vapply(syn_mutations_fitness_class_percent_per_species$mut_aa, 
                                                                function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)

# change fitness class letter to full name 
syn_mutations_fitness_class_percent_per_species$fitness_class = ifelse(syn_mutations_fitness_class_percent_per_species$fitness_class == 'B', 'beneficial', 
                                                                       ifelse(syn_mutations_fitness_class_percent_per_species$fitness_class == 'N', 'neutral', 
                                                                              ifelse(syn_mutations_fitness_class_percent_per_species$fitness_class == 'D', 'deleterious', 'lethal')))

# save
write.table(
  syn_mutations_fitness_class_percent_per_species, 
  '../syn_mutations_fitness_class_percent_per_species_p9.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)
  
################################################################################
  
  
  
  