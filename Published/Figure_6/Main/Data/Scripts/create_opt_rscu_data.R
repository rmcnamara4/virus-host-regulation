################################################################################
# script to create RSCU fold change and codon optimality table for human, 
# mosquito, and Dengue 2 strain 16681 
################################################################################

# load libraries
library(tidyverse)
library(reshape2)

# source helper functions
source('./helpers.R')

# load codon optimality and RSCU fold change data 
optimality_data = load_codon_optimality()
rscu_fc_data = load_rscu_fc()

# prepare the codon optimality data for merging 
optimality_data = optimality_data %>%
  mutate(
    mutcodon = .$codon, 
    .after = 1
  ) %>%
  rename(
    wtcodon = codon, 
    Human = human_csc, 
    Mosquito = mosquito_csc
  ) %>%
  select(wtcodon, mutcodon, Human, Mosquito) %>%
  melt(id = c('wtcodon', 'mutcodon')) %>%
  rename(
    host = variable, 
    csc = value
  )

# prepare the RSCU fold change data for merging 
rscu_fc_data = rscu_fc_data %>%
  mutate(
    mutcodon = .$codon, 
    .after = 1
  ) %>%
  rename(
    wtcodon = codon, 
    Human = human_rscu_fc, 
    Mosquito = mosquito_rscu_fc
  ) %>%
  select(wtcodon, mutcodon, Human, Mosquito) %>%
  filter(!(wtcodon %in% c('TAG', 'TGA', 'TAA'))) %>%
  melt(id = c('wtcodon', 'mutcodon')) %>%
  rename(
    host = variable, 
    rscu_fc = value
  )

# combine the two data.frames
opt_rscu_data = full_join(optimality_data, rscu_fc_data, by = c('wtcodon', 'mutcodon', 'host'))

# save 
write.table(
  opt_rscu_data, 
  '../opt_rscu_data.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################