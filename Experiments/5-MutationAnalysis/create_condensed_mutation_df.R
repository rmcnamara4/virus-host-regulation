################################################################################
# This script is to create the mutation data.frame that we need for all of the 
# plotting. 
################################################################################
# Load libraries 

library(dplyr)
library(reshape2)
library(hash)

################################################################################
# Source files 

source('./helper_functions.R')

################################################################################
# Load codon optimality, RSCU fold change, and mutation data 

optimality_data = load_codon_optimalities()
rscu_fc_data = load_rscu_fc()
mutation_data = load_mutation_data()

################################################################################
# Prepare the codon optimality data for merging 

optimality_data = optimality_data %>%
  mutate(
    mutcodon = .$codon, 
    .after = 1
  ) %>%
  select(codon, mutcodon, human_csc, mosquito_csc) %>%
  rename(
    wtcodon = codon, 
    Human = human_csc, 
    Mosquito = mosquito_csc
  ) %>%
  melt(id = c('wtcodon', 'mutcodon')) %>%
  rename(
    host = variable, 
    csc = value
  )

################################################################################
# Prepare the RSCU fold change data for merging

rscu_fc_data = rscu_fc_data %>%
  mutate(
    mutcodon = .$codon, 
    .after = 1
  ) %>%
  filter(!(codon %in% c('TAG', 'TGA', 'TAA'))) %>%
  select(codon, mutcodon, human_rscu_fc, mosquito_rscu_fc) %>%
  rename(
    wtcodon = codon, 
    Human = human_rscu_fc, 
    Mosquito = mosquito_rscu_fc
  ) %>%
  melt(id = c('wtcodon', 'mutcodon')) %>%
  rename(
    host = variable, 
    rscu_fc = value
  )

################################################################################
# Combine the two data.frames

opt_rscu_data = full_join(optimality_data, rscu_fc_data, by = c('wtcodon', 'mutcodon', 'host'))

################################################################################
# Prepare the mutation data and merge

mutation_data = mutation_data %>%
  select(host, wtcodon, mutcodon, wtRes, mutRes, wrel, wrel.ciLower, wrel.ciUpper, subTypes, muttype) %>%
  group_by(host, wtcodon, mutcodon, wtRes, mutRes) %>%
  summarize(
    mean_wrel = mean(wrel), 
    sd_wrel = sd(wrel), 
    mean_wrel_lower = mean(wrel.ciLower), 
    mean_wrel_upper = mean(wrel.ciUpper), 
    subTypes = subTypes, 
    muttype = muttype, 
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
  full_join(opt_rscu_data, by = c('host', 'mutcodon')) %>%
  select(-wtcodon.y) %>%
  rename(
    mutcodon_csc = csc, 
    mutcodon_rscu_fc = rscu_fc, 
    wtcodon = wtcodon.x
  ) %>%
  mutate(
    delta_csc = .$mutcodon_csc - .$wtcodon_csc, 
    delta_rscu_fc = .$mutcodon_rscu_fc - .$wtcodon_rscu_fc, 
    fitness_class = ifelse(.$mean_wrel_upper == 0, 'lethal', 
                           ifelse(.$mean_wrel_upper > 0 & .$mean_wrel_upper < 1, 'deleterious', 
                                  ifelse(.$mean_wrel_lower > 1, 'beneficial', 'neutral')))
  ) %>%
  rename(
    wt_aa = wtRes, 
    mut_aa = mutRes
  )

################################################################################
# Create a hash to convert the single letter amino acid residue to the three 
# letter abbreviation 

aa_conversion_hash = hash(
  c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'), 
  c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 
    'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'Stp')
)

mutation_data$wt_aa = vapply(mutation_data$wt_aa, function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)
mutation_data$mut_aa = vapply(mutation_data$mut_aa, function(x) aa_conversion_hash[[x]], character(1), USE.NAMES = FALSE)

################################################################################
# Write table for later use 

write.table(
  mutation_data, 
  './data/condensed_point_mutations_p9.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################







