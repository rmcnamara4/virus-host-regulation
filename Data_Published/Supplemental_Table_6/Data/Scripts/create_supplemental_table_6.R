################################################################################
# script to create supplemental table 6
################################################################################

# load libraries
library(tidyverse)

# load data 
optimalities = read.table(
  '../codon_optimalities.csv', 
  sep = ',', header = TRUE
) %>%
  dplyr::select(codon, human_csc, mosquito_csc)

syn_mrf_data = read.table(
  '../syn_mut_mrf_per_species_p9.csv', 
  sep = ',', header = TRUE
) %>%
  dplyr::select(host, mutcodon, mean_wrel) %>%
  rename(
    codon = mutcodon, 
    syn_mrf = mean_wrel
  )

rscu_fc_data = read.table(
  '../rscu_fc.csv', 
  sep = ',', header = TRUE
) %>%
  dplyr::select(codon, aa, human_rscu_fc, mosquito_rscu_fc)

# split data into different hosts 
mrf_human = syn_mrf_data %>%
  filter(host == 'Human') %>%
  dplyr::select(-host) %>%
  rename(
    human_syn_mrf = syn_mrf
  )

mrf_mosquito = syn_mrf_data %>%
  filter(host == 'Mosquito') %>%
  dplyr::select(-host) %>%
  rename(
    mosquito_syn_mrf = syn_mrf
  )

# merge all data 
final_data = Reduce(function(x, y) merge(x, y, by = 'codon'), list(mrf_human, mrf_mosquito, rscu_fc_data, optimalities))

# rearrange columns 
final_data = final_data[c(1, 4, 2:3, 5:8)]

# save 
write.table(
  final_data, 
  '../../supplemental_table_6.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################
