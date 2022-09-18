################################################################################
# script to create supplemental table 4
################################################################################

# load libraries 
library(tidyverse)

# load data 
human_opt = read.table(
  '../human_endo_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mean_endo_csc) %>%
  rename(
    human_csc = mean_endo_csc
  )

mosquito_opt = read.table(
  '../mosquito_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mosquito) %>%
  rename(
    mosquito_csc = mosquito
  )

viral_h_rscu = read.table(
  '../all_rscu_ratio_total_human.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(-c(realm, locus_tag, host, type, shape)) %>%
  filter(virus %in% c('Dengue virus 2', 'Zika virus', 'Chikungunya virus'))

viral_m_rscu = read.table(
  '../all_rscu_ratio_total_mosquito.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(-c(realm, locus_tag, host, type, shape)) %>%
  filter(virus %in% c('Dengue virus 2', 'Zika virus', 'Chikungunya virus'))

isolates_h_rscu = read.table(
  '../dengue_isolates_h-relative_rscu_fc.csv', 
  sep = ',', header = TRUE
) %>%
  group_by(species) %>%
  summarize_all(median)

isolates_m_rscu = read.table(
  '../dengue_isolates_m-relative_rscu_fc.csv', 
  sep = ',', header = TRUE
) %>%
  group_by(species) %>%
  summarize_all(median)

# transpose the data.frames
rownames(viral_h_rscu) = viral_h_rscu$virus
viral_h_rscu = as.data.frame(t(viral_h_rscu[-1]))

rownames(viral_m_rscu) = viral_m_rscu$virus
viral_m_rscu = as.data.frame(t(viral_m_rscu[-1]))

isolates_h_rscu = as.data.frame(isolates_h_rscu)
rownames(isolates_h_rscu) = isolates_h_rscu$species
isolates_h_rscu = as.data.frame(t(isolates_h_rscu[-1]))

isolates_m_rscu = as.data.frame(isolates_m_rscu)
rownames(isolates_m_rscu) = isolates_m_rscu$species
isolates_m_rscu = as.data.frame(t(isolates_m_rscu[-1]))

# rename columns 
colnames(viral_h_rscu) = c('Zika_h_rscu', 'Chikungunya_h_rscu', 'Dengue2_h_rscu')
colnames(viral_m_rscu) = c('Zika_m_rscu', 'Chikungunya_m_rscu', 'Dengue2_m_rscu')
colnames(isolates_h_rscu) = c('Dengue1_isolate_h_rscu', 'Dengue2_isolate_h_rscu', 'Dengue3_isolate_h_rscu', 'Dengue4_isolate_h_rscu')
colnames(isolates_m_rscu) = c('Dengue1_isolate_m_rscu', 'Dengue2_isolate_m_rscu', 'Dengue3_isolate_m_rscu', 'Dengue4_isolate_m_rscu')

# make codon column again 
viral_h_rscu$codon = rownames(viral_h_rscu)
isolates_h_rscu$codon = rownames(isolates_h_rscu) 
viral_m_rscu$codon = rownames(viral_m_rscu)
isolates_m_rscu$codon = rownames(isolates_m_rscu)

# merge all data 
data_list = list(viral_h_rscu, isolates_h_rscu, human_opt, viral_m_rscu, isolates_m_rscu, mosquito_opt) 
final_data = Reduce(function(x, y) merge(x, y, by = 'codon'), data_list)

# save 
write.table(
  final_data, 
  '../../supplemental_table_4.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################
