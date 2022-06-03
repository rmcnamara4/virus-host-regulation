################################################################################
# script to create virus rscu fold change relative to mosquito table
################################################################################

# load libraries 
library(tidyverse)

# load RSCU fold change data for mosquito 
rscu_fc_mosquito = read.table(
  '../../../../../Experiments/3-DownloadAllViruses/data/rscu_ratios/all_rscu_ratio_total_mosquito.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# load correlation data 
corr_data = read.table(
  '../../../../../Experiments/3-DownloadAllViruses/data/rscu_correlations/all_correlations_total_mosquito.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# match virus names 
iv = match(rscu_fc_mosquito$virus, corr_data$virus)

# combine tables 
final_table = rscu_fc_mosquito %>%
  mutate(
    corr = corr_data$corr[iv], 
    p.value = corr_data$p.value[iv]
  )

# write table 
write.table(
  final_table, 
  '../viral_rscu_fc_mosquito.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################
