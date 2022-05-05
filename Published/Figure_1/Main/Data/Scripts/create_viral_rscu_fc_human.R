################################################################################
# script to create virus rscu fold change relative to human table
################################################################################

# load libraries 
library(tidyverse)

# load RSCU fold change data for human 
rscu_fc_human = read.table(
  '../../../../../Experiments/3-DownloadAllViruses/data/rscu_ratios/all_rscu_ratio_total_human.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# load correlation data 
corr_data = read.table(
  '../../../../../Experiments/3-DownloadAllViruses/data/rscu_correlations/all_correlations_total_human.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# match virus names 
iv = match(rscu_fc_human$virus, corr_data$virus)

# combine tables 
final_table = rscu_fc_human %>%
  mutate(
    corr = corr_data$corr[iv], 
    p.value = corr_data$p.value[iv]
  )

# write table
write.table(
  final_table,
  '../viral_rscu_fc_human.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################