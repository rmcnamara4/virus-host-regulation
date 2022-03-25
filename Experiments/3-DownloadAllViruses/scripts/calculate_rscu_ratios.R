################################################################################
# Load libraries

library(dplyr)
library(readxl)
library(writexl)

################################################################################
# Source functions

source('../../Src/calc_rscu.R')
source('../../Src/synonymous_codons.R')

################################################################################
# Load rscu data 

viral_rscu = read.table(
  snakemake@input[['all_rscu']], 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

viral_rscu_total = read.table(
  snakemake@input[['all_rscu_total']], 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

################################################################################
# Load human and mosquito codon counts

human_counts = read_excel(
  '../0-Preprocessing/0.4-CreateSequenceStatsTables/data/human_hg38_stats.xlsx', 
  sheet = 'Codons'
)
human_counts = human_counts[-c(1:3)]

mosquito_counts = read_excel(
  '../0-Preprocessing/0.4-CreateSequenceStatsTables/data/mosquito_AaloF1_stats.xlsx', 
  sheet = 'Codons'
)
mosquito_counts = mosquito_counts[-c(1:3)]

################################################################################
# Get codon totals for human and mosquito 

human_counts_total = apply(human_counts, 2, sum)
mosquito_counts_total = apply(mosquito_counts, 2, sum)

################################################################################
# Calculate RSCU for human and mosquito 

human_rscu = calc_rscu(human_counts_total)
mosquito_rscu = calc_rscu(mosquito_counts_total)

################################################################################
# Calculate ratio for each 

rscu_ratio_human = viral_rscu
rscu_ratio_mosquito = viral_rscu
rscu_total_ratio_human = viral_rscu_total
rscu_total_ratio_mosquito = viral_rscu_total 

for (c in names(human_rscu)) {
  
  rscu_ratio_human[[c]] = log2(rscu_ratio_human[[c]] / human_rscu[[c]])
  rscu_ratio_mosquito[[c]] = log2(rscu_ratio_mosquito[[c]] / mosquito_rscu[[c]])
  
  rscu_total_ratio_human[[c]] = log2(rscu_total_ratio_human[[c]] / human_rscu[[c]])
  rscu_total_ratio_mosquito[[c]] = log2(rscu_total_ratio_mosquito[[c]] / mosquito_rscu[[c]])
  
}

################################################################################
# Write files 

write.table(
  rscu_ratio_human, 
  snakemake@output[['rscu_ratio_human']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

write.table(
  rscu_ratio_mosquito, 
  snakemake@output[['rscu_ratio_mosquito']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

write.table(
  rscu_total_ratio_human, 
  snakemake@output[['rscu_ratio_total_human']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

write.table(
  rscu_total_ratio_mosquito, 
  snakemake@output[['rscu_ratio_total_mosquito']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

################################################################################
