################################################################################
# script to create optimal and non-optimal frequency tables for human and 
# Dengue 2 isolates 
################################################################################

# load libraries 
library(tidyverse) 
library(readxl)

# load data 
human_codons = read_excel(
  '../../../../../Experiments/0-Preprocessing/0.4-CreateSequenceStatsTables/data/human_hg38_stats.xlsx', 
  sheet = 'Codons'
)

dengue_2_codons = read_excel(
  '../../../../../Experiments/0-Preprocessing/0.4-CreateSequenceStatsTables/data/dengue_2_isolates_stats.xlsx', 
  sheet = 'Codons'
)

human_opt = read.table(
  '../../Data/human_endo_csc.csv', 
  sep = ',', header = TRUE
)

# combine human_codons and dengue_2_codons 
codon_data = rbind(human_codons, dengue_2_codons)

# create table of optimal and non-optimal codons 
# optimal defined as 25% most optimal codons
# non-optimal defined as 25% most non-optimal codons 
quant = quantile(human_opt$mean_endo_csc)
codon_df = data.frame(
  optimal = human_opt[human_opt$mean_endo_csc >= quant[['75%']], ]$codon, 
  non_optimal = human_opt[human_opt$mean_endo_csc <= quant[['25%']], ]$codon
)

# calculate percentage of optimal and non-optimal codons for each gene/isolate
codon_data$optimal_perc = rowSums(codon_data[codon_df$optimal])
codon_data$non_optimal_perc = rowSums(codon_data[codon_df$non_optimal])

# calculate ratio of optimal to non-optimal 
codon_data$ratio = log2(codon_data$optimal_perc / codon_data$non_optimal_perc)

# save 
write.table(
  codon_data, 
  '../dengue_2_iso_and_human_codon_freq.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################