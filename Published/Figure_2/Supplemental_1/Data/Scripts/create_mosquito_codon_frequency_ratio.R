################################################################################
# script to create mosquito codon frequency table 
################################################################################

# load libraries
library(tidyverse)
library(readxl)

# load data 
mosquito_codons = read_excel(
  '../../../../../Experiments/0-Preprocessing/0.4-CreateSequenceStatsTables/data/mosquito_AaloF1_stats.xlsx', 
  sheet = 'Codons'
)

mosquito_opt = read.table(
  '../../Data/mosquito_csc.csv', 
  sep = ',', header = TRUE
)

# get quantiles of mosquito optimality 
quant = quantile(mosquito_opt$mosquito)

# define optimal and non-optimal codons
codon_class = list(
  optimal = mosquito_opt[mosquito_opt$mosquito >= quant[['75%']], ]$codon, 
  non_optimal = mosquito_opt[mosquito_opt$mosquito <= quant[['25%']], ]$codon
)

# calculate optimal and non-optimal percentage of each gene for mosquito 
mosquito_codons$optimal = rowSums(mosquito_codons[codon_class$optimal])
mosquito_codons$non_optimal = rowSums(mosquito_codons[codon_class$non_optimal])

# calculate ratio of optimal to non-optimal codon frequency for each gene for mosquito 
mosquito_codons$ratio = log2(mosquito_codons$optimal / mosquito_codons$non_optimal)

# save 
write.table(
  mosquito_codons,
  '../mosquito_codon_frequency_ratio.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################