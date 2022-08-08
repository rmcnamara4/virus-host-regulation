################################################################################
# script to create EGFP iCodon reporters codon frequency table 
################################################################################

# load libraries 
library(tidyverse)
library(Biostrings)

# load data 
egfp_reporters = read.table(
  '../../../../../Experiments/6-OptimalityReporters/egfp_icodon_reporters.csv', 
  sep = ',', header = TRUE
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

# calculate codon frequency of the reporters 
egfp_reporters_dna_string = DNAStringSet(egfp_reporters$seq)
egfp_reporters_codons = trinucleotideFrequency(egfp_reporters_dna_string, step = 3, as.prob = TRUE)

# column bind the metadata with the codon frequency for the reporters 
egfp_reporters = cbind(egfp_reporters, as.data.frame(egfp_reporters_codons))

# calculate optimal and non-optimal frequency for each reporter 
egfp_reporters$optimal = rowSums(egfp_reporters[codon_class$optimal])
egfp_reporters$non_optimal = rowSums(egfp_reporters[codon_class$non_optimal])

# calculate ratio of optimal to non-optimal for each reporter 
egfp_reporters$ratio = log2(egfp_reporters$optimal / egfp_reporters$non_optimal)

# save 
write.table(
  egfp_reporters, 
  '../egfp_reporters_codon_frequency_ratio.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################