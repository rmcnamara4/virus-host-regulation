################################################################################
# script to create 1nt out of frame reporters codon frequency table 
################################################################################

# load libraries 
library(tidyverse)
library(Biostrings)

# load data 
nt_out_of_frame_reporters = read.table(
  '../../../../../Experiments/6-OptimalityReporters/nt_out_of_frame_reporters.csv', 
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
nt_out_of_frame_reporters_dna_string = DNAStringSet(nt_out_of_frame_reporters$seq)
nt_out_of_frame_reporters_codons = trinucleotideFrequency(nt_out_of_frame_reporters_dna_string, 
                                                          step = 3, as.prob = TRUE)

# column bind the metadata with the codon frequency of the reporters 
nt_out_of_frame_reporters = cbind(nt_out_of_frame_reporters, 
                                  as.data.frame(nt_out_of_frame_reporters_codons))

# calculate optimal, non-optimal, and neutral frequency of each reporter
nt_out_of_frame_reporters$optimal = rowSums(nt_out_of_frame_reporters[codon_class$optimal])
nt_out_of_frame_reporters$non_optimal = rowSums(nt_out_of_frame_reporters[codon_class$non_optimal])
nt_out_of_frame_reporters$neutral =  1 - (nt_out_of_frame_reporters$optimal + nt_out_of_frame_reporters$non_optimal)

# calculate ratio of optimal to non-optimal for each reporter 
nt_out_of_frame_reporters$ratio = log2(nt_out_of_frame_reporters$optimal / nt_out_of_frame_reporters$non_optimal)

# save
write.table(
  nt_out_of_frame_reporters, 
  '../nt_out_of_frame_reporters_codon_frequency_ratio.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################
