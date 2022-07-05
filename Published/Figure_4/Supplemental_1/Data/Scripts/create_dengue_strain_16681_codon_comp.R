################################################################################
# script to create codon composition table for Dengue 2 strain 16681 
################################################################################

# load libraries 
library(tidyverse)
library(Biostrings)
library(seqinr)
library(reshape2)

# load Dengue strain 16681 sequence 
char_seq = read.fasta(
  '../../../../../Experiments/4-AnalyzeSingleCellSeq/data/fastas/dengue_2_16681.fa', 
  as.string = TRUE
)
char_seq = as.character(char_seq)

# convert sequence to a DNAString 
DNA_seq = DNAString(char_seq)

# calculate codon counts and codon frequency 
codon_counts = trinucleotideFrequency(DNA_seq, step = 3)
codon_freq = trinucleotideFrequency(DNA_seq, step = 3, as.prob = TRUE)

# melt both of the vectors 
codon_counts_melted = melt(codon_counts)
codon_freq_melted = melt(codon_freq)

# merge the two tables
codon_comp = merge(codon_counts_melted, codon_freq_melted, by = 'row.names')
names(codon_comp) = c('codon', 'count', 'freq')

# write table 
write.table(
  codon_comp, 
  '../dengue_strain_16681_codon_comp.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################