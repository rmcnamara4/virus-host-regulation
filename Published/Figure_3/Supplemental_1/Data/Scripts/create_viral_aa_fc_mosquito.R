################################################################################
# script to create Mosquito relative amino acid fold change data for mosquito-borne
# viruses 
################################################################################

# load libraries 
library(tidyverse)
library(Biostrings)
library(reshape2)

# source required files 
source('../../../../../Src/codon_to_aa.R')
source('../../../../../Src/synonymous_codons.R')

# load data 
mosquito_seq = read.table(
  '../../../../../Experiments/0-Preprocessing/0.3-CreateSequenceTables/data/mosquito_AaloF1_seq.csv', 
  sep = ',', header = TRUE
)

virus_codons = read.table(
  '../../../../../Experiments/3-DownloadAllViruses/data/codon_tables/all_codon_counts_total.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# get mosquito codons 
mosquito_dna_string = DNAStringSet(mosquito_seq$coding)
mosquito_codons = as.data.frame(trinucleotideFrequency(mosquito_dna_string), step = 3)

# melt the data 
mosquito_codons_melted = melt(mosquito_codons)
names(mosquito_codons_melted) = c('codon', 'count') 

# get the sum of each codon 
mosquito_codons_melted = mosquito_codons_melted %>%
  group_by(codon) %>%
  summarize_all(
    sum
  )

# make amino acid column 
mosquito_codons_melted$aa = codon_to_aa(mosquito_codons_melted$codon)

# get sum of each amino acid 
mosquito_aa_melted = mosquito_codons_melted %>%
  group_by(aa) %>%
  summarize(
    count = sum(count)
  )

# get frequency of each aa 
mosquito_aa_melted$count = mosquito_aa_melted$count / sum(mosquito_aa_melted$count)

# get frequency of viral amino acids 
virus_metadata = virus_codons[1:6]
virus_codons = virus_codons[-c(1:6)]

names(virus_codons) = codon_to_aa(names(virus_codons)) 

virus_aa = as.data.frame(t(apply(virus_codons, 1, function(x) tapply(x, colnames(virus_codons), sum))))
virus_aa = as.data.frame(t(apply(virus_aa, 1, function(x) x / sum(x))))

# get the log2 ratio of viral usage / mosquito usage 
aa_fc = as.data.frame(sapply(1:ncol(virus_aa), function(x) {
  virus_aa[[x]] = log2(virus_aa[[x]] / mosquito_aa_melted$count[x])
}))
names(aa_fc) = names(virus_aa) 

# add back metadata 
aa_fc = cbind(virus_metadata, aa_fc)

# remove stop amino acid 
aa_fc = aa_fc %>%
  select(-Stp)

# save 
write.table(
  aa_fc, 
  '../viral_aa_fc_mosquito.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################










