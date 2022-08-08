################################################################################
# script to create Dengue isolates rscu fold change relative to human table
################################################################################

# load libraries 
library(tidyverse)
library(readxl)
library(Biostrings)

# source required files
source('../../../../../Src/calc_rscu.R')
source('../../../../../Src/synonymous_codons.R')

# load data 
human_seq = read.table(
  '../../../../../Experiments/0-Preprocessing/0.3-CreateSequenceTables/data/human_hg38_seq.csv', 
  sep = ',', header = TRUE
)

files = list.files('../../../../../Experiments/0-Preprocessing/0.4-CreateSequenceStatsTables/data', 
                   pattern = 'dengue*', full.names = TRUE)
dengue_isolates = lapply(files, function(f) read_excel(f, sheet = 'RSCU'))
dengue_isolates = Reduce('rbind', dengue_isolates)
dengue_isolates = dengue_isolates[-c(2, 3)]

# calculate codon frequency of human sequences 
human_dna_string = DNAStringSet(human_seq$coding)
human_codons = trinucleotideFrequency(human_dna_string, step = 3)

# sum up codon counts for the whole genome 
human_codons_total = apply(human_codons, 2, sum)

# calculate RSCU 
human_rscu = calc_rscu(human_codons_total)

# calculate H-relative RSCU fold change for Dengue isolates 
for (i in 1:64) {
  dengue_isolates[[i + 1]] = log2(dengue_isolates[[i + 1]] / human_rscu[[i]])
}

# write table 
write.table(
  dengue_isolates, 
  '../dengue_isolates_h-relative_rscu_fc.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################