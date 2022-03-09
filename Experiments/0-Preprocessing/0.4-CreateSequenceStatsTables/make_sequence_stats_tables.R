################################################################################
# Set Java parameters to avoid overflow

options(java.parameters = "-Xmx8000m")

################################################################################
# Load libraries 

library(dplyr)
library(seqinr)
library(Biostrings)
library(tibble)
library(writexl)
library(stringr)

################################################################################
# Initiate file names

files = list.files('../0.3-CreateSequenceTables/data/') %>%
  str_replace('seq.csv', '')

input = paste0('../0.3-CreateSequenceTables/data/', files, 'seq.csv')
output = paste0('./data/', files, 'stats.xlsx')

################################################################################
# Source rscu function and synonymous codons

source('../../../Src/calc_rscu.R')
source('../../../Src/synonymous_codons.R')

################################################################################
# Write loop

sapply(1:length(input), function(x) {
  
  ################################################################################
  # Initiate names of amino acids
  
  amino_acids = aaa(AA_ALPHABET)
  amino_acids = amino_acids[!is.na(amino_acids)]
  
  ################################################################################
  # Load in data 
  
  data = read.table(
    input[x], 
    sep = ',', header = TRUE
  ) %>%
    select(coding, species, gene_ID)
  
  ################################################################################
  # Make the coding sequences a DNAStringSet
  
  data_seq = DNAStringSet(data$coding)
  
  ################################################################################
  # Get the overall nucleotides probability of the coding sequences
  
  nucleotides_all = oligonucleotideFrequency(data_seq, width = 1, as.prob = TRUE)
  
  ################################################################################
  # Get the first position nucleotide probability of the coding sequences
  
  nucleotides_first = oligonucleotideFrequency(data_seq, width = 1, step = 3, as.prob = TRUE)
  colnames(nucleotides_first) = paste0(colnames(nucleotides_first), '1')
  
  ################################################################################
  # Get the second position nucleotide probability of the coding sequences 
  
  nucleotides_second = oligonucleotideFrequency(subseq(data_seq, start = 2), width = 1, step = 3, 
                                                as.prob = TRUE)
  colnames(nucleotides_second) = paste0(colnames(nucleotides_second), '2')
  
  ################################################################################
  # Get the third position nucleotide probability of the coding sequences 
  
  nucleotides_third = oligonucleotideFrequency(subseq(data_seq, start = 3), width = 1, step = 3, 
                                               as.prob = TRUE)
  colnames(nucleotides_third) = paste0(colnames(nucleotides_third), '3')
  
  ################################################################################
  # Combine all nucleotide data into one data.frame
  
  nucleotides = cbind(nucleotides_all, nucleotides_first, nucleotides_second, nucleotides_third) %>%
    as.data.frame(.) %>%
    add_column(species = data$species, .before = 1) %>%
    add_column(gene_ID = data$gene_ID, .before = 2) %>%
    add_column(len = nchar(data$coding), .before = 3)
  
  ################################################################################
  # Get the codon frequency of the coding sequences 
  
  codons = trinucleotideFrequency(data_seq, step = 3, as.prob = TRUE) %>%
    as.data.frame(.) %>%
    add_column(species = data$species, .before = 1) %>%
    add_column(gene_ID = data$gene_ID, .before = 2) %>%
    add_column(len = nchar(data$coding), .before = 3)
  
  ################################################################################
  # Translate sequences into amino acids and get the probability of the 
  # amino acids in the coding sequences 
  
  aa = Biostrings::translate(data_seq) %>%
    alphabetFrequency(as.prob = TRUE)
  aa = aa[, colSums(aa) != 0]
  
  colnames(aa) = amino_acids
  
  aa = aa %>%
    as.data.frame(.) %>%
    add_column(species = data$species, .before = 1) %>%
    add_column(gene_ID = data$gene_ID, .before = 2) %>%
    add_column(len = nchar(data$coding), .before = 3)
  
  ################################################################################
  # Calculate the RSCU of each codon for each of the coding sequences 
  
  rscu = calc_rscu(codons[-1:-3]) %>%
    as.data.frame(.) %>%
    add_column(species = data$species, .before = 1) %>%
    add_column(gene_ID = data$gene_ID, .before = 2) %>%
    add_column(len = nchar(data$coding), .before = 3)
  
  ################################################################################
  # Write the data.frames to an Excel file as separate sheets 
  
  dataset_names = list(
    "Nucleotides" = nucleotides, 
    "Codons" = codons, 
    "Amino_Acids" = aa, 
    "RSCU" = rscu
  )
  
  write_xlsx(
    dataset_names, 
    output[x], 
    format_headers = F
  )
  
  ################################################################################
  
})
