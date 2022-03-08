################################################################################
# Load libraries 

library(dplyr)
library(seqinr)
library(Biostrings)
library(tibble)
library(stringr)

################################################################################
# Load human and mosquito fasta files

human_fa = readDNAStringSet('../0.1-GetHostCDSFastas/data/hg38.Ens_102.cds.fa')
mosquito_fa = readDNAStringSet('../0.1-GetHostCDSFastas/data/AaloF1.EnsGen_50.cds.fa')

################################################################################
# Define function to extract gene_ID

get_gene_ID = function(x) {
  str_extract(x, '\\|.*') %>%
    str_replace('\\|', '')
}

################################################################################
# Create data.frame of names and sequences 

human_seq_names = names(human_fa)
mosquito_seq_names = names(mosquito_fa)

human_seq = paste(human_fa)
mosquito_seq = paste(mosquito_fa)

human_df = data.frame(
  names = human_seq_names, 
  coding = human_seq
)
mosquito_df = data.frame(
  names = mosquito_seq_names, 
  coding = mosquito_seq
)

################################################################################
# Create final seq tables

# Human table
human_data = human_df %>%
  mutate(
    gene_ID = get_gene_ID(names) 
  ) %>%
  select(gene_ID, coding) %>%
  add_column(species = 'Human', .before = 1)

human_data = human_data[order(human_data$gene_ID, -nchar(human_data$coding)), ]
human_data = human_data[!duplicated(human_data$gene_ID), ]
human_data$gene_ID = str_replace(human_data$gene_ID, '\\..*', '')

# Remove genes that have letters other than A, C, G, T
iv = grep('[^ACGT]', human_data$coding)
if (length(iv) > 0) {
  human_data = human_data[-iv, ]
}

# Remove genes that aren't divisible by 3
iv = which(nchar(human_data$coding) %% 3 != 0)
if (length(iv) > 0) {
  human_data = human_data[-iv, ]
}

# Mosquito table
mosquito_data = mosquito_df %>%
  mutate(
    gene_ID = get_gene_ID(names) 
  ) %>%
  select(gene_ID, coding) %>%
  add_column(species = 'Mosquito', .before = 1)

mosquito_data = mosquito_data[order(mosquito_data$gene_ID, -nchar(mosquito_data$coding)), ]
mosquito_data = mosquito_data[!duplicated(mosquito_data$gene_ID), ]
mosquito_data$gene_ID = str_replace(mosquito_data$gene_ID, '\\..*', '')

# Remove genes that have letters other than A, C, G, T
iv = grep('[^ACGT]', mosquito_data$coding)
if (length(iv) > 0) {
  mosquito_data = mosquito_data[-iv, ]
}

# Remove genes that aren't divisible by 3
iv = which(nchar(mosquito_data$coding) %% 3 != 0)
if (length(iv) > 0) {
  mosquito_data = mosquito_data[-iv, ]
}

################################################################################
# Write tables 

write.table(
  human_data, 
  './data/human_hg38_seq.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

write.table(
  mosquito_data, 
  './data/mosquito_AaloF1_seq.csv', 
  sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################

