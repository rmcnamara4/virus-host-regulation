################################################################################
# script to create supplemental table 5
################################################################################

# load libraries
library(tidyverse)
library(readxl)

# load data 
edgeR_data = read.table(
  '../uninfected-high.csv', 
  sep = ',', header = TRUE
)

codon_data = read_excel(
  '../human_hg38_stats.xlsx', 
  sheet = 'Codons'
)

go_data = read.table(
  '../human_go_terms.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# remove columns from data.frames
codon_data = codon_data %>%
  select(-c(species, len))

edgeR_data = edgeR_data %>%
  select(-c('F', 'PValue'))

# merge all data 
final_data = merge(codon_data, edgeR_data, all.y = TRUE)
final_data = merge(go_data, final_data, all.y = TRUE)

# remove rows with NA that have duplicated entries 
iv = duplicated(final_data$gene_ID) & is.na(final_data$go_id)
final_data = final_data[!iv, ]

# add up/down regulated column 
final_data$dir = ifelse(
  final_data$padj < 0.01 & final_data$logFC > 1, 'upregulated', 
  ifelse(
    final_data$padj < 0.01 & final_data$logFC < -1, 'downregulated', 'none'
  )
)

# define denguenized and non-denguenized codons 
denguenized_codons = c('AGA', 'AAA', 'GGA', 'ATA', 'GCA', 'ACA', 'CCA', 'GAA')
non_denguenized_codons = c('CCC', 'CCG', 'CGG', 'CCT', 'TCC', 'GCG', 'CGC', 'CGA')

# create columns for denguenized/non-denguenized codon percentage and ratio 
final_data$denguenized_p = rowSums(final_data[denguenized_codons])
final_data$non_denguenized_p = rowSums(final_data[non_denguenized_codons])
final_data$ratio = log2(final_data$denguenized_p / final_data$non_denguenized_p)

# save 
write.table(
  final_data, 
  '../../supplemental_table_5.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################