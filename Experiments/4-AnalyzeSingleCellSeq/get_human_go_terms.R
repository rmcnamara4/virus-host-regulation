################################################################################
# script to create GO term table for human genome
################################################################################

# load libraries
library(tidyverse)
library(biomaRt)
library(stringr)

# create mart variable 
ensembl = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

# collect GO terms
go_terms = getBM(c('ensembl_gene_id', 'go_id', 'name_1006'), mart = ensembl, uniqueRows = TRUE)

# rename columns 
colnames(go_terms) = c('gene_ID', 'go_id', 'go_name')

# replace empty cells with NA 
go_terms[go_terms == ''] = NA

# save 
write.table(
  go_terms, 
  './data/go_terms/human_go_terms.tsv', 
  sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################
                 