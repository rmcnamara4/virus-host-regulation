################################################################################
# script to create supplemental table 1
################################################################################

# load libraries 
library(tidyverse)

# load data 
codon_counts = read.table(
  '../all_codon_counts_total.csv', 
  sep = ',', header = TRUE, quote = "\""
)

rscu_corr = read.table(
  '../all_correlations_total_human.csv', 
  sep = ',', header = TRUE, quote = "\""
)

# keep only viruses that have a human host 
iv = grepl('human', codon_counts$host)
codon_counts = codon_counts[iv, ]

# select specific columns from codon_counts 
codon_counts = codon_counts %>%
  select(virus, grep('^[AGCT]', colnames(.)))

# get codon composition in percentage 
codon_counts_p = as.data.frame(t(apply(codon_counts[-1], 1, function(x) {
  tot = sum(x)
  x = x / sum(x)
  return(x)
}))) %>%
  mutate(
    virus = codon_counts$virus, 
    .before = 1
  )

# merge with correlation data 
final_data = merge(rscu_corr, codon_counts_p, by = 'virus')

# get column index for codons 
iv = grep('^[ACGT]', colnames(final_data))

# rearrange columns 
final_data = final_data[, c(1:6, iv, 7:8)]

# remove columns 
final_data = final_data %>%
  select(-c(locus_tag, shape))

# fix host column 
final_data$host = str_replace(final_data$host, 'vertebrates: human', 'human:vertebrates')

# save 
write.table(
  final_data, 
  '../../supplemental_table_1.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################