################################################################################
# Load libraries

library(dplyr)

################################################################################
# Source helper functions

source('/n/projects/rm2498/Virus_Project/Viral_Seq_Analysis/Scripts/Functions/create_rscu_df_functions.R')

################################################################################
# Load in count and total count data 

counts = read.table(
  snakemake@input[['codon_counts']], 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

counts_total = read.table(
  snakemake@input[['codon_counts_total']], 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

################################################################################
# Get the codon data from both data.frames

iv = grep('^[ACGT]', names(counts))
codon_counts = counts[iv]
counts_info = counts[-iv]

iv = grep('^[ACGT]', names(counts_total))
codon_counts_total = counts_total[iv]
counts_total_info = counts_total[-iv]

################################################################################
# Calculate RSCU 

rscu = calc_rscu(codon_counts)
rscu_total = calc_rscu(codon_counts_total)

################################################################################
# rbind RSCUs to the other info 

rscu_df = cbind(counts_info, rscu)
rscu_total_df = cbind(counts_total_info, rscu_total)

################################################################################
# Write tables

write.table(
  rscu_df, 
  snakemake@output[['rscu']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

write.table(
  rscu_total_df, 
  snakemake@output[['rscu_total']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

################################################################################





