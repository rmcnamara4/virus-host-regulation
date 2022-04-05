################################################################################
# Load libraries

library(dplyr)

################################################################################
# Load in rscu files and rbind them all

rscu_files = lapply(
  snakemake@input[['rscu']],
  function(x) {
    read.table(x, header = TRUE, sep = ',', quote = "\"")
  }
)

rscu_df = do.call('rbind', rscu_files)

################################################################################
# Do the same thing for the rscu totals

rscu_total_files = lapply(
  snakemake@input[['rscu_total']],
  function(x) {
    read.table(x, header = TRUE, sep = ',', quote = "\"")
  }
)

rscu_total_df = do.call('rbind', rscu_total_files)

################################################################################
# Write the concatenated files

write.table(
  rscu_df,
  snakemake@output[['all_rscu']],
  sep = ',',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

write.table(
  rscu_total_df,
  snakemake@output[['all_rscu_total']],
  sep = ',',
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

################################################################################
