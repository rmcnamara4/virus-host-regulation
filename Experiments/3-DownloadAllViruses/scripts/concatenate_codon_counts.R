################################################################################
# Load libraries 

library(dplyr)

################################################################################
# Load in counts file and rbind all of them

counts_files = lapply(snakemake@input[['counts']], function(x) {
  read.table(x, header = TRUE, sep = ',', quote = "\"")
})

counts_df = do.call("rbind", counts_files)

################################################################################
# Do the same thing for the totals

counts_total_files = lapply(snakemake@input[['counts_total']], function(x) {
  read.table(x, header = TRUE, sep = ',', quote = "\"")
})

counts_total_df = do.call("rbind", counts_total_files)

################################################################################
# Write the concatenated files

write.table(
  counts_df, 
  snakemake@output[['all_counts']],
  col.names = TRUE, 
  row.names = FALSE, 
  sep = ',', 
  quote = FALSE
)

write.table(
  counts_total_df, 
  snakemake@output[['all_counts_total']],
  col.names = TRUE, 
  row.names = FALSE, 
  sep = ',', 
  quote = FALSE
)

################################################################################
