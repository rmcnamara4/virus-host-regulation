################################################################################
# Load libraries 

library(edgeR)
library(dplyr)

################################################################################
# Load count files

dengue_counts = read.table('./data/counts/counts_dengue.tsv', 
                           sep = '\t', header = TRUE)
zika_counts = read.table('./data/counts/counts_zika.tsv', 
                         sep = '\t', header = TRUE)

################################################################################
# Load metadata files

dengue_metadata = read.table('./data/metadata/cell_metadata_dengue.tsv', 
                             sep = '\t', header = TRUE)
zika_metadata = read.table('./data/metadata/cell_metadata_zika.tsv', 
                           sep = '\t', header = TRUE)

################################################################################
# Remove X from column names 

colnames(dengue_counts) = sub('^X', '', colnames(dengue_counts))
colnames(zika_counts) = sub('^X', '', colnames(zika_counts))

################################################################################
# Remove rows that aren't endogenous genes or the virus 

iv = grep('^ENSG', dengue_counts$EnsemblID)
dengue_counts = dengue_counts[iv, ]
write.table(dengue_counts, './data/counts/counts_dengue.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

iv = grep('^ENSG', zika_counts$EnsemblID)
zika_counts = zika_counts[iv, ]
write.table(zika_counts, './data/counts/counts_zika.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

################################################################################
# Add the counts of the Dengue and Zika genomes to the count files 

iv = match(colnames(dengue_counts)[-1], dengue_metadata$name)
dengue_counts = rbind(c('Dengue', dengue_metadata$numberDengueReads[iv]), dengue_counts)

iv = match(colnames(zika_counts)[-1], zika_metadata$name)
zika_counts = rbind(c('Zika', zika_metadata$numberZikaReads[iv]), zika_counts)

################################################################################
# Calculate Counts Per Million (CPM)

dengue_counts[-1] = apply(dengue_counts[-1], 2, as.numeric)
dengue_cpm = as.data.frame(cpm(dengue_counts[-1])) %>%
  mutate(
    EnsemblID = dengue_counts$EnsemblID, 
    .before = 1
  )

zika_counts[-1] = apply(zika_counts[-1], 2, as.numeric)
zika_cpm = as.data.frame(cpm(zika_counts[-1])) %>%
  mutate(
    EnsemblID = zika_counts$EnsemblID, 
    .before = 1
  )

################################################################################
# Write tables

write.table(dengue_cpm, './data/cpms/cpm_dengue.csv', sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(zika_cpm, './data/cpms/cpm_zika.csv', sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE)

################################################################################
