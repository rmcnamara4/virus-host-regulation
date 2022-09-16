################################################################################
# script to create supplemental table 2
################################################################################

# load libraries 
library(tidyverse)

# load data 
human_opt = read.table(
  '../human_endo_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, aa, mean_endo_csc) %>%
  rename(
    human_csc = mean_endo_csc
  )

mosquito_opt = read.table(
  '../mosquito_csc.csv', 
  sep = ',', header = TRUE
)

# append _csc to the columns in mosquito_opt 
colnames(mosquito_opt)[-1] = paste0(colnames(mosquito_opt)[-1], '_csc')

# merge both data.frames
final_data = merge(human_opt, mosquito_opt, by = 'codon')

# rearrange columns 
final_data = final_data[, c(1:2, 4:7, 3)]

# save 
write.table(
  final_data, 
  '../../supplemental_table_2.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################
