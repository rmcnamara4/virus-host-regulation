################################################################################
# Import libraries

library(dplyr)

################################################################################
# Load CSC data

csc_data = read.table(
  './data/human_csc_all_methods.csv', 
  sep = ',', header = TRUE
)

################################################################################
# Select only the first five columns (corresponding to codon, aa, and endo CSC)

csc_data = csc_data[1:5]

################################################################################
# Get average of endogenous CSCs

mean_endogenous_csc = rowMeans(csc_data[3:5])

################################################################################
# Add column for mean CSCs and rename columns 

csc_data = csc_data %>%
  mutate(
    mean_endo_csc = mean_endogenous_csc
  ) %>%
  rename(
    aa = Name, 
    endo_293T = X293T_endo, 
    endo_HeLa = HeLa_endo, 
    endo_RPE = RPE_endo
  )

################################################################################
# Write table

write.table(
  csc_data, 
  './data/human_endo_csc.csv', 
  sep = ',', 
  row.names = FALSE, 
  col.names = TRUE, 
  quote = FALSE
)

################################################################################