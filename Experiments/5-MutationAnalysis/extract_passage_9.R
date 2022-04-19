################################################################################
# Load libraries 

library(dplyr)

################################################################################
# Load point mutation data 

mutations = read.table('./data/point_mutations.csv', 
                       sep = ',', header = TRUE)

################################################################################
# Filter for passage 9 mutations in the CDS 

mutations_p9 = mutations %>%
  filter(passage == 'MF.9', muttype != 'UTR', reg != 'UTR-5', reg != 'UTR3')
mutations_p9 = mutations_p9[-1]

################################################################################
# Write table 

write.table(
  mutations_p9, 
  './data/point_mutations_p9.csv', 
  sep = ',', col.names = TRUE, row.names = FALSE, quote = FALSE
)

################################################################################