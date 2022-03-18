################################################################################
# Load libraries 

library(dplyr)
library(rio)
library(tibble)
library(stats)
library(writexl)

################################################################################
# Load and split data 

human_ratio = read.table(
  snakemake@input[['human_ratios']], 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

human_ratio = human_ratio %>%
  select(-c('TAG', 'TGA', 'TAA', 'ATG', 'TGG'))

mosquito_ratio = read.table(
  snakemake@input[['mosquito_ratios']], 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

mosquito_ratio = mosquito_ratio %>%
  select(-c('TAG', 'TGA', 'TAA', 'ATG', 'TGG'))

################################################################################
# Load codon optimalities 

optimalities = read.table('/n/projects/rm2498/Virus_Project/Data_Files/Basic_Files/codon_optimalities.csv', 
                          sep = ',', header = TRUE)

optimalities = optimalities %>%
  filter(!(codon %in% c('TAG', 'TGA', 'TAA', 'ATG', 'TGG')))

################################################################################
# Get indices for the codons and make sure the optimalities are in the same order

iv = grep('^[ACTG]', names(human_ratio))
ind = match(names(human_ratio)[iv], optimalities$codon)

optimalities = optimalities[ind, ]

################################################################################
# Create function to turn infinite values to NA

inf2NA = function(x) { x[is.infinite(x)] = NA; x }

################################################################################
# Get the correlation stats for humans

human_corr_stats = 
  sapply(1:nrow(human_ratio), function(x) {
    
  row = as.numeric(human_ratio[x, iv])
  correlation = cor.test(inf2NA(row), inf2NA(optimalities$human_csc), method = 'spearman')
  return(c(corr = correlation$estimate[[1]], p.value = correlation$p.value))
  
}) %>%
  t() %>%
  as.data.frame() 

################################################################################
# Get the correlation stats for mosquitos

mosquito_corr_stats = sapply(1:nrow(mosquito_ratio), function(x) {
  
  row = as.numeric(mosquito_ratio[x, iv])
  correlation = cor.test(inf2NA(row), inf2NA(optimalities$mosquito_csc), method = 'spearman')
  return(c(corr = correlation$estimate[[1]], p.value = correlation$p.value))
  
}) %>%
  t() %>%
  as.data.frame()

################################################################################
# Get rid of the codon columns and bind the new statistics to the data frames

human_corr = human_ratio[-iv] %>%
  cbind(human_corr_stats)

mosquito_corr = mosquito_ratio[-iv] %>%
  cbind(mosquito_corr_stats)

################################################################################
# Write files

write.table(
  human_corr, 
  snakemake@output[['human_corr']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

write.table(
  mosquito_corr, 
  snakemake@output[['mosquito_corr']], 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

################################################################################