################################################################################
# Load required files 

human_opt = read.table('../0-Preprocessing/0.5-GetHumanOptimality/data/human_endo_csc.csv', 
                       sep = ',', header = TRUE)
other_opts = read.table('../0-Preprocessing/0.6-OtherOptimalities/data/other_optimalities.csv', 
                          sep = ',', header = TRUE)
mosquito_opt = read.table('../1-CalculateMosquitoOptimality/data/mosquito_csc.csv', 
                          sep = ',', header = TRUE)

################################################################################
# Clean data.frames 

human_opt = human_opt[c(1, 6)]
names(human_opt) = c('codon', 'human_csc')

mosquito_opt = mosquito_opt[c(1, 5)]
names(mosquito_opt) = c('codon', 'mosquito_csc')

################################################################################
# Merge all optimalities 

df_list = list(human_opt, mosquito_opt, other_opts)
all_opt = Reduce(function(x, y) merge(x, y, by = 'codon'), df_list)

################################################################################
# Write file 

write.table(
  all_opt, 
  './data/codon_optimalities.csv', 
  sep = ',', row.names = F, col.names = T, quote = F
)

################################################################################

