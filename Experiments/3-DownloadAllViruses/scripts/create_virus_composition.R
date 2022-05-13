################################################################################
# Load libraries 

library(dplyr)
library(readr)
library(stringr)

################################################################################
# Define function to parse the GBFF file and extract the virus name and the 
# virus composition (DNA, RNA, etc.) and shape (linear, circular)

parse_gbff = function(gbff) {
  
  txt = read_lines(gbff) 
  
  virus = txt[grep('ORGANISM', txt)] %>%
    str_replace('\\s*ORGANISM\\s*', '') 
  
  locus = txt[grep('^LOCUS', txt)] %>%
    str_replace('bp', '  ') %>%
    str_replace('[[:alpha:]]{3}\\s\\d+', '  ') %>%
    str_split("\\s{2,}")
  
  type_and_shape = sapply(locus, function(x) {
    return(c(x[4], x[5]))
  }) %>%
    t %>%
    as.data.frame
  
  data = cbind(virus, type_and_shape)
  names(data) = c('virus', 'type', 'shape')
  
  return(data)
  
}

################################################################################
# Run the function on the downloaded GBFF file 

virus_info = parse_gbff('./viruses.gbff.gz')

################################################################################
# Write table 

write.table(
  virus_info, 
  './virus_composition.txt', 
  sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE
)

################################################################################


