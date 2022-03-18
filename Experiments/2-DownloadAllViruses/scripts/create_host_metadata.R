################################################################################
# Load libraries 

library(dplyr)
library(stringr)

################################################################################
# Load metadata

metadata = readLines(
  '../metadata.txt'
)
metadata = metadata[-c(1:2)]

################################################################################
# Get intervals of the realms of viruses and the realm names

intervals = c(grep('^\\S', metadata))
realm_names = metadata[intervals]

# Add a last interval to make the for loop easier 
intervals = c(intervals, length(metadata) + 1)

################################################################################
# Add the realm to the end of the string for each virus 

for (i in 2:length(intervals)) {
  
  metadata[(intervals[i - 1] + 1):(intervals[i] - 1)] = paste0(metadata[(intervals[i - 1] + 1):(intervals[i] - 1)], realm_names[i - 1])
  
}

################################################################################
# Get the entries that have a host specified 

iv = grep('algae|archaea|bacteria|eukaryotic algae|fungi|invertebrates|vertebrates|human|land plants|protozoa', 
     metadata)
host_metadata = metadata[iv]

################################################################################
# Split the strings on the tabs

host_metadata = str_split(host_metadata, '\t')

################################################################################
# Extract the virus name, host, and realm from the split strings

host_metadata_df = sapply(host_metadata, function(x) {
  return(c(x[1], x[9], x[12]))
}) %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(
    trimws
  ) %>%
  mutate_at(
    "V2", 
    str_replace_all, 
    ', ', 
    ':'
  )
names(host_metadata_df) = c('virus', 'host', 'realm')

################################################################################
# Write table

write.table(
  host_metadata_df, 
  '../host_metadata.txt', 
  sep = '\t', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

################################################################################
