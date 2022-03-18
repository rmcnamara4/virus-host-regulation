################################################################################
# Load data 

data = read.table('../complete_assembly.txt', 
                  sep = '\t', header = F, fill = TRUE, quote = "\"")

################################################################################
# Get which virus names are duplicated and keep only the reference one 

remove_duplicates = which(duplicated(data$V8))
keep = which(data[remove_duplicates, ]$V5 == 'reference genome')

mask = sample(TRUE, length(remove_duplicates), replace = TRUE)
mask[keep] = FALSE

data = data[-remove_duplicates[mask], ]

################################################################################
# Remove viruses that don't have a URL

iv = which(nchar(data$V20) < 20)
data = data[-iv, ]

################################################################################
# Write table

write.table(
  data, 
  '../complete_assembly.txt', 
  sep = '\t', 
  col.names = F, 
  row.names = F, 
  quote = F
)

################################################################################
