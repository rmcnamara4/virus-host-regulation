################################################################################
# script to define functions that are used in the analysis 
################################################################################
# function to load the mutation data that is used in the analysis 

load_mutation_data = function() {
  
  read.table('./data/point_mutations_p9.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to load the codon optimalities 

load_codon_optimalities = function() {
  
  read.table('../2-ConcatenateOptimalities/data/codon_optimalities.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to 