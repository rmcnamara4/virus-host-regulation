################################################################################
# script to define data loading functions 
################################################################################

# function to load mutation data 
load_p9_mutation_data = function() {
  
  read.table(
    '../point_mutations_p9.csv', 
    sep = ',', header = TRUE
  )
  
}

# function to load human optimality 
load_codon_optimality = function() {
  
  read.table(
    '../codon_optimalities.csv', 
    sep = ',', header = TRUE
  )
  
}

# function to load RSCU fold change data for Dengue 2 strain 16681
load_rscu_fc = function() {
  
  read.table(
    '../rscu_fc.csv', 
    sep = ',', header = TRUE
  )
  
}

# function to load combined optimality and RSCU fold change data 
load_opt_rscu_data = function() {
  
  read.table(
    '../opt_rscu_data.csv', 
    sep = ',', header = TRUE
  )
  
}

# function to load codon groups data 
load_codon_groups = function() {
  
  data = lapply(c('Human', 'Mosquito'), function(x) {
    readxl::read_excel(
      '../codon_groups.xlsx', 
      sheet = x
    )
  })
  
  names(data) = c('Human', 'Mosquito') 
  
  return(data)
  
}

# amino acid conversion hash 
aa_conversion_hash = hash::hash(
  c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'), 
  c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 
    'Ser', 'Thr', 'Trp', 'Tyr', 'Val', 'Stp')
)

# function to assign fitness class 
assign_fitness_class = function(table) {
  
  ifelse(table$mean_wrel_upper == 0, 'lethal', 
         ifelse(table$mean_wrel_upper > 0 & table$mean_wrel_upper < 1, 'deleterious', 
                ifelse(table$mean_wrel_lower > 1, 'beneficial', 'neutral')))
  
}


################################################################################