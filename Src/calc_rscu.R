calc_rscu = function(codon_counts_df) {
  
  if (class(codon_counts_df)[1] %in% c('data.frame', 'matrix', 'tibble')) flag = TRUE else flag = FALSE
  
  if (flag) codons = colnames(codon_counts_df) else codons = names(codon_counts_df)
  
  rscu = c()
  for (c in codons) {
    
    aa = names(synonymous_codons)[grep(c, synonymous_codons)]
    syn_codons = synonymous_codons[[aa]]
    n_codons = length(syn_codons)
    
    if (flag) {
      rscu_temp = n_codons * codon_counts_df[[c]] / rowSums(codon_counts_df[syn_codons])
    } else {
      rscu_temp = n_codons * codon_counts_df[[c]] / sum(codon_counts_df[syn_codons])
    }
    
    rscu = cbind(rscu, rscu_temp)
    
  }
  
  rscu = as.data.frame(rscu)
  names(rscu) = codons
  
  return(rscu)
  
}
