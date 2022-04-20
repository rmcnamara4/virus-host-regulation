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
# function to load the RSCU fold change data 

load_rscu_fc = function() {
  
  read.table('./data/rscu_fc.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to create Mean Relative Fitness plots per amino acid 

plot_mrf_per_aa = function(data, amino_acid, x_var, x_title) {
  
  temp = data %>%
    filter(aa == amino_acid) %>%
    ggplot(aes_string(x = x_var, y = 'mean_wrel')) +
    geom_point(aes(shape = set, fill = fitness_class, color = fitness_class), size = 4) +
    scale_fill_manual('Fitness Class', values = c('lethal' = 'black', 'deleterious' = 'purple',
                                                  'neutral' = 'darkgrey', 'beneficial' = 'yellow'), drop = FALSE) +
    scale_color_manual('Fitness Class', values = c('lethal' = 'black', 'deleterious' = 'purple', 
                                                   'neutral' = 'darkgrey', 'beneficial' = 'yellow'), drop = FALSE) +
    scale_shape_manual('Replicate', values = c(21, 24)) +
    stat_cor(aes(group = set), method = 'spearman', label.x = 0.1, label.y = c(2, 1.85)) +
    geom_text_repel(aes(label = codon), size = 4) +
    facet_grid(host ~ .) +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
    theme(panel.border = element_rect(color = 'black', fill = NA, size = 2), 
          axis.line = element_blank()) +
    labs(
      x = x_title, 
      y = 'Mean Relative Fitness', 
      title = amino_acid
    )
  
  if (x_var == 'csc') {
    temp = temp +
      coord_cartesian(ylim = c(0, 2), xlim = c(-0.19, 0.19))
  } else if (x_var == 'rscu_fc') {
    temp = temp + 
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
      coord_cartesian(ylim = c(0, 2), xlim = c(-3.1, 3.1))
  }
  
  print(temp)
  
}

################################################################################










