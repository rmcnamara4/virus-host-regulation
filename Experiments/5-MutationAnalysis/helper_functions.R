################################################################################
# SCRIPT TO DEFINE FUNCTIONS THAT ARE USED IN THE ANALYSIS 
################################################################################


################################################################################
# FUNCTIONS TO LOAD DATA 
################################################################################
# function to load the mutation data that is used to create the condensed version

load_mutation_data = function() {
  
  read.table('./data/point_mutations_p9.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to load the condensed mutation data that is used in the analysis 

load_condensed_mutation_data = function() {
  
  read.table('./data/condensed_point_mutations_p9.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to load the condensed synonymous mutation data that is used in the analysis 

load_condensed_synonymous_mutation_data = function() {
  
  read.table('./data/condensed_synonymous_point_mutations_p9.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to load the fitness class mutation percentages 

load_fitness_class_percentages_data = function() {
  
  read.table('./data/fitness_class_percentages_per_syn_mutation_p9.csv', 
             sep = ',', header = TRUE)
  
}

################################################################################
# function to load overall synonymous Mean Relative Fitness per species 

load_syn_mrf_per_species = function() {
  
  read.table('./data/syn_mut_mrf_per_species_p9.csv', 
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
# FUNCTIONS TO PLOT DATA
################################################################################
# function to create Mean Relative Fitness dot plots per amino acid 

plot_mrf_per_aa_dot = function(data, amino_acid, x_var, x_title) {
  
  temp = data %>%
    filter(mut_aa == amino_acid) %>%
    ggplot(aes_string(x = x_var, y = 'mean_wrel')) +
    geom_point(aes(shape = set, fill = fitness_class, color = fitness_class), size = 4) +
    scale_fill_manual('Fitness Class', values = c('lethal' = 'black', 'deleterious' = 'purple',
                                                  'neutral' = 'darkgrey', 'beneficial' = 'yellow'), drop = FALSE) +
    scale_color_manual('Fitness Class', values = c('lethal' = 'black', 'deleterious' = 'purple', 
                                                   'neutral' = 'darkgrey', 'beneficial' = 'yellow'), drop = FALSE) +
    scale_shape_manual('Replicate', values = c(21, 24)) +
    stat_cor(aes(group = set), method = 'spearman', label.x = 0.1, label.y = c(2, 1.85)) +
    geom_text_repel(aes(label = mutcodon), size = 4) +
    facet_grid(host ~ .) +
    theme_classic() +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
    theme(panel.border = element_rect(color = 'black', fill = NA, size = 2), 
          axis.line = element_blank()) +
    labs(
      x = x_title, 
      y = 'Mean Relative Fitness', 
      title = amino_acid
    )
  
  if (x_var == 'mutcodon_csc') {
    temp = temp +
      coord_cartesian(ylim = c(0, 2), xlim = c(-0.19, 0.19))
  } else if (x_var == 'mutcodon_rscu_fc') {
    temp = temp + 
      coord_cartesian(ylim = c(0, 2), xlim = c(-3.1, 3.1))
  }
  
  print(temp)
  
  return(temp)
  
}

################################################################################
# function to create Mean Relative Fitness dot plots for all codons together 

plot_mrf_all_dot = function(data, plot_vars, x_var, x_title) {
  
  temp = data %>%
    filter(host == plot_vars[1], set == plot_vars[2]) %>%
    ggplot(aes_string(x = x_var, y = 'mean_wrel')) +
    geom_point() +
    theme_classic() +
    geom_text_repel(aes(label = mutcodon)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
    labs(
      x = x_title, 
      y = 'Mean Relative Fitness', 
      title = paste0(plot_vars[1], ' ', plot_vars[2])
    )
  
  if (x_var == 'mutcodon_csc') {
    temp = temp + 
      stat_cor(method = 'spearman', label.x = -.25, label.y = 1.8) +
      coord_cartesian(xlim = c(-.25, .25), ylim = c(0.3, 2))
  } else if (x_var == 'mutcodon_rscu_fc') {
    temp = temp +
      stat_cor(method = 'spearman', label.x = -3, label.y = 1.8) +
      coord_cartesian(xlim = c(-3, 3), ylim = c(0.3, 2))
  }
  
  print(temp)
  
  return(temp)
  
}

################################################################################
# plot each replicates correlations for all amino acids as a dot plot 

plot_corr_all = function(data, host_var, color_title) {
  
  temp = data %>%
    dplyr::filter(host == host_var) %>%
    ggplot(aes(x = A, y = B)) +
    geom_point(aes(size = n_codons, color = range)) +
    scale_color_distiller(palette = 'YlOrBr', direction = 1) +
    geom_text_repel(aes(label = mut_aa), max.overlaps = 20) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(
      x = 'Rep A Spearman Cor', 
      y = 'Rep B Spearman Cor', 
      title = host_var, 
      color = color_title, 
      size = '# of Codons in AA'
    )
  
  print(temp)
  
  return(temp)
  
}

################################################################################
# plot barplots of synonymous mutations to an amino acid (split by codons) 

plot_mrf_per_aa_bar = function(data, amino_acid, fill_var, fill_title) {
  
  ordered_codons = read.table('../0-Preprocessing/0.5-GetHumanOptimality/data/ordered_human_codons.csv',
                              sep = ',', header = TRUE)$codon
  
  temp = data %>%
    filter(mut_aa == amino_acid) %>%
    ggplot(aes(x = factor(wtcodon, levels = ordered_codons))) +
    geom_bar(aes_string(y = 'mean_wrel', fill = fill_var), stat = "identity", color = 'black') +
    geom_errorbar(aes(ymin = mean_wrel, ymax = mean_wrel + sd_wrel), position = position_dodge(0.9), 
                  width = 0.5) +
    facet_grid(host ~ factor(mutcodon, levels = ordered_codons), scales = 'free_x', space = 'free_x') +
    geom_text(aes(label = n, y = mean_wrel + sd_wrel + 0.2), position = position_dodge(0.9)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
    theme_classic() +
    theme(panel.border = element_rect(color = 'black', fill = NA, size = 2), 
          axis.line = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = 'Mutated From', 
      y = 'Mean Relative Fitness', 
      title = amino_acid, 
      fill = fill_title
    )
  
  if (fill_var == 'delta_csc') {
    temp = temp +
      scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', 
                           midpoint = 0, limits = c(-0.3, 0.3))
  } else if (fill_var == 'delta_rscu_fc') {
    temp = temp + 
      scale_fill_gradient2(low = 'purple', high = 'darkgreen', mid = 'white', 
                           midpoint = 0, limits = c(-4.1, 4.1)) 
  }
  
  print(temp)
  
  return(temp)
  
}

################################################################################
# plot barplots of synonymous mutations to an amino acid (split by codons, filled 
# fitness class)

plot_status_per_aa_bar = function(data, amino_acid) {
  
  ordered_codons = read.table('../0-Preprocessing/0.5-GetHumanOptimality/data/ordered_human_codons.csv',
                              sep = ',', header = TRUE)$codon
  
  temp = data %>%
    filter(mut_aa == amino_acid) %>%
    ggplot(aes(x = factor(wtcodon, levels = ordered_codons))) +
    geom_bar(aes(y = perc, fill = fitness_class), stat = "identity", color = 'black') +
    scale_fill_manual(values = c('B' = 'gold2', 'N' = 'darkgrey', 'D' = 'purple', 'L' = 'black'), 
                      labels = c('B' = 'Beneficial', 'N' = 'Neutral', 'D' = 'Deleterious', 'L' = 'Lethal')) +
    facet_grid(host ~ factor(mutcodon, levels = ordered_codons), scales = 'free_x', space = 'free_x') +
    geom_text(aes(label = tot, y = 103)) +
    theme_classic() +
    theme(panel.border = element_rect(color = 'black', fill = NA, size = 2), 
          axis.line = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = 'Mutated From', 
      y = 'Percentage of Mutations', 
      fill = '', 
      title = amino_acid
    )
    
  print(temp)
  
  return(temp)
  
}




