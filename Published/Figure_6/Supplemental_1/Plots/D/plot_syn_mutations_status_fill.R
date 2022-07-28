################################################################################
# script to plot Figure 6 - figure supplement 1D
################################################################################

# load libraries
library(tidyverse)

# load data 
data = read.table(
  '../../Data/syn_mutations_fitness_class_percent_per_species_p9.csv', 
  sep = ',', header = TRUE
)

opt = read.table(
  '../../Data/codon_optimalities.csv', 
  sep = ',', header = TRUE
)

# get codons in increasing human optimality 
ordered_codons = opt[order(opt$human_csc), ]$codon

# plot function 
make_plot = function(table, amino_acid, title) {
  
  fig = table %>%
    filter(mut_aa == amino_acid) %>%
    ggplot(aes(x = factor(wtcodon, levels = ordered_codons))) + 
    geom_bar(aes(y = perc, fill = fitness_class), stat = "identity", color = 'black') + 
    scale_fill_manual(
      values = c('beneficial' = 'gold2', 'neutral' = 'darkgrey', 'deleterious' = 'purple', 'lethal' = 'black'), 
      labels = c('beneficial' = 'Beneficial', 'neutral' = 'Neutral', 'deleterious' = 'Deleterious', 'lethal' = 'Lethal')
    ) + 
    facet_grid(host ~ factor(mutcodon, levels = ordered_codons), scales = 'free_x', space = 'free_x') + 
    geom_text(aes(label = tot), y = 104, size = 5) + 
    coord_cartesian(ylim = c(0, 105)) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 2), 
      axis.text.x = element_text(size = 15, face = 'bold', angle = 45, hjust = 1), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1),
      title = element_text(size = 20, face = 'bold'), 
      legend.title = element_blank(), 
      legend.text = element_text(size = 15, face = 'bold'), 
      strip.background = element_rect(fill = NA, color = 'black', size = 2), 
      strip.text.x = element_text(face = 'bold', size = 12, margin = margin(3, 0, 3, 0, 'mm')), 
      strip.text.y = element_text(face = 'bold', size = 12, margin = margin(0, 3, 0, 3, 'mm'))
    ) + 
    labs(
      x = 'Mutated From', 
      y = 'Percentage of Mutations', 
      title = title
    )
  
}

# create plot 
fig = make_plot(data, 'Arg', 'Arginine')

# save 
ggsave(
  './syn_mutations_status_fill.pdf', 
  fig, 
  height = 8, 
  width = 12
)

################################################################################