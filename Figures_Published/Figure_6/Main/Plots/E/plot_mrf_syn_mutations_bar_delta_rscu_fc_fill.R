################################################################################
# script to plot Figure 6E
################################################################################

# load libraries
library(tidyverse)
library(scales)

# load data 
data = read.table(
  '../../Data/syn_mutations_per_species_p9.csv', 
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
    geom_bar(aes(y = mean_wrel, fill = delta_rscu_fc), stat = "identity", color = 'black') + 
    scale_fill_gradientn(colors = c('plum4', 'plum4', 'white', 'springgreen4', 'springgreen4'), 
                          values = rescale(c(-2, -0.75, 0, 0.75, 2), c(0, 1)), breaks = c(-2, 0, 2), 
                          limits = c(-4, 4)) + 
    geom_errorbar(aes(ymin = mean_wrel, ymax = mean_wrel + sd_wrel), position = position_dodge(0.9), 
                  width = 0.5) + 
    facet_grid(host ~ factor(mutcodon, levels = ordered_codons), scales = 'free_x', space = 'free_x') + 
    geom_text(aes(label = n, y = mean_wrel + sd_wrel + 0.2), position = position_dodge(0.9), size = 4) + 
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') + 
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 2), 
      axis.text.x = element_text(size = 15, face = 'bold', angle = 45, hjust = 1), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1),
      title = element_text(size = 20, face = 'bold'), 
      legend.title = element_text(size = 18, face = 'bold'), 
      legend.text = element_text(size = 15, face = 'bold'), 
      strip.background = element_rect(fill = NA, color = 'black', size = 2), 
      strip.text.x = element_text(face = 'bold', size = 12, margin = margin(3, 0, 3, 0, 'mm')), 
      strip.text.y = element_text(face = 'bold', size = 12, margin = margin(0, 3, 0, 3, 'mm'))
    ) + 
    labs(
      x = 'Mutated From', 
      y = 'log2 Mean Relative Fitness', 
      title = title,  
      fill = 'Delta RSCU FC'
    )
    
}

# plot 
fig = make_plot(data, 'Arg', 'Arginine')

# save
ggsave(
  './mrf_syn_mutations_bar_delta_rscu_fc_fill.pdf', 
  fig, 
  height = 8, 
  width = 12
)

################################################################################







