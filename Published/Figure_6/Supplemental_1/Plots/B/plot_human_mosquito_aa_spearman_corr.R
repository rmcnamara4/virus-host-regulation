################################################################################
# script to plot Figure 6 - figure supplement 1B
################################################################################

# load libraries
library(tidyverse)
library(ggrepel)

# load data 
mut_data = read.table(
  '../../Data/syn_mutations_mrf_per_rep_p9.csv', 
  sep = ',', header = TRUE
)

# get correlations against RSCU FC for each AA and each replicate 
rscu_fc_corr = mut_data %>%
  group_by(host, set, mut_aa) %>%
  summarize(
    correlation = cor(mutcodon_rscu_fc, mean_wrel, method = 'spearman'), 
    range = max(mutcodon_rscu_fc) - min(mutcodon_rscu_fc), 
    n_codons = n()
  )

# split the set variable into two columns, A and B
rscu_fc_corr = spread(rscu_fc_corr, set, correlation)

# plot function 
make_plot = function(data, host_var, color_title) {
  
  fig = data %>%
    filter(host == host_var) %>%
    ggplot(aes(x = A, y = B)) + 
    geom_point(aes(size = n_codons, color = range)) + 
    scale_color_distiller(palette = 'YlOrBr', direction = 1, 
                          limits = c(0, 5.5), breaks = c(0, 2.5, 5)) + 
    geom_text_repel(aes(label = mut_aa), max.overlaps = 20, size = 5) + 
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') + 
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 3), 
      axis.text.x = element_text(size = 15, face = 'bold'), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1),
      title = element_text(size = 20, face = 'bold'), 
      legend.title = element_text(size = 18, face = 'bold'), 
      legend.text = element_text(size = 15, face = 'bold')
    ) + 
    labs(
      x = 'Rep A Spearman Cor', 
      y = 'Rep B Spearman Cor', 
      title = host_var, 
      color = color_title, 
      size = '# of Synonymous Codons'
    )
  
} 


# create plot 
human_corr = make_plot(rscu_fc_corr, 'Human', 'AA RSCU FC Range')
mosquito_corr = make_plot(rscu_fc_corr, 'Mosquito', 'AA RSCU FC Range')

# save
ggsave(
  './human_aa_spearman_corr.pdf', 
  human_corr, 
  height = 8, 
  width = 12
)

ggsave(
  './mosquito_aa_spearman_corr.pdf', 
  mosquito_corr, 
  height = 8, 
  width = 12
)

################################################################################
