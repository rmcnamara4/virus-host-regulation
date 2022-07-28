################################################################################
# script to plot Figure 6 - figure supplement 1A
################################################################################

# load libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)

# load data 
data = read.table(
  '../../Data/syn_mutations_mrf_per_rep_p9.csv', 
  sep = ',', header = TRUE
)

# function for plotting 
make_plot = function(table, species, rep) {
  
  labeled_points = table %>%
    filter(host == species, set == rep, mut_aa == 'Arg') 
  
  fig = table %>%
    filter(host == species, set == rep) %>%
    ggplot(aes(x = mutcodon_rscu_fc, y = mean_wrel)) +
    geom_point(size = 3) + 
    scale_fill_gradient2(low = '#3852A3', high = '#E81F27', mid = 'white', 
                         midpoint = 0, limits = c(-.2, .2), breaks = c(-.2, 0, .2)) + 
    geom_label_repel(data = labeled_points, aes(label = mutcodon, fill = mutcodon_csc)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
    stat_cor(method = 'spearman', label.x = -3, label.y = 1.8, size = 4) + 
    coord_cartesian(xlim = c(-3, 2.5), ylim = c(.3, 2)) +
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
      x = 'log2 RSCU FC', 
      y = 'Mean Relative Fitness', 
      title = paste0(species, ' ', rep), 
      fill = paste0(species, ' ', 'CSC')
    )
  
}

# create plots 
human_b = make_plot(data, 'Human', 'B') 
mosquito_a = make_plot(data, 'Mosquito', 'A')

# save 
ggsave(
  './syn_mrf_human_b.pdf', 
  human_b, 
  height = 8, 
  width = 12
)

ggsave(
  './syn_mrf_mosquito_a.pdf', 
  mosquito_a, 
  height = 8, 
  width = 12
)

################################################################################