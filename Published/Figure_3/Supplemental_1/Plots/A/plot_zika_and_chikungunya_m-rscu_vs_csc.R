################################################################################
# script to plot Figure 3 - figure supplement 1A
################################################################################

# load libraries
library(tidyverse)
library(reshape2)
library(ggpubr)

# load data 
rscu_fc_mosquito = read.table(
  '../../Data/viral_rscu_fc_mosquito.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(-c(realm, locus_tag, host, type, shape, corr, p.value, TAG, TGA, TAA, ATG, TGG)) %>%
  filter(virus %in% c('Zika virus', 'Chikungunya virus')) 

mosquito_opt = read.table(
  '../../Data/codon_optimalities.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mosquito_csc)

# melt the rscu fc data 
melted_data = melt(rscu_fc_mosquito, id = 'virus')
names(melted_data) = c('virus', 'codon', 'rscu_fc')

# merge with optimality 
melted_data = merge(melted_data, mosquito_opt, by = 'codon')

# plot function 
make_plot = function(table, vir, col) {
  
  fig = melted_data %>%
    filter(virus == vir) %>%
    ggplot(aes(x = mosquito_csc, y = rscu_fc)) + 
    geom_point(color = col, size = 4) + 
    stat_cor(method = 'spearman', label.x = 0.01, label.y = 2.3, size = 7, color = col) + 
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') + 
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') + 
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 3), 
      axis.text.x = element_text(size = 15, face = 'bold'), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1),
      title = element_text(size = 20, face = 'bold') 
    ) + 
    coord_cartesian(ylim = c(-2.75, 2.75), xlim = c(-0.18, 0.18)) + 
    labs(
      x = 'Codon Stability Coefficient (CSC)', 
      y = 'M-relative RSCU', 
      title = vir
    )
  
  return(fig)
  
}

# make plots 
zika_fig = make_plot(melted_data, 'Zika virus', '#F961D5')
chikungunya_fig = make_plot(melted_data, 'Chikungunya virus', '#F8766D')

# save 
ggsave(
  './zika_m-rscu_vs_csc.pdf', 
  zika_fig, 
  height = 7, 
  width = 7
)

ggsave(
  './chikungunya_m-rscu_vs_csc.pdf', 
  chikungunya_fig, 
  height = 7, 
  width = 7
)

################################################################################