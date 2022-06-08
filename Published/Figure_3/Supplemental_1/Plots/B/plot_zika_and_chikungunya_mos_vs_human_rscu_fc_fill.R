################################################################################
# script to plot Figure 3 - figure supplement 1B
################################################################################

# load libraries 
library(tidyverse)
library(reshape2)
library(scales)
library(ggrepel)
library(ggpubr)

# load data 
opt = read.table(
  '../../Data/codon_optimalities.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, human_csc, mosquito_csc) %>%
  filter(!(codon %in% c('ATG', 'TGG', 'TAG', 'TAA', 'TGA')))

rscu_fc_mosquito = read.table(
  '../../Data/viral_rscu_fc_mosquito.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  filter(virus %in% c('Zika virus', 'Chikungunya virus')) %>%
  select(-c(realm, locus_tag, host, type, shape, corr, p.value, ATG, TGG, TAG, TAA, TGA))

rscu_fc_human = read.table(
  '../../Data/viral_rscu_fc_human.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  filter(virus %in% c('Zika virus', 'Chikungunya virus')) %>%
  select(-c(realm, locus_tag, host, type, shape, corr, p.value, ATG, TGG, TAG, TAA, TGA)) 

# melt rscu fc data
rscu_fc_mosquito = melt(rscu_fc_mosquito, id = 'virus')
names(rscu_fc_mosquito) = c('virus', 'codon', 'mosquito_fc')

rscu_fc_human = melt(rscu_fc_human, id = 'virus')
names(rscu_fc_human) = c('virus', 'codon', 'human_fc')

# merge rscu fc data 
rscu_fc = merge(rscu_fc_mosquito, rscu_fc_human, by = c('virus', 'codon'))

# merge optimality 
data = merge(opt, rscu_fc, by = 'codon')

# plot function 
make_plot = function(table, vir) {
  
  fig = table %>%
    filter(virus == vir) %>%
    ggplot() +
    geom_point(aes(x = human_csc - 0.0018, y = mosquito_csc, color = mosquito_fc), shape = '\u25D6', size = 4) + 
    geom_point(aes(x = human_csc + 0.0018, y = mosquito_csc, color = human_fc), shape = '\u25D7', size = 4) + 
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') + 
    geom_text_repel(aes(label = codon, x = human_csc, y = mosquito_csc)) +
    scale_color_gradientn(colors = c('plum4', 'plum4', 'white', 'springgreen4', 'springgreen4'), 
                          values = rescale(c(-2, -0.75, 0, 0.75, 2), c(0, 1)), breaks = c(-2, 0, 2), 
                          limits = c(-3, 3)) + 
    scale_y_continuous(limits = c(-.2, .2), breaks = c(-.2, -.1, 0, .1, .2)) +
    stat_cor(method = 'spearman', label.x = 0.025, label.y = -0.2, size = 7, aes(x = human_csc, y = mosquito_csc)) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 3), 
      axis.text.x = element_text(size = 15, face = 'bold'), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
      title = element_text(size = 20, face = 'bold'), 
      legend.key.size = unit(5, 'mm'), 
      legend.title = element_text(size = 10, face = 'bold', hjust = 0), 
      legend.text = element_text(size = 10, face = 'bold')
    ) + 
    labs(
      color = 'RSCU FC',
      x = 'Human CSC', 
      y = 'Mosquito CSC', 
      title = vir
    )
  
}

# make plots 
zika_fig = make_plot(data, 'Zika virus')
chikungunya_fig = make_plot(data, 'Chikungunya virus')

# save 
cairo_pdf(
  './zika_mos_vs_human_rscu_fc_fill.pdf', 
  family = 'Arial', 
  height = 7, 
  width = 10
)
zika_fig
dev.off()


cairo_pdf(
  './chikungunya_mos_vs_human_rscu_fc_fill.pdf', 
  family = 'Arial', 
  height = 7, 
  width = 10
)
chikungunya_fig
dev.off()

################################################################################