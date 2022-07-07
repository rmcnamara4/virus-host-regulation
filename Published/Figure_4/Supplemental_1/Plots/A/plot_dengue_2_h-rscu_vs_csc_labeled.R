################################################################################
# script to plot Figure 4 - supplemental figure 1A
################################################################################

# load libraries 
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggrepel)

# load data 
rscu_fc = read.table(
  '../../Data/rscu_fc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, human_rscu_fc) %>%
  filter(!(codon %in% c('TAG', 'TGA', 'TAA', 'ATG', 'TGG')))

human_opt = read.table(
  '../../Data/human_endo_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mean_endo_csc)

# create table of codons to be labeled 
rscu_fc$human_rscu_fc = round(rscu_fc$human_rscu_fc, digits = 2)

labeled_codons = rscu_fc %>%
  filter(human_rscu_fc >= 0.5 | human_rscu_fc <= -0.5)

# merge rscu fc with optimality 
data = merge(rscu_fc, human_opt, by = 'codon')
labeled_codons = merge(labeled_codons, human_opt, by = 'codon')

# plot 
fig = data %>%
  ggplot(aes(x = mean_endo_csc, y = human_rscu_fc)) + 
  geom_point(color = '#C29800', size = 4) + 
  stat_cor(method = 'spearman', label.x = 0.01, label.y = 2.3, size = 7, color = '#C29800') + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') + 
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') + 
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'red') + 
  geom_hline(yintercept = -0.5, linetype = 'dashed', color = 'red') +
  geom_label_repel(data = labeled_codons, aes(x = mean_endo_csc, y = human_rscu_fc, fill = mean_endo_csc, 
                                              label = codon), force_pull = 3) + 
  scale_fill_gradient2(low = '#3852A3', high = '#E81F27', mid = 'white', 
                       midpoint = 0, limits = c(-.2, .2), breaks = c(-.2, 0, .2)) + 
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    title = element_text(size = 20, face = 'bold'), 
    legend.key.size = unit(5, 'mm'), 
    legend.title = element_text(size = 15, face = 'bold', hjust = 0), 
    legend.text = element_text(size = 15, face = 'bold')
  ) + 
  coord_cartesian(ylim = c(-2.51, 2.51), xlim = c(-0.175, 0.175)) + 
  labs(
    x = 'Codon Stability Coefficient (CSC)', 
    y = 'H-relative RSCU', 
    fill = 'CSC',
    title = 'Dengue virus 2'
  )

# save
ggsave(
  './dengue_2_h-rscu_vs_csc_labeled.pdf', 
  fig,
  height = 8, 
  width = 8
)

################################################################################







