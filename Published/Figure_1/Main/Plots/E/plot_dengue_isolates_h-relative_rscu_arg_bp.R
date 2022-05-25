################################################################################
# script to plot Figure 1E
################################################################################

# load libraries 
library(tidyverse)
library(reshape2)

# load data 
dengue_isolates_fc = read.table(
  '../../Data/dengue_isolates_h-relative_rscu_fc.csv', 
  sep = ',', header = TRUE
)

# select only the codons that encode Arginine 
arg_codons = c('AGA', 'CGC', 'AGG', 'CGT', 'CGA', 'CGG')
dengue_isolates_fc_arg = dengue_isolates_fc[c('species', arg_codons)]

# get number of isolates for each serotype 
n_dengue_1 = sum(dengue_isolates_fc_arg$species == 'Dengue virus 1')
n_dengue_2 = sum(dengue_isolates_fc_arg$species == 'Dengue virus 2')
n_dengue_3 = sum(dengue_isolates_fc_arg$species == 'Dengue virus 3')
n_dengue_4 = sum(dengue_isolates_fc_arg$species == 'Dengue virus 4')

# melt data 
melted_data = melt(dengue_isolates_fc_arg, id = 'species')
names(melted_data) = c('species', 'codon', 'rscu_fc')

# plot 
plt = melted_data %>%
  ggplot(aes(x = factor(codon, levels = arg_codons), y = rscu_fc, fill = species)) + 
  geom_boxplot(outlier.shape = NA, size = .25) + 
  scale_fill_manual(
    values = c('Dengue virus 1' = '#009D45', 'Dengue virus 2' = '#CBA573', 'Dengue virus 3' = '#F2620F', 'Dengue virus 4' = '#32B3E9'), 
    labels = c('Dengue virus 1' = paste0('Dengue 1 Isolate\n(n = ', format(n_dengue_1, big.mark = ','), ')'), 
               'Dengue virus 2' = paste0('Dengue 2 Isolate\n(n = ', format(n_dengue_2, big.mark = ','), ')'), 
               'Dengue virus 3' = paste0('Dengue 3 Isolate\n(n = ', format(n_dengue_3, big.mark = ','), ')'), 
               'Dengue virus 4' = paste0('Dengue 4 Isolate\n(n = ', format(n_dengue_4, big.mark = ','), ')'))
  ) + 
  facet_grid(~ factor(codon, levels = arg_codons), scales = 'free_x') + 
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey', size = .5) + 
  coord_cartesian(ylim = c(-4.5, 2)) + 
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 1), 
    axis.line = element_blank(),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = 'bold'), 
    strip.text = element_text(face = 'bold', size = 12, margin = margin(3, 0, 3, 0, 'mm')), 
    strip.background = element_rect(size = 1, color = 'black', fill = NA), 
    title = element_text(size = 15, face = 'bold'), 
    legend.position = 'bottom'
  ) + 
  labs(
    title = 'Arginine',
    y = 'H-relative RSCU'
  )

# get p-value 
sig = aov(rscu_fc ~ codon, data = melted_data)

# save 
ggsave(
  './dengue_isolates_h-relative_rscu_arg_bp.pdf', 
  plt,
  height = 7, 
  width = 12
)

################################################################################