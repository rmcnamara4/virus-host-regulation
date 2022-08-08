################################################################################
# script to plot Figure 2 - figure supplement 1D
################################################################################

# load libraries
library(tidyverse)

# load data 
mosquito_codons = read.table(
  '../../Data/mosquito_codon_frequency_ratio.csv', 
  sep = ',', header = TRUE
) %>%
  select(gene_ID, ratio)

egfp_reporters_codons = read.table(
  '../../Data/egfp_reporters_codon_frequency_ratio.csv', 
  sep = ',', header = TRUE
) %>%
  select(name, ratio)

# get number of mosquito genes shown in plot 
n = nrow(mosquito_codons)

# define my colors 
my_colors = c(
  'EGFP' = 'dodgerblue2', 
  'EGFP_1' = '#E31A1C', 
  'EGFP_2' = 'green4', 
  'EGFP_3' = '#6A3D9A', 
  'EGFP_4' = '#FF7F00', 
  'EGFP_5' = 'brown',
  'EGFP_6' = 'skyblue2', 
  'EGFP_7' = '#FB9A99', 
  'EGFP_8' = 'palegreen2', 
  'EGFP_9' = '#CAB2D6', 
  'EGFP_10' = '#FDBF6F', 
  'EGFP_11' = 'deeppink1', 
  'EGFP_12' = 'maroon', 
  'EGFP_13' = 'orchid1'
)

# define order of EGFPs
my_order = c('EGFP', paste0('EGFP_', 1:13))

# plot 
plt = mosquito_codons %>%
  ggplot(aes(x = ratio)) + 
  geom_density(size = 1.5, aes(y = ..scaled..)) + 
  geom_vline(data = egfp_reporters_codons, aes(xintercept = ratio, color = factor(name, levels = my_order)), 
             linetype = 'dashed', size = 1) +
  geom_text(x = 4, y = 1, label = paste0('n = ', format(n, big.mark = ',')), size = 5, stat = "identity", check_overlap = TRUE) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = my_colors) +
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.title = element_blank(), 
    legend.text = element_text(size = 10, face = 'bold')
  ) +
  labs(
    x = 'log2(optimal / non-optimal)', 
    y = 'Density'
  ) 

# save
ggsave(
  './egfp_reporters_opt_ratio_density.pdf', 
  plt, 
  width = 8, 
  height = 6
)

################################################################################
