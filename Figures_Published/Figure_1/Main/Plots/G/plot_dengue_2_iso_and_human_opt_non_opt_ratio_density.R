################################################################################
# script to plot Figure 1G
################################################################################

# load libraries 
library(tidyverse)

# load data 
codon_freq = read.table(
  '../../Data/dengue_2_iso_and_human_codon_freq.csv', 
  sep = ',', header = TRUE
) %>%
  select(species, ratio)

# remove infinites
iv = is.infinite(codon_freq$ratio)
codon_freq = codon_freq[!iv, ]

# get human quantiles 
quant = quantile(codon_freq[codon_freq$species == 'Human', ]$ratio)
quant = as.numeric(quant[2:4])

# get n for human and Dengue isolates 
n_human = sum(codon_freq$species == 'Human')
n_dengue = sum(codon_freq$species == 'Dengue virus 2')

# plot 
plt = codon_freq %>%
  ggplot(aes(x = ratio, color = species)) + 
  geom_density(aes(y = ..scaled..), size = 1.5) + 
  scale_color_manual('Species', values = c('Dengue virus 2' = '#C29800', 'Human' = 'black'), 
                     labels = c('Dengue virus 2' = 'Dengue 2 Isolates', 'Human' = 'Human')) +
  geom_vline(xintercept = quant, linetype = 'dashed', color = 'red', size = 1.25) +
  geom_text(label = '25%', x = quant[1] - .1, y = 0.5, angle = 90, inherit.aes = FALSE, color = 'red', size = 5, 
            stat = "identity", check_overlap = TRUE) +
  geom_text(label = '50%', x = quant[2] - .1, y = 0.5, angle = 90, inherit.aes = FALSE, color = 'red', size = 5, 
            stat = "identity", check_overlap = TRUE) + 
  geom_text(label = '75%', x = quant[3] - .1, y = 0.5, angle = 90, inherit.aes = FALSE, color = 'red', size = 5, 
            stat = "identity", check_overlap = TRUE) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', size = 1.25) + 
  geom_text(label = paste0('n = ', format(n_human, big.mark = ',')), x = 3, y = 1.05, inherit.aes = FALSE, size = 5, 
            color = 'black', stat = "identity", check_overlap = TRUE) +
  geom_text(label = paste0('n = ', format(n_dengue, big.mark = ',')), x = 3, y = 1, inherit.aes = FALSE, size = 5, 
            color = '#C29800', stat = "identity", check_overlap = TRUE) + 
  scale_x_continuous(expand = c(0, 0), limits = c(-2.25, 3.5)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, .1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.title = element_text(size = 15, face = 'bold'), 
    legend.text = element_text(size = 12, face = 'bold')
  ) + 
  labs(
    x = 'log2(optimal / non-optimal)', 
    y = 'Density'
  )
  
# save 
ggsave(
  './dengue_2_iso_and_human_opt_non_opt_ratio_density.pdf', 
  plt, 
  height = 7, 
  width = 9
)
  
################################################################################