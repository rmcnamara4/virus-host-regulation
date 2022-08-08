################################################################################
# script to plot Figure 4 - supplemental figure 1B
################################################################################

# load libraries 
library(tidyverse)
library(scales)

# load data 
opt = read.table(
  '../../Data/human_endo_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mean_endo_csc)

codon_comp = read.table(
  '../../Data/dengue_strain_16681_codon_comp.csv', 
  sep = ',', header = TRUE
) %>%
  filter(!(codon %in% c('TAG', 'TGA', 'TAA')))

# filter for the top 16 and bottom 16 used codons 
codon_comp = codon_comp[order(-codon_comp$freq), ]

top_used = codon_comp[1:16, ] %>%
  mutate(
    group = 'Frequently Used'
  )

bottom_used = codon_comp[46:61, ] %>%
  mutate(
    group = 'Infrequently Used'
  )

data = rbind(top_used, bottom_used)

# merge data with opt
data = merge(data, opt, by = 'codon')

# get codons in decreasing frequency order 
ordered_codons = data[order(-data$freq), ]$codon

# plot 
fig = data %>%
  ggplot(aes(x = factor(codon, levels = ordered_codons), y = freq, fill = mean_endo_csc)) + 
  geom_bar(stat = "identity", color = 'black', size = .5) + 
  scale_fill_gradient2(low = '#3852A3', high = '#E81F27', mid = 'white', 
                       midpoint = 0, limits = c(-.2, .2), breaks = c(-.2, 0, .2)) + 
  facet_grid(~ factor(group, levels = c('Frequently Used', 'Infrequently Used')), 
             scales = 'free_x') + 
  scale_y_continuous(labels = percent, limits = c(0, 0.06)) + 
  geom_text(aes(label = codon, x = factor(codon, levels = ordered_codons), y = freq + 0.0025), 
            size = 6, angle = 90) + 
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 2), 
    axis.line = element_blank(),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.text = element_text(size = 10, face = 'bold'), 
    strip.text = element_text(face = 'bold', size = 12), 
    strip.background = element_rect(size = 2, color = 'black', fill = NA), 
    title = element_text(size = 15, face = 'bold')
  ) + 
  labs(
    y = 'Codon Frequency', 
    fill = 'Human CSC', 
    title = 'Dengue virus 2 strain 16681'
  )

# save
ggsave(
  './dengue_2_strain_16681_codon_frequency_labeled.pdf',
  fig, 
  height = 8, 
  width = 12
)
  
################################################################################
