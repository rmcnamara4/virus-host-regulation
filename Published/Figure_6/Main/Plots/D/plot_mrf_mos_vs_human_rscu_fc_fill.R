################################################################################
# script to plot Figure 6D
################################################################################

# load libraries
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(scales)


# load data 
data = read.table(
  '../../Data/syn_mutations_mrf_per_species_p9.csv', 
  sep = ',', header = TRUE
) %>%
  select(host, mutcodon, mean_wrel, mutcodon_rscu_fc)

# clean and rearrange data 
human_data = data %>%
  filter(host == 'Human') %>%
  rename(
    human_mean_wrel = mean_wrel, 
    human_rscu_fc = mutcodon_rscu_fc
  ) %>%
  select(-host)

mosquito_data = data %>%
  filter(host == 'Mosquito') %>%
  rename(
    mosquito_mean_wrel = mean_wrel, 
    mosquito_rscu_fc = mutcodon_rscu_fc
  ) %>%
  select(-host)

plot_data = merge(human_data, mosquito_data, by = 'mutcodon')

# plot 
fig = plot_data %>%
  ggplot() + 
  geom_point(aes(x = human_mean_wrel - 0.0035, y = mosquito_mean_wrel, color = mosquito_rscu_fc), shape = '\u25D6', size = 4) + 
  geom_point(aes(x = human_mean_wrel + 0.0035, y = mosquito_mean_wrel, color = human_rscu_fc), shape = '\u25D7', size = 4) + 
  scale_color_gradientn(colors = c('plum4', 'plum4', 'white', 'springgreen4', 'springgreen4'), 
                        values = rescale(c(-2, -0.75, 0, 0.75, 2), c(0, 1)), breaks = c(-2, 0, 2), 
                        limits = c(-3, 3)) + 
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  geom_text_repel(aes(label = mutcodon, x = human_mean_wrel, y = mosquito_mean_wrel)) + 
  stat_cor(method = 'spearman', aes(x = human_mean_wrel, y = mosquito_mean_wrel)) + 
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
    x = 'Human Mean Relative Fitness', 
    y = 'Mosquito Mean Relative Fitness', 
    color = 'RSCU FC'
  )

# save
cairo_pdf(
  './mrf_mos_vs_human_rscu_fc_fill.pdf', 
  family = 'Arial', 
  height = 8, 
  width = 12
)
fig 
dev.off()

################################################################################
  