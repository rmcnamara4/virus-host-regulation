################################################################################
# script to plot Figure 6F
################################################################################

# load libraries
library(tidyverse)
library(rstatix)
library(ggpubr)

# load data
data = read.table(
  '../../Data/syn_mutations_mrf_codon_groups_per_species_p9.csv', 
  sep = ',', header = TRUE
)

# calculate p-values between the codon group pairs
sig = data %>%
  group_by(host, group) %>%
  wilcox_test(mean_wrel ~ codon_class, paired = FALSE) %>%
  add_significance() %>%
  add_xy_position(fun = 'median_iqr', x = 'group')
sig$y.position = 1.45

# define function to get the sample size of the boxplots 
give.n = function(x) {
  return(c(y = 1.7, label = length(x))) 
}

# plot 
data$codon_class = factor(data$codon_class, 
                          levels = c('frequently_used', 'infrequently_used', 
                                     'preferentially_used', 'non_preferentially_used', 
                                     'denguenized', 'non_denguenized'))
fig = data %>%
  ggplot(aes(x = group, y = mean_wrel, fill = codon_class, group = codon_class)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual('Codon Group', values = c('frequently_used' = '#DD8505', 'infrequently_used' = '#B5025C', 
                                              'preferentially_used' = '#008A46', 'non_preferentially_used' = '#8A668A', 
                                              'denguenized' = '#CF3CD3', 'non_denguenized' = '#24C9C9'), 
                    labels = c('Frequently Used', 'Infrequently Used', 'Preferentially Used', 'Unpreferentially Used', 
                               'Denguenized', 'Non-Denguenized')) + 
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey') + 
  stat_pvalue_manual(sig, inherit.aes = FALSE, remove.bracket = TRUE) + 
  stat_summary(fun.data = give.n, geom = 'text', position = position_dodge(0.75)) + 
  facet_grid(host ~ .) + 
  theme_classic() + 
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    panel.border = element_rect(size = 2, color = 'black', fill = NA), 
    strip.background = element_rect(size = 2, color = 'black', fill = NA), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    title = element_text(size = 20, face = 'bold'), 
    legend.title = element_text(size = 18, face = 'bold'), 
    legend.text = element_text(size = 12, face = 'bold'), 
    strip.text = element_text(face = 'bold', size = 12, margin = margin(0, 3, 0, 3, 'mm'))
  ) + 
  coord_cartesian(ylim = c(.4, 1.7)) + 
  labs(
    y = 'Mean Relative Fitness', 
    title = 'Synonymous Mutations'
  )

# save
ggsave(
  './syn_mrf_codon_groups_box.pdf', 
  fig, 
  height = 8, 
  width = 12
)

################################################################################






