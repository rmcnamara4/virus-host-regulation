################################################################################
# script to plot Figure 1D
################################################################################

# load libraries 
library(tidyverse)
library(stringr)

# load data 
rscu_fc_human = read.table(
  '../../Data/viral_rscu_fc_human.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(virus, locus_tag, host, corr) %>%
  filter(grepl('human', host))

# get viruses to label 
labeled_vir = rscu_fc_human %>%
  filter(virus %in% c('Dengue virus 2', 'Zika virus', 'Chikungunya virus', 'Human alphaherpesvirus 1', 
                      'Human immunodeficiency virus 1', 'Influenza A virus (A/Puerto Rico/8/1934(H1N1))', 
                      'Severe acute respiratory syndrome coronavirus 2'))

# replace names
replacements = c('Dengue virus 2' = 'Dengue 2', 
                 'Zika virus' = 'Zika', 
                 'Chikungunya virus' = 'Chikungunya', 
                 'Human alphaherpesvirus 1' = 'Herpes', 
                 'Human immunodeficiency virus 1' = 'HIV', 
                 'Influenza A virus \\WA\\WPuerto Rico\\W8\\W1934\\WH1N1\\W\\W' = 'Flu', 
                 'Severe acute respiratory syndrome coronavirus 2' = 'Covid')
labeled_vir$virus = str_replace_all(string = labeled_vir$virus, pattern = replacements)


# plot 
plt = rscu_fc_human %>%
  ggplot(aes(x = corr)) + 
  geom_density(size = 1.5, aes(y = ..scaled..)) +
  geom_vline(data = labeled_vir, aes(xintercept = corr, color = virus), linetype = 'dashed', size = 1.1) + 
  geom_text(data = labeled_vir, aes(x = corr, label = virus, color = virus), y = .8, angle = 90, size = 5) +
  geom_text(x = 0.4, y = 1, label = paste0('n = ', nrow(rscu_fc_human)), size = 5) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual(values = c('Dengue 2' = '#C29800', 'Zika' = '#F961D5', 'Chikungunya' = '#F8766D', 
                                'Covid' = '#619CFF', 'Flu' = '#00BFC4', 'HIV' = '#9900AF', 'Herpes' = '#00BA38')) +
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.position = 'none'
  ) +
  labs(
    x = 'Spearman Correlation', 
    y = 'Density'
  )

# save
ggsave(
  './all_human_virus_corr_density.pdf', 
  plt, 
  height = 7, 
  width = 9
)
