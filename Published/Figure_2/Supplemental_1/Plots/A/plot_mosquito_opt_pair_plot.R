################################################################################
# script to plot Figure 2 - figure supplement 1A
################################################################################

# load libraries 
library(tidyverse)
library(GGally)

# load data 
mosquito_opt = read.table(
  '../../Data/mosquito_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(-c(codon, mosquito))

names(mosquito_opt) = c('DRB', 'Flavopiridol', 'Triptolide')

# make pair plot
pair_plot = ggpairs(mosquito_opt) + 
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = 'black', size = 2), 
        strip.text.x = element_text(face = 'bold', size = 13), 
        strip.background.x = element_rect(size = 1.5, color = 'black', fill = 'grey'), 
        strip.text.y = element_text(face = 'bold', size = 13), 
        strip.background.y = element_rect(size = 1.5, color = 'black', fill = 'grey'), 
        axis.text.x = element_text(size = 10, face = 'bold'), 
        axis.text.y = element_text(size = 10, face = 'bold'))

# save
ggsave(
  './mosquito_opt_pair_plot.pdf', 
  pair_plot, 
  width = 6, 
  height = 6
)

################################################################################