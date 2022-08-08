################################################################################
# script to plot Figure 1 - figure supplement 1C
################################################################################

# load libraries 
library(tidyverse)
library(reshape2)
library(scales)
library(grid)

# load data 
aa_fc_human = read.table(
  '../../Data/viral_aa_fc_human.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(-c(realm, locus_tag, host, type, shape)) %>%
  filter(virus %in% c('Dengue virus 2', 'Zika virus', 'Chikungunya virus')) 

# melt data 
melted_data = melt(aa_fc_human, id = 'virus')
names(melted_data) = c('virus', 'aa', 'aa_fc') 

# define order of amino acids 
ordered_aa = c('Trp', 'Met', 'Cys', 'His', 'Tyr', 'Asn', 'Phe', 'Ile', 'Asp', 'Gln', 'Thr', 'Arg',
               'Lys', 'Val', 'Pro', 'Gly', 'Ala', 'Glu', 'Ser', 'Leu')

# plot function 
plot_heatmap = function(table) {
  
  g = table %>%
    ggplot(aes(x = factor(aa, levels = ordered_aa), y = virus)) + 
    geom_tile(aes(fill = aa_fc), color = 'black', size = 0.05) + 
    facet_grid(. ~ factor(aa, levels = ordered_aa), scale = 'free_x', space = 'free_x', switch = 'x') +
    scale_fill_gradientn(colors = c('plum4', 'plum4', 'white', 'springgreen4', 'springgreen4'), 
                         values = rescale(c(-2, -1, 0, 1, 2), c(0, 1)), breaks = c(-2, 0, 2), 
                         limits = c(-2, 2)) +
    scale_y_discrete(labels = c('Zika virus' = 'Zika', 'Dengue virus 2' = 'Dengue 2',
                                'Chikungunya virus' = 'Chikungunya')) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(), 
      axis.text.y = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 13, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 13, face = 'bold', vjust = 1), 
      axis.ticks.x = element_blank(),
      title = element_text(size = 15, face = 'bold'), 
      panel.background = element_blank(), 
      strip.text = element_text(face = 'bold', size = 10), 
      strip.background = element_blank(), 
      legend.key.size = unit(5, 'mm'), 
      legend.title = element_text(size = 10, face = 'bold', hjust = 0), 
      legend.text = element_text(size = 10, face = 'bold')
    ) + 
    labs(
      x = 'Amino Acid', 
      y = 'Virus', 
      fill = 'H-relative\nAA Frequency',
      title = 'H-relative AA Frequency of common mosquito-borne viruses'
    )
  
  gt = ggplotGrob(g)
  
  ind = grep('null', gt$heights)
  gt$heights[ind] = unit(length(unique(table$virus)) + 0.2, 'null')
  
  gt$respect <- TRUE
  
  return(gt)
  
}

# make plot and save 
x = plot_heatmap(melted_data)
ggsave(
  './mosquito_borne_h-aa_hm.pdf', 
  height = 2, 
  width = 12
)
dev.off()

################################################################################
