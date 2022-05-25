################################################################################
# script to plot Figure 1C
################################################################################

# load libraries 
library(tidyverse)
library(reshape2)
library(scales)
library(grid)

# load data 
rscu_fc_human = read.table(
  '../../Data/viral_rscu_fc_human.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(-c(realm, locus_tag, host, type, shape, corr, p.value, TAG, TGA, TAA, ATG, TGG)) 

human_opt = read.table(
  '../../Data/human_endo_csc.csv', 
  sep = ',', header = TRUE
)

# filter for Zika, Dengue 2, and Chikungunya virus 
rscu_fc_human = rscu_fc_human %>%
  filter(virus %in% c('Dengue virus 2', 'Zika virus', 'Chikungunya virus'))

# melt rscu_fc_human
melted_data = melt(rscu_fc_human, id = 'virus')
names(melted_data) = c('virus', 'codon', 'rscu_fc')

# add amino acid column 
iv = match(melted_data$codon, human_opt$codon)
melted_data$aa = human_opt$aa[iv]

# get codons in increasing optimality 
ordered_codons = human_opt[order(human_opt$mean_endo_csc), ]$codon

# define order of amino acids 
ordered_aa = c('Cys', 'Tyr', 'His', 'Phe', 'Asn', 'Ile', 'Gln', 'Asp', 'Thr', 'Arg', 'Val', 'Lys',
               'Pro', 'Gly', 'Ala', 'Glu', 'Ser', 'Leu')

# plot function 
plot_heatmap = function(table) {
  
  g = melted_data %>%
    ggplot(aes(x = factor(codon, levels = ordered_codons), y = virus)) + 
    geom_tile(aes(fill = rscu_fc), color = 'black', size = 0.05) + 
    facet_grid(. ~ factor(aa, levels = ordered_aa), scale = 'free_x', space = 'free_x') +
    scale_fill_gradientn(colors = c('plum4', 'plum4', 'white', 'springgreen4', 'springgreen4'), 
                         values = rescale(c(-2, -0.75, 0, 0.75, 2), c(0, 1)), breaks = c(-2, 0, 2), 
                         limits = c(-3, 3)) +
    scale_y_discrete(labels = c('Zika virus' = 'Zika', 'Dengue virus 2' = 'Dengue 2',
                                'Chikungunya virus' = 'Chikungunya')) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, size = 10, face = 'bold', vjust = 0.5), 
      axis.text.y = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 13, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 13, face = 'bold', vjust = 1), 
      title = element_text(size = 15, face = 'bold'), 
      panel.background = element_blank(), 
      strip.text = element_text(face = 'bold', size = 10), 
      strip.background = element_rect(size = 0.5, color = 'black', fill = 'grey'), 
      legend.key.size = unit(5, 'mm'), 
      legend.title = element_text(size = 10, face = 'bold', hjust = 0), 
      legend.text = element_text(size = 10, face = 'bold')
    ) + 
    labs(
      x = 'Codon', 
      y = 'Virus', 
      fill = 'H-relative\nRSCU',
      title = 'H-relative RSCU of common mosquito-borne viruses'
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
  './mosquito_borne_h-rscu_hm.pdf', 
  x, 
  height = 3, 
  width = 12
)
dev.off()
  
  
  