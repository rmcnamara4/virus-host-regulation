################################################################################
# script to plot Figure 3D
################################################################################

# load libraries 
library(tidyverse)
library(reshape2)
library(scales)
library(grid)

# source required files
source('../../../../../Src/codon_to_aa.R')
source('../../../../../Src/synonymous_codons.R')

# load data 
dengue_isolates_fc = read.table(
  '../../Data/dengue_isolates_m-relative_rscu_fc.csv', 
  sep = ',', header = TRUE
) %>%
  select(-c(TAG, TGA, TAA, TGG, ATG))

mosquito_opt = read.table(
  '../../Data/codon_optimalities.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mosquito_csc)

# get median RSCU fold change for each codon for each serotype 
dengue_isolates_median_fc = dengue_isolates_fc %>%
  group_by(species) %>%
  summarize_all(
    median
  )

# melt data 
melted_data = melt(dengue_isolates_median_fc, id = 'species') 
names(melted_data) = c('species', 'codon', 'rscu_fc')

# add amino acid column 
melted_data$aa = codon_to_aa(melted_data$codon)

# get codons in increasing optimality 
ordered_codons = mosquito_opt[order(mosquito_opt$mosquito_csc), ]$codon

# define order of amion acids 
ordered_aa = c('Cys', 'His', 'Tyr', 'Phe', 'Gln', 'Asn', 'Pro', 'Ile', 'Arg', 'Asp', 'Thr', 'Gly', 'Lys', 
               'Val', 'Glu', 'Ala', 'Ser', 'Leu')

# plot function 
plot_heatmap = function(table) {
  
  g = table %>%
    ggplot(aes(x = factor(codon, levels = ordered_codons), y = species)) + 
    geom_tile(aes(fill = rscu_fc), color = 'black', size = 0.05) + 
    facet_grid(. ~ factor(aa, levels = ordered_aa), scale = 'free_x', space = 'free_x') +
    scale_fill_gradientn(colors = c('plum4', 'plum4', 'white', 'springgreen4', 'springgreen4'), 
                         values = rescale(c(-2, -0.75, 0, 0.75, 2), c(0, 1)), breaks = c(-2, 0, 2), 
                         limits = c(-3.3, 3.3)) +
    scale_y_discrete(labels = c('Dengue virus 1' = 'Dengue 1', 'Dengue virus 2' = 'Dengue 2', 
                                'Dengue virus 3' = 'Dengue 3', 'Dengue virus 4' = 'Dengue 4')) +
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
      y = 'Isolate', 
      fill = 'M-relative\nRSCU',
      title = 'M-relative RSCU of Dengue isolates'
    )
  
  gt = ggplotGrob(g)
  
  ind = grep('null', gt$heights)
  gt$heights[ind] = unit(length(unique(table$species)) + 0.2, 'null')
  
  gt$respect <- TRUE
  
  return(gt)
  
}

# make plot and save 
x = plot_heatmap(melted_data)
ggsave(
  './dengue_isolates_m-rscu_hm.pdf', 
  x, 
  height = 3, 
  width = 12
)
dev.off()

################################################################################
