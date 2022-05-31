################################################################################
# script to plot Figure 2B
################################################################################

# load libraries
library(tidyverse)
library(reshape2)
library(grid)

# source required files 
source('../../../../../Src/codon_to_aa.R')
source('../../../../../Src/synonymous_codons.R')

# load data 
mosquito_opt = read.table(
  '../../Data/mosquito_csc.csv', 
  sep = ',', header = TRUE
) 

all_opt = read.table(
  '../../Data/codon_optimalities.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, human_csc) %>%
  rename(
    human = human_csc
  )

# get codons in increasing optimality 
ordered_codons = mosquito_opt[order(mosquito_opt$mosquito), ]$codon

# merge data 
all_opt = merge(all_opt, mosquito_opt, by = 'codon') %>%
  select(-mosquito)

# melt data 
melted_data = melt(all_opt, id = 'codon')
names(melted_data) = c('codon', 'species', 'csc')

# make aa column 
melted_data$aa = codon_to_aa(melted_data$codon)

# distinguish between drug and human
melted_data$group = ifelse(melted_data$species == 'human', 'human', 'drug')

# define order of amino acids 
ordered_aa = c('Trp', 'Met', 'Cys', 'Tyr', 'His', 'Phe', 'Asn', 'Ile', 'Gln', 'Asp', 'Thr', 'Arg', 'Val', 'Lys',
               'Pro', 'Gly', 'Ala', 'Glu', 'Ser', 'Leu')

# plot function 
plot_heatmap = function(table) {
  
  g = table %>%
    ggplot(aes(x = factor(codon, levels = ordered_codons), y = species)) + 
    geom_tile(aes(fill = csc), color = 'black', size = 0.05) + 
    facet_grid(group ~ factor(aa, levels = ordered_aa), scale = 'free', space = 'free') +
    scale_fill_gradient2(low = '#3852A3', high = '#E81F27', mid = 'white', 
                         midpoint = 0, limits = c(-.2, .2), breaks = c(-.2, 0, .2)) +
    scale_y_discrete(labels = c('human' = 'Human', 'triptolide' = 'Triptolide', 
                                'drb' = 'DRB', 'flavopiridol' = 'Flavopiridol')) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, size = 10, face = 'bold', vjust = 0.5), 
      axis.text.y = element_text(size = 10, face = 'bold', angle = 20),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      title = element_blank(), 
      panel.background = element_blank(), 
      strip.text.x = element_text(face = 'bold', size = 10), 
      strip.background.x = element_rect(size = 0.5, color = 'black', fill = 'grey'),
      strip.background.y = element_blank(),
      strip.text.y = element_blank(),
      legend.key.size = unit(5, 'mm'), 
      legend.title = element_text(size = 10, face = 'bold', hjust = 0), 
      legend.text = element_text(size = 10, face = 'bold')
    ) + 
    labs(
      fill = 'CSC'
    )
  
  gt = ggplotGrob(g)
  
  ind = grep('null', gt$heights)
  gt$heights[ind[1]] = unit(3.2, 'null')
  gt$heights[ind[2]] = unit(1.2, 'null')
  
  gt$respect <- TRUE
  
  return(gt)
  
}

# plot and save
x = plot_heatmap(melted_data)
ggsave(
  './human_mos_drug_opt_hm.pdf', 
  x,
  height = 4, 
  width = 12
)
dev.off()

################################################################################
