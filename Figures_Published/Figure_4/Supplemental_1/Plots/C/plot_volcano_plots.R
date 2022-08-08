################################################################################
# script to plot Figure 4 - supplemental figure 1C
################################################################################

# load libraries
library(tidyverse)
library(stringr)
library(gridExtra)

# load data
files = list.files(
  path = '../../Data', 
  pattern = 'uninfected*', 
  full.names = TRUE
)

file_names = str_replace(files, '../../Data/uninfected-', '')
file_names = str_replace(file_names, '.csv', '')

data = lapply(files, function(x) read.table(x, sep = ',', header = TRUE))
names(data) = file_names

# function to plot 
make_vol_plot = function(table, infection) {
  
  table$c = ifelse(
    table$logFC > 1 & table$padj < 0.01, 'Upregulated', 
    ifelse(
      table$logFC < -1 & table$padj < 0.01, 'Downregulated', 'none'
    )
  )
  
  fig = table %>%
    ggplot(aes(x = logFC, y = -10 * log10(padj), color = c)) + 
    geom_point() + 
    geom_hline(yintercept = 20, linetype = 'dashed', color = 'grey') + 
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'grey') + 
    scale_color_manual(
      values = c('Upregulated' = 'orange', 
                 'Downregulated' = 'purple', 
                 'none' = 'black'), 
      breaks = c('Downregulated', 'Upregulated')
    ) + 
    coord_cartesian(xlim = c(-4, 4), ylim = c(0, 450)) + 
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 3), 
      axis.text.x = element_text(size = 15, face = 'bold'), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
      title = element_text(size = 20, face = 'bold'), 
      legend.title = element_blank(), 
      legend.text = element_text(size = 15, face = 'bold'), 
      legend.position = 'top'
    ) + 
    labs(
      x = 'log2(FC)', 
      y = '-10 * log10(padj)', 
      title = paste0(infection, ' Infection')
    )
  
}

# plot 
low_vol = make_vol_plot(data$low, 'Low')
med_vol = make_vol_plot(data$med, 'Medium')
high_vol = make_vol_plot(data$high, 'High')

# arrange in grid 
final_plot = grid.arrange(low_vol, med_vol, high_vol, nrow = 1)

# save 
ggsave(
  './volcano_plots.pdf', 
  final_plot,
  width = 15, 
  height = 5
)

################################################################################