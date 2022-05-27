################################################################################
# script to plot Figure 1 - figure supplement 1D 
################################################################################

# load libraries 
library(tidyverse)
library(stringr)
library(reshape2)

# load data
rscu_fc_human = read.table(
  '../../Data/viral_rscu_fc_human.csv', 
  sep = ',', header = TRUE, quote = "\""
) %>%
  select(virus, host, type, shape, corr) %>%
  filter(grepl('human', host))

# make sure human:vertebrates is counted the same as vertebrates:human
rscu_fc_human$host = str_replace_all(rscu_fc_human$host, 'vertebrates: human', 'human:vertebrates')

# get n values function
get_n = function(var) {
  
  n_df = table(rscu_fc_human[[var]]) %>%
    melt %>%
    rename(
      !!sym(var) := Var1, 
      n = value
    ) %>%
    mutate(
      x = 0.4, 
      y = seq(from = 1, by = -0.07, length.out = length(.[[var]]))
    )
  
  return(n_df)
  
} 

# get n values 
host_n = get_n('host')
host_n[c(1, 2)] = host_n[c(1, 2, 5, 3, 4), c(1, 2)]

rscu_fc_human$type = str_extract(rscu_fc_human$type, 'RNA|DNA')
type_n = get_n('type')

# plot function 
make_plot = function(table, col_var) {
  
  plt = subset(table, !is.na(get(col_var))) %>%
    ggplot(aes(x = corr)) + 
    geom_density(size = 1.5, aes_string(y = '..scaled..', color = col_var)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_classic() + 
    theme(
      panel.background = element_rect(fill = NA, color = 'black', size = 3), 
      axis.text.x = element_text(size = 15, face = 'bold'), 
      axis.text.y = element_text(size = 15, face = 'bold'), 
      axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
      axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
      legend.title = element_text(size = 15, face = 'bold', ), 
      legend.text = element_text(size = 12, face = 'bold'), 
      legend.position = 'top'
    ) +
    labs(
      x = 'Spearman Correlation', 
      y = 'Density'
    )
  
}

# plot
host_plot = make_plot(rscu_fc_human, 'host')
host_plot = host_plot + 
  scale_color_manual(
    'Host',
    values = c('human' = '#750D37', 'human:invertebrates' = '#B3DEC1', 
               'human:vertebrates' = '#EFAAC4', 'human:invertebrates:vertebrates' = '#D66853', 
               'human:protozoa' = '#6B717E'), 
    labels = c('human' = 'Human', 'human:invertebrates' = 'Human & Invertebrates',
               'human:vertebrates' = 'Human & Vertebrates',
               'human:invertebrates:vertebrates' = 'Human & Invertebrates & Vertebrates', 
               'human:protozoa' = 'Human & Protozoa')
  ) +
  guides(color = guide_legend(nrow = 2, title.position = 'top')) +
  geom_text(data = host_n, aes(label = paste0('n = ', n), color = host, x = x, y = y), size = 5, 
            stat = 'identity', check_overlap = TRUE, show.legend = FALSE) 

type_plot = make_plot(rscu_fc_human, 'type')
type_plot = type_plot + 
  scale_color_manual(
    'Type', 
    values = c('DNA' = '#9A7AA0', 'RNA' = '#78A1BB')
  ) + 
  guides(color = guide_legend(title.position = 'top')) + 
  geom_text(data = type_n, aes(label = paste0('n = ', n), color = type, x = x, y = y), size = 5, 
            stat = 'identity', check_overlap = TRUE, show.legend = FALSE)

# save
ggsave(
  './all_human_virus_host_corr_density.pdf', 
  host_plot, 
  height = 8, 
  width = 9
)

ggsave(
  './all_human_virus_type_corr_density.pdf', 
  type_plot, 
  height = 8, 
  width = 9
)

################################################################################
