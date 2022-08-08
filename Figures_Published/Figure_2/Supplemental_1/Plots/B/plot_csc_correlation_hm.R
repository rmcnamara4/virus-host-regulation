################################################################################
# script to plot Figure 2 - figure supplement 1B
################################################################################

# load libraries 
library(tidyverse)
library(stringr)
library(Hmisc)
library(reshape2)
library(RColorBrewer)

# define functions
# function to flatten correlation matrix 
flatten_cormat <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut],
    p = pmat[ut]
  )
}

# function to get upper triangle of correlation matrix 
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# function to order correlation matrix using correlation as distance 
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# load data 
opts = read.table(
  '../../Data/codon_optimalities.csv', 
  sep = ',', header = TRUE
) %>%
  select(-codon)

names(opts) = str_replace(names(opts), '_csc', '')

# get correlation matrix between the different optimalities 
res = rcorr(as.matrix(opts), type = 'spearman')

cor_vals = res$r
p_vals = res$P

# reorder correlation matrix 
cor_vals = reorder_cormat(cor_vals)

# get upper triangle of correlation matrix 
cor_vals = get_upper_tri(cor_vals)

# get upper triangle of p-value matrix
p_vals = get_upper_tri(p_vals)

# melt the correlation matrix and p-values
cor_vals_melted = melt(cor_vals, na.rm = TRUE)
names(cor_vals_melted)[3] = 'corr'

p_vals_melted = melt(p_vals, na.rm = TRUE)
names(p_vals_melted)[3] = 'p_val'

# merge correlations and p-values 
all_data = full_join(cor_vals_melted, p_vals_melted, by = c('Var1', 'Var2'))

# add p-value column

# replace p-values with stars
all_data$p_val = symnum(all_data$p_val, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', ' '), abbr.colnames = FALSE)
all_data$p_val[all_data$p_val == '?'] = '***'

# remove na rows 
iv = apply(is.na(all_data), 1, function(x) any(x))
all_data = all_data[!iv, ]
all_data$p_val = as.character(all_data$p_val)

# plot 
plt = all_data %>%
  ggplot(aes(x = factor(Var2), y = factor(Var1))) +
  geom_tile(color = 'black', aes(fill = cut(corr, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)))) +
  scale_fill_manual('Spearman\nCorrelation', values = brewer.pal(5, 'Greens')) +
  scale_x_discrete(labels = str_to_title) + 
  scale_y_discrete(labels = str_to_title) +
  theme_classic() +
  coord_fixed() +
  geom_text(aes(x = Var2, y = Var1, label = p_val), color = 'white', size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, face = 'bold'), 
    axis.text.y = element_text(size = 10, face = 'bold'), 
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(), 
    legend.key.size = unit(5, 'mm'), 
    legend.title = element_text(size = 10, face = 'bold', hjust = 0), 
    legend.text = element_text(size = 10, face = 'bold')) +
  guides(fill = guide_legend(reverse=TRUE))

# save 
ggsave(
  './csc_correlation_hm.pdf', 
  plt, 
  width = 6, 
  height = 6
)

################################################################################
