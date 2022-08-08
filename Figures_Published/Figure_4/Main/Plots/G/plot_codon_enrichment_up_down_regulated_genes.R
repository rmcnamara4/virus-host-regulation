################################################################################
# script to plot Figure 4G
################################################################################

# load libraries
library(tidyverse) 
library(readxl)
library(ggpubr)

# load data 
diff_expression = read.table(
  '../../Data/uninfected-high.csv', 
  sep = ',', header = TRUE
) %>%
  select(gene_ID, logFC, padj)

human_codons = read_excel(
  '../../Data/human_hg38_stats.xlsx', 
  sheet = 'Codons'
) %>%
  select(-c(species, len, TAG, TAA, TGA))

human_opt = read.table(
  '../../Data/human_endo_csc.csv', 
  sep = ',', header = TRUE
) %>%
  select(codon, mean_endo_csc)

# filter for significantly upregulated and downregulated genes
diff_expression$group = ifelse(diff_expression$logFC > 1 & diff_expression$padj < 0.01, 'upregulated', 
                               ifelse(diff_expression$logFC < -1 & diff_expression$padj < 0.01, 'downregulated', 'none'))
diff_expression = diff_expression %>%
  filter(group != 'none') %>%
  select(gene_ID, group)

# merge with human codons
all_data = merge(human_codons, diff_expression, by = 'gene_ID') 
  
# get median codon composition for upregulated and downregulated genes
all_data_median = all_data %>%
  select(-gene_ID) %>%
  group_by(group) %>%
  summarize_all(
    median
  )

# get the fold change between the median frequency of the upregulated genes and the median frequency 
# of the downregulated genes 
all_data_fc = log2(as.numeric(all_data_median[2, -1]) / as.numeric(all_data_median[1, -1]))

# create data.frame with codon, fold change, and optimality 
fc_df = data.frame(
  codon = names(all_data_median)[-1], 
  enrichment = all_data_fc
)

fc_df = merge(fc_df, human_opt, by = 'codon')
names(fc_df)[3] = 'opt'

# get codons ordered in increasing optimality 
ordered_codons = fc_df[order(fc_df$opt), ]$codon

# plot 
fig = fc_df %>%
  ggplot(aes(x = factor(codon, levels = ordered_codons), y = enrichment, fill = opt)) + 
  geom_bar(color = 'black', stat = 'identity') + 
  scale_fill_gradient2('Human CSC', low = '#3852A3', high = '#E81F27', mid = 'white', 
                       midpoint = 0, limits = c(-.2, .2), breaks = c(-.2, 0, .2)) + 
  scale_x_discrete(expand = c(0.013, 0.013)) + 
  scale_y_continuous(limits = c(-0.63, 0.63), breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) + 
  stat_cor(aes(x = opt, y = enrichment), method = 'spearman', label.x = 40, label.y = 0.6, size = 5, 
           inherit.aes = FALSE) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(angle = 90, size = 10, face = 'bold', vjust = 0.5), 
    axis.text.y = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 13, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 13, face = 'bold', vjust = 1), 
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    legend.title = element_text(size = 13, face = 'bold', hjust = 0), 
    legend.text = element_text(size = 10, face = 'bold')
  ) + 
  labs(
    x = 'Codon', 
    y = 'log2 Fold Change of Codon Frequency'
  )

# save 
ggsave(
  './codon_enrichment_up_down_regulated_genes.pdf', 
  fig, 
  height = 4, 
  width = 10
)

################################################################################
