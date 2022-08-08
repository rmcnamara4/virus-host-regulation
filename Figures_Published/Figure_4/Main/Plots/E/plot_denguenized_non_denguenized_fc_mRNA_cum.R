################################################################################
# script to plot Figure 4E
################################################################################

# load libraries
library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)

# load data 
codon_groups = read_excel(
  '../../Data/codon_groups.xlsx', 
  sheet = 'Human'
)

diff_expression = read.table(
  '../../Data/uninfected-high.csv', 
  sep = ',', header = TRUE
) %>%
  select(gene_ID, logFC) 

human_codons = read_excel(
  '../../Data/human_hg38_stats.xlsx', 
  sheet = 'Codons'
) %>%
  select(-c(species, len)) %>%
  as.data.frame
rownames(human_codons) = human_codons$gene_ID
human_codons = human_codons[-1]

# get the sum of Denguenized codons for each human gene
denguenized_codons = codon_groups$denguenized[!is.na(codon_groups$denguenized)]
denguenized_sum = rowSums(human_codons[denguenized_codons])

# get the sum of Not Denguenized codons for each human gene
non_denguenized_codons = codon_groups$non_denguenized[!is.na(codon_groups$non_denguenized)]
non_denguenized_sum = rowSums(human_codons[non_denguenized_codons])

# keep only genes that have differential expression information 
denguenized_sum = denguenized_sum[names(denguenized_sum) %in% diff_expression$gene_ID]
non_denguenized_sum = non_denguenized_sum[names(non_denguenized_sum) %in% diff_expression$gene_ID]

# get the 250 genes that contain the highest frequency of codons from each group 
denguenized_sum = denguenized_sum[order(-denguenized_sum)][1:250]
non_denguenized_sum = non_denguenized_sum[order(-non_denguenized_sum)][1:250]

# get the mRNA fold change data for each group of 250 genes 
denguenized_fc = diff_expression[diff_expression$gene_ID %in% names(denguenized_sum), ] %>%
  mutate(
    group = 'Denguenized'
  )

non_denguenized_fc = diff_expression[diff_expression$gene_ID %in% names(non_denguenized_sum), ] %>%
  mutate(
    group = 'Not Denguenized'
  )

# rbind Denguenized and Not Denguenized data 
dengue_fc = rbind(denguenized_fc, non_denguenized_fc)

# isolate all other genes not included in the groups 
gene_ids_remove = unique(dengue_fc$gene_ID)
all_other_fc = diff_expression[!(diff_expression$gene_ID %in% gene_ids_remove), ] %>%
  mutate(
    group = 'All'
  )

# get 250 groups of 250 random genes from all other genes and add to the group data 
for (i in 1:250) {
  
  iv = sample(1:nrow(all_other_fc), 250, replace = FALSE)
  hold = all_other_fc[iv, ] %>%
    mutate(
      group = paste0('group_', i)
    )
  if (i == 1) {
    all_genes_fc = rbind(dengue_fc, hold) 
  } else {
    all_genes_fc = rbind(all_genes_fc, hold)
  }
  
}

# add all other genes to the data.frame
all_genes_fc = rbind(all_genes_fc, all_other_fc) 

# add a column for the color 
all_genes_fc$c = ifelse(all_genes_fc$group == 'Denguenized', 'Denguenized', 
                        ifelse(all_genes_fc$group == 'Not Denguenized', 'Not Denguenized', 
                               ifelse(all_genes_fc$group == 'All', 'All', 'other')))

# get count for each group 
n_non_denguenized = nrow(all_genes_fc[all_genes_fc$c == 'Not Denguenized', ])
n_denguenized = nrow(all_genes_fc[all_genes_fc$c == 'Denguenized', ])
n_all = nrow(all_genes_fc[all_genes_fc$c == 'All', ])

count_df = data.frame(
  group = c('Not Denguenized', 'Denguenized', 'All'), 
  count = c(n_non_denguenized, n_denguenized, n_all), 
  x.pos = c(1, 1, 1), 
  y.pos = c(0.35, 0.25, 0.15)
)

# calculate significance between Denguenized, Not Denguenized, and All
sig = all_genes_fc %>%
  filter(c != 'other') %>%
  wilcox_test(logFC ~ c, paired = FALSE)

# plot 
levs = c(paste0('group_', 1:250), 'All', 'Not Denguenized', 'Denguenized')
all_genes_fc$group = factor(all_genes_fc$group, levels = levs)

fig = all_genes_fc %>%
  ggplot(aes(x = logFC, group = group, color = c, alpha = c)) +
  geom_line(aes(y = ..y.., size = c), stat = 'ecdf') +
  scale_alpha_manual(values = c('Not Denguenized' = 1, 
                                'Denguenized' = 1, 
                                'other' = 0.1, 
                                'All' = 1)) +
  scale_size_manual(values = c('Not Denguenized' = 1.5, 
                               'Denguenized' = 1.5, 
                               'other' = 0.5, 
                               'All' = 1.5)) +
  scale_color_manual(values = c('Not Denguenized' = '#24C9C9', 
                                'Denguenized' = '#CF3CD3',
                                'other' = 'grey', 
                                'All' = 'black')) +
  geom_text(data = count_df, aes(label = paste0(group, '\n', '(n = ', format(count, big.mark = ','), ')'), x = x.pos, y = y.pos), size = 5, 
            color = c('#24C9C9', '#CF3CD3', 'black'), inherit.aes = FALSE, check_overlap = TRUE) +
  geom_text(label = 'Random\n(n = 250 samplings of 250 genes)', x = 1, y = 0.05, color = 'grey', size = 5, 
            inherit.aes = FALSE, check_overlap = TRUE) +  
  geom_text(data = sig, aes(label = paste0(group1, '-', group2, ' = ', p.adj)), x = -1.25, y = c(1, 0.95, 0.9), size = 5, 
            inherit.aes = FALSE, check_overlap = TRUE) +
  theme_classic() +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  theme(    
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.position = 'none'
  ) + 
  labs(
    x = 'log2 Fold Change of mRNA Level', 
    y = 'Fraction'
  )

# save 
ggsave(
  './denguenized_non_denguenized_fc_mRNA_cum.pdf', 
  fig, 
  height = 7, 
  width = 8
)

################################################################################