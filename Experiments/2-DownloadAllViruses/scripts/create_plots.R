################################################################################
# Load libraries

library(dplyr)
library(ggplot2)
library(rio)

################################################################################
# Load data and split

human_corr = read.table(
  '/n/projects/rm2498/Virus_Project/Viral_Seq_Analysis/RSCU_Correlations/all_correlations_total_human.csv', 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

mosquito_corr = read.table(
  '/n/projects/rm2498/Virus_Project/Viral_Seq_Analysis/RSCU_Correlations/all_correlations_total_mosquito.csv', 
  sep = ',', 
  header = TRUE, 
  quote = "\""
)

################################################################################

iv = grep('human', human_corr$host)
test = human_corr[iv, ]

test$host = str_replace_all(test$host, 'vertebrates: human', 'human:vertebrates')

test %>%
  ggplot(aes(x = corr, y = -log10(p.value), color = host)) +
  geom_point() +
  geom_hline(yintercept = 2, linetype = 'dashed', color = 'grey') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  theme_classic()

test %>%
  dplyr::filter(host == 'human') %>%
  ggplot(aes(x = corr)) +
  geom_histogram(aes(fill = host, y = (..count..) / sum(..count..)), bins = 45, position = "dodge2") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual('Host', labels = c('Human', 'Human & Invertebrates', 'Human & Invertebrates & Vertebrates',
                              'Human & Protozoa', 'Human & Vertebrates'),
         values = c('#246EB9', '#4CB944', '#F5EE9E', '#FDFFFC', '#F06543')) +
  theme_classic()

specific_viruses = test %>%
  dplyr::filter(
    locus_tag == 'DV1' | 
    grepl('Zika', virus) |
    grepl('CHIKV', locus_tag) |
    grepl('HHV1', locus_tag) | 
    grepl('GU280', locus_tag) | 
    locus_tag == 'FLUAV'
  )

test %>%
  ggplot(aes(x = corr)) +
  geom_density(aes(y = ..scaled..)) +
  geom_vline(data = specific_viruses, aes(xintercept = corr, color = virus), linetype = 'dashed') +
  theme_classic() 

iv = grep('invertebrates', mosquito_corr$host)
test2 = mosquito_corr[iv, ]

test2 %>%
  ggplot(aes(x = corr, y = -log10(p.value), color = type)) +
  geom_point() +
  geom_hline(yintercept = 2, linetype = 'dashed', color = 'grey') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  theme_classic()

test2 %>%
  ggplot(aes(x = corr, fill =type)) +
  geom_histogram(color = 'black', aes(y = (..count..) / sum(..count..))) +
  theme_classic()
