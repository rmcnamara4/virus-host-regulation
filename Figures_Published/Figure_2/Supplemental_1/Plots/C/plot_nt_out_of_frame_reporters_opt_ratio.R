################################################################################
# script to plot Figure 2 - figure supplement 1C
################################################################################

# load libraries 
library(tidyverse)

# load data 
mosquito_codons = read.table(
  '../../Data/mosquito_codon_frequency_ratio.csv', 
  sep = ',', header = TRUE
) %>%
  select(gene_ID, ratio)

nt_out_of_frame_reporters_codons = read.table(
  '../../Data/nt_out_of_frame_reporters_codon_frequency_ratio.csv', 
  sep = ',', header = TRUE
) %>%
  select(frame, optimal, non_optimal, neutral, ratio)

# get number of mosquito genes shown in plot 
n = nrow(mosquito_codons)

# plot density 
plt_density = mosquito_codons %>%
  ggplot(aes(x = ratio)) + 
  geom_density(size = 1.5, aes(y = ..scaled..)) + 
  geom_vline(data = nt_out_of_frame_reporters_codons, aes(xintercept = ratio, color = factor(frame, levels = c('first', 'second'))), 
             linetype = 'dashed', size = 1) + 
  geom_text(x = 4, y = 1, label = paste0('n = ', format(n, big.mark = ',')), size = 5, stat = "identity", check_overlap = TRUE) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_color_manual('Frame', values = c('first' = 'red', 'second' = '#3852A3'), labels = str_to_title) + 
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = NA, color = 'black', size = 3), 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.title = element_text(size = 10, face = 'bold'), 
    legend.text = element_text(size = 10, face = 'bold')
  ) + 
  labs(
    x = 'log2(optimal / non-optimal)', 
    y = 'Density'
  )

# plot bar plot 
coeff = 50

nt_out_of_frame_reporters_codons$x = c(1, 3.5)

plt_bar = nt_out_of_frame_reporters_codons %>%
  ggplot() + 
  geom_bar(aes(x = x - 0.5, y = ratio, fill = factor(frame, levels = c('first', 'second'))), stat = "identity", width = 1, 
           color = 'black', size = 1.5) + 
  geom_bar(aes(x = x + 0.5, y = 100 * neutral / coeff), fill = 'darkgrey', stat = "identity", width = 1, 
           color = 'black', size = 1.5) + 
  scale_y_continuous(
    name = 'log2(optimal / non-optimal)', 
    limits = c(-1.1, 1.1),
    sec.axis = sec_axis(~.*coeff, name = '% neutral codons', breaks = c(0, 25, 50))
  ) +
  scale_x_continuous('Frame', breaks = c(1, 3.5), labels = c('First', 'Second'), 
                     sec.axis = sec_axis(~.)) +
  scale_fill_manual(values = c('first' = 'red', 'second' = '#3852A3')) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1.5) +
  geom_text(aes(x = x - 0.5, label = round(ratio, 2), y = ratio / 2), size = 7, color = 'white', stat = "identity", check_overlap = TRUE) + 
  geom_text(aes(x = x + 0.5, label = paste0(round(neutral * 100, 2), '%'), y = neutral * 100 / (coeff * 2)), 
            size = 7, color = 'white', stat = "identity", check_overlap = TRUE) +
  theme_classic() + 
  theme( 
    axis.text.x = element_text(size = 15, face = 'bold'), 
    axis.text.y = element_text(size = 15, face = 'bold'), 
    axis.text.y.right = element_text(color = 'darkgrey'),
    axis.text.x.top = element_blank(),
    axis.title.y.right = element_text(color = 'darkgrey'),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_line(color = 'darkgrey'),
    axis.line.y.right = element_line(color = 'darkgrey', size = 2),
    axis.line.x.top = element_line(color = 'black', size = 2),
    axis.line.x.bottom = element_line(color = 'black', size = 2), 
    axis.line.y.left = element_line(color = 'black', size = 2),
    axis.title.x = element_text(size = 18, face = 'bold', vjust = -1), 
    axis.title.y = element_text(size = 18, face = 'bold', vjust = 1), 
    legend.position = 'none'
  )

# save
ggsave(
  './nt_out_of_frame_reporters_opt_ratio_density.pdf', 
  plt_density, 
  width = 8, 
  height = 6
)

ggsave(
  './nt_out_of_frame_reporters_opt_ratio_bar.pdf', 
  plt_bar, 
  width = 7, 
  height = 6
)
  
################################################################################










