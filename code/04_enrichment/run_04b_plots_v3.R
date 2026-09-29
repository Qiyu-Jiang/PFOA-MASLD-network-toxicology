# run_04b_plots_v3.R -- KEGG barplot, fill = -log10(p.adjust) (v2.5, Figure 6)
# Change vs v2: colour scale now -log10(p.adjust) so the saturated top-10 bars
# (p.adjust 5.7e-24 .. 4.0e-11, all below the old 5e-11 legend floor) become
# distinguishable. Same top-15 pathways, counts and canvas; data read as-is
# (data/enrichment/run04_KEGG.csv), zero recomputation of the enrichment itself.
suppressMessages(library(ggplot2))
k <- read.csv('data/enrichment/run04_KEGG.csv')
k <- k[order(k\.adjust), ][1:15, ]
k\ <- -log10(k\.adjust)
k\ <- factor(k\, levels = rev(as.character(k\)))
p <- ggplot(k, aes(x = Count, y = Description, fill = neg)) +
  geom_col(width = 0.72, colour = 'grey30', linewidth = 0.2) +
  scale_fill_gradient(low = '#2166AC', high = '#B2182B',
                      name = expression(-log[10](p.adjust))) +
  labs(title = 'KEGG pathway enrichment (top 15)', x = 'Count', y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 13, face = 'bold'),
        axis.text.y = element_text(size = 9.5, colour = 'black'),
        axis.title.x = element_text(size = 11),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))
ggsave('docs/figures/Figure6.png', p, width = 2090/300, height = 1540/300, dpi = 300)
