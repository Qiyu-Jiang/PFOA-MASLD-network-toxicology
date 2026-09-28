# run_04b_plots_v2.R - Figure 5 rebuild: y-axis ordered by p.adjust, colour = -log10(adj. P)
# Data source: run04_GO_BP.csv (same file as Table 3 / SI enrichment tables -> zero recomputation)
suppressMessages(library(ggplot2))

BASE <- "C:/Users/user/workspace/PFOA_NAFLD_project"
cands <- c(file.path(BASE, "04_analysis", "run04_GO_BP.csv"),
           file.path(BASE, "\u8f6c\u6295\u5305_v2.2_20260926", "06_GitHub_Repository",
                     "PFOA-MASLD-network-toxicology", "data", "enrichment", "run04_GO_BP.csv"))
CSV <- cands[file.exists(cands)][1]
cat("csv:", CSV, "\n")

df <- read.csv(CSV, check.names = FALSE)
df <- df[order(df$p.adjust), ]
df <- df[1:12, ]
cat("top12 terms:\n"); print(df$Description)
stopifnot(df$Description[1] == "response to nutrient levels")
stopifnot(df$Description[4] == "response to molecule of bacterial origin")

gr <- sapply(strsplit(as.character(df$GeneRatio), "/"),
             function(x) as.numeric(x[1]) / as.numeric(x[2]))
df$GeneRatioNum <- gr
df$Description  <- factor(df$Description, levels = rev(df$Description))
df$nl <- -log10(df$p.adjust)

p <- ggplot(df, aes(x = GeneRatioNum, y = Description)) +
  geom_point(aes(size = Count, colour = nl)) +
  scale_colour_gradient(low = "#2166ac", high = "#b2182b",
                        name = expression(-log[10]("adj. P"))) +
  scale_size_continuous(range = c(3, 9), name = "Count") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = "GO-BP enrichment (PFOA x NAFLD intersection, n=479)",
       x = "GeneRatio", y = NULL) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 11, colour = "black"))

out <- file.path(BASE, "05_figures", "run04_GO_BP_dotplot_v2.png")
ggsave(out, p, width = 9.5, height = 7, dpi = 220)
cat("saved:", out, file.size(out), "\n")
cat("FIG5V2-DONE\n")
