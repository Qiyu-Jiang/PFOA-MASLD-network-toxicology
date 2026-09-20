# run_04b_plots.R - plots only (GO-BP dotplot + KEGG barplot), recomputed from existing results (fast cached run)
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2); library(enrichplot)
})
BASE <- "C:/Users/user/workspace/PFOA_NAFLD_project"
ANA  <- file.path(BASE, "04_analysis"); FIG <- file.path(BASE, "05_figures")

map <- read.csv(file.path(ANA, "run04_symbol2entrez.csv"))
cat("entrez genes:", nrow(map), "\n")

bp <- enrichGO(gene = as.character(map$ENTREZID), OrgDb = org.Hs.eg.db, ont = "BP",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
p1 <- dotplot(bp, showCategory = 12, title = "GO-BP enrichment (PFOA x NAFLD intersection, n=479)")
ggsave(file.path(FIG, "run04_GO_BP_dotplot.png"), p1, width = 9.5, height = 7, dpi = 220)
cat("saved GO-BP dotplot\n")

kegg <- enrichKEGG(gene = as.character(map$ENTREZID), organism = "hsa",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
p2 <- barplot(kegg, showCategory = 15, title = "KEGG pathway enrichment (top 15)")
ggsave(file.path(FIG, "run04_KEGG_barplot.png"), p2, width = 9.5, height = 7, dpi = 220)
cat("saved KEGG barplot\n")

cat("RUN04B-DONE\n")
