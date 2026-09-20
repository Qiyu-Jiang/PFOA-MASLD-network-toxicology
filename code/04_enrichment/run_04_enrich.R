# ============================================================
# run_04_enrich.R - four-way enrichment of the 479 intersection genes (GO BP/CC/MF + KEGG + Hallmark)
# outputs: 04_analysis/*.csv, 05_figures/*.png
# ============================================================
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
  library(enrichplot); library(msigdbr)
})

BASE <- "C:/Users/user/workspace/PFOA_NAFLD_project"
ANA  <- file.path(BASE, "04_analysis"); FIG <- file.path(BASE, "05_figures")

genes <- readLines(file.path(ANA, "intersection_genes_479.txt"))
genes <- unique(trimws(genes)); genes <- genes[genes != ""]
cat("input genes:", length(genes), "\n")

map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("mapped to Entrez:", nrow(map), "/", length(genes), "\n")
write.csv(map, file.path(ANA, "run04_symbol2entrez.csv"), row.names = FALSE)

run_go <- function(ont, tag) {
  r <- tryCatch(enrichGO(gene = map$ENTREZID, OrgDb = org.Hs.eg.db, ont = ont,
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                         readable = TRUE), error = function(e) { cat(tag, "FAIL:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
    df <- as.data.frame(r)
    write.csv(df, file.path(ANA, sprintf("run04_GO_%s.csv", tag)), row.names = FALSE)
    cat(sprintf("GO-%s: %d terms | top: %s\n", tag, nrow(df), paste(head(df$Description, 5), collapse = " / ")))
  } else cat("GO-", tag, ": no significant terms (or failed)\n", sep = "")
  r
}
bp <- run_go("BP", "BP"); cc <- run_go("CC", "CC"); mf <- run_go("MF", "MF")

kegg <- tryCatch(enrichKEGG(gene = as.character(map$ENTREZID), organism = "hsa",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05), error = function(e) { cat("KEGG FAIL:", conditionMessage(e), "\n"); NULL })
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  kdf <- as.data.frame(setReadable(kegg, org.Hs.eg.db, keyType = "ENTREZID"))
  write.csv(kdf, file.path(ANA, "run04_KEGG.csv"), row.names = FALSE)
  cat(sprintf("KEGG: %d pathways | top: %s\n", nrow(kdf), paste(head(kdf$Description, 5), collapse = " / ")))
} else cat("KEGG: no significant pathways (or failed)\n")

# Hallmark (compatible with old and new msigdbr APIs)
hs <- tryCatch({
  msigdbr(species = "Homo sapiens", collection = "H")
}, error = function(e) {
  msigdbr(species = "Homo sapiens", category = "H")
})
cat("Hallmark gene sets:", length(unique(hs$gs_name)), "\n")
hal <- enricher(gene = map$ENTREZID, TERM2GENE = hs[, c("gs_name", "entrez_gene")],
                pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
if (!is.null(hal) && nrow(as.data.frame(hal)) > 0) {
  hdf <- as.data.frame(hal)
  write.csv(hdf, file.path(ANA, "run04_Hallmark.csv"), row.names = FALSE)
  cat(sprintf("Hallmark: %d sets | top: %s\n", nrow(hdf), paste(head(hdf$ID, 6), collapse = " / ")))
} else cat("Hallmark: no significant sets\n")

# ---- plots ----
if (!is.null(bp) && nrow(as.data.frame(bp)) > 0) {
  p1 <- dotplot(bp, showCategory = 12, title = "GO-BP enrichment (PFOA x NAFLD intersection, n=479)")
  ggsave(file.path(FIG, "run04_GO_BP_dotplot.png"), p1, width = 9, height = 7, dpi = 200)
}
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  p2 <- barplot(kegg, showCategory = 15, title = "KEGG pathway enrichment (top 15)")
  ggsave(file.path(FIG, "run04_KEGG_barplot.png"), p2, width = 9, height = 7, dpi = 200)
}
if (exists("hal") && !is.null(hal) && nrow(as.data.frame(hal)) > 0) {
  p3 <- dotplot(hal, showCategory = 12, title = "Hallmark enrichment (top 12)")
  ggsave(file.path(FIG, "run04_Hallmark_dotplot.png"), p3, width = 9, height = 7, dpi = 200)
}
cat("figures saved to 05_figures\n")
cat("RUN04-DONE\n")
