# GSE126848 validation analysis: hub & signature gene expression across NAFLD groups
# Rscript analysis_gse126848.R

base <- "C:\\Users\\user\\workspace\\PFOA_NAFLD_project"
counts_path <- file.path(base, "02_data_raw", "cohort_GSE126848", "GSE126848_Gene_counts_raw.txt.gz")
sm_path <- file.path(base, "02_data_raw", "cohort_GSE126848", "GSE126848_series_matrix.txt.gz")
out_dir <- file.path(base, "04_analysis", "cohort_GSE126848")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(org.Hs.eg.db))

counts <- read.delim(gzfile(counts_path), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cat("counts dim:", dim(counts), "\n")

ens <- sub("\\..*$", "", rownames(counts))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep_sym <- ifelse(is.na(sym), ens, sym)
# collapse duplicates: keep row with max mean count
ord <- order(-rowMeans(counts))
counts <- counts[ord, ]
keep_sym <- keep_sym[ord]
dup <- duplicated(keep_sym)
counts <- counts[!dup, ]
keep_sym <- keep_sym[!dup]
rownames(counts) <- keep_sym
cat("after symbol collapse:", dim(counts), "\n")

libsize <- colSums(counts)
cpm <- sweep(counts, 2, libsize, "/") * 1e6
logcpm <- log2(cpm + 1)

# sample groups from series matrix
sm <- readLines(gzfile(sm_path))
getrow <- function(prefix) {
  l <- grep(paste0("^", prefix, "\t"), sm, value = TRUE)[1]
  if (is.na(l)) return(character(0))
  x <- sub(paste0("^", prefix, "\t"), "", l)
  s <- strsplit(x, "\t")[[1]]
  gsub('"', "", s)
}
titles <- getrow("!Sample_title")
descs <- getrow("!Sample_description")
gsms <- getrow("!Sample_geo_accession")
mapping <- data.frame(id = gsub("^0+", "", descs), group = sub("_[0-9]+$", "", titles), title = titles,
                      gsm = gsms, stringsAsFactors = FALSE)
cat("groups:\n"); print(table(mapping$group))

cols_norm <- gsub("^0+", "", colnames(counts))
grp <- mapping$group[match(cols_norm, mapping$id)]
cat("unmatched samples:", sum(is.na(grp)), "\n")
names(grp) <- colnames(counts)

# target genes
TARGETS <- c("TNF", "IL6", "TP53", "ALB", "IL1B", "INS", "AKT1", "PPARG", "EGFR", "BCL2",
             "EPHA2", "FMO1", "MYC", "PHLDA1", "TMPO", "TNFSF10", "TRIB1", "WNT5A")
SIG5 <- c("EPHA2", "FMO1", "MYC", "PHLDA1", "WNT5A")

present <- TARGETS[TARGETS %in% rownames(logcpm)]
missing <- setdiff(TARGETS, rownames(logcpm))
cat("present genes:", length(present), "| missing:", paste(missing, collapse = ","), "\n")

gsel <- logcpm[present, , drop = FALSE]

# group definitions for tests
mask_nash <- grp == "NASH"
mask_nw <- grp %in% c("Normal-weight")
mask_nafl <- grp == "NAFL"

tt <- function(v, g1, g0) {
  x <- v[g1]; y <- v[g0]
  if (length(x) < 2 || length(y) < 2) return(c(NA, NA))
  t <- tryCatch(t.test(x, y), error = function(e) NULL)
  if (is.null(t)) return(c(NA, NA))
  c(t$p.value, mean(x) - mean(y))
}

res <- data.frame()
for (g in present) {
  v <- gsel[g, ]
  a <- tt(v, mask_nash, mask_nw)
  b <- tt(v, mask_nafl, mask_nw)
  res <- rbind(res, data.frame(
    gene = g,
    mean_NASH = mean(v[mask_nash]), mean_NAFL = mean(v[mask_nafl]), mean_NW = mean(v[mask_nw]),
    diff_NASH_NW = a[2], p_NASH_NW = a[1],
    diff_NAFL_NW = b[2], p_NAFL_NW = b[1],
    stringsAsFactors = FALSE))
}
res$p_NASH_NW_BH <- p.adjust(res$p_NASH_NW, method = "BH")

# signature score (z-scored relative to Normal-weight group)
sig_present <- SIG5[SIG5 %in% rownames(logcpm)]
if (length(sig_present) >= 4) {
  z <- t(scale(t(logcpm[sig_present, , drop = FALSE])))
  z <- z[, !is.na(grp)]
  score <- colMeans(z, na.rm = TRUE)
  grp2 <- grp[!is.na(grp)]
  # AUC (mann-whitney) NASH vs NW
  x <- score[grp2 == "NASH"]; y <- score[grp2 == "Normal-weight"]
  if (length(x) >= 2 && length(y) >= 2) {
    r <- rank(c(x, y))
    auc <- (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
    tt2 <- t.test(x, y)
    cat(sprintf("Signature score: NASH vs Normal-weight AUC = %.3f, p = %.4g\n", auc, tt2$p.value))
    writeLines(sprintf("signature score NASH vs Normal-weight AUC=%.3f p=%.4g genes=%d", auc, tt2$p.value, length(sig_present)),
               file.path(out_dir, "signature_score_summary.txt"))

    # score boxplot + per-sample values
    png(file.path(out_dir, "gse126848_score_boxplot.png"), width = 2600, height = 900, res = 260)
    par(mar = c(5, 4, 4, 1))
    boxplot(score ~ factor(grp2, levels = c("Normal-weight", "NAFL", "NASH"), labels = c("NW", "NAFL", "NASH")),
            main = sprintf("Signature score (AUC = %.3f, p = %.3f)", auc, tt2$p.value),
            ylab = "signature score (mean z)", xlab = "", col = c("#bcd8cf", "#e9f1ee", "#f2b8ad"),
            cex.axis = 0.9)
    dev.off()
    gsm_map <- mapping$gsm[match(gsub("^0+", "", names(score)), mapping$id)]
    write.csv(data.frame(gsm = gsm_map, sample = names(score), group = grp2, score = as.numeric(score)),
              file.path(out_dir, "gse126848_score_per_sample.csv"), row.names = FALSE)
  }
}

write.csv(res, file.path(out_dir, "gse126848_gene_stats.csv"), row.names = FALSE)

# figure 1: signature genes (2 x 3)
png(file.path(out_dir, "gse126848_signature_boxplots.png"), width = 2600, height = 1050, res = 260)
par(mfrow = c(2, 3), mar = c(5, 4, 3, 1))
for (g in sig_present) {
  v <- gsel[g, ]
  boxplot(v ~ factor(grp, levels = c("Normal-weight", "NAFL", "NASH"), labels = c("NW", "NAFL", "NASH")),
          main = g, ylab = "log2 CPM", xlab = "", col = c("#bcd8cf", "#e9f1ee", "#f2b8ad"),
          cex.axis = 0.85)
}
dev.off()

# figure 2: hub genes (2 x 5)
hub <- c("TNF", "IL6", "IL1B", "TP53", "ALB", "INS", "AKT1", "PPARG", "EGFR", "BCL2")
hub_p <- hub[hub %in% rownames(logcpm)]
png(file.path(out_dir, "gse126848_hub_boxplots.png"), width = 2600, height = 1000, res = 260)
par(mfrow = c(2, 5), mar = c(5, 4, 3, 1))
for (g in hub_p) {
  v <- gsel[g, ]
  boxplot(v ~ factor(grp, levels = c("Normal-weight", "NAFL", "NASH"), labels = c("NW", "NAFL", "NASH")),
          main = g, ylab = "log2 CPM", xlab = "", col = c("#bcd8cf", "#e9f1ee", "#f2b8ad"),
          cex.axis = 0.85)
}
dev.off()

cat("top results (NASH vs NW):\n")
show <- res[order(res$p_NASH_NW), c("gene", "diff_NASH_NW", "p_NASH_NW", "p_NASH_NW_BH")]
print(head(show, 18))
cat("DONE\n")
