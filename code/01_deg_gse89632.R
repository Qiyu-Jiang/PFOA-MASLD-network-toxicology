# ============================================================
# run_01_deg_geo_v2.R - differential-expression analysis for GSE89632 / GSE48452 (v2)
# Usage: Rscript run_01_deg_geo_v2.R
# Group definitions (derived from sample titles):
#   GSE89632: liver_HC_*=Control(n=24), liver_SS_*=NAFL(n=20), liver_NASH_*=NASH(n=19)
#   GSE48452: "liver, H/C/S/N, ID" - H+C pooled as Control, S=NAFL, N=NASH
#             (C and H semantics verified against the official GEO series description)
# ============================================================

BASE  <- "C:/Users/user/workspace/PFOA_NAFLD_project"
RAW   <- file.path(BASE, "02_data_raw")
OUT_A <- file.path(BASE, "04_analysis")
OUT_F <- file.path(BASE, "05_figures")
dir.create(OUT_A, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_F, showWarnings = FALSE, recursive = TRUE)

suppressMessages({ library(GEOquery); library(limma); library(ggplot2) })

## ---------- utility functions ----------
get_symbol_vec <- function(fdata) {
  cn  <- colnames(fdata)
  sym <- grep("^(Symbol|GENE_SYMBOL|Gene Symbol)$", cn, ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sym)) {
    s <- as.character(fdata[[sym]])
  } else {
    ga <- grep("gene.assignment", cn, ignore.case = TRUE, value = TRUE)[1]
    if (is.na(ga)) stop("no symbol column; columns: ", paste(cn, collapse = " | "))
    s <- sub("^([^ ]+) //.*$", "\\1", as.character(fdata[[ga]]))
  }
  s[is.na(s)] <- ""
  s[s %in% c("", "---")] <- NA
  s
}

volcano <- function(tt, title, fname) {
  tt$sig <- with(tt, ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                            ifelse(logFC > 0, "Up", "Down"), "NS"))
  p <- ggplot(tt, aes(logFC, -log10(adj.P.Val), color = sig)) +
    geom_point(size = 1.1, alpha = .6) +
    scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey75")) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    labs(title = title, x = "log2FC", y = "-log10(adj.P)") +
    theme_bw()
  ggsave(fname, p, width = 7, height = 5.5, dpi = 300)
}

load_eset <- function(gse_id, local_gz) {
  eset <- NULL
  if (file.exists(local_gz)) {
    eset <- tryCatch(getGEO(filename = local_gz), error = function(e) {
      cat("local parse failed:", conditionMessage(e), "\n"); NULL
    })
  }
  if (is.null(eset)) {
    cat("falling back to online getGEO for", gse_id, "\n")
    gse <- getGEO(gse_id, GSEMatrix = TRUE)
    eset <- gse[[1]]
  }
  eset
}

run_dataset <- function(gse_id, local_gz, code_fun, code_map, contrasts) {
  cat("=====", gse_id, "=====\n")
  eset <- load_eset(gse_id, local_gz)
  expr <- exprs(eset); pd <- pData(eset)
  cat("samples:", ncol(expr), " probes:", nrow(expr), "\n")
  if (max(expr, na.rm = TRUE) > 200) { expr <- log2(expr + 1); cat("log2 transform applied\n") }
  codes <- code_fun(as.character(pd$title))
  grp   <- unname(code_map[codes])
  if (any(is.na(grp))) {
    cat("UNMAPPED samples:", sum(is.na(grp)), "\n")
    print(table(codes[is.na(grp)]))
    grp[is.na(grp)] <- "UNK"
  }
  cat("group counts:\n"); print(table(grp))
  keep <- grp != "UNK"
  expr <- expr[, keep, drop = FALSE]; grp <- factor(grp[keep])
  keep_rows <- complete.cases(expr)
  cat("dropping rows with NA:", sum(!keep_rows), "\n")
  expr <- expr[keep_rows, , drop = FALSE]
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  fit <- lmFit(expr, design)
  fit2 <- eBayes(contrasts.fit(fit, contrasts))
  sym <- get_symbol_vec(fData(eset))
  for (i in seq_along(contrasts)) {
    nm <- names(contrasts)[i]
    tt <- topTable(fit2, number = Inf, coef = nm, adjust.method = "BH")
    tt$Symbol <- sym[rownames(tt)]
    write.csv(tt, file.path(OUT_A, sprintf("%s_%s_topTable.csv", gse_id, nm)),
              row.names = TRUE)
    deg <- subset(tt, adj.P.Val < 0.05 & abs(logFC) > 1)
    cat(sprintf("%s (%s): DEG=%d up=%d down=%d\n", nm, unname(contrasts[i]),
                nrow(deg), sum(deg$logFC > 0), sum(deg$logFC < 0)))
    volcano(tt, paste(gse_id, nm), file.path(OUT_F, sprintf("%s_%s_volcano.png", gse_id, nm)))
  }
  invisible(NULL)
}

## ---------- GSE89632 ----------
code89632 <- function(titles) {
  vapply(titles, function(t) {
    parts <- strsplit(t, "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2) parts[2] else "UNK"
  }, character(1), USE.NAMES = FALSE)
}
map89632 <- c(HC = "Control", SS = "NAFL", NASH = "NASH")
c89632 <- c("NASH - Control", "NAFL - Control", "NASH - NAFL")
names(c89632) <- c("NASH_vs_Control", "NAFL_vs_Control", "NASH_vs_NAFL")
run_dataset("GSE89632", file.path(RAW, "GSE89632_series_matrix.txt.gz"),
            code89632, map89632, c89632)

## ---------- GSE48452 ----------
code48452 <- function(titles) {
  vapply(titles, function(t) {
    parts <- strsplit(t, ",", fixed = TRUE)[[1]]
    if (length(parts) >= 2) trimws(parts[2]) else "UNK"
  }, character(1), USE.NAMES = FALSE)
}
map48452 <- c(H = "Control", C = "Control", S = "NAFL", N = "NASH")
c48452 <- c("NASH - Control", "NAFL - Control")
names(c48452) <- c("NASH_vs_Control", "NAFL_vs_Control")
run_dataset("GSE48452", file.path(RAW, "GSE48452_series_matrix.txt.gz"),
            code48452, map48452, c48452)

writeLines(capture.output(sessionInfo()), file.path(OUT_A, "run01_sessionInfo.txt"))
cat("RUN-ALL-DONE\n")
