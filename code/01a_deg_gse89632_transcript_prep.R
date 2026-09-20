# ============================================================
# run_01_deg_geo.R - v2 (statically reviewed revision, 2026-09-07)
#
# Purpose
#   limma differential-expression analysis for GSE89632 (discovery) and GSE48452 (staging/validation).
#   Probe-level output (multiple probes per gene are not collapsed); csv files carry Symbol for later gene-level aggregation.
#
# Usage
#   Rscript run_01_deg_geo.R   (non-interactive; all paths are absolute constants, no setwd)
#
# Outputs
#   04_analysis/<GSE>_<contrast>_topTable.csv : full limma result for each contrast
#                                               (ProbeID, Symbol, logFC, ..., adj.P.Val)
#   04_analysis/deg_summary.csv               : DEG counts per contrast
#   05_figures/<GSE>_<contrast>_volcano.png   : volcano plot for each contrast
#
# Key assumptions (manually reviewed)
#   1) GSE89632 / GPL14951 (Illumina HumanHT-12 V4): group codes are encoded in the sample title
#      prefix of pData(eset)$title: liver_HC_* = healthy control HC (n=24),
#      liver_SS_* = simple steatosis SS (n=20), liver_NASH_* = NASH (n=19); 63 samples in total.
#      Note: group labels for this dataset are NOT in characteristics_ch1.
#   2) GSE48452 / GPL11532 (Affymetrix Human Gene 1.0 ST): titles look like
#      "liver, H, A1359-50"; codes H(27)/C(14)/S(14)/N(18); 73 samples in total.
#      [verified] Semantics confirmed against the official GEO series description: C=control(14),
#      H=healthy obese(27), S=steatosis(14), N=NASH(18).
#      Pooling H+C as Control (n=41) is justified; if an obesity-matched control is needed later,
#      only the mapping table inside parse_gse48452_group() must be changed.
#   3) DEG thresholds: |log2FC| > 1 and BH-adjusted adj.P.Val < 0.05; probe-level output,
#      no collapsing of multiple probes per gene (csv carries Symbol for later aggregation).
#   4) Scale heuristic: maximum expression > 200 is treated as linear intensity -> log2(expr + 1)
#      (triggered by Illumina raw intensities of GSE89632; GSE48452 is already on log2 scale).
#   5) Data loading: local series matrix (.txt.gz) under 02_data_raw takes precedence;
#      if the file is missing/truncated (parse failure or failed sample/symbol checks), fall back to online download.
# ============================================================

## dependency check (non-interactive: stop with installation guidance when a package is missing)
for (pkg in c("GEOquery", "limma", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("missing dependency: ", pkg,
         "\nPlease run first: install.packages(c('BiocManager','ggplot2')); ",
         "BiocManager::install(c('GEOquery','limma'))\n")
  }
}
suppressMessages({
  library(GEOquery); library(limma); library(ggplot2)
})

## ---------- paths and global parameters (absolute constants; do not use setwd/sys.frame tricks) ----------
BASE_DIR <- "C:/Users/user/workspace/PFOA_NAFLD_project"
RAW_DIR  <- file.path(BASE_DIR, "02_data_raw")
OUT_A    <- file.path(BASE_DIR, "04_analysis")
OUT_F    <- file.path(BASE_DIR, "05_figures")
dir.create(OUT_A, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_F, showWarnings = FALSE, recursive = TRUE)

LFC_CUT         <- 1      # DEG: |log2FC| threshold
PADJ_CUT        <- 0.05   # DEG: BH-adjusted P threshold
LOG2_LINEAR_MAX <- 200    # max expr above this -> treated as linear intensity -> log2(expr+1)

## ---------- probe symbol annotation: cross-platform generic detection ----------
## GPL14951: fData contains a "Symbol" column directly;
## GPL11532: symbols are in the "gene assignment" column, format "SYMBOL // description // location"; take the first token.
clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x %in% c("---", "NA", "N/A", "-")] <- NA_character_
  x
}
probe_symbols <- function(fdat) {
  cn   <- colnames(fdat)
  norm <- tolower(gsub("[._ -]", "", cn))
  sym <- NULL; method <- NA_character_
  hit <- which(norm %in% c("symbol", "genesymbol"))[1]
  if (!is.na(hit)) {
    sym <- clean_symbol(fdat[[hit]]); method <- cn[hit]
  } else {
    hit <- which(norm == "geneassignment")[1]
    if (!is.na(hit)) {
      first_tok <- vapply(strsplit(as.character(fdat[[hit]]), "//", fixed = TRUE),
                          function(p) if (length(p) >= 1) p[1] else "", character(1))
      sym <- clean_symbol(first_tok)
      method <- paste0(cn[hit], " (first token)")
    }
  }
  if (is.null(sym)) {
    cat("  [warn] no symbol column found in fData; available columns: ", paste(cn, collapse = ", "), "\n")
    return(NULL)
  }
  attr(sym, "method") <- method
  sym
}

## ---------- scale check (Illumina linear intensity -> log2; applied uniformly to both datasets) ----------
ensure_log2 <- function(expr, label) {
  m <- suppressWarnings(max(expr, na.rm = TRUE))
  if (!is.finite(m)) stop(label, ": expression matrix has no finite values (all NA?)")
  if (m > LOG2_LINEAR_MAX) {
    cat("  [scale] max =", format(m, digits = 6), ">", LOG2_LINEAR_MAX,
        "-> treated as linear intensity, applying log2(expr + 1)\n")
    log2(expr + 1)
  } else {
    cat("  [scale] max =", format(m, digits = 6), "<=", LOG2_LINEAR_MAX,
        "-> treated as already log2 scale, no transform\n")
    expr
  }
}

## ---------- group parsing, implemented separately per dataset (input = title vector, output = group vector) ----------
## GSE89632: titles look like liver_HC_xxx / liver_SS_xxx / liver_NASH_xxx
parse_gse89632_group <- function(title) {
  tv <- toupper(as.character(title))
  out <- rep(NA_character_, length(tv))
  out[grepl("LIVER_HC",   tv, fixed = TRUE)] <- "HC"
  out[grepl("LIVER_SS",   tv, fixed = TRUE)] <- "SS"
  out[grepl("LIVER_NASH", tv, fixed = TRUE)] <- "NASH"
  out
}
## GSE48452: titles look like "liver, H, A1359-50"; codes H/C/S/N
## [note] semantic assumption: H+C -> Control, S -> NAFL, N -> NASH (independently verified; if refuted, change one map line)
parse_gse48452_group <- function(title) {
  map <- c(H = "Control", C = "Control", S = "NAFL", N = "NASH")
  parts <- strsplit(toupper(as.character(title)), ",", fixed = TRUE)
  vapply(parts, function(p) {
    p <- trimws(p)
    hit <- p[p %in% names(map)]
    if (length(hit) >= 1) unname(map[[hit[1]]]) else NA_character_
  }, character(1))
}

## ---------- eset validation (sample count + symbol annotation coverage; triggers fallback for truncated local files) ----------
check_eset_valid <- function(eset, expected_n, label) {
  if (is.null(eset)) return(FALSE)
  ok_n <- (ncol(exprs(eset)) == expected_n)
  if (!ok_n) cat("  [check]", label, ": sample count", ncol(exprs(eset)),
                "!= expected", expected_n, "\n")
  sym <- probe_symbols(fData(eset))
  ok_sym <- !is.null(sym) && mean(!is.na(sym)) >= 0.5
  if (!ok_sym) cat("  [check]", label, ": probe symbol annotation missing or coverage < 50%\n")
  ok_n && ok_sym
}

## ---------- data loading: local series matrix first, fall back online on failure/invalid file ----------
load_eset <- function(gse_id, platform, local_file, expected_n) {
  eset <- NULL
  if (!is.null(local_file) && file.exists(local_file)) {
    cat("trying local file:", local_file, "\n")
    eset <- tryCatch({
      e <- getGEO(filename = local_file, GSEMatrix = TRUE)
      if (check_eset_valid(e, expected_n, paste0(gse_id, "/local"))) e else NULL
    }, error = function(e) {
      cat("  [warn] local read failed: ", conditionMessage(e), "\n  -> falling back to online download\n")
      NULL
    })
  } else {
    cat("local file does not exist (",
        ifelse(is.null(local_file), "<unspecified>", local_file), ") -> downloading online directly\n", sep = "")
  }
  if (is.null(eset)) {
    cat("downloading", gse_id, "online (GSEMatrix, getGPL=FALSE)...\n")
    lst <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
    idx <- which(grepl(platform, names(lst)))
    if (length(idx) == 0L) idx <- which.max(vapply(lst, function(e) ncol(exprs(e)), integer(1)))
    e <- lst[[idx[1L]]]
    if (!check_eset_valid(e, expected_n, paste0(gse_id, "/online"))) {
      cat("  [warn] embedded annotation insufficient -> re-download with getGPL=TRUE and merge full GPL annotation\n")
      lst2 <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)
      idx2 <- which(grepl(platform, names(lst2))); if (length(idx2) == 0L) idx2 <- 1L
      e <- lst2[[idx2[1L]]]
    }
    if (!check_eset_valid(e, expected_n, paste0(gse_id, "/online-final")))
      stop(gse_id, " online data validation failed (abnormal sample count or symbol annotation); please inspect manually")
    eset <- e
  }
  eset
}

## ---------- single-dataset full limma workflow ----------
## contrasts_def: named character vector, e.g. c(NASH_vs_HC = "NASH-HC", ...)
run_dataset <- function(gse_id, platform, local_file, expected_n,
                        parse_group, group_levels, expected_counts,
                        contrasts_def, note = "") {
  cat("\n==== ", gse_id, " ====\n", sep = "")
  if (nzchar(note)) cat("  assumption/note:", note, "\n")

  eset <- load_eset(gse_id, platform, local_file, expected_n)

  ## GSM alignment: anchor on pData rownames (GSM); expression columns aligned strictly to pData rows
  pd   <- pData(eset)
  expr <- exprs(eset)
  common <- intersect(rownames(pd), colnames(expr))
  if (length(common) < ncol(expr))
    cat("  [align] dropping", ncol(expr) - length(common), "expression columns missing pData\n")
  expr <- expr[, common, drop = FALSE]
  pd   <- pd[common, , drop = FALSE]
  cat("  samples:", ncol(expr), "  probes:", nrow(expr), "\n")
  if (nrow(expr) < 20000)
    cat("  [warn] probe count <20000; if this is a truncated local file please check (GPL14951~47k, GPL11532~33k)\n")

  ## group parsing + per-group sample-count assertions (prevent silent misclassification)
  grp_vec <- parse_group(pd$title)
  names(grp_vec) <- rownames(pd)
  tbl <- table(grp_vec, useNA = "ifany")
  cat("  group parsing result:\n"); print(tbl)
  if (anyNA(grp_vec)) {
    bad <- pd$title[is.na(grp_vec)]
    stop(gse_id, ": ", length(bad), " samples could not be parsed from title, examples:\n",
         paste(head(bad, 5), collapse = "\n"))
  }
  for (g in names(expected_counts)) {
    if (sum(grp_vec == g) != expected_counts[[g]])
      stop(gse_id, ": group ", g, " sample count ", sum(grp_vec == g),
           " != expected ", expected_counts[[g]], "; please check the title parsing rules")
  }

  ## scale + NA logging
  expr <- ensure_log2(expr, gse_id)
  if (anyNA(expr))
    cat("  [warn] expression matrix contains", sum(is.na(expr)),
        "NAs; fitted coefficients for the corresponding probes will be NA (no automatic imputation)\n")

  ## symbol annotation; probes without symbols are discarded before fitting (probe-level output, no multi-probe collapse)
  fd  <- fData(eset)
  sym <- probe_symbols(fd)
  if (is.null(sym)) stop(gse_id, ": could not obtain probe symbol annotation, aborting")
  sym_method <- attr(sym, "method")
  sym <- sym[match(rownames(expr), rownames(fd))]
  keep <- !is.na(sym)
  cat("  [annotation] symbol source:", sym_method,
      "; discarding", sum(!keep), "probes without symbols, keeping", sum(keep), "\n")
  expr <- expr[keep, , drop = FALSE]
  sym  <- sym[keep]; names(sym) <- rownames(expr)

  ## design matrix and contrasts (explicit makeContrasts construction, compatible across limma versions)
  grp <- factor(grp_vec, levels = group_levels)
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- gsub("^grp", "", colnames(design))
  stopifnot(identical(colnames(design), group_levels))
  fit <- lmFit(expr, design)
  cm  <- makeContrasts(contrasts = unname(contrasts_def), levels = design)
  colnames(cm) <- names(contrasts_def)
  fit2 <- eBayes(contrasts.fit(fit, cm))

  ## per contrast: topTable csv (ProbeID+Symbol as first two columns) + volcano plot
  summary_rows <- list()
  for (cn in colnames(cm)) {
    tt <- topTable(fit2, coef = cn, number = Inf,
                   adjust.method = "BH", sort.by = "P")
    tt$ProbeID <- rownames(tt)
    tt$Symbol  <- unname(sym[tt$ProbeID])
    tt <- tt[, c("ProbeID", "Symbol",
                 setdiff(colnames(tt), c("ProbeID", "Symbol"))), drop = FALSE]
    csv_name <- paste0(gse_id, "_", cn, "_topTable.csv")
    write.csv(tt, file.path(OUT_A, csv_name),
              row.names = FALSE, fileEncoding = "UTF-8")

    deg  <- subset(tt, !is.na(adj.P.Val) & adj.P.Val < PADJ_CUT & abs(logFC) > LFC_CUT)
    n_up <- sum(deg$logFC > 0); n_dn <- sum(deg$logFC < 0)
    cat("  contrast", cn, ": probes", nrow(tt), ", DEG", nrow(deg),
        "(up", n_up, "/ down", n_dn, ") ->", csv_name, "\n")
    summary_rows[[length(summary_rows) + 1L]] <-
      data.frame(dataset = gse_id, contrast = cn, n_probe = nrow(tt),
                 n_deg = nrow(deg), n_up = n_up, n_down = n_dn)

    ## volcano plot (horizontal line = adj.P threshold; vertical lines = +/-logFC threshold)
    vol <- tt
    vol$neglog10P <- -log10(pmax(vol$adj.P.Val, 1e-30))  # prevent adj.P~0 -> Inf breaking the axis
    vol$sig <- ifelse(vol$adj.P.Val < PADJ_CUT & abs(vol$logFC) > LFC_CUT,
                      ifelse(vol$logFC > 0, "Up", "Down"), "NS")
    lv <- intersect(c("Up", "Down", "NS"), unique(vol$sig))
    vol$sig <- factor(vol$sig, levels = lv)
    p <- ggplot(vol, aes(x = logFC, y = neglog10P, color = sig)) +
      geom_point(size = 1.1, alpha = 0.6) +
      geom_hline(yintercept = -log10(PADJ_CUT), linetype = 2, colour = "grey35") +
      geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, colour = "grey35") +
      scale_color_manual(values = c(Up = "#d62728", Down = "#1f77b4",
                                    NS = "grey75")[lv], name = NULL) +
      labs(title = paste0(gse_id, ": ", cn),
           x = "log2 fold change", y = "-log10(adjusted P)") +
      theme_bw()
    png_name <- paste0(gse_id, "_", cn, "_volcano.png")
    ggsave(file.path(OUT_F, png_name), p, width = 7, height = 5.5, dpi = 300)
    cat("  volcano plot ->", png_name, "\n")
  }
  do.call(rbind, summary_rows)
}

## ---------- main workflow: two datasets ----------
res1 <- run_dataset(
  gse_id          = "GSE89632",
  platform        = "GPL14951",
  local_file      = file.path(RAW_DIR, "GSE89632_series_matrix.txt.gz"),
  expected_n      = 63L,
  parse_group     = parse_gse89632_group,
  group_levels    = c("HC", "SS", "NASH"),
  expected_counts = c(HC = 24L, SS = 20L, NASH = 19L),
  contrasts_def   = c(NASH_vs_HC = "NASH-HC",   # primary contrast
                      SS_vs_HC   = "SS-HC",
                      NASH_vs_SS = "NASH-SS"),
  note = "discovery cohort; group codes encoded in title prefixes liver_HC_ / liver_SS_ / liver_NASH_")

res2 <- run_dataset(
  gse_id          = "GSE48452",
  platform        = "GPL11532",
  local_file      = file.path(RAW_DIR, "GSE48452_series_matrix.txt.gz"),
  expected_n      = 73L,
  parse_group     = parse_gse48452_group,
  group_levels    = c("Control", "NAFL", "NASH"),
  expected_counts = c(Control = 41L, NAFL = 14L, NASH = 18L),
  contrasts_def   = c(NASH_vs_Control = "NASH-Control",
                      NAFL_vs_Control = "NAFL-Control"),
  note = paste("staging/validation cohort; title codes H/C/S/N.",
               "Semantic assumption (verified): H+C->Control(27+14=41), S->NAFL(14), N->NASH(18)"))

summary_all <- rbind(res1, res2)
write.csv(summary_all, file.path(OUT_A, "deg_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n==== DEG summary ====\n"); print(summary_all)
cat("\nDone. csv -> 04_analysis/, figures -> 05_figures/. Next: merge target pools in 03_data_processed.\n")
