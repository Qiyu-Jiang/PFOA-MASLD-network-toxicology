# NHANES 2015-2018 PFOA - liver enzyme/FLI association analysis (revision round 2)
# Rscript nhanes_analysis.R
suppressPackageStartupMessages({ library(survey); library(foreign) })

base <- "C:\\Users\\user\\workspace\\PFOA_NAFLD_project"
nha <- "C:\\Users\\user\\workspace\\nhanes"
out_dir <- file.path(base, "04_analysis", "nhanes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_cyc <- function(cyc) {
  dem <- read.xport(file.path(nha, paste0("DEMO_", cyc, ".XPT")))
  pfas <- read.xport(file.path(nha, paste0("PFAS_", cyc, ".XPT")))
  bio <- read.xport(file.path(nha, paste0("BIOPRO_", cyc, ".XPT")))
  bmx <- read.xport(file.path(nha, paste0("BMX_", cyc, ".XPT")))
  cot <- read.xport(file.path(nha, paste0("COT_", cyc, ".XPT")))
  d <- Reduce(function(a, b) merge(a, b, by = "SEQN", all = FALSE), list(
    dem[, c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "SDMVPSU", "SDMVSTRA")],
    pfas[, c("SEQN", "LBXNFOA", "LBXBFOA", "WTSB2YR")],
    bio[, c("SEQN", "LBXSATSI", "LBXSASSI", "LBXSGTSI", "LBXSTR")],
    bmx[, c("SEQN", "BMXBMI", "BMXWAIST")],
    cot[, c("SEQN", "LBXCOT")]
  ))
  d$cycle <- cyc
  d
}

df <- rbind(read_cyc("I"), read_cyc("J"))
cat("merged rows:", nrow(df), "\n")

df$pfoa <- df$LBXNFOA + df$LBXBFOA
df$FLI_y <- 0.953 * log(df$LBXSTR) + 0.139 * df$BMXBMI + 0.718 * log(df$LBXSGTSI) + 0.053 * df$BMXWAIST - 15.745
df$FLI <- 100 * exp(df$FLI_y) / (1 + exp(df$FLI_y))

d2 <- subset(df,
             RIDAGEYR >= 20 &
             !is.na(pfoa) & pfoa > 0 &
             !is.na(LBXSATSI) & LBXSATSI > 0 &
             !is.na(LBXSASSI) & LBXSASSI > 0 &
             !is.na(LBXSGTSI) & LBXSGTSI > 0 &
             !is.na(BMXBMI) & !is.na(LBXCOT) & !is.na(WTSB2YR) &
             !is.na(SDMVPSU) & !is.na(SDMVSTRA))
cat("complete-case adults:", nrow(d2), "\n")
cat("FLI component missingness (waist, TG):", sum(is.na(d2$BMXWAIST)), sum(is.na(d2$LBXSTR)), "\n")

d2$w <- d2$WTSB2YR / 2
d2$strata <- interaction(d2$cycle, d2$SDMVSTRA, drop = TRUE)
d2$psu <- interaction(d2$cycle, d2$SDMVPSU, drop = TRUE)
d2$log2pfoa <- log2(d2$pfoa)

des <- svydesign(ids = ~psu, strata = ~strata, weights = ~w, nest = TRUE, data = d2)

fit_outcome <- function(f, label, scale = c("log", "linear")) {
  scale <- match.arg(scale)
  m <- svyglm(f, design = des)
  co <- summary(m)$coefficients
  b <- co["log2pfoa", 1]; se <- co["log2pfoa", 2]; pv <- co["log2pfoa", 4]
  n_model <- tryCatch(length(m$residuals), error = function(e) NA_integer_)
  if (scale == "log") {
    est <- (2^b - 1) * 100
    lo <- (2^(b - 1.96 * se) - 1) * 100
    hi <- (2^(b + 1.96 * se) - 1) * 100
    unit <- "percent"
    cat(sprintf("%-4s per doubling PFOA: %+0.2f%% (95%%CI %+0.2f..%+0.2f), p=%.4g, n=%d\n", label, est, lo, hi, pv, n_model))
  } else {
    est <- b
    lo <- b - 1.96 * se
    hi <- b + 1.96 * se
    unit <- "points"
    cat(sprintf("%-4s per doubling PFOA: %+0.3f points (95%%CI %+0.3f..%+0.3f), p=%.4g, n=%d\n", label, est, lo, hi, pv, n_model))
  }
  data.frame(outcome = label, scale = scale, unit = unit, est = est, ci_lo = lo, ci_hi = hi, p = pv, n = n_model)
}

covs <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH3) + BMXBMI + log1p(LBXCOT) + factor(cycle)"
res <- rbind(
  fit_outcome(as.formula(paste("log(LBXSATSI) ~ log2pfoa +", covs)), "ALT", "log"),
  fit_outcome(as.formula(paste("log(LBXSASSI) ~ log2pfoa +", covs)), "AST", "log"),
  fit_outcome(as.formula(paste("log(LBXSGTSI) ~ log2pfoa +", covs)), "GGT", "log"),
  fit_outcome(as.formula(paste("FLI ~ log2pfoa +", covs)), "FLI", "linear")
)
print(res)
write.csv(res, file.path(out_dir, "nhanes_results.csv"), row.names = FALSE)

# ALT quartile trend (secondary)
q <- cut(d2$pfoa, breaks = quantile(d2$pfoa, probs = 0:4 / 4), include.lowest = TRUE, labels = 1:4)
d2$pfoa_q <- q
des2 <- update(des, pfoa_q = q)
m_q <- svyglm(log(LBXSATSI) ~ pfoa_q + RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH3) + BMXBMI + log1p(LBXCOT) + factor(cycle), design = des2)
co <- summary(m_q)$coefficients
q_rows <- data.frame()
for (lv in c("2", "3", "4")) {
  nm <- paste0("pfoa_q", lv)
  if (nm %in% rownames(co)) {
    b <- co[nm, 1]; se <- co[nm, 2]; pv <- co[nm, 4]
    q_rows <- rbind(q_rows, data.frame(quartile = paste0("Q", lv, " vs Q1"),
                                       pct_change = (2^b - 1) * 100, ci_lo = (2^(b - 1.96 * se) - 1) * 100,
                                       ci_hi = (2^(b + 1.96 * se) - 1) * 100, p = pv))
  }
}
print(q_rows)
write.csv(q_rows, file.path(out_dir, "nhanes_alt_quartiles.csv"), row.names = FALSE)

write.csv(res, file.path(out_dir, "nhanes_results.csv"), row.names = FALSE)

keep <- c("SEQN", "cycle", "SDMVPSU", "SDMVSTRA", "w", "RIDAGEYR", "RIAGENDR", "RIDRETH3",
          "BMXBMI", "BMXWAIST", "LBXCOT", "LBXNFOA", "LBXBFOA", "pfoa",
          "LBXSATSI", "LBXSASSI", "LBXSGTSI", "LBXSTR", "FLI")
write.csv(d2[, keep], file.path(out_dir, "nhanes_analysis_dataset.csv"), row.names = FALSE)
cat("saved outputs to", out_dir, "\n")
cat("NHANES-R DONE\n")
