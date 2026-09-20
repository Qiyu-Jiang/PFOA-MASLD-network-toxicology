# NHANES sensitivity analysis: additional adjustment for heavy alcohol use,
# diabetes, and lipid-lowering medication (revision round 2)
# Rscript nhanes_sensitivity.R
suppressPackageStartupMessages({ library(survey); library(foreign) })

base <- "C:\\Users\\user\\workspace\\PFOA_NAFLD_project"
nha <- "C:\\Users\\user\\workspace\\nhanes"
out_dir <- file.path(base, "04_analysis", "nhanes")

d <- read.csv(file.path(out_dir, "nhanes_analysis_dataset.csv"))
read1 <- function(f) read.xport(file.path(nha, f))

# ---------- alcohol ----------
a1 <- read1("ALQ_I.XPT"); a2 <- read1("ALQ_J.XPT")
cat("ALQ_I vars:", paste(grep("^ALQ", names(a1), value = TRUE), collapse = ","), "\n")
cat("ALQ_J vars:", paste(grep("^ALQ", names(a2), value = TRUE), collapse = ","), "\n")
alqv <- if ("ALQ151" %in% names(a1) && "ALQ151" %in% names(a2)) "ALQ151" else "ALQ130"
mk_heavy <- function(a) {
  v <- a[[alqv]]
  h <- if (alqv == "ALQ151") {
    ifelse(is.na(v), NA, ifelse(v == 1, 1, ifelse(v %in% c(7, 9), NA, 0)))
  } else {
    ifelse(is.na(v), NA, ifelse(v >= 2, 1, 0))
  }
  data.frame(SEQN = a$SEQN, heavy_alc = h)
}
heavy <- rbind(mk_heavy(a1), mk_heavy(a2))

# ---------- diabetes ----------
di1 <- read1("DIQ_I.XPT"); di2 <- read1("DIQ_J.XPT")
cat("DIQ vars sample:", paste(head(grep("^DIQ", names(di1), value = TRUE), 8), collapse = ","), "\n")
mk_diab <- function(a) {
  v <- a$DIQ010
  data.frame(SEQN = a$SEQN, diab = ifelse(is.na(v), NA, ifelse(v == 1, 1, ifelse(v %in% c(7, 9), NA, 0))))
}
dia <- rbind(mk_diab(di1), mk_diab(di2))

# ---------- lipid-lowering medication ----------
r1 <- read1("RXQ_RX_I.XPT"); r2 <- read1("RXQ_RX_J.XPT")
rxcol <- if ("RXDDRUG" %in% names(r1)) "RXDDRUG" else grep("DRUG", names(r1), value = TRUE)[1]
rxcol2 <- if ("RXDDRUG" %in% names(r2)) "RXDDRUG" else rxcol
cat("RX col:", rxcol, "/", rxcol2, "\n")
mk_rx <- function(a, col) data.frame(SEQN = a$SEQN, drug = toupper(as.character(a[[col]])))
rx <- rbind(mk_rx(r1, rxcol), mk_rx(r2, rxcol2))
patterns <- c("STATIN", "EZETIMIBE", "FIBRATE", "GEMFIBROZIL", "FENOFIBR", "CHOLESTYRAMINE",
              "COLESTIPOL", "COLESEVELAM", "NIACIN", "OMEGA-3", "ICOSAPENT", "EVOLOCUMAB",
              "ALIROCUMAB", "BEMPEDOIC", "INCLISIRAN", "LIPITOR", "CRESTOR", "ZOCOR",
              "PRAVACHOL", "MEVACOR", "LESCOL", "LIVALO", "ZETIA", "TRICOR", "LOPID",
              "WELCHOL", "NIASPAN", "LOVAZA", "VASCEPA", "REPATHA", "PRALUENT", "NEXLETOL", "LEQVIO")
rx$lipid <- grepl(paste(patterns, collapse = "|"), rx$drug)
lipid <- aggregate(lipid ~ SEQN, data = rx, FUN = any)

d <- merge(d, heavy, by = "SEQN", all.x = TRUE)
d <- merge(d, dia, by = "SEQN", all.x = TRUE)
d <- merge(d, lipid, by = "SEQN", all.x = TRUE)
d$lipid[is.na(d$lipid)] <- FALSE
cat("heavy alcohol yes:", sum(d$heavy_alc == 1, na.rm = TRUE),
    "| diabetes yes:", sum(d$diab == 1, na.rm = TRUE),
    "| lipid-med users:", sum(d$lipid, na.rm = TRUE), "\n")

d2 <- d[!is.na(d$heavy_alc) & !is.na(d$diab), ]
d2$strata <- interaction(d2$cycle, d2$SDMVSTRA, drop = TRUE)
d2$psu <- interaction(d2$cycle, d2$SDMVPSU, drop = TRUE)
d2$log2pfoa <- log2(d2$pfoa)
des2 <- svydesign(ids = ~psu, strata = ~strata, weights = ~w, nest = TRUE, data = d2)
cat("sensitivity analysis n:", nrow(d2), "\n")

covs <- "RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH3) + BMXBMI + log1p(LBXCOT) + factor(cycle)"

fit2 <- function(f, label, scale, tag) {
  m <- svyglm(f, design = des2)
  co <- summary(m)$coefficients
  b <- co["log2pfoa", 1]; se <- co["log2pfoa", 2]; pv <- co["log2pfoa", 4]
  n_model <- tryCatch(length(m$residuals), error = function(e) NA_integer_)
  if (scale == "log") {
    est <- (2^b - 1) * 100; lo <- (2^(b - 1.96 * se) - 1) * 100; hi <- (2^(b + 1.96 * se) - 1) * 100
    unit <- "percent"
  } else {
    est <- b; lo <- b - 1.96 * se; hi <- b + 1.96 * se; unit <- "points"
  }
  cat(sprintf("%-4s [%s]: %+0.3f %s (95%%CI %+0.3f..%+0.3f), p=%.4g, n=%d\n", label, tag, est, unit, lo, hi, pv, n_model))
  data.frame(analysis = tag, outcome = label, scale = scale, unit = unit, est = est,
             ci_lo = lo, ci_hi = hi, p = pv, n = n_model)
}

res2 <- rbind(
  fit2(as.formula(paste("log(LBXSATSI) ~ log2pfoa +", covs)), "ALT", "log", "primary_same_subset"),
  fit2(as.formula(paste("log(LBXSATSI) ~ log2pfoa + heavy_alc + diab + lipid +", covs)), "ALT", "log", "plus_alcohol_diabetes_lipidmeds"),
  fit2(as.formula(paste("log(LBXSASSI) ~ log2pfoa +", covs)), "AST", "log", "primary_same_subset"),
  fit2(as.formula(paste("log(LBXSASSI) ~ log2pfoa + heavy_alc + diab + lipid +", covs)), "AST", "log", "plus_alcohol_diabetes_lipidmeds"),
  fit2(as.formula(paste("log(LBXSGTSI) ~ log2pfoa +", covs)), "GGT", "log", "primary_same_subset"),
  fit2(as.formula(paste("log(LBXSGTSI) ~ log2pfoa + heavy_alc + diab + lipid +", covs)), "GGT", "log", "plus_alcohol_diabetes_lipidmeds"),
  fit2(as.formula(paste("FLI ~ log2pfoa +", covs)), "FLI", "linear", "primary_same_subset"),
  fit2(as.formula(paste("FLI ~ log2pfoa + heavy_alc + diab + lipid +", covs)), "FLI", "linear", "plus_alcohol_diabetes_lipidmeds")
)
write.csv(res2, file.path(out_dir, "nhanes_sensitivity.csv"), row.names = FALSE)
cat("saved nhanes_sensitivity.csv\n")
cat("NHANES-SENS DONE\n")
