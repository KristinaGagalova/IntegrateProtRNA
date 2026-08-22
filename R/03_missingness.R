#!/usr/bin/env Rscript
## =============================================================================
## 03_missingness.R -- diagnose and handle proteomics missing values
##
## Missingness in MS proteomics is a MIXTURE and must be treated as one:
##   MNAR  left-censored, intensity-dependent (below the detection limit)
##   MAR   stochastic (MS/MS undersampling, ID-transfer failure)
##
## The rule used here to separate them, per protein:
##   an entry missing in ALL replicates of a design cell, while the protein is
##   observed in some other cell   -> MNAR (the protein really is low/absent)
##   an entry missing in only SOME replicates of its cell -> MAR (stochastic)
##
## THE CRITICAL RULE: imputed values are written for Universe A only (factor
## models, clustering, pathway scores). Universe B analyses -- concordance,
## archetypes, kinetics, delay -- must use OBSERVED values only. Imputing
## left-censored entries with a low constant and then correlating RNA against
## protein manufactures anti-correlation and produces false
## "post-transcriptional repression" calls. The mask from 01 enforces this.
##
## Three schemes are written so every headline result can be repeated under
## each of them; that sensitivity panel is a standard reviewer request.
##
## Usage: Rscript R/03_missingness.R [config/config.yaml]
## =============================================================================

source("R/utils.R")
need_pkgs("impute")

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- load_config(if (length(args)) args[1] else "config/config.yaml")
qcdir  <- file.path(cfg$paths$results, "qc")
outdir <- ensure_dir(cfg$paths$results, "qc")
figdir <- ensure_dir(cfg$paths$figures)
set.seed(cfg$project$seed)

prot <- read_mat(file.path(qcdir, "protein_norm.tsv"))
cd   <- read.delim(file.path(qcdir, "protein_coldata.tsv"), stringsAsFactors = FALSE)
cd$cell <- factor(cd$cell, levels = cell_levels(cfg))
prot <- prot[, cd$sample_id, drop = FALSE]
M    <- cfg$missingness

log_step("protein matrix ", nrow(prot), "x", ncol(prot),
         sprintf("; missing %.1f%%", 100 * mean(is.na(prot))))

## ----------------------------------------------------------- 1. diagnose ---
mean_int  <- rowMeans(prot, na.rm = TRUE)
miss_rate <- rowMeans(is.na(prot))
ok  <- is.finite(mean_int)
rho <- suppressWarnings(cor(mean_int[ok], miss_rate[ok], method = "spearman"))
log_step("Spearman(mean intensity, missing rate) = ", round(rho, 3),
         if (rho < -0.3) "  -> MNAR-dominant" else "  -> MAR-dominant")

## Per-entry MNAR/MAR labels (see header for the rule).
cells        <- levels(cd$cell)
cell_idx     <- lapply(cells, function(cl) which(cd$cell == cl))
n_obs_cell   <- vapply(cell_idx, function(ii)
  rowSums(!is.na(prot[, ii, drop = FALSE])), numeric(nrow(prot)))
colnames(n_obs_cell) <- cells

is_na  <- is.na(prot)
lab    <- matrix("obs", nrow(prot), ncol(prot), dimnames = dimnames(prot))
for (k in seq_along(cells)) {
  ii   <- cell_idx[[k]]
  gone <- n_obs_cell[, k] == 0                       # whole cell missing
  lab[gone,  ii] <- "MNAR"
  lab[!gone, ii][is_na[!gone, ii]] <- "MAR"
}
lab[!is_na] <- "obs"
tab <- table(factor(lab, c("obs", "MAR", "MNAR")))
log_step("entries: observed=", tab[["obs"]], " MAR=", tab[["MAR"]],
         " MNAR=", tab[["MNAR"]])

## ---------------------------------------------------------- 2. imputers ----
## Down-shifted normal (Perseus convention): per SAMPLE, draw from a normal
## centred `downshift_sd` SDs below the observed mean with `width` * SD.
impute_downshift <- function(m, shift = 1.8, width = 0.3) {
  out <- m
  for (j in seq_len(ncol(m))) {
    v  <- m[, j]; nas <- is.na(v)
    if (!any(nas)) next
    mu <- mean(v, na.rm = TRUE); sdv <- sd(v, na.rm = TRUE)
    out[nas, j] <- rnorm(sum(nas), mu - shift * sdv, width * sdv)
  }
  out
}

## QRILC if imputeLCMD is available, otherwise fall back to the down-shifted
## normal, which targets the same left-censored region.
impute_mnar <- function(m) {
  if (M$mnar_method == "QRILC" && has_pkg("imputeLCMD")) {
    log_step("MNAR imputation: QRILC (imputeLCMD)")
    return(imputeLCMD::impute.QRILC(m)[[1]])
  }
  if (M$mnar_method == "MinProb" && has_pkg("imputeLCMD")) {
    log_step("MNAR imputation: MinProb (imputeLCMD)")
    return(imputeLCMD::impute.MinProb(m))
  }
  log_step("MNAR imputation: down-shifted normal",
           if (M$mnar_method != "downshift") "  [imputeLCMD unavailable, fallback]" else "")
  impute_downshift(m, M$downshift_sd, M$downshift_width)
}

impute_mar <- function(m) {
  if (M$mar_method == "missForest" && has_pkg("missForest")) {
    log_step("MAR imputation: missForest")
    return(missForest::missForest(m)$ximp)
  }
  log_step("MAR imputation: kNN (k=10)")
  ## impute.knn works on rows=genes; it needs some observed values per row,
  ## which the 01 validity filter already guarantees.
  suppressWarnings(impute::impute.knn(as.matrix(m), k = 10,
                                      rowmax = 0.95, colmax = 0.95)$data)
}

## -------------------------------------------------------- 3. mixed scheme --
## Impute the two mechanisms separately, then splice by label. Using kNN on
## MNAR entries would borrow from samples where the protein IS present and
## erase the very on/off biology the validity filter was designed to keep.
mar_full  <- impute_mar(prot)
mnar_full <- impute_mnar(prot)
mixed <- prot
mixed[lab == "MAR"]  <- mar_full[lab == "MAR"]
mixed[lab == "MNAR"] <- mnar_full[lab == "MNAR"]
stopifnot(!anyNA(mixed))

## ------------------------------------------------- 4. sensitivity schemes --
schemes <- list()
if ("mixed" %in% M$sensitivity_schemes)      schemes$mixed <- mixed
if ("downshift" %in% M$sensitivity_schemes)  schemes$downshift <-
  impute_downshift(prot, M$downshift_sd, M$downshift_width)
if ("complete_case" %in% M$sensitivity_schemes) {
  cc <- prot[complete.cases(prot), , drop = FALSE]
  log_step("complete-case retains ", nrow(cc), " / ", nrow(prot), " proteins")
  schemes$complete_case <- cc
}
if ("knn" %in% M$sensitivity_schemes)        schemes$knn <- mar_full

for (nm in names(schemes))
  write_mat(schemes[[nm]], file.path(outdir, paste0("protein_imputed_", nm, ".tsv")),
            "group_id")
write_mat(lab, file.path(outdir, "protein_missing_labels.tsv"), "group_id")

## ------------------------------------------------------------- 5. figure ---
open_pdf(file.path(figdir, "fig01c_missingness.pdf"), 11, 7.5)
op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.5, 1))

plot(mean_int[ok], miss_rate[ok], pch = 16, cex = 0.35, col = "#00000044",
     xlab = "mean log2 intensity", ylab = "missing rate",
     main = sprintf("MNAR diagnostic\nSpearman rho = %.2f", rho))
lines(lowess(mean_int[ok], miss_rate[ok], f = 0.4), col = "#D1495B", lwd = 2)

barplot(tab[c("MAR", "MNAR")], col = c("#7FB2E5", "#D1495B"),
        ylab = "missing entries", main = "Missingness mechanism split")

d_obs <- density(prot[!is_na])
plot(d_obs, lwd = 2, main = "Observed vs imputed distribution",
     xlab = "log2 intensity", ylim = c(0, max(d_obs$y) * 1.35))
lines(density(mixed[lab == "MNAR"]), col = "#D1495B", lwd = 2)
lines(density(mixed[lab == "MAR"]),  col = "#7FB2E5", lwd = 2)
legend("topright", bty = "n", cex = 0.8, lwd = 2,
       col = c("black", "#D1495B", "#7FB2E5"),
       legend = c("observed", "MNAR imputed", "MAR imputed"))

## PCA under each scheme: if the sample layout moves between schemes, the
## imputation is driving the structure and must be reported as such.
for (nm in names(schemes)) {
  m <- schemes[[nm]]
  m <- m[order(row_var(m), decreasing = TRUE, na.last = NA)[seq_len(min(500, nrow(m)))], ,
         drop = FALSE]
  p <- prcomp(t(m[complete.cases(m), , drop = FALSE]), scale. = FALSE)
  v <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  plot(p$x[, 1], p$x[, 2], cex = 1.5,
       col = c("#2C6FBB", "#D1495B")[as.integer(factor(cd$treatment))],
       pch = c(15, 16, 17, 18)[as.integer(factor(cd$time))],
       xlab = sprintf("PC1 (%.1f%%)", v[1]), ylab = sprintf("PC2 (%.1f%%)", v[2]),
       main = paste0("PCA: ", nm))
}
par(op); dev.off()

write_tsv(data.frame(
  metric = c("spearman_intensity_vs_missing", "n_obs", "n_MAR", "n_MNAR",
             "missing_rate", "complete_case_proteins"),
  value = c(round(rho, 4), tab[["obs"]], tab[["MAR"]], tab[["MNAR"]],
            round(mean(is_na), 4),
            if (!is.null(schemes$complete_case)) nrow(schemes$complete_case) else NA)),
  file.path(outdir, "missingness_summary.tsv"))

log_step("03_missingness done; schemes written: ", paste(names(schemes), collapse = ", "))
