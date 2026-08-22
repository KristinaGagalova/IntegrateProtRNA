#!/usr/bin/env Rscript
## =============================================================================
## 05_concordance_archetypes.R -- cross-omics concordance and regulatory
##                                archetype classification (Universe B)
##
## THE CENTRAL IDEA
## ----------------
## Naive quadrant classification ("RNA up, protein flat => post-transcriptional
## repression") is wrong: a gene whose protein simply has a long half-life
## shows exactly that pattern with no regulation whatsoever. Separating the two
## is the reason for running a time course at all.
##
## Write everything as fold-change relative to the time-matched control. If the
## control is at steady state and the treatment changes ONLY transcript level,
## the protein must follow
##
##     dp/dt = kd * ( r(t) - p(t) ),      p(0) = 1                        (M0)
##
## with r(t) = 2^logFC_RNA(t), p(t) = 2^logFC_protein(t), and kd the protein
## degradation rate as the single free parameter. Post-transcriptional
## regulation is then a DEPARTURE from M0, which we test as a nested model:
##
##     logFC_protein(t) = a * log2 p(t; kd),   a free                     (M1)
##
## a is the amplitude-scaling factor. a = 1 is M0. a < 1 means the protein
## response is attenuated relative to what its mRNA predicts (buffering,
## accelerated degradation, translational repression); a > 1 means amplified
## (stabilisation, increased translational efficiency); a < 0 means the protein
## moves opposite to its mRNA. Both models are fitted by inverse-variance
## weighted least squares using the moderated SEs from script 04, and compared
## by a 1-df likelihood-ratio test. With n = 3, weighting by each gene's own
## measurement error is what stops noise being read as regulation.
##
## IDENTIFIABILITY -- stated explicitly rather than hidden
## ------------------------------------------------------
## A half-life much longer than the time course cannot be estimated from the
## time course. When the best fit runs to the edge of the identifiable range,
## attenuation (a < 1) and slow turnover (large 1/kd) are confounded: both
## flatten the predicted trajectory. Those genes are reported as
## "kinetics_limited" rather than being silently called lag or buffering. The
## fix in a real study is an informative prior on kd -- published half-lives,
## or a pulse-SILAC arm -- not a longer grid search.
##
## Also computed:
##   - per-timepoint cross-omics concordance, disattenuated for measurement
##     error (the raw correlation is biased toward zero by n=3 noise)
##   - total-least-squares slope (both axes carry error, so OLS is wrong)
##   - validation against the simulator ground truth when present
##
## Usage: Rscript R/05_concordance_archetypes.R [config/config.yaml]
## =============================================================================

source("R/utils.R")

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- load_config(if (length(args)) args[1] else "config/config.yaml")
dedir  <- file.path(cfg$paths$results, "de")
outdir <- ensure_dir(cfg$paths$results, "concordance")
figdir <- ensure_dir(cfg$paths$figures)
set.seed(cfg$project$seed)

tps   <- cfg$design$timepoints
CC    <- cfg$concordance
FDR   <- CC$sig_fdr
MINLF <- CC$min_abs_lfc

## ------------------------------------------------------------------ load ----
Rl <- read_mat(file.path(dedir, "rna_lfc.tsv"));  Rs <- read_mat(file.path(dedir, "rna_se.tsv"))
Pl <- read_mat(file.path(dedir, "prot_lfc.tsv")); Ps <- read_mat(file.path(dedir, "prot_se.tsv"))
Rf <- read.delim(file.path(dedir, "rna_Ftest_any.tsv"),  stringsAsFactors = FALSE)
Pf <- read.delim(file.path(dedir, "prot_Ftest_any.tsv"), stringsAsFactors = FALSE)

g  <- intersect(rownames(Rl), rownames(Pl))
Rl <- Rl[g, , drop = FALSE]; Rs <- Rs[g, , drop = FALSE]
Pl <- Pl[g, , drop = FALSE]; Ps <- Ps[g, , drop = FALSE]
rownames(Rf) <- Rf$feature_id; rownames(Pf) <- Pf$feature_id
r_any <- Rf[g, "FDR"]; p_any <- Pf[g, "FDR"]
log_step("Universe B genes with both layers: ", length(g))

## Centre on t0: the kinetic model starts from p(0) = r(0) = 1, and a non-zero
## baseline contrast is a batch artefact, not biology.
Rl <- Rl - Rl[, 1]; Pl <- Pl - Pl[, 1]

## ================================================== 1. per-timepoint stats ==
tls_slope <- function(x, y) {          # orthogonal regression: both axes noisy
  ok <- is.finite(x) & is.finite(y)
  v  <- svd(scale(cbind(x[ok], y[ok]), TRUE, FALSE))$v
  -v[1, 2] / v[2, 2]
}

## Reliability = (observed variance - mean error variance) / observed variance;
## the corrected correlation divides by sqrt(rel_x * rel_y) (Spearman 1904).
disattenuate <- function(x, y, sex, sey) {
  ok <- is.finite(x) & is.finite(y)
  rx <- max(1e-6, (var(x[ok]) - mean(sex[ok]^2)) / var(x[ok]))
  ry <- max(1e-6, (var(y[ok]) - mean(sey[ok]^2)) / var(y[ok]))
  c(rel_rna = rx, rel_prot = ry,
    r_corrected = min(1, cor(x[ok], y[ok]) / sqrt(rx * ry)))
}

conc <- do.call(rbind, lapply(seq_along(tps), function(k) {
  x <- Rl[, k]; y <- Pl[, k]
  ## t0 is the centring anchor: both vectors are identically 0 by
  ## construction, so the correlation there is undefined, not zero.
  if (sd(x) < 1e-12 || sd(y) < 1e-12)
    return(data.frame(time = tps[k], pearson = NA_real_, spearman = NA_real_,
                      ci_lo = NA_real_, ci_hi = NA_real_, tls_slope = NA_real_,
                      rel_rna = NA_real_, rel_prot = NA_real_,
                      pearson_disattenuated = NA_real_))
  bs <- replicate(min(CC$n_boot, 2000),
                  { i <- sample(length(x), replace = TRUE); cor(x[i], y[i]) })
  d <- disattenuate(x, y, Rs[, k], Ps[, k])
  data.frame(time = tps[k], pearson = cor(x, y),
             spearman = cor(x, y, method = "spearman"),
             ci_lo = quantile(bs, .025), ci_hi = quantile(bs, .975),
             tls_slope = tls_slope(x, y),
             rel_rna = d[["rel_rna"]], rel_prot = d[["rel_prot"]],
             pearson_disattenuated = d[["r_corrected"]])
}))
rownames(conc) <- NULL
cat("\n=== cross-omics concordance of logFC, per timepoint ===\n")
print(round(conc, 3))
write_tsv(conc, file.path(outdir, "concordance_by_timepoint.tsv"))

## ========================================== 2. nested kinetic model fitting =
grid   <- seq(0, max(tps), by = 0.25)
gi     <- match(tps, grid)
r_grid <- t(apply(Rl, 1, function(v) approx(tps, v, grid, rule = 2)$y))
r_lin  <- 2^r_grid

## Identifiability bound: half-lives beyond the time course cannot be
## estimated from it. The grid extends past the bound so the optimiser is not
## artificially truncated, but fits landing above it are flagged.
th_bound <- max(tps)
th_grid  <- exp(seq(log(0.5), log(6 * th_bound), length.out = 80))
kd_grid  <- log(2) / th_grid

## OPTIONAL, and the single most effective way to shrink the unidentifiable
## region: an informative prior on each protein's half-life from an external
## source (published turnover datasets, or a pulse-SILAC arm of your own).
## Supply concordance$half_life_prior_file as a TSV with columns
##   gene_id, t_half_h [, sd_log]
## and the grid search gains a log-normal penalty, so a gene is only allowed to
## claim a very long half-life if its own data really demand it. Without this,
## a free kd can absorb almost any amount of attenuation.
prior_pen <- matrix(0, nrow(Pl), 1)
hp_file <- CC$half_life_prior_file
if (!is.null(hp_file) && nzchar(hp_file) && file.exists(hp_file)) {
  hp <- read.delim(hp_file, stringsAsFactors = FALSE)
  sd_log <- if ("sd_log" %in% names(hp)) hp$sd_log else rep(0.7, nrow(hp))
  mu  <- log(hp$t_half_h[match(g, hp$gene_id)])
  sdl <- sd_log[match(g, hp$gene_id)]
  sdl[is.na(mu)] <- Inf                     # no prior for unmatched genes
  mu[is.na(mu)]  <- 0
  log_step("half-life prior applied to ", sum(is.finite(sdl)), " / ", length(g),
           " genes")
  ## Penalty matrix: genes x half-life grid, a log-normal prior on t_half.
  prior_pen <- sapply(log(th_grid), function(lt) ((lt - mu) / sdl)^2)
  prior_pen[!is.finite(prior_pen)] <- 0
} else {
  prior_pen <- matrix(0, length(g), length(th_grid))
  if (!is.null(hp_file) && nzchar(hp_file))
    log_step("half-life prior file not found, proceeding uninformed: ", hp_file)
}

W    <- 1 / pmax(Ps, 0.05)^2                     # inverse-variance weights
SST  <- rowSums(W * Pl^2)
n    <- nrow(Pl); J <- length(kd_grid)
ssr0 <- matrix(NA_real_, n, J)   # M0: a fixed at 1
ssr1 <- matrix(NA_real_, n, J)   # M1: a profiled out analytically
ahat <- matrix(NA_real_, n, J)
preds <- vector("list", J)

log_step("fitting nested kinetic models over ", J, " half-lives")
for (j in seq_len(J)) {
  kd <- kd_grid[j]
  ## ks = kd makes the steady state equal r, and P0 = 1 is the shared
  ## pre-treatment state: exactly dp/dt = kd*(r - p).
  pl <- log2(pmax(integrate_protein_mat(grid, r_lin, ks = kd, kd = kd,
                                        P0 = rep(1, n))[, gi, drop = FALSE], 1e-9))
  preds[[j]]  <- pl
  ssr0[, j]   <- rowSums(W * (pl - Pl)^2) + prior_pen[, j]
  ## Weighted least squares through the origin: a_hat has a closed form, so
  ## profiling it out costs nothing.
  num <- rowSums(W * pl * Pl); den <- rowSums(W * pl^2)
  a   <- ifelse(den > 1e-12, num / den, NA_real_)
  ahat[, j] <- a
  ssr1[, j] <- SST - ifelse(den > 1e-12, num^2 / den, 0) + prior_pen[, j]
}

j0 <- max.col(-ssr0, ties.method = "first")
j1 <- max.col(-ssr1, ties.method = "first")
ix <- function(j) cbind(seq_len(n), j)
S0 <- ssr0[ix(j0)]; S1 <- ssr1[ix(j1)]
t_half_M0 <- th_grid[j0]; t_half_M1 <- th_grid[j1]
a_hat <- ahat[ix(j1)]
pred0 <- t(vapply(seq_len(n), function(i) preds[[j0[i]]][i, ], numeric(length(tps))))
colnames(pred0) <- colnames(Pl)

## Likelihood-ratio test of M1 (amplitude free) against M0 (amplitude = 1).
## Under Gaussian errors with known weights, the LRT statistic is the drop in
## weighted SSR, on 1 df.
lrt  <- pmax(S0 - S1, 0)
p_lrt <- pchisq(lrt, 1, lower.tail = FALSE)
fdr_lrt <- p.adjust(p_lrt, "BH")

## Absolute goodness of fit of M0, reported alongside (df = T - anchor - kd).
p_gof <- pchisq(S0, max(1, length(tps) - 2), lower.tail = FALSE)

## ======================================================= 3. classification ==
max_abs   <- function(m) apply(abs(m[, -1, drop = FALSE]), 1, max)
rna_resp  <- (r_any < FDR) & (max_abs(Rl) > MINLF)
prot_resp <- (p_any < FDR) & (max_abs(Pl) > MINLF)
regulated <- rna_resp & (fdr_lrt < FDR)
identifiable <- t_half_M0 <= th_bound
slow      <- t_half_M0 > 0.5 * th_bound

archetype <- rep("unchanged", n)
archetype[!rna_resp &  prot_resp] <- "protein_only"
archetype[rna_resp & !regulated &  identifiable & !slow] <- "concordant"
archetype[rna_resp & !regulated &  identifiable &  slow] <- "kinetic_lag"
archetype[rna_resp & !regulated & !identifiable]         <- "kinetics_limited"
archetype[regulated & a_hat >= 1]            <- "amplified"
archetype[regulated & a_hat < 1 & a_hat > 0] <- "buffered"
archetype[regulated & a_hat <= 0]            <- "anticorrelated"

res <- data.frame(
  gene_id = g, archetype = archetype,
  rna_responds = rna_resp, prot_responds = prot_resp, regulated = regulated,
  amplitude_a = round(a_hat, 3),
  t_half_M0_h = round(t_half_M0, 2), t_half_M1_h = round(t_half_M1, 2),
  identifiable = identifiable,
  lrt = round(lrt, 3), lrt_p = p_lrt, lrt_fdr = fdr_lrt,
  gof_M0_p = p_gof, rna_any_fdr = r_any, prot_any_fdr = p_any,
  stringsAsFactors = FALSE)
res <- cbind(res,
  setNames(as.data.frame(round(Rl, 3)),    paste0("lfc_rna_",   colnames(Rl))),
  setNames(as.data.frame(round(Pl, 3)),    paste0("lfc_prot_",  colnames(Pl))),
  setNames(as.data.frame(round(pred0, 3)), paste0("pred_prot_", colnames(Pl))))
write_tsv(res, file.path(outdir, "archetypes.tsv"))

cat("\n=== inferred regulatory archetypes ===\n"); print(table(archetype))
cat("\n=== median implied half-life (h), M0 fit ===\n")
print(round(tapply(t_half_M0, archetype, median), 1))
cat("\n=== median amplitude factor a ===\n")
print(round(tapply(a_hat, archetype, median), 2))

## ==================================================== 4. truth validation ===
truth_path <- file.path(cfg$paths$data_sim, "truth.tsv")
if (file.exists(truth_path)) {
  tr <- read.delim(truth_path, stringsAsFactors = FALSE); rownames(tr) <- tr$gene_id
  tt <- tr[g, ]
  cat("\n=== TRUE archetype (rows) vs INFERRED (cols) ===\n")
  cm <- table(true = tt$archetype, inferred = archetype)
  print(cm)

  reg_true  <- tt$archetype %in% c("buffered", "mrna_only", "anticorrelated")
  lag_true  <- tt$archetype == "lag_only"
  conc_true <- tt$archetype == "concordant"
  called_reg <- archetype %in% c("buffered", "amplified", "anticorrelated")

  cat(sprintf(
    "\nregulated genes called regulated (recall)          : %.1f%% (n=%d)",
    100 * mean(called_reg[reg_true]), sum(reg_true)))
  cat(sprintf(
    "\ntrue lag_only wrongly called regulated (the trap)  : %.1f%% (n=%d)",
    100 * mean(called_reg[lag_true]), sum(lag_true)))
  cat(sprintf(
    "\ntrue concordant wrongly called regulated (FPR)     : %.1f%% (n=%d)",
    100 * mean(called_reg[conc_true]), sum(conc_true)))
  cat(sprintf(
    "\nunchanged genes wrongly called regulated (FPR)     : %.1f%% (n=%d)\n",
    100 * mean(called_reg[tt$archetype == "unchanged"]),
    sum(tt$archetype == "unchanged")))
  cat(sprintf("genes parked in kinetics_limited (unidentifiable)   : %d (%.1f%%)\n",
              sum(archetype == "kinetics_limited"),
              100 * mean(archetype == "kinetics_limited")))

  ok <- rna_resp & is.finite(tt$t_half_h)
  if (sum(ok) > 20)
    cat(sprintf("Spearman(implied t_half, true t_half), responders   : %.2f (n=%d)\n",
                cor(t_half_M0[ok], tt$t_half_h[ok], method = "spearman"), sum(ok)))
  ## The amplitude factor should recover the true attenuation directly.
  bt <- tt$archetype == "buffered" & rna_resp
  if (sum(bt) > 10)
    cat(sprintf("median a for true buffered / true concordant       : %.2f / %.2f\n",
                median(a_hat[bt]), median(a_hat[conc_true & rna_resp])))
  write_tsv(as.data.frame.matrix(cm), file.path(outdir, "truth_confusion.tsv"))
}

## ============================================================== 5. figures ==
pal <- c(unchanged = "#BFBFBF", concordant = "#2C6FBB", kinetic_lag = "#7FB2E5",
         kinetics_limited = "#C9CBA3", buffered = "#E5A11F",
         amplified = "#4C9F70", protein_only = "#8E7CC3",
         anticorrelated = "#D1495B")
acol <- unname(pal[archetype]); acol[is.na(acol)] <- "#999999"
tcol <- adjustcolor(acol, alpha.f = 0.55)

open_pdf(file.path(figdir, "fig03_concordance.pdf"), 13, 9)
op <- par(mfrow = c(3, 4), mar = c(4.2, 4.2, 3.2, 1))

for (k in seq_along(tps)) {                      # (a) concordance scatters
  if (is.na(conc$pearson[k])) next
  plot(Rl[, k], Pl[, k], pch = 16, cex = 0.4, col = tcol,
       xlab = "RNA log2FC", ylab = "protein log2FC",
       main = sprintf("t%g   r=%.2f  (corrected %.2f)", tps[k],
                      conc$pearson[k], conc$pearson_disattenuated[k]))
  abline(0, 1, lty = 3, col = "grey50"); abline(h = 0, v = 0, col = "grey85")
  abline(0, conc$tls_slope[k], col = "#D1495B", lwd = 2)
}

cc <- conc[!is.na(conc$pearson), , drop = FALSE]  # (b) coupling over time
plot(cc$time, cc$pearson, type = "b", pch = 16, lwd = 2, col = "#2C6FBB",
     ylim = range(0, cc$ci_lo, cc$ci_hi, cc$pearson_disattenuated, na.rm = TRUE),
     xlab = "time (h)", ylab = "cross-omics correlation",
     main = "Global RNA-protein coupling")
arrows(cc$time, cc$ci_lo, cc$time, cc$ci_hi, angle = 90, code = 3,
       length = 0.03, col = "#2C6FBB")
lines(cc$time, cc$pearson_disattenuated, type = "b", pch = 17, lty = 2,
      lwd = 2, col = "#D1495B")
legend("topleft", bty = "n", cex = 0.75, lwd = 2, pch = c(16, 17),
       col = c("#2C6FBB", "#D1495B"),
       legend = c("observed", "measurement-error corrected"))

tb <- sort(table(archetype), decreasing = TRUE)   # (c) composition
barplot(tb, las = 2, cex.names = 0.6, col = unname(pal[names(tb)]),
        ylab = "genes", main = "Regulatory archetypes")

## (d) the decision surface: amplitude factor vs evidence against a = 1
plot(pmax(pmin(a_hat, 3), -2), -log10(pmax(p_lrt, 1e-12)), pch = 16, cex = 0.45,
     col = tcol, xlab = "amplitude factor a  (1 = RNA-predicted)",
     ylab = "-log10 p (LRT vs a = 1)", main = "Regulation test")
abline(v = 1, lty = 3); abline(h = -log10(FDR), lty = 2, col = "grey40")

## (e) identifiability: where the kinetic fit runs out of information
plot(pmin(t_half_M0, 6 * th_bound), pmax(pmin(a_hat, 3), -2), log = "x",
     pch = 16, cex = 0.45, col = tcol,
     xlab = "implied half-life (h), M0", ylab = "amplitude factor a",
     main = "Attenuation vs turnover\n(right of the line is unidentifiable)")
abline(v = th_bound, lwd = 2, col = "#D1495B"); abline(h = 1, lty = 3)

if (file.exists(truth_path)) {                    # (f) half-life recovery
  ok <- rna_resp
  plot(tt$t_half_h[ok], t_half_M0[ok], log = "xy", pch = 16, cex = 0.5,
       col = tcol[ok], xlab = "true half-life (h)", ylab = "implied half-life (h)",
       main = "Half-life recovery")
  abline(0, 1, col = "#D1495B", lwd = 2); abline(h = th_bound, lty = 2)
}

for (a in c("concordant", "kinetic_lag", "buffered", "anticorrelated",
            "protein_only")) {                    # (g) exemplars
  i <- which(archetype == a)
  if (!length(i)) { plot.new(); title(paste0(a, ": none")); next }
  i  <- i[which.max(abs(Rl[i, ncol(Rl)]) + abs(Pl[i, ncol(Pl)]))]
  yl <- range(0, Rl[i, ], Pl[i, ], pred0[i, ])
  plot(tps, Rl[i, ], type = "b", pch = 16, lwd = 2, col = "#2C6FBB", ylim = yl,
       xlab = "time (h)", ylab = "log2FC",
       main = sprintf("%s\n%s  (a=%.2f, t1/2=%.0fh)", a, g[i], a_hat[i], t_half_M0[i]))
  lines(tps, Pl[i, ], type = "b", pch = 17, lwd = 2, col = "#D1495B")
  lines(tps, pred0[i, ], lty = 2, lwd = 2, col = "grey40")
  abline(h = 0, col = "grey85")
  legend("topleft", bty = "n", cex = 0.6, lwd = 2, pch = c(16, 17, NA),
         lty = c(1, 1, 2), col = c("#2C6FBB", "#D1495B", "grey40"),
         legend = c("RNA", "protein", "RNA-predicted (M0)"))
}
par(op); dev.off()
log_step("05_concordance_archetypes done -> ", outdir)
