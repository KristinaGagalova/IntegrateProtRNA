#!/usr/bin/env Rscript
## =============================================================================
## 06_integration_caseB.R -- cross-omics latent integration for UNMATCHED
##                           samples (Case B)
##
## WHY THIS SCRIPT EXISTS, AND WHY IT IS NOT JUST "RUN MOFA"
## --------------------------------------------------------
## MOFA+/MEFISTO, DIABLO and O2PLS all decompose blocks that share their ROWS.
## In Case B the RNA and protein samples are different biological material, so
## there are no shared rows and no information links the two blocks at the
## sample level. Feeding 48 rows (24 RNA + 24 protein) to MOFA with NAs in the
## off-blocks does not integrate anything -- the factors are unidentifiable
## across views and whatever comes out is an artefact.
##
## What IS shared is the DESIGN CELL: treatment x timepoint, 8 of them. So:
##
##   * the integrable signal is BETWEEN-cell covariance -- how the two layers
##     co-vary across treatment and time. That is the biology of interest.
##   * WITHIN-cell covariance (co-variation across replicates) is NOT estimable
##     in Case B, at all, by any method. State this in the manuscript.
##
## Three complementary analyses, all dependency-free (base R):
##
##  1. PLS on design-cell means (8 pseudo-samples per block). Gives joint
##     components whose scores are factor TRAJECTORIES over time per treatment.
##  2. Replicate bootstrap. Resample the 3 replicates within each cell, rebuild
##     the cell means, refit. With n = 3 the point estimate alone is
##     meaningless; every trajectory gets a bootstrap CI.
##  3. Permutation null. Shuffle the design-cell labels of one block and refit.
##     With only 8 pseudo-samples and thousands of features, a PLS will ALWAYS
##     find an apparently strong component; the permutation null is the only
##     thing that says whether it means anything. This is not optional.
##
## Also provided: an exchangeable replicate-PAIRING bootstrap, which randomly
## pairs RNA replicate i with protein replicate j inside each cell to build 24
## pseudo-matched rows. Because the pairing is random, the within-cell
## covariance it induces is zero in expectation, so it estimates exactly the
## between-cell signal -- and it lets you run any sample-matched tool
## (mixOmics/DIABLO) unchanged, with stability-selection frequencies across
## pairings instead of one arbitrary result.
##
## Usage: Rscript R/06_integration_caseB.R [config/config.yaml]
## =============================================================================

source("R/utils.R")

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- load_config(if (length(args)) args[1] else "config/config.yaml")
qcdir  <- file.path(cfg$paths$results, "qc")
outdir <- ensure_dir(cfg$paths$results, "integration")
figdir <- ensure_dir(cfg$paths$figures)
set.seed(cfg$project$seed)

IN   <- cfg$integration
tps  <- cfg$design$timepoints
trts <- cfg$design$treatments
cells <- cell_levels(cfg)

rna  <- read_mat(file.path(qcdir, "rna_vst.tsv"))
prot <- read_mat(file.path(qcdir, "protein_imputed_mixed.tsv"))
cd_r <- read.delim(file.path(qcdir, "rna_coldata.tsv"), stringsAsFactors = FALSE)
cd_p <- read.delim(file.path(qcdir, "protein_coldata.tsv"), stringsAsFactors = FALSE)
cd_r$cell <- factor(cd_r$cell, levels = cells); cd_p$cell <- factor(cd_p$cell, levels = cells)
rna <- rna[, cd_r$sample_id, drop = FALSE]; prot <- prot[, cd_p$sample_id, drop = FALSE]

## --------------------------------------------------------- feature choice ---
hv_r <- top_variable(rna,  IN$n_hvg_rna)
hv_p <- top_variable(prot, IN$n_hvp_prot)
log_step("features: ", length(hv_r), " RNA, ", length(hv_p), " protein")

## Per-block scaling to unit TOTAL variance. Without it the block with more
## features and larger dynamic range (usually RNA) dominates every component.
prep_block <- function(m) {
  m <- t(m)                                   # cells x features
  m <- scale(m, center = TRUE, scale = FALSE)
  m / sqrt(sum(m^2) / nrow(m))
}

## ------------------------------------------------------------- PLS (NIPALS) --
## Plain PLS2. Implemented here rather than pulled from a package so the
## bootstrap and permutation loops stay fast and the algorithm is auditable.
pls2 <- function(X, Y, ncomp = 2, iters = 200, tol = 1e-9) {
  X0 <- X; Y0 <- Y
  n  <- nrow(X)
  Tt <- Uu <- matrix(0, n, ncomp)
  Wt <- Pp <- matrix(0, ncol(X), ncomp); Cc <- Qy <- matrix(0, ncol(Y), ncomp)
  vx <- vy <- rxy <- numeric(ncomp)
  ssx <- sum(X0^2); ssy <- sum(Y0^2)
  for (k in seq_len(ncomp)) {
    u <- Y[, which.max(apply(Y, 2, var))]
    for (it in seq_len(iters)) {
      w <- crossprod(X, u); w <- w / sqrt(sum(w^2))
      t <- X %*% w
      cc <- crossprod(Y, t); cc <- cc / sqrt(sum(cc^2))
      u_new <- Y %*% cc
      if (sum((u_new - u)^2) < tol) { u <- u_new; break }
      u <- u_new
    }
    p <- crossprod(X, t) / drop(crossprod(t))
    q <- crossprod(Y, u) / drop(crossprod(u))
    vx[k]  <- sum((t %*% t(p))^2) / ssx
    vy[k]  <- sum((u %*% t(q))^2) / ssy
    rxy[k] <- cor(drop(t), drop(u))
    ## qy regresses Y on the X-score t (not on u): this is the form needed to
    ## predict Y for a design cell the model has not seen.
    qy <- crossprod(Y, t) / drop(crossprod(t))
    Tt[, k] <- t; Uu[, k] <- u; Wt[, k] <- w; Cc[, k] <- cc
    Pp[, k] <- p; Qy[, k] <- qy
    X <- X - t %*% t(p); Y <- Y - t %*% t(qy)
  }
  list(t = Tt, u = Uu, w = Wt, c = Cc, p = Pp, qy = Qy,
       var_x = vx, var_y = vy, cor = rxy)
}

## Prediction for new rows. B = W (P'W)^-1 Qy' is never materialised: with
## thousands of features that matrix is millions of entries and dominates the
## whole runtime. Associating left-to-right keeps every intermediate at
## ncomp width.
pls_predict <- function(f, xnew, ncomp) {
  W <- f$w[,  seq_len(ncomp), drop = FALSE]
  P <- f$p[,  seq_len(ncomp), drop = FALSE]
  Q <- f$qy[, seq_len(ncomp), drop = FALSE]
  ((xnew %*% W) %*% solve(crossprod(P, W))) %*% t(Q)
}

## Exact row-space reduction.
##
## With 8 design cells both blocks have rank <= 8, so X = U D V' with V
## orthonormal. PLS on X is equivalent to PLS on A = U D (8 x 8), and because V
## is orthonormal all the Euclidean distances that PRESS and TSS are built from
## are preserved exactly. So this is not an approximation -- Q2 computed on A
## equals Q2 computed on the full matrix, at a few thousandths of the cost.
## That is what makes a proper permutation null affordable here.
row_space <- function(M) {
  s <- svd(M)
  k <- sum(s$d > max(dim(M)) * .Machine$double.eps * s$d[1])
  s$u[, seq_len(k), drop = FALSE] %*% diag(s$d[seq_len(k)], k, k)
}

## Leave-one-design-cell-out cross-validation.
##
## THIS is the valid test statistic for Case B, not cor(t, u). With 8
## pseudo-samples and thousands of features a PLS component can be steered onto
## almost any target, so the in-sample cross-block correlation approaches 1
## even for PERMUTED data -- the permutation median printed below demonstrates
## it directly. Q2 asks the one question that cannot be gamed: given the RNA
## profile of a design cell the model has never seen, can it predict that
## cell's protein profile? Centring is recomputed inside each fold, so nothing
## leaks from the held-out cell.
q2_loo <- function(Xr, Yr, ncomp) {
  n <- nrow(Xr)
  press <- numeric(ncomp); tss <- 0
  for (i in seq_len(n)) {
    Xtr <- Xr[-i, , drop = FALSE]; Ytr <- Yr[-i, , drop = FALSE]
    mx  <- colMeans(Xtr); my <- colMeans(Ytr)
    f   <- pls2(sweep(Xtr, 2, mx), sweep(Ytr, 2, my), ncomp = ncomp)
    xte <- matrix(Xr[i, ] - mx, nrow = 1)
    for (k in seq_len(ncomp)) {
      yhat <- drop(pls_predict(f, xte, k)) + my
      press[k] <- press[k] + sum((Yr[i, ] - yhat)^2)
    }
    tss <- tss + sum((Yr[i, ] - my)^2)
  }
  1 - press / tss
}

## Sign of a PLS component is arbitrary; align to a reference before averaging
## bootstrap replicates or the mean collapses toward zero.
align_sign <- function(m, ref) {
  s <- sign(colSums(m * ref)); s[s == 0] <- 1
  sweep(m, 2, s, "*")
}

## ================================================= 1. point-estimate fit ====
Xc <- prep_block(cell_means(rna[hv_r, ], cd_r))
Yc <- prep_block(cell_means(prot[hv_p, ], cd_p))
stopifnot(identical(rownames(Xc), rownames(Yc)))
K  <- min(IN$mofa_factors, nrow(Xc) - 2)
fit <- pls2(Xc, Yc, ncomp = K)
log_step("PLS on ", nrow(Xc), " design-cell pseudo-samples, ", K, " components")

cellmeta <- data.frame(cell = rownames(Xc), stringsAsFactors = FALSE)
cellmeta$treatment <- sub("_t.*$", "", cellmeta$cell)
cellmeta$time <- as.numeric(sub("^.*_t", "", cellmeta$cell))

## ================================================= 2. replicate bootstrap ===
B <- IN$n_boot_cells
log_step("replicate bootstrap, B = ", B)
boot_t <- array(NA_real_, c(nrow(Xc), K, B))
boot_u <- array(NA_real_, c(nrow(Xc), K, B))
boot_vx <- boot_vy <- boot_r <- matrix(NA_real_, B, K)
idx_r <- lapply(cells, function(cl) which(cd_r$cell == cl))
idx_p <- lapply(cells, function(cl) which(cd_p$cell == cl))

for (b in seq_len(B)) {
  ir <- unlist(lapply(idx_r, function(ii) sample(ii, length(ii), TRUE)))
  ip <- unlist(lapply(idx_p, function(ii) sample(ii, length(ii), TRUE)))
  Xb <- prep_block(cell_means(rna[hv_r, ir, drop = FALSE], cd_r[ir, ]))
  Yb <- prep_block(cell_means(prot[hv_p, ip, drop = FALSE], cd_p[ip, ]))
  fb <- pls2(Xb, Yb, ncomp = K)
  boot_t[, , b] <- align_sign(fb$t, fit$t)
  boot_u[, , b] <- align_sign(fb$u, fit$u)
  boot_vx[b, ] <- fb$var_x; boot_vy[b, ] <- fb$var_y; boot_r[b, ] <- abs(fb$cor)
}

## ================================================== 3. permutation null =====
## Shuffle which design cell each protein profile belongs to. This destroys
## the cross-block correspondence while preserving each block's own structure,
## so it is the correct null for "is the shared component real".
Kq    <- min(K, nrow(Xc) - 3)
Ax    <- row_space(Xc); Ay <- row_space(Yc)   # exact, see row_space() above
q2    <- q2_loo(Ax, Ay, Kq)
nperm <- IN$n_perm
log_step("permutation null on Q2, nperm = ", nperm)
perm_q2 <- matrix(NA_real_, nperm, Kq)
perm_r  <- matrix(NA_real_, nperm, K)
for (b in seq_len(nperm)) {
  o <- sample(nrow(Ay))
  perm_q2[b, ] <- q2_loo(Ax, Ay[o, , drop = FALSE], Kq)
  perm_r[b, ]  <- abs(pls2(Ax, Ay[o, , drop = FALSE], ncomp = K)$cor)
}
p_q2  <- vapply(seq_len(Kq), function(k)
  (1 + sum(perm_q2[, k] >= q2[k])) / (nperm + 1), 0)
p_cor <- vapply(seq_len(K), function(k)
  (1 + sum(perm_r[, k] >= abs(fit$cor[k]))) / (nperm + 1), 0)

pad <- function(v, n) c(v, rep(NA_real_, n - length(v)))
comp <- data.frame(
  n_components        = seq_len(K),
  var_rna             = round(fit$var_x, 4),
  var_prot            = round(fit$var_y, 4),
  cor_blocks_INSAMPLE = round(fit$cor, 3),
  cor_perm_median     = round(apply(perm_r, 2, median), 3),
  p_cor               = signif(p_cor, 3),
  Q2_cumulative       = pad(round(q2, 4), K),
  Q2_perm_q95         = pad(round(apply(perm_q2, 2, quantile, .95), 4), K),
  p_Q2                = pad(signif(p_q2, 3), K))
cat("\n=== joint components (design-cell level) ===\n"); print(comp)
cat("\nREAD p_Q2, NOT cor_blocks_INSAMPLE. With 8 pseudo-samples and thousands\n",
    "of features a PLS reaches cor ~ ", round(median(perm_r[, 1]), 2),
    " even on PERMUTED data, so the in-sample\n",
    "correlation carries no evidence. Leave-one-cell-out Q2 cannot be gamed.\n",
    sep = "")
write_tsv(comp, file.path(outdir, "joint_components.tsv"))

## Component trajectories, with bootstrap CIs, in the form that becomes the
## Case B analogue of a MEFISTO factor-trajectory plot.
traj <- do.call(rbind, lapply(seq_len(K), function(k) {
  data.frame(cellmeta, component = k,
             rna_score  = fit$t[, k],
             rna_lo = apply(boot_t[, k, ], 1, quantile, .025, na.rm = TRUE),
             rna_hi = apply(boot_t[, k, ], 1, quantile, .975, na.rm = TRUE),
             prot_score = fit$u[, k],
             prot_lo = apply(boot_u[, k, ], 1, quantile, .025, na.rm = TRUE),
             prot_hi = apply(boot_u[, k, ], 1, quantile, .975, na.rm = TRUE))
}))
write_tsv(traj, file.path(outdir, "component_trajectories.tsv"))
Wout <- fit$w; rownames(Wout) <- hv_r; colnames(Wout) <- paste0("C", seq_len(K))
Cout <- fit$c; rownames(Cout) <- hv_p; colnames(Cout) <- paste0("C", seq_len(K))
write_mat(Wout, file.path(outdir, "loadings_rna.tsv"),     "gene_id")
write_mat(Cout, file.path(outdir, "loadings_protein.tsv"), "group_id")

## ======================================= 4. replicate-pairing bootstrap =====
## Random within-cell pairing gives 24 pseudo-matched rows. The induced
## within-cell cross-block covariance is zero in expectation, so this estimates
## the between-cell signal -- and it lets sample-matched tools run unchanged.
## Reported as SELECTION FREQUENCY across pairings, never one arbitrary run.
Bp <- IN$rep_pairing_B
log_step("replicate-pairing bootstrap, B = ", Bp)
sel_r <- setNames(numeric(length(hv_r)), hv_r)
sel_p <- setNames(numeric(length(hv_p)), hv_p)
topn <- 100
for (b in seq_len(Bp)) {
  ir <- unlist(lapply(idx_r, function(ii) sample(ii)))
  ip <- unlist(lapply(idx_p, function(ii) sample(ii)))
  Xb <- scale(t(rna[hv_r, ir, drop = FALSE]), TRUE, TRUE)
  Yb <- scale(t(prot[hv_p, ip, drop = FALSE]), TRUE, TRUE)
  Xb[!is.finite(Xb)] <- 0; Yb[!is.finite(Yb)] <- 0
  fb <- pls2(Xb / sqrt(sum(Xb^2) / nrow(Xb)), Yb / sqrt(sum(Yb^2) / nrow(Yb)), 1)
  sel_r[order(abs(fb$w[, 1]), decreasing = TRUE)[seq_len(topn)]] <-
    sel_r[order(abs(fb$w[, 1]), decreasing = TRUE)[seq_len(topn)]] + 1
  sel_p[order(abs(fb$c[, 1]), decreasing = TRUE)[seq_len(topn)]] <-
    sel_p[order(abs(fb$c[, 1]), decreasing = TRUE)[seq_len(topn)]] + 1
}
stab <- rbind(
  data.frame(layer = "rna",     feature = names(sel_r), freq = sel_r / Bp),
  data.frame(layer = "protein", feature = names(sel_p), freq = sel_p / Bp))
stab <- stab[order(-stab$freq), ]
write_tsv(stab, file.path(outdir, "pairing_stability_selection.tsv"))
cat("\n=== stability selection (top 10, freq = fraction of pairings in top",
    topn, ") ===\n")
print(head(stab, 10), row.names = FALSE)

## ================================= 5. optional MOFA / MEFISTO on pseudo-cells
if (has_pkg("MOFA2")) {
  log_step("MOFA2 present: fitting on the 8 design-cell pseudo-samples")
  message("  NOTE: 8 samples is far below MOFA's comfortable regime. This is a ",
          "cross-check on the PLS result, not the primary analysis.")
  ok <- try({
    mo <- MOFA2::create_mofa(list(rna = t(Xc), protein = t(Yc)))
    dopt <- MOFA2::get_default_data_options(mo)
    mopt <- MOFA2::get_default_model_options(mo); mopt$num_factors <- min(3, K)
    topt <- MOFA2::get_default_training_options(mo)
    topt$seed <- cfg$project$seed; topt$verbose <- FALSE
    ## MEFISTO: a Gaussian-process prior on the factors over time, so factors
    ## are smooth functions of the time covariate rather than free per cell.
    MOFA2::covariates(mo) <- matrix(cellmeta$time, nrow = 1,
                                    dimnames = list("time", rownames(Xc)))
    mo <- MOFA2::prepare_mofa(mo, data_options = dopt, model_options = mopt,
                              training_options = topt)
    mo <- MOFA2::run_mofa(mo, use_basilisk = TRUE, save_data = FALSE)
    vexp <- MOFA2::get_variance_explained(mo)$r2_per_factor[[1]]
    write_tsv(data.frame(factor = rownames(vexp), vexp),
              file.path(outdir, "mofa_variance_explained.tsv"))
    print(round(vexp, 2))
  }, silent = TRUE)
  if (inherits(ok, "try-error"))
    log_step("MOFA2 run failed (expected with 8 samples); PLS result stands: ",
             conditionMessage(attr(ok, "condition")))
} else {
  log_step("MOFA2 not installed -- skipping the MEFISTO cross-check. ",
           "Install with: Rscript env/install_r_deps.R")
}

## ================================================================= figure ===
open_pdf(file.path(figdir, "fig04_integration_caseB.pdf"), 12, 8)
op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.5, 1))

vm <- rbind(RNA = fit$var_x, Protein = fit$var_y)
bp <- barplot(vm, beside = TRUE, names.arg = paste0("C", seq_len(K)),
              col = c("#2C6FBB", "#D1495B"), ylab = "variance explained",
              main = "Variance captured per block")
legend("topright", bty = "n", fill = c("#2C6FBB", "#D1495B"),
       legend = c("RNA", "Protein"), cex = 0.8)

## Observed cross-block correlation against its permutation null. This panel
## is the honest answer to "is the integration real".
yl <- range(0, q2, perm_q2, na.rm = TRUE)
plot(seq_len(Kq), q2, type = "b", pch = 16, lwd = 2, col = "#2C6FBB", ylim = yl,
     xlab = "number of components", ylab = expression(Q^2~"(leave-one-cell-out)"),
     main = "Predictive shared signal vs null", xaxt = "n")
axis(1, seq_len(Kq)); abline(h = 0, lty = 2, col = "grey50")
for (k in seq_len(Kq)) {
  segments(k, quantile(perm_q2[, k], .05), k, quantile(perm_q2[, k], .95),
           col = "grey60", lwd = 6, lend = 1)
  text(k, yl[1] + 0.04 * diff(yl), sprintf("p=%.3f", p_q2[k]), cex = 0.65)
}
points(seq_len(Kq), q2, pch = 16, col = "#2C6FBB", cex = 1.4)
legend("bottomright", bty = "n", cex = 0.7, pch = c(16, NA), lwd = c(NA, 6),
       col = c("#2C6FBB", "grey60"),
       legend = c("observed Q2", "permutation 5-95%"))

## Component trajectories with replicate-bootstrap CIs, per treatment.
for (k in seq_len(min(2, K))) {
  for (blk in c("rna", "prot")) {
    sc <- traj[traj$component == k, ]
    yl <- range(sc[, paste0(blk, c("_lo", "_hi"))])
    plot(NA, xlim = range(tps), ylim = yl, xlab = "time (h)",
         ylab = sprintf("component %d score", k),
         main = sprintf("C%d %s trajectory", k,
                        ifelse(blk == "rna", "RNA", "protein")))
    for (i in seq_along(trts)) {
      s <- sc[sc$treatment == trts[i], ]
      s <- s[order(s$time), ]
      cl <- c("#2C6FBB", "#D1495B")[i]
      polygon(c(s$time, rev(s$time)),
              c(s[[paste0(blk, "_lo")]], rev(s[[paste0(blk, "_hi")]])),
              col = adjustcolor(cl, alpha.f = 0.18), border = NA)
      lines(s$time, s[[paste0(blk, "_score")]], type = "b", pch = 16, lwd = 2, col = cl)
    }
    abline(h = 0, col = "grey85")
    legend("topleft", bty = "n", cex = 0.75, lwd = 2, col = c("#2C6FBB", "#D1495B"),
           legend = trts)
  }
}
par(op); dev.off()
log_step("06_integration_caseB done -> ", outdir)
