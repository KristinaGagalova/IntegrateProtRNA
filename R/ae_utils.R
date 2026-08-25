## ae_utils.R -- condition-aligned autoencoder, extracted from
## dimensionality_reduction_wheat.qmd Sec.3 so a second consumer
## (analysis/enrichment/enrichment_wheat.qmd) does not duplicate it.
## adam_init/adam_update/cond_ae_fit/pred_cond_ae/fold_scaler/
## permute_Y_by_cell are verbatim copies of that notebook's own
## implementation -- same audited code, not a reimplementation.

adam_init <- function(par) {
  list(m = lapply(par, function(p) array(0, dim(p))),
       v = lapply(par, function(p) array(0, dim(p))),
       t = 0)
}

adam_update <- function(par, grad, state, lr = 0.01,
                        b1 = 0.9, b2 = 0.999, eps = 1e-8) {
  state$t <- state$t + 1
  for (nm in names(par)) {
    state$m[[nm]] <- b1 * state$m[[nm]] + (1 - b1) * grad[[nm]]
    state$v[[nm]] <- b2 * state$v[[nm]] + (1 - b2) * grad[[nm]]^2
    mhat <- state$m[[nm]] / (1 - b1^state$t)
    vhat <- state$v[[nm]] / (1 - b2^state$t)
    par[[nm]] <- par[[nm]] - lr * mhat / (sqrt(vhat) + eps)
  }
  list(par = par, state = state)
}

#' Condition-aligned autoencoder (linear), after `rnaprot/unpaired.py`.
#' Two encoders, two decoders, tied by CONDITION CENTROIDS -- never by sample
#' pairing. See dimensionality_reduction_wheat.qmd Sec.3 (fn-condition-ae) for
#' the full derivation and the four loss terms' rationale.
cond_ae_fit <- function(R, P, cond_r, cond_p, k = 3,
                        w = c(0.5, 0.5, 2.0, 0.5),
                        epochs = 3000, lr = 0.01, seed = 1) {

  set.seed(seed)
  nr <- nrow(R); pr <- ncol(R)
  np <- nrow(P); pp <- ncol(P)

  init <- function(a, b) matrix(rnorm(a * b, 0, 1 / sqrt(a)), a, b)
  par <- list(Wr = init(pr, k), Ar = init(k, pr),
              Wp = init(pp, k), Ap = init(k, pp))

  shared <- intersect(unique(cond_r), unique(cond_p))
  idx_r <- lapply(shared, function(cc) which(cond_r == cc))
  idx_p <- lapply(shared, function(cc) which(cond_p == cc))
  nc <- length(shared)

  st <- adam_init(par); hist <- numeric(epochs)

  for (e in seq_len(epochs)) {
    z_r <- R %*% par$Wr; z_p <- P %*% par$Wp
    Rh  <- z_r %*% par$Ar; Ph <- z_p %*% par$Ap

    Er <- Rh - R; Ep <- Ph - P
    L1 <- sum(Er^2) / (nr * pr); L2 <- sum(Ep^2) / (np * pp)

    dRh <- w[1] * 2 * Er / (nr * pr)
    dPh <- w[2] * 2 * Ep / (np * pp)

    g <- list(Ar = t(z_r) %*% dRh, Ap = t(z_p) %*% dPh)
    dz_r <- dRh %*% t(par$Ar)
    dz_p <- dPh %*% t(par$Ap)

    L3 <- 0; L4 <- 0
    for (i in seq_len(nc)) {
      ri <- idx_r[[i]]; pi <- idx_p[[i]]

      mzr  <- colMeans(z_r[ri, , drop = FALSE])
      mzp  <- colMeans(z_p[pi, , drop = FALSE])
      pred <- matrix(mzr, 1) %*% par$Ap
      obs  <- colMeans(P[pi, , drop = FALSE])

      e3 <- pred - matrix(obs, 1); L3 <- L3 + sum(e3^2) / pp
      e4 <- matrix(mzr - mzp, 1);  L4 <- L4 + sum(e4^2) / k

      c3 <- w[3] * 2 / (nc * pp)
      g$Ap <- g$Ap + c3 * matrix(mzr, ncol = 1) %*% e3
      d_mzr <- c3 * e3 %*% t(par$Ap)

      c4 <- w[4] * 2 / (nc * k)
      d_mzr <- d_mzr + c4 * e4
      d_mzp <- -c4 * e4

      dz_r[ri, ] <- dz_r[ri, ] + matrix(rep(d_mzr / length(ri), each = length(ri)),
                                        length(ri))
      dz_p[pi, ] <- dz_p[pi, ] + matrix(rep(d_mzp / length(pi), each = length(pi)),
                                        length(pi))
    }

    g$Wr <- t(R) %*% dz_r
    g$Wp <- t(P) %*% dz_p

    hist[e] <- w[1]*L1 + w[2]*L2 + w[3]*L3/nc + w[4]*L4/nc
    up <- adam_update(par, g, st, lr = lr); par <- up$par; st <- up$state
  }

  list(par = par, k = k, loss = hist,
       predict_protein = function(Rn) (Rn %*% par$Wr) %*% par$Ap,
       encode_r = function(Rn) Rn %*% par$Wr)
}

pred_cond_ae <- function(Xtr, Ytr, Xte, cond_tr, k = 3, epochs = 600, seed = 1, ...) {
  m <- cond_ae_fit(Xtr, Ytr, cond_tr, cond_tr, k = k, epochs = epochs, seed = seed)
  m$predict_protein(Xte)
}

#' Fold-safe centre+scale: fit mean + RMS scale from TRAINING rows only.
#' See dimensionality_reduction_wheat.qmd Sec.3/Sec.7 for the CV-centering
#' leak this specifically guards against for a homogeneous, no-intercept
#' linear predictor like cond_ae_fit's `predict_protein`.
fold_scaler <- function(train_mat) {
  mu <- colMeans(train_mat)
  s  <- sqrt(sum(scale(train_mat, center = TRUE, scale = FALSE)^2) / nrow(train_mat))
  function(M) sweep(M, 2, mu, "-") / s
}

#' Reshuffle which cell's protein data is matched to which cell's RNA data.
permute_Y_by_cell <- function(Y, cell_of) {
  cells <- unique(cell_of)
  perm  <- setNames(sample(cells), cells)
  Yp <- Y
  for (cl in cells) {
    src <- which(cell_of == perm[[cl]])
    dst <- which(cell_of == cl)
    Yp[dst, ] <- Y[src[sample(length(src))], , drop = FALSE]
  }
  Yp
}

## ---------------------------------------------------------------------
## New in this project (not in dimensionality_reduction_wheat.qmd): a
## permutation-null Q2 for the condition-aligned AE at an arbitrary k,
## and a replicate-resampling bootstrap of its per-gene DECODER loading
## (Ap) -- the AE analogue of pls_boot_loadings() (R/pls_utils.R).
## ---------------------------------------------------------------------

#' LOCO Q2 + permutation null for the condition-aligned AE, condition-mean
#' level (each training row is its own singleton "cell" -- the same scheme
#' dimensionality_reduction_wheat.qmd's D4 uses for its own headline Q2/
#' p_perm numbers; D9 in that notebook notes this is a degenerate use of the
#' architecture's namesake centroid-alignment mechanism, which needs
#' replicate-level data to be exercised as designed -- not extended to a
#' permutation null here; flagged as a caveat in the enrichment writeup).
cond_ae_q2_perm <- function(X_raw, Y_raw, k = 1, seeds = 1:3, epochs = 3000,
                            nperm = 49, perm_epochs = 600, seed_base = 4000) {
  n <- nrow(X_raw)

  loco <- function(Yr, ep, sds) {
    press <- 0
    for (i in seq_len(n)) {
      tr <- setdiff(seq_len(n), i)
      sx <- fold_scaler(X_raw[tr, , drop = FALSE]); sy <- fold_scaler(Yr[tr, , drop = FALSE])
      Xtr <- sx(X_raw[tr, , drop = FALSE]); Xte <- sx(X_raw[i, , drop = FALSE])
      Ytr <- sy(Yr[tr, , drop = FALSE]);   Yte <- sy(Yr[i, , drop = FALSE])
      pr <- rowMeans(vapply(sds, function(s)
        as.vector(pred_cond_ae(Xtr, Ytr, Xte, rownames(Xtr), k = k, epochs = ep, seed = s)),
        numeric(ncol(Ytr))))
      press <- press + sum((as.vector(Yte) - pr)^2)
    }
    tss <- sum(vapply(seq_len(n), function(i) {
      tr <- setdiff(seq_len(n), i)
      sy <- fold_scaler(Yr[tr, , drop = FALSE])
      Yte <- sy(Yr[i, , drop = FALSE]); Ytr <- sy(Yr[tr, , drop = FALSE])
      sum(sweep(Yte, 2, colMeans(Ytr))^2)
    }, numeric(1)))
    1 - press / tss
  }

  q2_obs <- loco(Y_raw, epochs, seeds)
  null <- vapply(seq_len(nperm), function(b) {
    set.seed(seed_base + b)
    Yp <- permute_Y_by_cell(Y_raw, rownames(Y_raw))
    loco(Yp, perm_epochs, 1)
  }, numeric(1))
  p <- (1 + sum(null >= q2_obs)) / (nperm + 1)
  data.frame(Q2 = round(q2_obs, 4),
            null_p95 = round(unname(quantile(null, .95)), 4),
            p_perm = signif(p, 3))
}

#' Replicate-resampling bootstrap of the condition-aligned AE's protein
#' decoder loading (Ap), condition-mean level. k = 1 is deliberate, not a
#' default: at k >= 2 the latent dimensions carry no ordering or identity
#' constraint (unlike PLS's sequentially-deflated components), so which
#' resampled dimension corresponds to which point-estimate dimension is not
#' well defined. k = 1 sidesteps that alignment problem entirely, leaving
#' only a sign ambiguity, which align_sign() (R/pls_utils.R) already handles.
#'
#' Each bootstrap draw gets its OWN fold_scaler fit from its own resampled
#' cell means (never the point estimate's scaler) -- the same convention
#' pls_boot_loadings() and 06_integration_caseB.R's own bootstrap use.
cond_ae_boot_loadings <- function(rna_mat, prot_mat, cd_r, cd_p, cells, hv_r, hv_p,
                                  k = 1, B = 100, epochs = 1000, seed = 1) {
  set.seed(seed)
  idx_r <- lapply(cells, function(cl) which(cd_r$cell == cl))
  idx_p <- lapply(cells, function(cl) which(cd_p$cell == cl))

  Xc <- t(cell_means(rna_mat[hv_r, , drop = FALSE],  cd_r))
  Yc <- t(cell_means(prot_mat[hv_p, , drop = FALSE], cd_p))
  sx <- fold_scaler(Xc); sy <- fold_scaler(Yc)
  Xcs <- sx(Xc); Ycs <- sy(Yc)

  fit <- cond_ae_fit(Xcs, Ycs, rownames(Xcs), rownames(Ycs), k = k, epochs = epochs, seed = seed)
  point <- as.vector(fit$par$Ap)

  boot_ap <- matrix(NA_real_, length(hv_p), B)
  for (b in seq_len(B)) {
    ir <- unlist(lapply(idx_r, function(ii) sample(ii, length(ii), TRUE)))
    ip <- unlist(lapply(idx_p, function(ii) sample(ii, length(ii), TRUE)))
    Xb <- t(cell_means(rna_mat[hv_r, ir, drop = FALSE],  cd_r[ir, , drop = FALSE]))
    Yb <- t(cell_means(prot_mat[hv_p, ip, drop = FALSE], cd_p[ip, , drop = FALSE]))
    sxb <- fold_scaler(Xb); syb <- fold_scaler(Yb)
    fb <- cond_ae_fit(sxb(Xb), syb(Yb), rownames(Xb), rownames(Yb),
                      k = k, epochs = epochs, seed = seed + b)
    boot_ap[, b] <- as.vector(fb$par$Ap)
  }

  ## align_sign()'s `m * ref` requires ref to be a plain vector, not a
  ## single-column matrix -- R's `*` treats two dim'd objects of mismatched
  ## shape as non-conformable (it will not recycle a p x 1 matrix down a
  ## p x B one), whereas a plain length-p vector recycles correctly. `point`
  ## is already a plain vector (as.vector() above), so pass it as-is.
  boot_aligned <- align_sign(boot_ap, point)
  pt_sign <- sign(point); pt_sign[pt_sign == 0] <- 1
  stability <- rowMeans(sign(boot_aligned) == pt_sign, na.rm = TRUE)

  list(protein = data.frame(feature = hv_p, loading = point,
                            boot_mean = rowMeans(boot_aligned, na.rm = TRUE),
                            boot_sd = apply(boot_aligned, 1, sd, na.rm = TRUE),
                            stability = stability),
       fit = fit)
}
