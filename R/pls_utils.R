## =============================================================================
## pls_utils.R -- PLS2 / leave-one-cell-out Q2 helpers for Case B integration
##
## Extracted from 06_integration_caseB.R so the same audited implementation is
## shared between that script and the analysis notebooks (analysis/*/*.qmd)
## rather than duplicated. See 06_integration_caseB.R's header for why Q2
## against a permutation null -- not the in-sample cross-block correlation --
## is the valid statistic for Case B design-cell integration.
## =============================================================================

## Per-block scaling to unit TOTAL variance. Without it the block with more
## features and larger dynamic range (usually RNA) dominates every component.
prep_block <- function(m) {
  m <- t(m)                                   # cells x features
  m <- scale(m, center = TRUE, scale = FALSE)
  m / sqrt(sum(m^2) / nrow(m))
}

## Plain PLS2 (NIPALS). Implemented here rather than pulled from a package so
## the bootstrap and permutation loops stay fast and the algorithm is auditable.
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
## With K design cells both blocks have rank <= K, so X = U D V' with V
## orthonormal. PLS on X is equivalent to PLS on A = U D (K x K), and because V
## is orthonormal all the Euclidean distances that PRESS and TSS are built from
## are preserved exactly. So this is not an approximation -- Q2 computed on A
## equals Q2 computed on the full matrix, at a few thousandths of the cost.
## That is what makes a proper permutation null affordable at this scale.
row_space <- function(M) {
  s <- svd(M)
  k <- sum(s$d > max(dim(M)) * .Machine$double.eps * s$d[1])
  s$u[, seq_len(k), drop = FALSE] %*% diag(s$d[seq_len(k)], k, k)
}

## Leave-one-design-cell-out cross-validation.
##
## THIS is the valid test statistic for Case B, not cor(t, u). With few
## pseudo-samples and thousands of features a PLS component can be steered
## onto almost any target, so the in-sample cross-block correlation approaches
## 1 even for PERMUTED data. Q2 asks the one question that cannot be gamed:
## given the RNA profile of a design cell the model has never seen, can it
## predict that cell's protein profile? Centring is recomputed inside each
## fold, so nothing leaks from the held-out cell.
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

## Permutation-null p-value for each component's leave-one-cell-out Q2.
##
## Q2 alone cannot be trusted as "validated" at face value here: with as few
## as 8 design-cell pseudo-samples, q2_loo() can still look deceptively good
## by chance for a component with no real cross-block signal (the same reason
## 06_integration_caseB.R and dimensionality_reduction_wheat.qmd's D7 both
## build a permutation null for Q2 rather than reading it directly). This is
## that check, factored out so any PLS-derived per-gene statistic -- not just
## the ones in 06_integration_caseB.R -- can be gated on it before use.
q2_perm_pvalue <- function(Xc, Yc, ncomp, nperm = 300, seed = 1) {
  set.seed(seed)
  Ax <- row_space(Xc); Ay <- row_space(Yc)
  Kq <- min(ncomp, nrow(Ax) - 3)
  q2_obs <- q2_loo(Ax, Ay, Kq)
  perm_q2 <- matrix(NA_real_, nperm, Kq)
  for (b in seq_len(nperm)) {
    o <- sample(nrow(Ay))
    perm_q2[b, ] <- q2_loo(Ax, Ay[o, , drop = FALSE], Kq)
  }
  p <- vapply(seq_len(Kq), function(k)
    (1 + sum(perm_q2[, k] >= q2_obs[k])) / (nperm + 1), 0)
  data.frame(component = seq_len(Kq), Q2 = round(q2_obs, 4),
            Q2_perm_q95 = round(apply(perm_q2, 2, quantile, .95), 4),
            p_Q2 = signif(p, 3))
}

## Replicate-resampling bootstrap of PLS LOADINGS (not just scores).
##
## 06_integration_caseB.R's replicate bootstrap (its step 2) already resamples
## replicates within each design cell, rebuilds cell means, and refits -- but
## it only tracks the bootstrap distribution of SCORES (t, u), for component
## trajectory CIs. It does not track LOADINGS (w, c), so no per-gene stability
## statistic has existed anywhere in this project until now.
##
## This matters because dimensionality_reduction_wheat.qmd's D7 established,
## on the real wheat data, which resampling scheme is trustworthy: condition-
## mean aggregation with THIS replicate-within-cell resampling (Q2 = 0.309,
## p = 0.015 for Cadenza) versus a FIXED individual-replicate pairing (an
## unrepeatable, apparently cherry-picked draw -- see enrichment_wheat.qmd
## E4 discussion) versus a fully randomised cross-block replicate pairing
## (decisively negative, does not beat its own null). Bootstrapping loadings
## under the FIRST scheme is therefore the one way to get a per-gene PLS
## statistic that inherits a validated resampling procedure rather than an
## invalidated one.
##
## Returns, for each feature in each block, the bootstrap mean loading and a
## STABILITY score = the fraction of bootstrap draws whose loading sign
## agrees with the point estimate's sign (1 = always agrees, 0.5 = coin flip).
## A gene with a large point-estimate loading but low stability is exactly
## the failure mode this guards against: at n = 8 pseudo-samples one gene's
## loading can be large by chance, and only resampling reveals that.
pls_boot_loadings <- function(rna_mat, prot_mat, cd_r, cd_p, cells, hv_r, hv_p,
                              ncomp, B = 200, seed = 1) {
  set.seed(seed)
  idx_r <- lapply(cells, function(cl) which(cd_r$cell == cl))
  idx_p <- lapply(cells, function(cl) which(cd_p$cell == cl))

  Xc <- prep_block(cell_means(rna_mat[hv_r, , drop = FALSE],  cd_r))
  Yc <- prep_block(cell_means(prot_mat[hv_p, , drop = FALSE], cd_p))
  stopifnot(identical(rownames(Xc), rownames(Yc)))
  fit <- pls2(Xc, Yc, ncomp = ncomp)

  boot_w <- array(NA_real_, c(length(hv_r), ncomp, B))
  boot_c <- array(NA_real_, c(length(hv_p), ncomp, B))
  for (b in seq_len(B)) {
    ir <- unlist(lapply(idx_r, function(ii) sample(ii, length(ii), TRUE)))
    ip <- unlist(lapply(idx_p, function(ii) sample(ii, length(ii), TRUE)))
    Xb <- prep_block(cell_means(rna_mat[hv_r, ir, drop = FALSE],  cd_r[ir, ]))
    Yb <- prep_block(cell_means(prot_mat[hv_p, ip, drop = FALSE], cd_p[ip, ]))
    fb <- pls2(Xb, Yb, ncomp = ncomp)
    boot_w[, , b] <- align_sign(fb$w, fit$w)
    boot_c[, , b] <- align_sign(fb$c, fit$c)
  }

  summarise_block <- function(point, boot, feat_names) {
    do.call(rbind, lapply(seq_len(ncomp), function(k) {
      bk <- boot[, k, ]                      # features x B
      pt_sign <- sign(point[, k]); pt_sign[pt_sign == 0] <- 1
      stability <- rowMeans(sign(bk) == pt_sign, na.rm = TRUE)
      data.frame(feature = feat_names, component = k,
                loading = point[, k],
                boot_mean = rowMeans(bk, na.rm = TRUE),
                boot_sd = apply(bk, 1, sd, na.rm = TRUE),
                stability = stability)
    }))
  }

  list(
    rna     = summarise_block(fit$w, boot_w, hv_r),
    protein = summarise_block(fit$c, boot_c, hv_p),
    fit     = fit
  )
}
