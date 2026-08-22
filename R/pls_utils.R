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
