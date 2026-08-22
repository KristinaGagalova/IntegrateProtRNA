#!/usr/bin/env Rscript
## =============================================================================
## 07_bayesian_kinetics.R -- hierarchical Bayesian upgrade of the kinetic null
##                            model in 05_concordance_archetypes.R
##
## WHY THIS EXISTS
## ----------------
## 05's kinetic null / amplitude-LRT fit is a per-gene profile-likelihood
## point estimate: each gene identifies its own half-life and amplitude from 4
## timepoints and n=3 replicates alone, and calls a gene "regulated" with an
## arbitrary FDR cutoff on a per-gene LRT. That is why a large fraction of
## genes end up in the "kinetics_limited" bucket -- individually, a gene with
## a long half-life and a gene with genuine attenuation can produce nearly
## identical trajectories, and 4 points cannot always tell them apart.
##
## The fix that does not require more data per gene: LET GENES BORROW STRENGTH
## FROM EACH OTHER on the AMPLITUDE parameter specifically. Most genes' protein
## response should track their mRNA response one-for-one (a=1); a minority are
## genuinely attenuated/amplified/reversed. That is exactly the assumption a
## horseshoe prior encodes, and it is the model proposed as "4A" in
## docs/planning.md:
##
##     delta_g = a_g - 1      ~ Normal(0, tau0^2 * tau^2 * lambda_g^2)  (local
##                                                    horseshoe scale mixture)
##     lambda_g, tau          ~ half-Cauchy                            (heavy
##                                                    tails: most delta_g get
##                                                    shrunk to ~0, a few are
##                                                    let through unshrunk)
##
## Likelihood, identical in form to 05's weighted SSR (moderated SEs from
## limma are treated as known/fixed -- a pragmatic simplification, since they
## are already pooled, well-calibrated estimates from thousands of genes):
##
##     logFC_prot_g(t) ~ Normal( a_g * log2 p(t; kd_g),  se_prot_g(t)^2 )
##
## where p(t; kd) solves dp/dt = kd*(r(t) - p), p(0) = 1, via the SAME
## integrate_protein_mat() used by the simulator and by 05.
##
## WHY kd IS FIXED, NOT JOINTLY SAMPLED (an earlier version of this script
## tried that, and the self-test below is what caught the problem)
## -----------------------------------------------------------------
## A long half-life and a strongly attenuated amplitude produce nearly the
## same predicted trajectory -- that confounding IS the "kinetics_limited"
## problem this script exists to address. Jointly sampling log(half-life) and
## delta with independent single-site Metropolis proposals turns that
## confounding into a ridge-shaped posterior that component-wise MCMC mixes
## along very slowly: proposing a change to one parameter without a matched
## compensating change to the other steps off the ridge and is rejected. In
## a controlled test (40 simulated `buffered` genes, 6000 iterations) the
## joint sampler's delta visibly drifted toward 0 over the run instead of
## stabilising near the true attenuation, while 05's per-gene grid search
## (which profiles the amplitude out in closed form at every half-life on
## the grid, i.e. searches the *joint* optimum directly rather than by MCMC)
## recovered the right answer without difficulty. So: each gene's half-life
## is fixed here at ITS OWN grid-search joint-optimum (the same value
## 05 reports as `t_half_M1_h`), and only the amplitude is given a hierarchical
## prior. This sacrifices joint uncertainty propagation between the two
## parameters (a real limitation, stated once here rather than re-argued
## throughout) in exchange for a sampler that is simple enough to verify is
## doing the right thing given the time available for this analysis. A fully
## joint fit would need a sampler that proposes correlated moves along the
## ridge (e.g. NUTS, which is exactly why Stan/NumPyro are the standard tool
## for this class of model -- see docs/planning.md 4A) -- worth revisiting if
## joint kinetic/amplitude uncertainty turns out to matter for a specific
## downstream claim.
##
## INFERENCE: Metropolis-within-Gibbs, vectorised across genes
## -------------------------------------------------------------
## Conditional on the hyperparameters (tau, lambda), genes are independent in
## both prior and likelihood. That means the per-gene Metropolis step for
## delta_g can be PROPOSED AND ACCEPTED SIMULTANEOUSLY for every gene in one
## vectorised pass (integrate_protein_mat() is already vectorised over genes).
## The horseshoe hyperparameters get exact conjugate Gibbs updates (scale-
## mixture representation: Makalic & Schmidt 2016).
##
## No Stan/JAGS/C++ toolchain is required -- this runs with base R + the
## Bioconductor packages already in the pipeline, so it is renv-lockable
## without exotic build dependencies.
## =============================================================================

## here::here() locates the project root by walking up from the current
## working directory looking for a marker file (IntegrateProtRNA.Rproj,
## DESCRIPTION, renv.lock, .git -- this project has all four), so these two
## source() calls resolve correctly regardless of HOW this file is loaded:
## `Rscript R/07_bayesian_kinetics.R` from the project root (CWD already
## there), source()d from another R/*.R script (also CWD = project root, by
## convention), or source()d from inside a knitr/quarto chunk in
## analysis/*/*.qmd (CWD = that notebook's own directory). A working-
## directory-relative "R/utils.R" only covers the first two cases; a
## sys.frame(1)$ofile trick covers none of them reliably inside a knitr
## chunk, since there is no ofile in that evaluation frame.
source(here::here("R", "utils.R"))
source(here::here("R", "pls_utils.R"))

#' Fit the hierarchical horseshoe amplitude model at a fixed per-gene kd.
#'
#' @param lfc_rna,lfc_prot  genes x length(tps) matrices, CENTRED ON t0 (i.e.
#'   column 1 is all zero), same convention as 05_concordance_archetypes.R.
#' @param se_prot           genes x length(tps) moderated SE matrix (protein
#'   layer only -- the RNA trajectory is used deterministically, as in 05).
#' @param tps               numeric vector of sampled timepoints.
#' @param t_half            optional per-gene half-life vector (hours), fixed
#'   during MCMC. If NULL, computed internally via the same joint
#'   grid-search-with-profiled-amplitude procedure 05 uses for `t_half_M1_h`
#'   (recommended: pass 05's own `t_half_M1_h` column directly when available,
#'   so the two analyses are anchored to literally the same kinetic fit).
#' @param p0_frac  EXPECTED fraction of genes under genuine post-transcriptional
#'   regulation. The vanilla horseshoe's default half-Cauchy(0,1) scale on tau
#'   assumes high sparsity (a handful of non-null coefficients among many); if
#'   the true non-null fraction is well above ~10%, that default over-shrinks
#'   real signal along with the noise. Piironen & Vehtari (2017) recommend
#'   scaling tau's prior to tau0 = p0/(n - p0), p0 = p0_frac * n_genes --
#'   implemented below by running the standard scale-mixture Gibbs sampler on
#'   delta/tau0 rather than on delta directly. Set from domain expectation, or
#'   from the fraction 05_concordance_archetypes.R already called non-null.
#' @return list(summary = data.frame per gene, hyper = data.frame of
#'   hyperparameter posterior draws, accept = acceptance-rate diagnostic)
#'
#' CALLER RESPONSIBILITY: pre-filter to RNA-RESPONSIVE genes (the same
#' `rna_resp` gate 05_concordance_archetypes.R applies before classification:
#' RNA any-effect FDR below threshold AND max |logFC| above a floor). The
#' amplitude parameter `a` is only meaningful when the RNA trajectory carries
#' real signal -- for a gene with a flat RNA response, `pl` (the kinetic
#' prediction) is near zero at every timepoint, so (1+delta)*pl stays near
#' zero for ANY delta, and the likelihood has no information to constrain it.
#' Found in dev testing on real data: without this filter, a small fraction of
#' flat-RNA/large-protein-change genes drove delta to the numerical clip bound
#' even with the horseshoe prior in place, because the (mis-specified, for
#' those genes) likelihood still nominally preferred an extreme value. Genes
#' excluded this way are not lost -- they are the `protein_only` archetype and
#' belong in a separate, simpler test (protein FDR low, RNA FDR high), not
#' this model.
fit_hierarchical_kinetics <- function(lfc_rna, se_rna, lfc_prot, se_prot, tps,
                                       t_half = NULL, init_delta = NULL,
                                       p0_frac = 0.2,
                                       n_iter = 2000, warmup = 800,
                                       step_delta = 0.25, delta_bound = 10,
                                       seed = 1, verbose = TRUE) {
  set.seed(seed)
  g <- rownames(lfc_prot)
  n <- length(g)
  stopifnot(identical(dim(lfc_rna), dim(lfc_prot)), ncol(lfc_prot) == length(tps))

  grid <- seq(0, max(tps), by = 0.25)
  gi   <- match(tps, grid); stopifnot(!anyNA(gi))
  r_grid <- t(apply(lfc_rna, 1, function(v) approx(tps, v, grid, rule = 2)$y))
  r_lin  <- 2^r_grid
  W <- 1 / pmax(se_prot, 0.05)^2
  W[, 1] <- 0                                   # t0 is the centring anchor, not data

  ## ---- fixed kinetics: joint grid search with profiled-out amplitude -------
  ## Same procedure as 05's M1 fit: for each candidate half-life, the WLS-
  ## optimal amplitude has a closed form, so scanning the grid finds the joint
  ## (kd, a) optimum directly -- no MCMC needed for this part, and it is what
  ## anchors the horseshoe model at the right kinetics (see header note above).
  if (is.null(t_half) || is.null(init_delta)) {
    if (verbose) log_step("07: joint (kd, a) grid search over ", n, " genes")
    th_grid <- exp(seq(log(0.5), log(6 * max(tps)), length.out = 60))
    kd_grid <- log(2) / th_grid
    best_ssr <- rep(Inf, n); best_j <- rep(1L, n); best_a <- rep(1, n)
    for (j in seq_along(kd_grid)) {
      pl  <- log2(pmax(integrate_protein_mat(grid, r_lin, ks = kd_grid[j], kd = kd_grid[j],
                                             P0 = rep(1, n))[, gi, drop = FALSE], 1e-9))
      num <- rowSums(W * pl * lfc_prot); den <- rowSums(W * pl^2)
      a   <- ifelse(den > 1e-8, num / den, 1)
      ssr <- rowSums(W * (a * pl - lfc_prot)^2)
      better <- ssr < best_ssr
      best_ssr[better] <- ssr[better]; best_j[better] <- j; best_a[better] <- a[better]
    }
    ## Amplitude is not an identified quantity for a gene with no RNA signal
    ## (a * 0 = 0 for any a: den = sum(W*pl^2) ~= 0 makes the closed-form `a`
    ## numerically unstable). Start those genes at the neutral a=1 so the
    ## horseshoe prior -- correctly -- dominates their posterior instead of an
    ## MCMC chain anchored to a numerical artefact.
    rna_amplitude <- rowSums(W * lfc_rna^2)
    weak_rna <- rna_amplitude < quantile(rna_amplitude, 0.15)
    best_a[weak_rna] <- 1
    t_half      <- if (is.null(t_half))      th_grid[best_j]              else t_half
    init_delta  <- if (is.null(init_delta))
      pmax(pmin(best_a - 1, delta_bound), -delta_bound) else init_delta
  }
  kd_fixed <- log(2) / t_half

  ## Kinetic prediction at the FIXED kd, computed once -- this is the entire
  ## reason fixing kd removes the ridge: `pl` below never changes during MCMC.
  pl <- log2(pmax(integrate_protein_mat(grid, r_lin, ks = kd_fixed, kd = kd_fixed,
                                        P0 = rep(1, n))[, gi, drop = FALSE], 1e-9))

  ## ---- state ----------------------------------------------------------------
  ## tau0 rescales the horseshoe's implicit sparsity assumption to p0_frac
  ## (see the p0_frac doc above). The Gibbs conjugate updates for lambda2/tau2
  ## below are the textbook scheme applied to the RESCALED coefficient
  ## delta/tau0; the MH likelihood step then uses prior variance
  ## tau0^2 * tau2 * lambda2 for the unscaled delta.
  p0   <- max(1, min(n - 1, p0_frac * n))
  tau0 <- p0 / (n - p0)
  delta   <- init_delta
  lambda2 <- rep(1, n)
  tau2    <- 1
  cur_ssr <- rowSums(W * ((1 + delta) * pl - lfc_prot)^2)

  n_keep <- n_iter - warmup
  keep_delta <- matrix(NA_real_, n, n_keep)
  hyper <- data.frame(iter = seq_len(n_iter), tau = NA_real_, accept_delta = NA_real_)

  if (verbose) log_step("07: MCMC (fixed kd, horseshoe amplitude), ",
                        n_iter, " iterations (", warmup, " warmup), ", n, " genes")
  for (it in seq_len(n_iter)) {
    ## -- amplitude step, vectorised across all genes simultaneously ---------
    ## delta_bound is a pure numerical safety valve, not a scientific prior:
    ## real data (unlike the simulator) occasionally has a gene with a near-
    ## degenerate weighted likelihood (e.g. a near-zero moderated SE from a
    ## heavily-imputed low-abundance protein), and without a bound such a gene
    ## can random-walk to an implausible amplitude (|a| in the tens) before its
    ## own horseshoe lambda_g has inflated enough to pull it back. The bound is
    ## set generously (delta in [-10, 10], i.e. a in [-9, 11]) so it essentially
    ## never engages for a well-behaved gene.
    prop_delta <- pmax(pmin(delta + rnorm(n, 0, step_delta), delta_bound), -delta_bound)
    prop_ssr   <- rowSums(W * ((1 + prop_delta) * pl - lfc_prot)^2)
    prior_var  <- pmax(tau0^2 * tau2 * lambda2, 1e-8)
    log_prior_cur  <- -0.5 * delta^2      / prior_var
    log_prior_prop <- -0.5 * prop_delta^2 / prior_var
    log_ratio <- -0.5 * (prop_ssr - cur_ssr) + (log_prior_prop - log_prior_cur)
    acc <- log(runif(n)) < log_ratio
    delta[acc] <- prop_delta[acc]; cur_ssr[acc] <- prop_ssr[acc]

    ## -- horseshoe scale mixture, exact Gibbs conjugate updates -------------
    ## (Makalic & Schmidt 2016), applied to delta_resc = delta / tau0 -- a
    ## standard-scale-mixture coefficient under a half-Cauchy(0,1) prior on
    ## tau, which makes these updates exactly conjugate regardless of p0_frac.
    delta_resc <- delta / tau0
    nu_g    <- 1 / rgamma(n, 1, rate = 1 + 1 / lambda2)
    lambda2 <- 1 / rgamma(n, 1, rate = 1 / nu_g + delta_resc^2 / (2 * tau2))
    xi      <- 1 / rgamma(1, 1, rate = 1 + 1 / tau2)
    tau2    <- 1 / rgamma(1, (n + 1) / 2, rate = 1 / xi + 0.5 * sum(delta_resc^2 / lambda2))

    hyper$tau[it] <- sqrt(tau2); hyper$accept_delta[it] <- mean(acc)
    if (it > warmup) keep_delta[, it - warmup] <- delta
    if (verbose && it %% max(1, n_iter %/% 10) == 0)
      log_step("  iter ", it, "/", n_iter, "  accept=", round(mean(acc), 2),
               " tau_eff=", round(tau0 * sqrt(tau2), 3))
  }

  qtl <- function(m, p) apply(m, 1, quantile, probs = p, na.rm = TRUE)
  summary_df <- data.frame(
    gene_id = g, t_half_h = t_half,
    a_mean = 1 + rowMeans(keep_delta),
    a_lo   = 1 + qtl(keep_delta, .025),
    a_hi   = 1 + qtl(keep_delta, .975),
    ## a horseshoe-native "is this gene regulated" call: the 95% posterior
    ## interval for a excludes 1 (the RNA-only-kinetics null) -- replaces
    ## 05's per-gene LRT + arbitrary FDR cutoff with a hierarchically shrunk
    ## interval that borrows strength across all genes fit together.
    regulated_95 = (1 + qtl(keep_delta, .025) > 1) | (1 + qtl(keep_delta, .975) < 1),
    ## Flags genes whose posterior sits at the delta_bound numerical safety
    ## valve (see the argument doc) -- their a_mean/a_lo/a_hi should be read
    ## as "at least this extreme", not as a precise point estimate. Expect a
    ## small fraction (a couple of percent) even after RNA-responsive gating:
    ## real data legitimately contains a few genes with very strong apparent
    ## amplification/attenuation.
    at_bound = abs(rowMeans(keep_delta)) > 0.9 * delta_bound,
    stringsAsFactors = FALSE)

  list(summary = summary_df, hyper = hyper, trace_delta = keep_delta,
       accept_delta = mean(hyper$accept_delta[-seq_len(warmup)]))
}

## ---------------------------------------------------------------------------
## Standalone self-test: fit on the simulator's ground truth and check the
## shrinkage recovers the same buffered-vs-concordant separation as 05, with
## honest uncertainty intervals. Run directly:  Rscript R/07_bayesian_kinetics.R
## ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  cfg <- load_config()
  dedir <- file.path(cfg$paths$results, "de")
  Rl <- read_mat(file.path(dedir, "rna_lfc.tsv"));  Rs <- read_mat(file.path(dedir, "rna_se.tsv"))
  Pl <- read_mat(file.path(dedir, "prot_lfc.tsv")); Ps <- read_mat(file.path(dedir, "prot_se.tsv"))
  arch <- read.delim(file.path(cfg$paths$results, "concordance", "archetypes.tsv"),
                     stringsAsFactors = FALSE)
  rownames(arch) <- arch$gene_id
  g  <- intersect(rownames(Rl), rownames(Pl))
  set.seed(1); g <- sample(g, min(800, length(g)))
  Rl <- Rl[g, ] - Rl[g, 1]; Pl <- Pl[g, ] - Pl[g, 1]
  Rs <- Rs[g, ]; Ps <- Ps[g, ]
  tps <- cfg$design$timepoints

  ## Anchor at 05's own joint-optimal half-life for these genes (t_half_M1_h)
  ## rather than recomputing it -- keeps the two analyses on the same kinetics.
  fit <- fit_hierarchical_kinetics(Rl, Rs, Pl, Ps, tps,
                                   t_half = arch[g, "t_half_M1_h"], n_iter = 2500, warmup = 800, seed = 1)
  cat("\nacceptance rate (delta):", round(fit$accept_delta, 2), "\n")

  truth <- read.delim(file.path(cfg$paths$data_sim, "truth.tsv"), stringsAsFactors = FALSE)
  rownames(truth) <- truth$gene_id
  tt <- truth[g, ]
  cat("\nmedian posterior a (this script) by TRUE archetype:\n")
  print(round(tapply(fit$summary$a_mean, tt$archetype, median), 2))
  cat("\nmedian amplitude_a (05, profile-likelihood) by TRUE archetype, same genes:\n")
  print(round(tapply(arch[g, "amplitude_a"], tt$archetype, median), 2))
  cat("\nfraction posterior-regulated (95% CI excludes a=1) by TRUE archetype:\n")
  print(round(tapply(fit$summary$regulated_95, tt$archetype, mean), 2))
  cat("\nfraction 05 called regulated (LRT FDR<0.05) by TRUE archetype:\n")
  print(round(tapply(arch[g, "regulated"], tt$archetype, mean), 2))
}
