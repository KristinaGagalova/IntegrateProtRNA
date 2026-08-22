#!/usr/bin/env Rscript
## =============================================================================
## 00b_check_simulation.R -- verify the simulator produces the biology it claims
##
## Archetype amplitudes are sign-randomised, so raw means cancel. Everything
## here is SIGN-ALIGNED to the direction of the RNA response (or, for
## protein_only genes, to the direction of the protein response).
##
## Run this after any change to 00_simulate_data.R. If these tables stop
## looking like the archetype definitions, the downstream validation is
## measuring the wrong thing.
## =============================================================================

source("R/utils.R")
cfg <- load_config()
tps <- cfg$design$timepoints
tr  <- read.delim(file.path(cfg$paths$data_sim, "truth.tsv"),
                  stringsAsFactors = FALSE)

rcols <- paste0("true_lfc_rna_t",  tps)
pcols <- paste0("true_lfc_prot_t", tps)

## Align each gene to the sign of its own largest-magnitude response, so that
## "up" and "down" instances of the same archetype reinforce instead of cancel.
ref <- ifelse(tr$archetype == "protein_only",
              tr[[pcols[length(pcols)]]], tr[[rcols[length(rcols)]]])
s   <- ifelse(ref >= 0, 1, -1)
s[tr$archetype == "unchanged"] <- 1

R <- as.matrix(tr[, rcols]) * s
P <- as.matrix(tr[, pcols]) * s
colnames(R) <- colnames(P) <- paste0("t", tps)

fmt <- function(m) round(aggregate(m, list(archetype = tr$archetype), mean)[, -1], 2)
tab <- data.frame(archetype = sort(unique(tr$archetype)))
cat("\n=== sign-aligned mean RNA log2FC (Trt vs Ctrl) ===\n")
print(cbind(tab, fmt(R)))
cat("\n=== sign-aligned mean PROTEIN log2FC (Trt vs Ctrl) ===\n")
print(cbind(tab, fmt(P)))

## Attenuation: protein amplitude relative to RNA amplitude, at the last
## timepoint. Should be ~1 for concordant, small-positive for buffered,
## ~0 for mrna_only, negative for anticorrelated.
last <- length(tps)
ratio <- P[, last] / ifelse(abs(R[, last]) < 0.25, NA, R[, last])
cat("\n=== median protein/RNA amplitude ratio at t", tps[last], " ===\n", sep = "")
print(round(tapply(ratio, tr$archetype, median, na.rm = TRUE), 2))

## Empirical delay: shift (h) maximising correlation between the interpolated
## RNA and protein trajectories. lag_only genes must dominate here.
delay_of <- function(r, p) {
  if (max(abs(r)) < 0.3) return(NA_real_)
  g  <- seq(0, max(tps), by = 1)
  ri <- approx(tps, r, g)$y; pi_ <- approx(tps, p, g)$y
  d  <- 0:round(max(tps) * 0.75)
  sc <- vapply(d, function(k) {
    n <- length(g) - k
    if (n < 4) return(-Inf)
    suppressWarnings(cor(ri[1:n], pi_[(1 + k):length(g)]))
  }, numeric(1))
  d[which.max(sc)]
}
set.seed(1)
idx <- sample(nrow(tr), min(1200, nrow(tr)))
dl  <- vapply(idx, function(i) delay_of(R[i, ], P[i, ]), numeric(1))
cat("\n=== median empirical RNA->protein delay (h), responders only ===\n")
print(round(tapply(dl, tr$archetype[idx], median, na.rm = TRUE), 1))

cat("\n=== median protein half-life (h) ===\n")
print(round(tapply(tr$t_half_h, tr$archetype, median), 1))

## Hard assertions: the simulator is broken if any of these fail.
mu <- function(a, m, tt) mean(m[tr$archetype == a, tt])
stopifnot(
  abs(mu("unchanged", R, last)) < 0.05,
  abs(mu("unchanged", P, last)) < 0.05,
  mu("concordant",     R, last) > 1.0,
  mu("concordant",     P, last) > 0.5,
  mu("mrna_only",      R, last) > 1.0,
  abs(mu("mrna_only",  P, last)) < 0.55,
  abs(mu("protein_only", R, last)) < 0.05,
  mu("protein_only",   P, last) > 0.3,
  mu("anticorrelated", R, last) > 1.0,
  mu("anticorrelated", P, last) < -0.5,
  mu("buffered",       R, last) > 1.0,
  mu("buffered",       P, last) > 0.30,
  mu("buffered",       P, last) < 0.75 * mu("concordant", P, last)
)
cat("\nAll simulator assertions passed.\n")
