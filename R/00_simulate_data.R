#!/usr/bin/env Rscript
## =============================================================================
## 00_simulate_data.R -- Case B ground-truth simulator
##
## Generates TWO INDEPENDENT sets of 24 samples (RNA and protein) from one
## shared latent biology, so the pipeline can be validated end to end before
## it ever touches real data.
##
## Generative model, per gene g and condition c:
##     log2FC_R(t) = A_time * s(t)  +  A_trt * s(t) * 1[c = Trt]
##     R_gc(t)     = base_g * 2^log2FC_R(t)
##     dP/dt       = ks_gc * R_gc(t) - kd_gc * P_gc(t)
##     ks_gc       = ks_g * exp(ds_g * 1[c = Trt])
##     kd_gc       = kd_g * exp(dd_g * 1[c = Trt])
##
## Archetypes are produced by MECHANISM, not by post-hoc labelling:
##   buffered      -> dd > 0        (treatment accelerates degradation)
##   mrna_only     -> ds = -A_trt   (translation cancels the mRNA change)
##   protein_only  -> A_trt = 0, dd != 0
##   lag_only      -> long protein half-life; discordance is PURE KINETICS
## The lag_only class is the important one: it is what naive quadrant
## classification wrongly calls post-transcriptional regulation.
##
## Usage:  Rscript R/00_simulate_data.R [config/config.yaml]
## =============================================================================

source("R/utils.R")
args <- commandArgs(trailingOnly = TRUE)
cfg  <- load_config(if (length(args)) args[1] else "config/config.yaml")
set.seed(cfg$project$seed)

S      <- cfg$simulate
tps    <- cfg$design$timepoints
G      <- S$n_genes
outdir <- ensure_dir(cfg$paths$data_sim)
log_step("simulating ", G, " genes; timepoints = ", paste(tps, collapse = ","))

## ------------------------------------------------------- 1. gene identity ---
gene_id   <- sprintf("ENSG%08d", seq_len(G))
log2_base <- pmin(pmax(rnorm(G, 7, 2.2), 1), 18)   # realistic dynamic range
base_expr <- 2^log2_base

## ------------------------------------------------- 2. archetype assignment --
props <- unlist(S$archetype_props)
stopifnot(abs(sum(props) - 1) < 1e-8)
archetype <- sample(names(props), G, replace = TRUE, prob = props)

## Protein half-life (h) controls how much discordance is pure kinetics.
t_half <- ifelse(archetype == "lag_only",
                 exp(rnorm(G, log(90), 0.35)),    # slow turnover -> long lag
                 exp(rnorm(G, log(12), 0.75)))    # median ~12 h
t_half <- pmin(pmax(t_half, 1.5), 400)
kd     <- log(2) / t_half
ks     <- kd                 # sets P0 == R0; the MS response factor absorbs
                             # the true absolute scale (see note in 08)

## Response shape: half monotone-saturating, half impulse.
shape_type <- sample(c("mono", "impulse"), G, TRUE, c(0.5, 0.5))
tau_on     <- exp(rnorm(G, log(5), 0.5))
tau_off    <- exp(rnorm(G, log(40), 0.5))

shape_fun <- function(t, type, ton, toff) {
  if (type == "mono") {
    1 - exp(-t / ton)
  } else {
    y <- (1 - exp(-t / ton)) * exp(-t / toff)
    y / max(max(y), 1e-9)
  }
}

## A_time: a shared time effect present in BOTH arms (media change, handling).
## Including it is what makes the treatment:time interaction test meaningful.
A_time <- ifelse(runif(G) < 0.15, rnorm(G, 0, 0.8), 0)

## Treatment-specific amplitude and rate-constant shifts, set per archetype.
A_trt <- ds <- dd <- numeric(G)
sgn   <- sample(c(-1, 1), G, TRUE)
amp   <- abs(rnorm(G, 1.6, 0.5)) + 0.6
is_   <- function(a) archetype == a

A_trt[is_("concordant")] <- (sgn * amp)[is_("concordant")]

## buffered: mRNA moves, protein amplitude attenuated to a fraction f
f_buf <- runif(G, 0.30, 0.55)
A_trt[is_("buffered")] <- amp[is_("buffered")]
dd[is_("buffered")]    <- (log(2) * amp * (1 - f_buf))[is_("buffered")]

## mrna_only: a translation-rate shift cancels the mRNA change.
## NOTE: ds is a CONSTANT rate shift while the mRNA change A_trt*s(t) is
## time-varying, so cancellation can only ever be approximate -- exact at the
## timepoints where s(t) equals its mean, residual elsewhere. This is
## biologically honest (a fixed change in translational efficiency cannot
## perfectly track a dynamic transcript) and it is why mrna_only genes retain
## a small residual protein signal rather than being perfectly flat.
## Scaling by the mean shape over the SAMPLED timepoints centres that residual.
A_trt[is_("mrna_only")] <- (sgn * amp)[is_("mrna_only")]

## protein_only: flat mRNA, protein driven purely by a degradation shift
delta_p <- sgn * (abs(rnorm(G, 1.3, 0.4)) + 0.5)
dd[is_("protein_only")] <- (-log(2) * delta_p)[is_("protein_only")]

## anticorrelated: mRNA up, protein down (over-compensating degradation)
A_trt[is_("anticorrelated")] <- amp[is_("anticorrelated")]
dd[is_("anticorrelated")] <-
  (log(2) * (amp + abs(rnorm(G, 1.0, 0.3))))[is_("anticorrelated")]

## lag_only: fully concordant mechanism, but slow protein turnover
A_trt[is_("lag_only")] <- (sgn * amp)[is_("lag_only")]

## --------------------------------------------- 3. latent trajectories -------
grid <- seq(0, max(tps), by = 0.25)
gi   <- match(tps, grid); stopifnot(!anyNA(gi))

shape_grid <- t(vapply(seq_len(G), function(g)
  shape_fun(grid, shape_type[g], tau_on[g], tau_off[g]), numeric(length(grid))))

## Deferred from section 2: the mrna_only translation shift needs the mean
## response shape over the sampled timepoints (see the note there).
shp_bar <- rowMeans(shape_grid[, gi, drop = FALSE])
ds[is_("mrna_only")] <- (-log(2) * sgn * amp * shp_bar)[is_("mrna_only")]

R_ctrl <- base_expr * 2^(A_time * shape_grid)
R_trt  <- base_expr * 2^((A_time + A_trt) * shape_grid)
rownames(R_ctrl) <- rownames(R_trt) <- gene_id

P_ctrl <- integrate_protein_mat(grid, R_ctrl, ks, kd)
## Both arms start from the same pre-treatment steady state at t = 0.
P_trt  <- integrate_protein_mat(grid, R_trt, ks * exp(ds), kd * exp(dd),
                                P0 = ks * R_trt[, 1] / kd)

true_lfc_rna  <- log2(R_trt[, gi, drop = FALSE] / R_ctrl[, gi, drop = FALSE])
true_lfc_prot <- log2(P_trt[, gi, drop = FALSE] / P_ctrl[, gi, drop = FALSE])
colnames(true_lfc_rna) <- colnames(true_lfc_prot) <- sprintf("t%g", tps)

## ----------------------------------------------- 4. RNA-seq observations ----
cd_rna <- make_coldata(cfg, "rna")
R_true_cell <- cbind(R_ctrl[, gi, drop = FALSE], R_trt[, gi, drop = FALSE])
colnames(R_true_cell) <- cell_levels(cfg)

lib_factor <- exp(rnorm(nrow(cd_rna), 0, 0.18))
mu_rna <- R_true_cell[, as.character(cd_rna$cell), drop = FALSE]
mu_rna <- sweep(mu_rna, 2, lib_factor, "*")
mu_rna <- mu_rna * (S$rna_lib_size / mean(colSums(mu_rna)))

disp   <- S$rna_dispersion_a0 + S$rna_dispersion_a1 / rowMeans(mu_rna)
counts <- matrix(rnbinom(length(mu_rna), mu = as.vector(mu_rna),
                         size = rep(1 / disp, times = ncol(mu_rna))),
                 nrow = G, dimnames = list(gene_id, cd_rna$sample_id))

## ------------------------------------------ 5. protein-group construction ---
has_prot  <- runif(G) > S$frac_rna_only     # Tier 2 (RNA-only) genes excluded
prot_gene <- gene_id[has_prot]
np        <- length(prot_gene)
log_step("genes with protein evidence: ", np, " / ", G)

map <- data.frame(gene_id = prot_gene,
                  group   = sprintf("PG%05d", seq_len(np)),
                  frac    = 1, stringsAsFactors = FALSE)

## (a) one gene -> several protein groups (isoform / PTM-form splitting)
n_multi <- round(S$frac_gene_multi_group * np)
if (n_multi > 0) {
  gsel  <- sample(prot_gene, n_multi)
  extra <- data.frame(gene_id = gsel,
                      group   = sprintf("PG%05db", match(gsel, prot_gene)),
                      frac    = 0.45, stringsAsFactors = FALSE)
  map$frac[map$gene_id %in% gsel] <- 0.55
  map <- rbind(map, extra)
}

## (b) several genes -> one protein group (shared peptides, paralogue families)
n_shared <- round(S$frac_shared_peptide_groups * np)
if (n_shared > 0) {
  fam <- split(sample(prot_gene, n_shared * 2), rep(seq_len(n_shared), each = 2))
  for (k in seq_along(fam)) {
    map$group[map$gene_id %in% fam[[k]]] <- sprintf("PGSHARED%04d", k)
  }
}

## Latent protein level per group = weighted sum over its member genes.
groups      <- sort(unique(map$group))
P_true_cell <- cbind(P_ctrl[, gi, drop = FALSE], P_trt[, gi, drop = FALSE])
colnames(P_true_cell) <- cell_levels(cfg)

Pg <- matrix(0, length(groups), ncol(P_true_cell),
             dimnames = list(groups, colnames(P_true_cell)))
for (i in seq_len(nrow(map))) {
  Pg[map$group[i], ] <- Pg[map$group[i], ] +
    map$frac[i] * P_true_cell[map$gene_id[i], ]
}

## ------------------------------------------ 6. proteomics observations ------
cd_prot   <- make_coldata(cfg, "prot")          # 24 DIFFERENT samples: Case B
ms_offset <- rnorm(length(groups), 0, 1.5)      # ionisation efficiency

int <- log2(Pg[, as.character(cd_prot$cell), drop = FALSE] + 1e-6) + ms_offset
int <- int + rnorm(length(int), 0, S$prot_noise_sd)
int <- sweep(int, 2, rnorm(ncol(int), 0, 0.12), "+")   # per-run loading drift
dimnames(int) <- list(groups, cd_prot$sample_id)

## Missingness = MNAR (intensity-dependent censoring) + MAR (stochastic).
mid    <- as.numeric(quantile(int, S$prot_mnar_midpoint_q, na.rm = TRUE))
p_mnar <- expit(-(int - mid) / S$prot_mnar_slope)
miss   <- (matrix(runif(length(int)), nrow(int)) < p_mnar) |
          (matrix(runif(length(int)), nrow(int)) < S$prot_mar_rate)
int[miss] <- NA_real_
log_step("protein missing rate: ", sprintf("%.1f%%", 100 * mean(is.na(int))))

## ------------------------------------ 7. MaxQuant-style metadata + noise ----
acc_of_gene <- setNames(sprintf("P%05d", seq_len(G)), gene_id)
grp_acc <- vapply(groups, function(gr)
  paste(acc_of_gene[map$gene_id[map$group == gr]], collapse = ";"), character(1))

## A few groups carry an accession with no gene mapping at all (Tier 3).
n_orphan <- max(3L, as.integer(round(0.01 * length(groups))))
orphan   <- sample(groups, n_orphan)
grp_acc[orphan] <- paste0(grp_acc[orphan], ";Q9ORPH", seq_len(n_orphan))

pmeta <- data.frame(
  group_id                = groups,
  majority_protein_ids    = unname(grp_acc[groups]),
  unique_peptides         = pmax(1, rpois(length(groups), 5)),
  reverse                 = "",
  potential_contaminant   = "",
  only_identified_by_site = "",
  stringsAsFactors = FALSE)

## Decoy / contaminant / site-only rows that script 01 must remove.
n_junk   <- as.integer(round(0.03 * length(groups)))
junk_id  <- sprintf("JUNK%04d", seq_len(n_junk))
junk_int <- matrix(rnorm(n_junk * ncol(int), mean(int, na.rm = TRUE), 2),
                   n_junk, ncol(int), dimnames = list(junk_id, colnames(int)))
junk_int[matrix(runif(length(junk_int)) < 0.4, n_junk)] <- NA
int   <- rbind(int, junk_int)
pmeta <- rbind(pmeta, data.frame(
  group_id                = junk_id,
  majority_protein_ids    = paste0("REV__P", seq_len(n_junk)),
  unique_peptides         = pmax(1, rpois(n_junk, 2)),
  reverse                 = rep(c("+", "", ""),  length.out = n_junk),
  potential_contaminant   = rep(c("", "+", ""),  length.out = n_junk),
  only_identified_by_site = rep(c("", "", "+"),  length.out = n_junk),
  stringsAsFactors = FALSE))

## UniProt-style raw mapping table, including unmappable accessions.
id_map_raw <- rbind(
  data.frame(accession = unname(acc_of_gene), gene_id = names(acc_of_gene),
             status = "reviewed", stringsAsFactors = FALSE),
  data.frame(accession = sprintf("Q9ORPH%d", seq_len(n_orphan)),
             gene_id = NA_character_, status = "unmapped",
             stringsAsFactors = FALSE))

## ------------------------------------------------------------ 8. truth -----
truth <- data.frame(
  gene_id = gene_id, archetype = archetype, log2_base = log2_base,
  t_half_h = t_half, kd = kd, ks = ks, A_time = A_time, A_trt = A_trt,
  delta_log_ks = ds, delta_log_kd = dd, shape = shape_type,
  has_protein = has_prot, stringsAsFactors = FALSE)
truth <- cbind(truth,
  setNames(as.data.frame(true_lfc_rna),  paste0("true_lfc_rna_",  colnames(true_lfc_rna))),
  setNames(as.data.frame(true_lfc_prot), paste0("true_lfc_prot_", colnames(true_lfc_prot))))

## ------------------------------------------------------------- 9. write ----
write_mat(counts, file.path(outdir, "rna_counts.tsv"), "gene_id")
write_tsv(cd_rna,     file.path(outdir, "rna_coldata.tsv"))
write_mat(int,        file.path(outdir, "protein_intensities.tsv"), "group_id")
write_tsv(cd_prot,    file.path(outdir, "protein_coldata.tsv"))
write_tsv(pmeta,      file.path(outdir, "protein_meta.tsv"))
write_tsv(id_map_raw, file.path(outdir, "id_map_raw.tsv"))
write_tsv(truth,      file.path(outdir, "truth.tsv"))
write_tsv(map,        file.path(outdir, "truth_gene_group_map.tsv"))

log_step("wrote simulated data to ", outdir)
print(table(truth$archetype))
cat("\nRNA counts:     ", nrow(counts), "x", ncol(counts),
    "\nProtein groups: ", nrow(int), "x", ncol(int),
    "\nShared sample IDs between layers:",
    length(intersect(colnames(counts), colnames(int))),
    "  <- must be 0 for Case B\n")
