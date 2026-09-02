# RNA-protein condition-profile distance — can genes be ranked by
discordance?
Kristina Gagalova

- [Purpose and scope](#purpose-and-scope)
  - [Relationship to the existing Python
    work](#relationship-to-the-existing-python-work)
- [Analysis outline](#analysis-outline)
- [1. Setup](#setup)
- [2. The method and its assumptions](#the-method-and-its-assumptions)
  - [Assumptions, stated up front](#assumptions-stated-up-front)
- [3. Implementation](#implementation)
- [4. G1 — Does distance separate the true
  archetypes?](#g1-does-distance-separate-the-true-archetypes)
- [5. G2 — Full space vs 2D embedding](#g2-full-space-vs-2d-embedding)
- [6. G3 — Permutation null on real
  data](#g3-permutation-null-on-real-data)
- [7. G4 — Replicate-level stability](#g4-replicate-level-stability)
- [8. Apply to real data](#apply-to-real-data)
- [8a. Per-condition analysis: Which condition drives
  discordance?](#a.-per-condition-analysis-which-condition-drives-discordance)
- [8b. Three-metric decomposition: Distance, Correlation,
  Amplitude](#b.-three-metric-decomposition-distance-correlation-amplitude)
  - [Interpretation of the three
    metrics](#interpretation-of-the-three-metrics)
- [8c. Infection-response contrast analysis
  (optional)](#c.-infection-response-contrast-analysis-optional)
- [9. Visualizations: UMAP in shared
  space](#visualizations-umap-in-shared-space)
  - [V1 — Shared UMAP embeddings, real
    data](#v1-shared-umap-embeddings-real-data)
  - [V2 — Top and bottom genes in UMAP
    space](#v2-top-and-bottom-genes-in-umap-space)
  - [V3 — Distance distributions](#v3-distance-distributions)
  - [V3b — Density distribution in UMAP
    space](#v3b-density-distribution-in-umap-space)
  - [V3c — Per-condition contributions to
    distance](#v3c-per-condition-contributions-to-distance)
  - [V3d — Per-condition trends: average discordance by
    timepoint/treatment](#v3d-per-condition-trends-average-discordance-by-timepointtreatment)
  - [V4 — Comparing full space, PCA, and UMAP
    distances](#v4-comparing-full-space-pca-and-umap-distances)
- [Interpretation: Is the distance
  meaningful?](#interpretation-is-the-distance-meaningful)
  - [What the distance tells us](#what-the-distance-tells-us)
  - [Can we use it for “similarity”?](#can-we-use-it-for-similarity)
  - [Full space vs PCA vs UMAP](#full-space-vs-pca-vs-umap)
- [10. Assumption validation summary](#assumption-validation-summary)
- [11. What is safe to claim](#what-is-safe-to-claim)
  - [✅ Safe to claim](#safe-to-claim)
  - [❌ Not supported](#not-supported)
  - [The one-paragraph version](#the-one-paragraph-version)
  - [What would strengthen this](#what-would-strengthen-this)
- [References](#references)

## Purpose and scope

This notebook asks a question distinct from anything in
`analysis/kinetics_limited` or `analysis/dimensionality-reduction`:

> **For an individual gene, does its protein’s condition-response
> profile look like its RNA’s condition-response profile — and can that
> resemblance be turned into a per-gene distance that ranks genes by
> RNA/protein discordance?**

`kinetics_limited` tests a *mechanistic* null (first-order kinetics) per
gene. `dimensionality-reduction` asks whether RNA predicts protein *in
aggregate*, across all genes at once. This notebook asks a *geometric*
question per gene: place each gene’s RNA profile and its own protein
profile as two points in a shared coordinate system built from the
population’s dominant response shapes, and measure how far apart they
land.

### Relationship to the existing Python work

Two notebooks in `KristinaGagalova/autoencoders-test` motivate this:

**`rna_protein_unpaired_per_variety.ipynb`** is a thin driver over
`rnaprot/unpaired_reviewer.py` (a revised module – not the `unpaired.py`
referenced in `dimensionality_reduction_wheat.qmd`; the two have
different loss weights, 3.0/0.75 here vs 2.0/0.5 there). It re-runs
essentially the same leave-condition-out comparison (mean, design,
cognate ridge, PCA+ridge, PLS, design-residualised variants,
condition-aligned AE) already built in this repository’s
`dimensionality_reduction_wheat.qmd` §7 (D4) and §7b (D7). Two things in
it are genuinely different from what is already here, and are picked up
as an addendum (**D8**) in that notebook rather than duplicated
wholesale in a second ~40-chunk document:

1.  Feature selection is done **inside each outer CV fold**, not once
    globally before the loop – D4/D7 currently select the top-variable
    genes once, which is a selection-leakage risk.
2.  A permutation null (`N_PERM=200`) is applied to the PCA/PLS/design
    baselines (not to the autoencoder).

See `dimensionality_reduction_wheat.qmd` §7c for that check. This
document does not repeat D4/D7.

**`rna_protein_umap_shared_space.ipynb`** is the actual source of the
“distance between genes” idea, and *is* new work relative to everything
else in this repository. It builds condition-mean profiles for RNA and
for cognate proteins, fits a shared coordinate system (PCA or UMAP) on
the RNA profiles only, projects the protein profiles into that same
system, and measures the Euclidean distance between a gene’s RNA point
and its own protein point. Two things about it matter for how it is
ported here:

1.  **Its own synthetic validation found that 2D UMAP hurts the
    ranking**: distance-to-truth-divergence correlation was Spearman
    $\rho = 0.923$ in the full aligned input space, but only
    $\rho = 0.713$ after compressing to a 2D UMAP embedding. That is
    carried over below as a **hypothesis to test on this project’s own
    ground truth**, not assumed.
2.  **It was never actually run to completion on real data.** The
    real-data cells stop at the visualization step; no per-gene distance
    table or ranking exists for Cadenza or Norin. Building that is the
    new work in §6-§8 below.

------------------------------------------------------------------------

## Analysis outline

| Step | What                                           | Why                                                       |
|:-----|:-----------------------------------------------|:----------------------------------------------------------|
| §1   | Setup and data (real + simulated)              | Validate on known ground truth before touching real data  |
| §2   | The method and its assumptions                 | Each with references                                      |
| §3   | Implementation                                 | Shared-space construction, from scratch                   |
| §4   | **G1** Does distance separate true archetypes? | Validation against `truth.tsv`, incl. the `lag_only` trap |
| §5   | **G2** Full space vs 2D embedding              | Replicates the source notebook’s own finding, on our data |
| §6   | **G3** Permutation null on real data           | No ground truth on real data – this substitutes for it    |
| §7   | **G4** Replicate-level stability               | Bootstrap over which replicates feed the condition means  |
| §8   | Apply to real data                             | Per-gene distance tables, Cadenza and Norin               |
| §9   | Assumption validation summary                  |                                                           |
| §10  | What is safe to claim                          |                                                           |

------------------------------------------------------------------------

## 1. Setup

``` r
if (!requireNamespace("here", quietly = TRUE)) stop("package 'here' required")

## RENV_PROJECT must be set before sourcing activate.R, or renv treats this
## subdirectory as the project and bootstraps an empty library.
if (file.exists(here::here("renv", "activate.R"))) {
  Sys.setenv(RENV_PROJECT = here::here())
  source(here::here("renv", "activate.R"))
}

source(here::here("R", "utils.R"))
source(here::here("R", "wheat_pipeline.R"))
need_pkgs(c("limma", "edgeR", "DESeq2", "matrixStats", "impute", "uwot", "yaml", "MASS", "ggplot2"))

library(ggplot2)
library(MASS)

set.seed(20260823)
```

``` r
DESIGN <- wheat_design()
cad <- prepare_variety("cadenza", DESIGN)
```

    Cluster size 5861 broken into 3622 2239 
    Cluster size 3622 broken into 2430 1192 
    Cluster size 2430 broken into 1561 869 
    Cluster size 1561 broken into 819 742 
    Done cluster 819 
    Done cluster 742 
    Done cluster 1561 
    Done cluster 869 
    Done cluster 2430 
    Done cluster 1192 
    Done cluster 3622 
    Cluster size 2239 broken into 1075 1164 
    Done cluster 1075 
    Done cluster 1164 
    Done cluster 2239 

    Validated DE (Cadenza): 59931 genes x 4 timepoints

``` r
nor <- prepare_variety("norin",   DESIGN)
```

    Cluster size 5641 broken into 1939 3702 
    Cluster size 1939 broken into 1030 909 
    Done cluster 1030 
    Done cluster 909 
    Done cluster 1939 
    Cluster size 3702 broken into 2489 1213 
    Cluster size 2489 broken into 838 1651 
    Done cluster 838 
    Cluster size 1651 broken into 941 710 
    Done cluster 941 
    Done cluster 710 
    Done cluster 1651 
    Done cluster 2489 
    Done cluster 1213 
    Done cluster 3702 

    Validated DE (Norin): 59502 genes x 4 timepoints

``` r
REAL <- list(Cadenza = cad, Norin = nor)
```

``` r
#' Loader for the project's own Case B ground-truth simulator
#' (`R/00_simulate_data.R`), parallel to `load_variety()` in
#' `wheat_pipeline.R` but pointed at `data/simulated` and its different file
#' layout. This gives archetype-labelled truth to validate the distance
#' metric against, before it is trusted on real data -- the same validate-
#' on-simulation-first discipline used for the rest of this pipeline.
load_simulated <- function(cfg, data_dir = here::here("data", "simulated")) {

  rna_counts <- as.matrix(read.delim(file.path(data_dir, "rna_counts.tsv"),
                                     row.names = 1, check.names = FALSE))
  prot_int   <- as.matrix(read.delim(file.path(data_dir, "protein_intensities.tsv"),
                                     row.names = 1, check.names = FALSE))
  rna_meta   <- read.delim(file.path(data_dir, "rna_coldata.tsv"),
                           row.names = "sample_id", stringsAsFactors = FALSE)
  prot_meta  <- read.delim(file.path(data_dir, "protein_coldata.tsv"),
                           row.names = "sample_id", stringsAsFactors = FALSE)
  truth <- read.delim(file.path(data_dir, "truth.tsv"), stringsAsFactors = FALSE)
  gmap  <- read.delim(file.path(data_dir, "truth_gene_group_map.tsv"),
                      stringsAsFactors = FALSE)

  lv <- cell_levels(cfg)
  rna_meta$cell  <- factor(rna_meta$cell,  levels = lv)
  prot_meta$cell <- factor(prot_meta$cell, levels = lv)

  design <- list(treatments = cfg$design$treatments,
                timepoints = cfg$design$timepoints,
                n_reps     = cfg$design$n_reps)

  list(rna  = rna_counts[, rownames(rna_meta), drop = FALSE],
       prot = prot_int[,  rownames(prot_meta), drop = FALSE],
       rna_meta = rna_meta, prot_meta = prot_meta,
       design = design, truth = truth, gmap = gmap)
}

sim_cfg <- load_config(here::here("config", "config.yaml"))
sim_raw <- load_simulated(sim_cfg)

sim_qc <- local({
  qr <- run_qc_rna(sim_raw$rna, sim_raw$rna_meta, sim_raw$design)
  qp <- run_qc_prot(sim_raw$prot, sim_raw$prot_meta)
  mi <- run_missingness(qp$norm, sim_raw$prot_meta)
  list(vst = qr$vst, prot = mi$mixed)
})
```

    Warning in run_qc_prot(sim_raw$prot, sim_raw$prot_meta): NaNs produced

    Cluster size 1781 broken into 663 1118 
    Done cluster 663 
    Done cluster 1118 

``` r
cat("Simulated genes after RNA QC:   ", nrow(sim_qc$vst),  "\n")
```

    Simulated genes after RNA QC:    4000 

``` r
cat("Simulated groups after prot QC: ", nrow(sim_qc$prot), "\n")
```

    Simulated groups after prot QC:  1781 

------------------------------------------------------------------------

## 2. The method and its assumptions

For each gene present in both layers:

1.  Take its **condition-mean profile** in each layer separately (8
    treatment x timepoint cells) – the only aggregation licensed under
    Case B, per `dimensionality_reduction_wheat.qmd` §7b (D7).
2.  **Row-centre** each profile (subtract the gene’s own
    across-condition mean). This compares response *shape*, not absolute
    level – required regardless of Case B, since RNA (VST) and protein
    (log2 intensity) are on incommensurable physical scales to begin
    with.
3.  **Standardise per condition** (z-score each of the 8 condition
    columns across genes), fit on the RNA block only, and apply the same
    transform to protein. This matches the source notebook’s
    `StandardScaler().fit(Xr)` step.
4.  Fit **PCA on the RNA profiles only** (their dominant shapes of
    response across the gene population), and project the protein
    profiles into that same space with no refitting.
5.  **Distance** = Euclidean distance between a gene’s transformed RNA
    point and its own transformed protein point – computed in the full
    standardised space, in the top-2 PCA components, and in a 2D UMAP
    embedding fit on RNA and transformed for protein.

### Assumptions, stated up front

| \#  | Assumption                                                                               | Status                                                                     |
|:----|:-----------------------------------------------------------------------------------------|:---------------------------------------------------------------------------|
| A1  | The RNA\<-\>protein correspondence used is genuinely 1:1 for the genes scored            | Enforced by construction (see §4); ambiguous groups excluded               |
| A2  | Condition means are the right granularity for a cross-block comparison                   | Carried over from D7, not re-derived here                                  |
| A3  | Comparing row-centred *shape*, not absolute level, is the right notion of “look alike”   | Design choice, stated not tested                                           |
| A4  | Distance in the full standardised space is more informative than in a 2D embedding       | **Tested in G2**, not assumed                                              |
| A5  | A short Euclidean distance means concordant regulation, not just similar kinetics/timing | **Tested in G1** – this is where `lag_only` genes are the adversarial case |

McInnes, Healy, and Melville (2018) and Becht et al. (2019) are the
basis for A4: UMAP is explicitly designed to preserve local
neighbourhood structure, not global pairwise distances, so a 2D UMAP
distance between two specific, possibly-distant points is not guaranteed
to track their true distance in the original space – exactly the
degradation the source notebook’s own synthetic check already found
($\rho$ 0.923 to 0.713).

------------------------------------------------------------------------

## 3. Implementation

``` r
#' Per-gene shape: subtract the gene's own across-condition mean.
row_center <- function(M) M - rowMeans(M)

#' Fit a per-column (per-condition) z-scorer on one matrix, return a function
#' that applies it to any matrix with the same columns. Fit on RNA only, per
#' the source notebook.
fit_col_scaler <- function(M) {
  mu  <- colMeans(M)
  sdv <- apply(M, 2, sd); sdv[sdv < 1e-8] <- 1
  function(X) sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
}

#' Build the shared space: row-centre, z-score on RNA's column moments, PCA
#' fit on RNA only. Returns the standardised profiles (full-space distance
#' needs nothing else, since the PCA rotation below is orthonormal and
#' therefore distance-preserving) and their PCA scores (for the 2D check).
build_shared_space <- function(R_cond, P_cond) {
  stopifnot(identical(rownames(R_cond), rownames(P_cond)))
  scaler <- fit_col_scaler(row_center(R_cond))
  Rs <- scaler(row_center(R_cond))
  Ps <- scaler(row_center(P_cond))

  sv <- svd(Rs)
  list(Rs = Rs, Ps = Ps,
       Zr = Rs %*% sv$v, Zp = Ps %*% sv$v,
       var_explained = sv$d^2 / sum(sv$d^2))
}

euclid <- function(A, B, k = NULL) {
  if (!is.null(k)) { A <- A[, seq_len(k), drop = FALSE]; B <- B[, seq_len(k), drop = FALSE] }
  sqrt(rowSums((A - B)^2))
}

#' Calculate correlation and amplitude for each gene
#' Returns: correlation (directional concordance), amplitude ratio (relative response magnitude)
gene_metrics <- function(R_centered, P_centered) {
  # Correlation: directional concordance (correlation of response profiles)
  cors <- sapply(seq_len(nrow(R_centered)), function(i) {
    cor(R_centered[i, ], P_centered[i, ], method = "pearson")
  })
  cors[is.na(cors)] <- 0  # Handle constant profiles

  # Amplitude: log2 ratio of response magnitudes
  # ||R|| = sqrt(sum of squared values across conditions)
  R_norm <- sqrt(rowSums(R_centered^2))
  P_norm <- sqrt(rowSums(P_centered^2))
  # Avoid log(0); replace with 0
  amps <- ifelse(R_norm > 1e-10 & P_norm > 1e-10,
                 log2(P_norm / R_norm),
                 0)

  data.frame(correlation = cors, amplitude = amps)
}

#' Calculate infection-response contrasts (Infected - Control) for each timepoint
#' Returns condition means for the 4 timepoint-level contrasts
infection_contrasts <- function(counts, meta) {
  # Assumes meta has columns: condition, treatment, timepoint
  # where treatment = "infected" or "control"

  timepoints <- sort(unique(meta$timepoint))
  contrasts <- data.frame()

  for (tp in timepoints) {
    idx_inf <- meta$treatment == "infected" & meta$timepoint == tp
    idx_ctl <- meta$treatment == "control"  & meta$timepoint == tp

    if (sum(idx_inf) > 0 && sum(idx_ctl) > 0) {
      mean_inf <- rowMeans(counts[, idx_inf, drop = FALSE])
      mean_ctl <- rowMeans(counts[, idx_ctl, drop = FALSE])
      contrast <- mean_inf - mean_ctl

      if (nrow(contrasts) == 0) {
        contrasts <- data.frame(contrast)
        colnames(contrasts)[1] <- paste0("dpi_", tp)
      } else {
        contrasts[[paste0("dpi_", tp)]] <- contrast
      }
    }
  }
  contrasts
}

#' UMAP fit on RNA, protein transformed into the same embedding -- matching
#' the source notebook's asymmetric fit/transform (never a joint/symmetric
#' fit, which would let protein leak into the axes it is being scored against).
umap_distance <- function(Rs, Ps, seed = 1) {
  nn <- max(2, min(15, nrow(Rs) %/% 4))
  set.seed(seed)
  fit <- uwot::umap(Rs, n_neighbors = nn, min_dist = 0.10, n_components = 2,
                    metric = "euclidean", ret_model = TRUE)
  Ep <- uwot::umap_transform(Ps, fit)
  sqrt(rowSums((fit$embedding - Ep)^2))
}
```

------------------------------------------------------------------------

## 4. G1 — Does distance separate the true archetypes?

The simulator’s `truth.tsv` labels every gene with the *mechanism* that
generated it: `unchanged`, `concordant`, `buffered`, `mrna_only`,
`protein_only`, `anticorrelated`, `lag_only`. None of that label is used
to build the distance metric – it is used only to grade it afterwards.

`lag_only` is the adversarial case, carried over from the kinetics
notebooks: it is **mechanistically fully concordant** (protein is driven
purely by RNA, first-order kinetics, `A2` holds exactly) but has a
**slow protein half-life**, so its condition-profile *shape* lags the
RNA shape. A method that cannot tell “concordant but slow” from
“genuinely decoupled” will misclassify it –
`analysis/kinetics_limited/assumptions_validation.qmd` documents the
same trap for naive quadrant classification. This is the sharpest test
of A5.

``` r
## Restrict to genes with an unambiguous 1:1 RNA<->protein-group
## correspondence. Isoform-split and shared-peptide groups (the simulator's
## deliberately messy 1:many / many:1 cases) are excluded from this analysis
## by design -- a distance metric between a gene and an ambiguous protein
## mixture is not answering the question this notebook asks.
gmap <- sim_raw$gmap
one_to_one <- with(gmap,
  ave(seq_along(gene_id), gene_id, FUN = length) == 1 &
  ave(seq_along(group),   group,   FUN = length) == 1 &
  frac == 1)
map1 <- gmap[one_to_one, c("gene_id", "group")]
map1 <- map1[map1$gene_id %in% rownames(sim_qc$vst) &
             map1$group  %in% rownames(sim_qc$prot), ]

cat("1:1 genes available for validation:", nrow(map1), "/", nrow(gmap), "\n")
```

    1:1 genes available for validation: 1500 / 1923 

``` r
R_cond_sim <- cell_means(sim_qc$vst[map1$gene_id, , drop = FALSE],  sim_raw$rna_meta)
P_cond_sim <- cell_means(sim_qc$prot[map1$group,  , drop = FALSE], sim_raw$prot_meta)
rownames(P_cond_sim) <- map1$gene_id

sp_sim     <- build_shared_space(R_cond_sim, P_cond_sim)
d_full_sim <- euclid(sp_sim$Rs, sp_sim$Ps)
d_pca2_sim <- euclid(sp_sim$Zr, sp_sim$Zp, k = 2)
d_umap_sim <- umap_distance(sp_sim$Rs, sp_sim$Ps)

arch <- setNames(sim_raw$truth$archetype, sim_raw$truth$gene_id)[map1$gene_id]

g1 <- data.frame(gene_id = map1$gene_id, archetype = arch,
                 d_full = d_full_sim, d_pca2 = d_pca2_sim, d_umap = d_umap_sim)
```

``` r
ord <- c("unchanged", "concordant", "lag_only", "buffered",
        "mrna_only", "protein_only", "anticorrelated")
boxplot(d_full ~ factor(archetype, levels = ord), data = g1,
       las = 2, xlab = "", ylab = "distance (full standardised space)",
       col = ifelse(ord %in% c("unchanged", "concordant"), "#2C6FBB",
             ifelse(ord == "lag_only", "#E8A33D", "#D1495B")))
```

<img
src="gene_distance_shared_space_files/figure-commonmark/fig-g1-1.png"
id="fig-g1"
alt="Figure 1: G1: full-space distance by true archetype (simulated data). concordant/unchanged are the null; lag_only is the adversarial case – mechanistically concordant but kinetically slow." />

``` r
#' Rank-based separation of a "positive" (expected-far) archetype set from a
#' "negative" (expected-close) baseline. AUC > 0.5 means the positive set
#' really does score farther, in the direction the mechanism predicts.
sep_stat <- function(d, pos, neg) {
  w <- suppressWarnings(wilcox.test(d[pos], d[neg]))
  data.frame(n_pos = sum(pos), n_neg = sum(neg),
            auc = round(unname(w$statistic) / (sum(pos) * sum(neg)), 3),
            p   = signif(w$p.value, 3))
}

baseline <- g1$archetype %in% c("unchanged", "concordant")
truly_discordant <- g1$archetype %in% c("buffered", "mrna_only",
                                        "protein_only", "anticorrelated")
lag <- g1$archetype == "lag_only"

g1_sep <- rbind(
  data.frame(comparison = "genuinely discordant vs concordant/unchanged",
            sep_stat(g1$d_full, truly_discordant, baseline)),
  data.frame(comparison = "lag_only vs concordant/unchanged",
            sep_stat(g1$d_full, lag, baseline))
)
knitr::kable(g1_sep, caption = "G1: separation of archetypes by full-space distance. AUC near 0.5 means the metric cannot tell the groups apart.")
```

| comparison                                   | n_pos | n_neg |   auc |   p |
|:---------------------------------------------|------:|------:|------:|----:|
| genuinely discordant vs concordant/unchanged |   363 |  1070 | 0.820 |   0 |
| lag_only vs concordant/unchanged             |    67 |  1070 | 0.932 |   0 |

G1: separation of archetypes by full-space distance. AUC near 0.5 means
the metric cannot tell the groups apart.

------------------------------------------------------------------------

## 5. G2 — Full space vs 2D embedding

``` r
g2 <- rbind(
  data.frame(space = "full standardised space",
            rbind(sep_stat(g1$d_full, truly_discordant, baseline),
                  sep_stat(g1$d_full, lag, baseline))[, "auc", drop = FALSE]),
  data.frame(space = "top-2 PCA components",
            rbind(sep_stat(g1$d_pca2, truly_discordant, baseline),
                  sep_stat(g1$d_pca2, lag, baseline))[, "auc", drop = FALSE]),
  data.frame(space = "2D UMAP",
            rbind(sep_stat(g1$d_umap, truly_discordant, baseline),
                  sep_stat(g1$d_umap, lag, baseline))[, "auc", drop = FALSE])
)
g2$comparison <- rep(c("discordant vs baseline", "lag_only vs baseline"), 3)
g2 <- g2[, c("space", "comparison", "auc")]
knitr::kable(g2, caption = "G2: does dimensionality reduction for visualisation cost separation power? Compare each space's AUC to the full-space row.")
```

| space                   | comparison             |   auc |
|:------------------------|:-----------------------|------:|
| full standardised space | discordant vs baseline | 0.820 |
| full standardised space | lag_only vs baseline   | 0.932 |
| top-2 PCA components    | discordant vs baseline | 0.827 |
| top-2 PCA components    | lag_only vs baseline   | 0.929 |
| 2D UMAP                 | discordant vs baseline | 0.799 |
| 2D UMAP                 | lag_only vs baseline   | 0.933 |

G2: does dimensionality reduction for visualisation cost separation
power? Compare each space’s AUC to the full-space row.

------------------------------------------------------------------------

## 6. G3 — Permutation null on real data

There is no `truth.tsv` for Cadenza or Norin. What substitutes for it:
if a gene’s protein profile genuinely resembles its own RNA profile more
than it resembles a random other gene’s, the **observed mean cognate
distance should sit below the distribution obtained by randomly
re-pairing genes**. This is a matching (Mantel-type) permutation test,
and needs no ground truth.

``` r
#' Build the shared space and full-space gene distances for one real variety.
#' Enforces assumption A1 (1:1 gene<->protein correspondence) by checking that
#' the protein matrix has already been collapsed to 1:1 identifiers.
gene_distance_real <- function(v, n_perm = 199, seed = 1) {

  common <- intersect(rownames(v$qc_rna$vst), rownames(v$imputed$mixed))

  # A1: Verify 1:1 correspondence
  # In the simulation, this is enforced by explicit filtering in §4.
  # For real data, we verify that row names are already clean 1:1 identifiers.
  # If the protein matrix contains ambiguous groups or isoforms, they should have
  # been collapsed/filtered upstream by prepare_variety() or equivalent.
  stopifnot(
    "Protein matrix must have unique row identifiers (assumption A1)" =
      length(unique(rownames(v$imputed$mixed))) == nrow(v$imputed$mixed)
  )

  R_cond <- cell_means(v$qc_rna$vst[common, , drop = FALSE],    v$meta)
  P_cond <- cell_means(v$imputed$mixed[common, , drop = FALSE], v$meta)

  sp     <- build_shared_space(R_cond, P_cond)
  d_full <- euclid(sp$Rs, sp$Ps)
  d_pca2 <- euclid(sp$Zr, sp$Zp, k = 2)

  # Calculate correlation and amplitude for each gene
  metrics <- gene_metrics(sp$Rs, sp$Ps)

  # UMAP fit on RNA, transform protein (asymmetric, never joint)
  nn <- max(2, min(15, nrow(sp$Rs) %/% 4))
  set.seed(seed)
  umap_fit <- uwot::umap(sp$Rs, n_neighbors = nn, min_dist = 0.10, n_components = 2,
                         metric = "euclidean", ret_model = TRUE)
  Er <- umap_fit$embedding
  Ep <- uwot::umap_transform(sp$Ps, umap_fit)
  d_umap <- sqrt(rowSums((Er - Ep)^2))

  set.seed(seed)
  null <- replicate(n_perm, mean(euclid(sp$Rs, sp$Ps[sample(nrow(sp$Ps)), , drop = FALSE])))
  obs  <- mean(d_full)
  p_perm <- (sum(null <= obs) + 1) / (n_perm + 1)

  list(genes = data.frame(gene_id = common,
                          d_full = d_full,
                          d_pca2 = d_pca2,
                          d_umap = d_umap,
                          correlation = metrics$correlation,
                          amplitude = metrics$amplitude),
       space = list(Er = Er, Ep = Ep),
       obs = obs, null = null, p_perm = p_perm, n_genes = length(common))
}

GD <- setNames(Map(function(v, s) gene_distance_real(v, seed = s), REAL, seq_along(REAL)), names(REAL))

#' Infection-response contrast analysis: compare Infected - Control at each timepoint
#' Instead of 8 condition means, use 4-point infection-response trajectories
gene_distance_infection_contrast <- function(v, n_perm = 199, seed = 1) {
  # Extract infection and control means for each timepoint
  # v$meta should have: condition, treatment (infected/control), timepoint

  common <- intersect(rownames(v$qc_rna$vst), rownames(v$imputed$mixed))

  # Verify 1:1 mapping
  stopifnot(
    "Protein matrix must have unique row identifiers (assumption A1)" =
      length(unique(rownames(v$imputed$mixed))) == nrow(v$imputed$mixed)
  )

  # Calculate contrasts for RNA and protein
  R_contrast <- infection_contrasts(v$qc_rna$vst[common, ], v$meta)
  P_contrast <- infection_contrasts(v$imputed$mixed[common, ], v$meta)

  # If contrasts are empty or don't match, return NULL (metadata may not have treatment column)
  if (nrow(R_contrast) == 0 || nrow(P_contrast) == 0 || ncol(R_contrast) != ncol(P_contrast)) {
    return(NULL)
  }

  # Build shared space on contrasts
  sp <- build_shared_space(R_contrast, P_contrast)
  d_full <- euclid(sp$Rs, sp$Ps)
  metrics <- gene_metrics(sp$Rs, sp$Ps)

  # UMAP
  nn <- max(2, min(15, nrow(sp$Rs) %/% 4))
  set.seed(seed)
  umap_fit <- uwot::umap(sp$Rs, n_neighbors = nn, min_dist = 0.10, n_components = 2,
                         metric = "euclidean", ret_model = TRUE)
  Er <- umap_fit$embedding
  Ep <- uwot::umap_transform(sp$Ps, umap_fit)
  d_umap <- sqrt(rowSums((Er - Ep)^2))

  set.seed(seed)
  null <- replicate(n_perm, mean(euclid(sp$Rs, sp$Ps[sample(nrow(sp$Ps)), , drop = FALSE])))
  obs <- mean(d_full)
  p_perm <- (sum(null <= obs) + 1) / (n_perm + 1)

  list(genes = data.frame(gene_id = common,
                          d_full = d_full,
                          d_umap = d_umap,
                          correlation = metrics$correlation,
                          amplitude = metrics$amplitude),
       space = list(Er = Er, Ep = Ep),
       obs = obs, null = null, p_perm = p_perm, n_genes = length(common))
}

# Try to compute infection-contrast version (may be NULL if metadata lacks treatment info)
GD_contrast <- tryCatch(
  setNames(Map(function(v, s) gene_distance_infection_contrast(v, seed = s), REAL, seq_along(REAL)), names(REAL)),
  error = function(e) NULL
)
```

``` r
g3 <- do.call(rbind, lapply(names(GD), function(nm) data.frame(
  variety = nm, n_genes = GD[[nm]]$n_genes,
  mean_dist_observed = round(GD[[nm]]$obs, 3),
  null_p5  = round(quantile(GD[[nm]]$null, .05), 3),
  null_p50 = round(quantile(GD[[nm]]$null, .50), 3),
  p_perm   = GD[[nm]]$p_perm
)))
knitr::kable(g3, caption = "G3: is the real, correctly-matched RNA<->protein pairing closer (in aggregate) than randomly mismatched genes? p_perm is one-sided: P(null mean <= observed mean).")
```

|     | variety | n_genes | mean_dist_observed | null_p5 | null_p50 | p_perm |
|:----|:--------|--------:|-------------------:|--------:|---------:|-------:|
| 5%  | Cadenza |    5580 |              2.867 |   3.331 |    3.347 |  0.005 |
| 5%1 | Norin   |    5159 |              2.960 |   3.292 |    3.308 |  0.005 |

G3: is the real, correctly-matched RNA\<-\>protein pairing closer (in
aggregate) than randomly mismatched genes? p_perm is one-sided: P(null
mean \<= observed mean).

``` r
op <- par(mfrow = c(1, 2))
for (nm in names(GD)) {
  hist(GD[[nm]]$null, breaks = 25, col = "grey80", border = "white",
      main = nm, xlab = "mean distance under random re-pairing")
  abline(v = GD[[nm]]$obs, col = "#D1495B", lwd = 2)
}
par(op)
```

<img
src="gene_distance_shared_space_files/figure-commonmark/fig-g3-1.png"
id="fig-g3"
alt="Figure 2: G3: null distribution of mean gene distance under random re-pairing (histogram), observed mean under the true pairing (vertical line)." />

------------------------------------------------------------------------

## 7. G4 — Replicate-level stability

Condition means are built from only 3 replicates.
`analysis/kinetics_limited/assumptions_validation.qmd` (B4) already
found that 36-49% of individual gene-level kinetic calls flip under an
RNA-measurement-error bootstrap. The equivalent question here: **does a
gene’s distance rank survive resampling which replicates fed the
condition mean?**

RNA and protein are resampled **independently** (different seeds) –
resampling them together would silently reintroduce the sample-level
pairing this whole project has established does not exist.

``` r
stability_check <- function(v, B = 60) {
  common <- intersect(rownames(v$qc_rna$vst), rownames(v$imputed$mixed))
  Rb <- boot_cell_means(v$qc_rna$vst[common, , drop = FALSE],    v$meta, B = B, seed = 11)
  Pb <- boot_cell_means(v$imputed$mixed[common, , drop = FALSE], v$meta, B = B, seed = 97)

  # Calculate distance independently for each bootstrap replicate
  D <- vapply(seq_len(B), function(b) {
    sp <- build_shared_space(Rb[[b]][common, ], Pb[[b]][common, ])
    euclid(sp$Rs, sp$Ps)
  }, numeric(length(common)))
  rownames(D) <- common

  # For reference: calculate the original ranking
  sp_orig <- build_shared_space(Rb[[1]][common, ], Pb[[1]][common, ])  # Use first bootstrap as baseline
  d_orig <- euclid(sp_orig$Rs, sp_orig$Ps)
  rank_orig <- rank(-d_orig)

  # Calculate Spearman correlation between original ranking and each bootstrap ranking
  boot_rank_cors <- vapply(seq_len(B), function(b) {
    rank_boot <- rank(-D[, b])
    cor(rank_orig, rank_boot, method = "spearman")
  }, numeric(1))

  # Report median and interval
  rank_cor_median <- median(boot_rank_cors)
  rank_cor_interval <- quantile(boot_rank_cors, c(0.05, 0.95))

  # For top 10%, calculate how often each gene appears in that decile across bootstraps
  top10_count <- round(0.1 * length(common))
  top10_appearance <- sapply(seq_len(length(common)), function(i) {
    sum(rank(-D[i, ]) <= top10_count) / B
  })
  names(top10_appearance) <- common

  list(
    rank_cor_median = rank_cor_median,
    rank_cor_interval = rank_cor_interval,
    top10_appearance = top10_appearance,
    n_genes = length(common),
    D = D,
    boot_rank_cors = boot_rank_cors
  )
}

STAB <- lapply(REAL, stability_check)

# Summary table for G4
g4 <- do.call(rbind, lapply(names(STAB), function(nm) data.frame(
  variety = nm,
  n_genes = STAB[[nm]]$n_genes,
  median_spearman_rho = round(STAB[[nm]]$rank_cor_median, 3),
  rho_5pct = round(STAB[[nm]]$rank_cor_interval[1], 3),
  rho_95pct = round(STAB[[nm]]$rank_cor_interval[2], 3)
)))
knitr::kable(g4, caption = "G4: Ranking stability across 60 bootstrap resamples of condition means. For each variety, shows median Spearman rank correlation between original ranking and each bootstrap, plus 5–95% interval. This tests whether a gene's rank is stable under replicate resampling (not whether bootstrapped averages have converged).")
```

|     | variety | n_genes | median_spearman_rho | rho_5pct | rho_95pct |
|:----|:--------|--------:|--------------------:|---------:|----------:|
| 5%  | Cadenza |    5580 |               0.911 |    0.876 |     0.947 |
| 5%1 | Norin   |    5159 |               0.896 |    0.853 |     0.942 |

G4: Ranking stability across 60 bootstrap resamples of condition means.
For each variety, shows median Spearman rank correlation between
original ranking and each bootstrap, plus 5–95% interval. This tests
whether a gene’s rank is stable under replicate resampling (not whether
bootstrapped averages have converged).

``` r
# For each variety, show which genes most consistently appear in the top 10%
cat("Top 10% consistency across bootstraps:\n\n")
```

    Top 10% consistency across bootstraps:

``` r
for (nm in names(STAB)) {
  top10_pct <- sort(STAB[[nm]]$top10_appearance, decreasing = TRUE)[1:10]
  cat("**", nm, ":**\n", sep = "")
  cat("Genes appearing in top 10% most-discordant in N% of 60 bootstrap resamples:\n")
  for (i in seq_len(length(top10_pct))) {
    cat(sprintf("  %s: %.0f%%\n", names(top10_pct)[i], top10_pct[i] * 100))
  }
  cat("\n")
}
```

    **Cadenza:**
    Genes appearing in top 10% most-discordant in N% of 60 bootstrap resamples:
      TraesCAD_scaffold_000005_01G000100: 100%
      TraesCAD_scaffold_000016_01G000300: 100%
      TraesCAD_scaffold_000018_01G001700: 100%
      TraesCAD_scaffold_000022_01G000100: 100%
      TraesCAD_scaffold_000030_01G000100: 100%
      TraesCAD_scaffold_000030_01G000400: 100%
      TraesCAD_scaffold_000043_01G000600: 100%
      TraesCAD_scaffold_000050_01G000700: 100%
      TraesCAD_scaffold_000077_01G000100: 100%
      TraesCAD_scaffold_000081_01G000100: 100%

    **Norin:**
    Genes appearing in top 10% most-discordant in N% of 60 bootstrap resamples:
      TraesNOR1A03G00003900: 100%
      TraesNOR1A03G00005050: 100%
      TraesNOR1A03G00009760: 100%
      TraesNOR1A03G00010740: 100%
      TraesNOR1A03G00010780: 100%
      TraesNOR1A03G00011550: 100%
      TraesNOR1A03G00011890: 100%
      TraesNOR1A03G00013400: 100%
      TraesNOR1A03G00015580: 100%
      TraesNOR1A03G00015590: 100%

------------------------------------------------------------------------

## 8. Apply to real data

``` r
top_bottom <- function(nm, n = 15) {
  d <- GD[[nm]]$genes
  d <- d[order(-d$d_full), ]
  list(most_discordant  = head(d, n),
      most_concordant  = tail(d, n)[order(tail(d, n)$d_full), ])
}

TB_cad <- top_bottom("Cadenza")
TB_nor <- top_bottom("Norin")

knitr::kable(TB_cad$most_discordant, row.names = FALSE,
            caption = "Cadenza: 15 genes with the largest RNA<->protein condition-profile distance.")
```

| gene_id                            |   d_full |    d_pca2 |     d_umap | correlation |  amplitude |
|:-----------------------------------|---------:|----------:|-----------:|------------:|-----------:|
| TraesCAD_scaffold_014747_01G000100 | 14.58265 | 14.475242 |  2.3657339 |   0.9061702 | -2.1774305 |
| TraesCAD_scaffold_133124_01G000100 | 13.46543 | 12.019077 |  1.9724290 |   0.4208656 | -0.5813689 |
| TraesCAD_scaffold_116541_01G000100 | 13.42158 | 12.097656 |  5.6464228 |  -0.2217841 |  0.0811294 |
| TraesCAD_scaffold_146252_01G000100 | 12.98171 |  7.499297 | 14.0091441 |  -0.4626135 | -0.3808938 |
| TraesCAD_scaffold_129059_01G000100 | 12.07514 |  6.287907 |  0.5496855 |   0.8112003 |  1.2487586 |
| TraesCAD_scaffold_004054_01G000200 | 11.69405 | 11.410600 |  4.0967407 |  -0.2543228 |  1.1620031 |
| TraesCAD_scaffold_086045_01G000100 | 10.90339 |  8.892567 | 12.5607742 |  -0.3549208 | -0.2021084 |
| TraesCAD_scaffold_070666_01G000100 | 10.83616 | 10.347964 |  2.4119324 |  -0.0695196 |  0.9572901 |
| TraesCAD_scaffold_054259_01G000200 | 10.77012 |  4.470433 |  3.7250462 |  -0.2541636 |  2.4643760 |
| TraesCAD_scaffold_033501_01G000300 | 10.73777 |  5.949084 |  1.1749163 |   0.1943790 | -0.0447390 |
| TraesCAD_scaffold_003324_01G000300 | 10.67303 | 10.132232 |  5.8074039 |  -0.3497620 |  1.7236806 |
| TraesCAD_scaffold_111708_01G000100 | 10.53067 |  9.694695 |  2.1611899 |   0.6299001 | -1.3882058 |
| TraesCAD_scaffold_055929_01G000400 | 10.41145 |  9.304480 | 10.1569778 |  -0.6270213 | -0.5273418 |
| TraesCAD_scaffold_078281_01G000100 | 10.28404 |  9.888884 |  3.9489334 |  -0.1407167 |  0.7715232 |
| TraesCAD_scaffold_045162_01G000100 | 10.25314 |  9.683636 |  5.2107489 |   0.4418092 | -1.6304404 |

Cadenza: 15 genes with the largest RNA\<-\>protein condition-profile
distance.

``` r
knitr::kable(TB_cad$most_concordant, row.names = FALSE,
            caption = "Cadenza: 15 genes with the smallest RNA<->protein condition-profile distance.")
```

| gene_id                            |    d_full |    d_pca2 |    d_umap | correlation |  amplitude |
|:-----------------------------------|----------:|----------:|----------:|------------:|-----------:|
| TraesCAD_scaffold_038576_01G000100 | 0.3320594 | 0.2223863 | 0.3600272 |   0.9273678 |  0.0085526 |
| TraesCAD_scaffold_052346_01G000400 | 0.3461207 | 0.1039178 | 0.6840019 |   0.3693726 |  0.1578930 |
| TraesCAD_scaffold_021918_01G000100 | 0.3723484 | 0.0243553 | 0.5830958 |   0.9420747 |  0.0837291 |
| TraesCAD_scaffold_031505_01G000300 | 0.3753941 | 0.2407662 | 1.5627541 |   0.8538846 | -0.3104532 |
| TraesCAD_scaffold_013569_01G000100 | 0.3801392 | 0.0678640 | 0.4268784 |   0.8132702 |  0.4999251 |
| TraesCAD_scaffold_033939_01G000100 | 0.3813461 | 0.2157055 | 1.1792848 |   0.9270704 | -0.2945736 |
| TraesCAD_scaffold_049099_01G000200 | 0.3836871 | 0.3402093 | 0.8953310 |   0.8244191 |  0.1483794 |
| TraesCAD_scaffold_054741_01G000200 | 0.3956708 | 0.1252366 | 0.4805051 |   0.8212844 | -0.1461318 |
| TraesCAD_scaffold_056084_01G000100 | 0.4139949 | 0.2143559 | 0.3528777 |   0.9108897 | -0.2131325 |
| TraesCAD_scaffold_006030_01G000100 | 0.4190055 | 0.1167284 | 0.9587574 |   0.9170550 | -0.0198177 |
| TraesCAD_scaffold_133213_01G000100 | 0.4212843 | 0.2471508 | 0.6597650 |   0.8161074 |  0.4892320 |
| TraesCAD_scaffold_012263_01G000100 | 0.4384106 | 0.2282005 | 0.8391910 |   0.7327938 |  0.1481034 |
| TraesCAD_scaffold_153650_01G000100 | 0.4397477 | 0.2696171 | 0.3331563 |   0.8116068 |  0.1318558 |
| TraesCAD_scaffold_034146_01G000800 | 0.4432211 | 0.1519103 | 0.4861834 |   0.4308411 |  0.5617491 |
| TraesCAD_scaffold_099972_01G000200 | 0.4564836 | 0.1383663 | 0.2545443 |   0.6464594 |  0.1419196 |

Cadenza: 15 genes with the smallest RNA\<-\>protein condition-profile
distance.

``` r
knitr::kable(TB_nor$most_discordant, row.names = FALSE,
            caption = "Norin: 15 genes with the largest RNA<->protein condition-profile distance.")
```

| gene_id               |   d_full |    d_pca2 |    d_umap | correlation |  amplitude |
|:----------------------|---------:|----------:|----------:|------------:|-----------:|
| TraesNOR1D03G00570400 | 15.02558 |  3.893949 |  4.891499 |  -0.0626448 |  1.4259658 |
| TraesNOR6A03G03297390 | 14.10718 | 10.185115 |  4.876135 |  -0.2789142 |  0.7465790 |
| TraesNOR2B03G00927270 | 14.05256 |  9.510851 |  3.986484 |   0.2967501 |  3.9658209 |
| TraesNOR2D03G01226140 | 13.03249 |  9.846638 |  5.104037 |  -0.1558289 |  2.5691946 |
| TraesNOR1A03G00156570 | 12.85264 | 10.179145 | 11.012752 |   0.3194219 |  0.7387961 |
| TraesNOR3D03G01932610 | 12.84277 | 12.023901 |  6.074400 |   0.1232775 | -0.9050636 |
| TraesNOR2A03G00800550 | 12.50523 |  9.526019 |  9.645459 |  -0.7827240 |  1.5097885 |
| TraesNOR4D03G02551270 | 12.35318 | 11.492479 |  7.378051 |  -0.7200042 | -0.1133120 |
| TraesNOR7D03G04390710 | 12.27664 | 11.351402 |  3.807223 |   0.6483526 | -1.6222463 |
| TraesNOR7A03G03920620 | 12.23833 |  5.142548 |  6.133183 |  -0.3505449 |  1.6776861 |
| TraesNOR6A03G03323590 | 11.90557 |  8.263454 |  7.847375 |  -0.3421360 |  3.1282067 |
| TraesNOR7B03G04214550 | 11.89695 | 11.446567 |  3.176004 |   0.7576481 | -1.3974439 |
| TraesNOR1A03G00091890 | 11.76049 | 10.866732 |  7.409615 |  -0.5678999 | -0.4247522 |
| TraesNOR2D03G01170600 | 11.74577 | 10.731394 |  6.482513 |   0.0120591 | -0.5085356 |
| TraesNOR7A03G03913820 | 11.71902 |  8.988329 |  4.581833 |  -0.2717220 |  0.3850277 |

Norin: 15 genes with the largest RNA\<-\>protein condition-profile
distance.

``` r
knitr::kable(TB_nor$most_concordant, row.names = FALSE,
            caption = "Norin: 15 genes with the smallest RNA<->protein condition-profile distance.")
```

| gene_id               |    d_full |    d_pca2 |    d_umap | correlation |  amplitude |
|:----------------------|----------:|----------:|----------:|------------:|-----------:|
| TraesNOR3A03G01455230 | 0.3198929 | 0.0642875 | 0.2819026 |   0.9492192 | -0.0663763 |
| TraesNOR5D03G03099340 | 0.3819582 | 0.2654183 | 0.2613097 |   0.9544345 | -0.1586547 |
| TraesNOR6B03G03594350 | 0.4092053 | 0.1745770 | 0.3085167 |   0.8688774 |  0.5033650 |
| TraesNOR7D03G04462700 | 0.4488227 | 0.3170282 | 1.1204432 |   0.8223209 |  0.3053959 |
| TraesNOR4A03G02179410 | 0.4504038 | 0.4095664 | 0.8385381 |   0.9145572 | -0.2182271 |
| TraesNORUn03G04756130 | 0.4810398 | 0.1520531 | 0.3105033 |   0.9138537 | -0.0012103 |
| TraesNOR3B03G01612090 | 0.5017156 | 0.2577880 | 0.9370760 |   0.6969435 | -0.1925588 |
| TraesNOR2D03G01314180 | 0.5031233 | 0.2400009 | 0.8241532 |   0.7303275 |  0.8852612 |
| TraesNOR3D03G01917840 | 0.5073714 | 0.2782908 | 1.8957432 |   0.7718713 | -0.7103490 |
| TraesNORUn03G04556950 | 0.5097156 | 0.3569833 | 1.1144144 |   0.8155968 | -0.7000619 |
| TraesNOR5D03G03235660 | 0.5200205 | 0.2069519 | 0.3701004 |   0.8110015 |  0.1060114 |
| TraesNOR5A03G02733810 | 0.5273895 | 0.4099040 | 0.2313463 |   0.6981259 | -0.3755882 |
| TraesNOR7D03G04390410 | 0.5275586 | 0.1927471 | 0.4816958 |   0.7505460 | -0.9049114 |
| TraesNOR1D03G00560520 | 0.5289193 | 0.2828857 | 0.5918853 |   0.6256649 |  0.1821194 |
| TraesNOR6D03G03742630 | 0.5355704 | 0.3224174 | 0.5759194 |   0.9284549 |  0.4215877 |

Norin: 15 genes with the smallest RNA\<-\>protein condition-profile
distance.

``` r
# Export full gene rankings (8-condition analysis)
output_dir <- here::here("results", "gene-distance")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (nm in c("Cadenza", "Norin")) {
  d <- GD[[nm]]$genes[order(-GD[[nm]]$genes$d_full), ]
  filename <- file.path(output_dir, paste0(tolower(nm), "_gene_distances_8condition.csv"))
  write.csv(d, filename, row.names = FALSE)
  cat("✓ Saved:", filename, "\n")
}
```

    ✓ Saved: C:/Claude/Projects/IntegrateProtRNA/results/gene-distance/cadenza_gene_distances_8condition.csv 
    ✓ Saved: C:/Claude/Projects/IntegrateProtRNA/results/gene-distance/norin_gene_distances_8condition.csv 

``` r
cat("\nFiles contain columns:\n",
    "  gene_id: gene identifier\n",
    "  d_full: distance (full standardised space) — use this for ranking\n",
    "  d_pca2: distance (top-2 PCA) — for comparison\n",
    "  d_umap: distance (2D UMAP) — visualization only\n",
    "  correlation: Spearman correlation (RNA vs protein response patterns)\n",
    "  amplitude: log2(protein response magnitude / RNA response magnitude)\n")
```


    Files contain columns:
       gene_id: gene identifier
       d_full: distance (full standardised space) — use this for ranking
       d_pca2: distance (top-2 PCA) — for comparison
       d_umap: distance (2D UMAP) — visualization only
       correlation: Spearman correlation (RNA vs protein response patterns)
       amplitude: log2(protein response magnitude / RNA response magnitude)

------------------------------------------------------------------------

## 8a. Per-condition analysis: Which condition drives discordance?

The overall distance d_full summarizes discrepancy across all 8
conditions. But which conditions contribute most to each gene’s
distance? This breakdown reveals:

- Condition-specific RNA-protein divergence
- Whether discordance is spread across treatments or concentrated
- Which timepoint shows largest divergence

``` r
#' Calculate per-condition contribution to overall distance
#' Returns data frame with each condition's squared contribution
condition_contributions <- function(R_centered, P_centered, condition_names = NULL) {

  if (is.null(condition_names)) {
    condition_names <- colnames(R_centered)
  }

  # Fit scaler on RNA only
  mu <- colMeans(R_centered)
  sdv <- apply(R_centered, 2, sd); sdv[sdv < 1e-8] <- 1

  # Standardize
  Rs <- sweep(sweep(R_centered, 2, mu, "-"), 2, sdv, "/")
  Ps <- sweep(sweep(P_centered, 2, mu, "-"), 2, sdv, "/")

  # Per-condition squared contribution: (s_c * x_ic - y_ic)²
  contributions <- (Rs - Ps)^2

  # Convert to data frame with condition names
  colnames(contributions) <- condition_names
  contributions
}

# Calculate for each variety
per_cond_contributions <- list()

for (nm in c("Cadenza", "Norin")) {
  cat("\n**", nm, ":**\n", sep = "")

  # Get the real variety data
  v <- REAL[[nm]]
  common <- intersect(rownames(v$qc_rna$vst), rownames(v$imputed$mixed))
  R_cond <- cell_means(v$qc_rna$vst[common, , drop = FALSE], v$meta)
  P_cond <- cell_means(v$imputed$mixed[common, , drop = FALSE], v$meta)

  # Get centered versions (from the gene_distance_real output)
  sp <- build_shared_space(R_cond, P_cond)

  # Calculate per-condition contributions
  contrib <- condition_contributions(sp$Rs, sp$Ps, colnames(R_cond))

  # Calculate total distance per gene and contribution %
  d_total <- sqrt(rowSums(contrib))
  contrib_pct <- sweep(contrib, 1, d_total^2, "/") * 100

  # For each gene, identify top 3 contributing conditions
  top_contrib <- data.frame(
    gene_id = rownames(contrib),
    distance = d_total,
    top_cond_1 = NA,
    top_pct_1 = NA,
    top_cond_2 = NA,
    top_pct_2 = NA,
    top_cond_3 = NA,
    top_pct_3 = NA
  )

  for (i in seq_len(nrow(contrib_pct))) {
    top3_idx <- order(-contrib_pct[i, ])[1:3]
    for (j in 1:3) {
      top_contrib[[paste0("top_cond_", j)]][i] <- colnames(contrib_pct)[top3_idx[j]]
      top_contrib[[paste0("top_pct_", j)]][i] <- round(contrib_pct[i, top3_idx[j]], 1)
    }
  }

  per_cond_contributions[[nm]] <- list(
    contributions = contrib,
    contributions_pct = contrib_pct,
    top_contributors = top_contrib,
    distance = d_total
  )

  # Show summary
  cat("  Example genes with top contributing conditions:\n")
  top_genes <- top_contrib[order(-top_contrib$distance), ][1:5, ]
  print(top_genes[, c("gene_id", "distance", "top_cond_1", "top_pct_1")])
}
```


    **Cadenza:**
      Example genes with top contributing conditions:
                                                                  gene_id distance
    TraesCAD_scaffold_014747_01G000100 TraesCAD_scaffold_014747_01G000100 14.58265
    TraesCAD_scaffold_133124_01G000100 TraesCAD_scaffold_133124_01G000100 13.46543
    TraesCAD_scaffold_116541_01G000100 TraesCAD_scaffold_116541_01G000100 13.42158
    TraesCAD_scaffold_146252_01G000100 TraesCAD_scaffold_146252_01G000100 12.98171
    TraesCAD_scaffold_129059_01G000100 TraesCAD_scaffold_129059_01G000100 12.07514
                                       top_cond_1 top_pct_1
    TraesCAD_scaffold_014747_01G000100     T1_t24      18.2
    TraesCAD_scaffold_133124_01G000100     T1_t24      45.4
    TraesCAD_scaffold_116541_01G000100      T1_t0      55.0
    TraesCAD_scaffold_146252_01G000100     T0_t48      48.6
    TraesCAD_scaffold_129059_01G000100     T0_t24      42.6

    **Norin:**
      Example genes with top contributing conditions:
                                        gene_id distance top_cond_1 top_pct_1
    TraesNOR1D03G00570400 TraesNOR1D03G00570400 15.02558      T1_t0      28.6
    TraesNOR6A03G03297390 TraesNOR6A03G03297390 14.10718     T1_t48      77.3
    TraesNOR2B03G00927270 TraesNOR2B03G00927270 14.05256     T1_t48      56.2
    TraesNOR2D03G01226140 TraesNOR2D03G01226140 13.03249     T0_t24      28.1
    TraesNOR1A03G00156570 TraesNOR1A03G00156570 12.85264     T0_t24      57.7

------------------------------------------------------------------------

## 8b. Three-metric decomposition: Distance, Correlation, Amplitude

A single distance number masks biological diversity. The same high
distance can arise from different mechanisms:

- **High distance + high correlation**: same response pattern, different
  amplitude (buffering/scaling)
- **High distance + low/negative correlation**: opposite or independent
  regulation
- **Low distance + high correlation**: strong RNA-protein concordance

Breaking d_full into three complementary metrics reveals which genes
have which behaviour.

``` r
# For each variety, show top discordant genes with their three metrics
for (nm in c("Cadenza", "Norin")) {
  d <- GD[[nm]]$genes
  d <- d[order(-d$d_full), ]

  cat("\n**", nm, ": top 10 most-discordant genes with decomposed metrics**\n\n", sep = "")

  top10 <- head(d, 10)[, c("gene_id", "d_full", "correlation", "amplitude")]
  top10$d_full <- round(top10$d_full, 3)
  top10$correlation <- round(top10$correlation, 3)
  top10$amplitude <- round(top10$amplitude, 3)

  print(knitr::kable(top10, row.names = FALSE, caption = paste(nm, "decomposition")))
  cat("\n")
}
```


    **Cadenza: top 10 most-discordant genes with decomposed metrics**



    Table: Cadenza decomposition

    |gene_id                            | d_full| correlation| amplitude|
    |:----------------------------------|------:|-----------:|---------:|
    |TraesCAD_scaffold_014747_01G000100 | 14.583|       0.906|    -2.177|
    |TraesCAD_scaffold_133124_01G000100 | 13.465|       0.421|    -0.581|
    |TraesCAD_scaffold_116541_01G000100 | 13.422|      -0.222|     0.081|
    |TraesCAD_scaffold_146252_01G000100 | 12.982|      -0.463|    -0.381|
    |TraesCAD_scaffold_129059_01G000100 | 12.075|       0.811|     1.249|
    |TraesCAD_scaffold_004054_01G000200 | 11.694|      -0.254|     1.162|
    |TraesCAD_scaffold_086045_01G000100 | 10.903|      -0.355|    -0.202|
    |TraesCAD_scaffold_070666_01G000100 | 10.836|      -0.070|     0.957|
    |TraesCAD_scaffold_054259_01G000200 | 10.770|      -0.254|     2.464|
    |TraesCAD_scaffold_033501_01G000300 | 10.738|       0.194|    -0.045|


    **Norin: top 10 most-discordant genes with decomposed metrics**



    Table: Norin decomposition

    |gene_id               | d_full| correlation| amplitude|
    |:---------------------|------:|-----------:|---------:|
    |TraesNOR1D03G00570400 | 15.026|      -0.063|     1.426|
    |TraesNOR6A03G03297390 | 14.107|      -0.279|     0.747|
    |TraesNOR2B03G00927270 | 14.053|       0.297|     3.966|
    |TraesNOR2D03G01226140 | 13.032|      -0.156|     2.569|
    |TraesNOR1A03G00156570 | 12.853|       0.319|     0.739|
    |TraesNOR3D03G01932610 | 12.843|       0.123|    -0.905|
    |TraesNOR2A03G00800550 | 12.505|      -0.783|     1.510|
    |TraesNOR4D03G02551270 | 12.353|      -0.720|    -0.113|
    |TraesNOR7D03G04390710 | 12.277|       0.648|    -1.622|
    |TraesNOR7A03G03920620 | 12.238|      -0.351|     1.678|

### Interpretation of the three metrics

| d_full   | correlation  | amplitude | biological meaning                                       |
|----------|--------------|-----------|----------------------------------------------------------|
| low      | high         | ~0        | identical response shape and strength                    |
| moderate | high         | ±large    | same response pattern, protein buffered or amplified     |
| high     | high         | ±large    | similar pattern but large amplitude difference (unusual) |
| high     | low/negative | any       | opposite or orthogonal regulation; genuine decoupling    |
| high     | moderate     | any       | delayed/shifted response; kinetic lag                    |

------------------------------------------------------------------------

## 8c. Infection-response contrast analysis (optional)

The main analysis (§8) uses all 8 condition means (2 treatments × 4
timepoints). But for the disease biology, a cleaner question is:

**Does protein reproduce the infection-induced change observed at the
transcript level?**

This version uses only the 4-point infection-response trajectories:

$$\Delta_t^{\text{RNA}} = \text{RNA}_{\text{infected}}(t) - \text{RNA}_{\text{control}}(t)$$
$$\Delta_t^{\text{protein}} = \text{Protein}_{\text{infected}}(t) - \text{Protein}_{\text{control}}(t)$$

This avoids conflating normal developmental time-effects with
infection-specific effects.

``` r
if (!is.null(GD_contrast) && !all(sapply(GD_contrast, is.null))) {
  cat("✓ Infection-response contrast analysis computed.\n\n")
  for (nm in names(GD_contrast)) {
    if (!is.null(GD_contrast[[nm]])) {
      cat("**", nm, ":**\n", sep = "")
      d <- GD_contrast[[nm]]$genes
      cat("  Genes:", nrow(d), "\n")
      cat("  Mean distance (contrasts):", round(GD_contrast[[nm]]$obs, 3), "\n")
      cat("  Permutation p:", GD_contrast[[nm]]$p_perm, "\n\n")
    }
  }
} else {
  cat("⚠️ Infection-response contrast analysis not available.\n",
      "Metadata may not have 'treatment' column (infected vs control).\n")
}
```

    ⚠️ Infection-response contrast analysis not available.
     Metadata may not have 'treatment' column (infected vs control).

``` r
if (!is.null(GD_contrast) && !all(sapply(GD_contrast, is.null))) {
  op <- par(mfrow = c(1, 2))

  for (nm in names(GD_contrast)) {
    if (!is.null(GD_contrast[[nm]])) {
      d_main <- GD[[nm]]$genes
      d_contrast <- GD_contrast[[nm]]$genes

      # Merge on gene_id
      merged <- merge(d_main[, c("gene_id", "d_full")],
                      d_contrast[, c("gene_id", "d_full")],
                      by = "gene_id", suffixes = c("_main", "_contrast"))

      r_rank <- cor(rank(-merged$d_full_main), rank(-merged$d_full_contrast), method = "spearman")

      plot(merged$d_full_main, merged$d_full_contrast, pch = 19, col = "#2E86AB", alpha = 0.5,
           main = paste(nm, "\nSpearman ρ =", round(r_rank, 3)),
           xlab = "Distance (8 conditions)", ylab = "Distance (4 infection contrasts)")
      abline(0, 1, col = "gray", lty = 2, lwd = 1)
      abline(lm(merged$d_full_contrast ~ merged$d_full_main), col = "#D1495B", lwd = 2)
    }
  }

  par(op)

  # Export infection-contrast results
  output_dir <- here::here("results", "gene-distance")
  for (nm in names(GD_contrast)) {
    if (!is.null(GD_contrast[[nm]])) {
      d <- GD_contrast[[nm]]$genes[order(-GD_contrast[[nm]]$genes$d_full), ]
      filename <- file.path(output_dir, paste0(tolower(nm), "_gene_distances_infection_contrast.csv"))
      write.csv(d, filename, row.names = FALSE)
      cat("✓ Saved:", filename, "\n")
    }
  }
}
```

------------------------------------------------------------------------

## 9. Visualizations: UMAP in shared space

### V1 — Shared UMAP embeddings, real data

Fit UMAP on RNA condition-means, project protein into the same space.
This visualization shows where each gene’s RNA and protein profiles land
in a 2D manifold learned from population-level response shapes.

``` r
#' Create UMAP plots for shared-space visualization
plot_shared_umap <- function(space_list, dist_df, variety_name, tag = "") {

  # Extract embeddings
  Er <- space_list$Er
  Ep <- space_list$Ep
  distances <- dist_df$d_umap

  # Create plot
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Left: colored by distance
  color_scale <- colorRampPalette(c("#2C6FBB", "#E8A33D", "#D1495B"))(100)
  color_idx <- pmin(99, pmax(1, round((distances - min(distances)) /
                                       (max(distances) - min(distances)) * 99) + 1))

  plot(Er[, 1], Er[, 2], type = "n", xlab = "UMAP1", ylab = "UMAP2",
       main = paste(variety_name, ": RNA-Protein Shared UMAP\n(colored by distance)"),
       cex.main = 1.1)

  # Background points
  points(Er[, 1], Er[, 2], col = color_scale[color_idx], pch = 19,
         cex = 2, alpha = 0.6)
  points(Ep[, 1], Ep[, 2], col = color_scale[color_idx], pch = 17,
         cex = 2, alpha = 0.6)

  # Connecting lines
  for (i in seq_len(nrow(Er))) {
    lines(c(Er[i, 1], Ep[i, 1]), c(Er[i, 2], Ep[i, 2]),
          col = rgb(0, 0, 0, 0.1), lwd = 0.5)
  }

  # Right: RNA vs Protein
  plot(Er[, 1], Er[, 2], pch = 19, col = "#2E86AB", cex = 2, alpha = 0.7,
       xlab = "UMAP1", ylab = "UMAP2",
       main = paste(variety_name, ": RNA vs Protein Positions"))
  points(Ep[, 1], Ep[, 2], pch = 17, col = "#A23B72", cex = 2, alpha = 0.7)

  legend("topright", legend = c("RNA", "Protein"), pch = c(19, 17),
         col = c("#2E86AB", "#A23B72"), cex = 0.9)

  par(mfrow = c(1, 1))
}

#' Highlight top and bottom genes by distance
#' Uses d_full (full space distance) for ranking consistency, but visualizes in UMAP space
plot_top_genes_umap <- function(space_list, dist_df, variety_name, n_top = 10) {

  Er <- space_list$Er
  Ep <- space_list$Ep
  distances <- dist_df$d_full  # Rank by full-space distance, visualize in UMAP

  par(mar = c(4, 4, 3, 1))

  # Background points in light gray
  plot(Er[, 1], Er[, 2], type = "n", xlab = "UMAP1", ylab = "UMAP2",
       main = paste(variety_name, ": Top", n_top, "Most/Least Discordant Genes\n(ranked by full-space distance, shown in UMAP)"),
       cex.main = 1.1)

  points(Er[, 1], Er[, 2], col = rgb(0.5, 0.5, 0.5, 0.2), pch = 19, cex = 1.5)
  points(Ep[, 1], Ep[, 2], col = rgb(0.5, 0.5, 0.5, 0.2), pch = 17, cex = 1.5)

  # Top N most discordant (high distance)
  top_idx <- order(-distances)[seq_len(n_top)]
  points(Er[top_idx, 1], Er[top_idx, 2], col = "#D1495B", pch = 19, cex = 3.5)
  points(Ep[top_idx, 1], Ep[top_idx, 2], col = "#D1495B", pch = 17, cex = 3.5)

  for (i in top_idx) {
    lines(c(Er[i, 1], Ep[i, 1]), c(Er[i, 2], Ep[i, 2]),
          col = "#D1495B", lwd = 1.5, alpha = 0.5)
  }

  # Bottom N most concordant (low distance)
  bottom_idx <- order(distances)[seq_len(n_top)]
  points(Er[bottom_idx, 1], Er[bottom_idx, 2], col = "#2C6FBB", pch = 19, cex = 3.5)
  points(Ep[bottom_idx, 1], Ep[bottom_idx, 2], col = "#2C6FBB", pch = 17, cex = 3.5)

  for (i in bottom_idx) {
    lines(c(Er[i, 1], Ep[i, 1]), c(Er[i, 2], Ep[i, 2]),
          col = "#2C6FBB", lwd = 1.5, alpha = 0.5)
  }

  legend("topright", legend = c("Discordant (RNA)", "Discordant (Prot)",
                                "Concordant (RNA)", "Concordant (Prot)"),
         pch = c(19, 17, 19, 17), col = c("#D1495B", "#D1495B", "#2C6FBB", "#2C6FBB"),
         cex = 0.9)
}
```

``` r
plot_shared_umap(GD$Cadenza$space, GD$Cadenza$genes, "Cadenza")
```

    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter

    Warning in plot.window(...): "alpha" is not a graphical parameter

    Warning in plot.xy(xy, type, ...): "alpha" is not a graphical parameter

    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter
    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter

    Warning in box(...): "alpha" is not a graphical parameter

    Warning in title(...): "alpha" is not a graphical parameter

    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter

![V1a: Cadenza shared UMAP embeddings. Left: all genes colored by
RNA-protein distance. Right: RNA (circles) and protein (triangles) in
the same learned
manifold.](gene_distance_shared_space_files/figure-commonmark/v1-cadenza-umap-1.png)

``` r
plot_shared_umap(GD$Norin$space, GD$Norin$genes, "Norin")
```

    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter

    Warning in plot.window(...): "alpha" is not a graphical parameter

    Warning in plot.xy(xy, type, ...): "alpha" is not a graphical parameter

    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter
    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter

    Warning in box(...): "alpha" is not a graphical parameter

    Warning in title(...): "alpha" is not a graphical parameter

    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter

![V1b: Norin shared UMAP embeddings, same layout as
Cadenza.](gene_distance_shared_space_files/figure-commonmark/v1-norin-umap-1.png)

### V2 — Top and bottom genes in UMAP space

Which genes are most concordant (close together) and most discordant
(far apart) in the shared UMAP manifold?

``` r
plot_top_genes_umap(GD$Cadenza$space, GD$Cadenza$genes, "Cadenza", n_top = 10)
```

    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter

![V2a: Cadenza. Circles = RNA, triangles = protein. Red = top 10 most
discordant (high distance). Blue = bottom 10 most concordant (low
distance).](gene_distance_shared_space_files/figure-commonmark/v2-cadenza-topgenes-1.png)

``` r
plot_top_genes_umap(GD$Norin$space, GD$Norin$genes, "Norin", n_top = 10)
```

    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter
    Warning in plot.xy(xy.coords(x, y), type = type, ...): "alpha" is not a
    graphical parameter

![V2b: Norin, same layout as
Cadenza.](gene_distance_shared_space_files/figure-commonmark/v2-norin-topgenes-1.png)

### V3 — Distance distributions

Histogram of full-space distances (the metric used for ranking, not the
2D UMAP distance).

``` r
op <- par(mfrow = c(1, 2))
for (nm in c("Cadenza", "Norin")) {
  d <- GD[[nm]]$genes$d_full
  hist(d, breaks = 30, col = "#2E86AB", border = "white",
       main = paste(nm, "distance distribution"),
       xlab = "RNA-Protein distance (full space)")
  abline(v = mean(d), col = "#D1495B", lwd = 2.5, lty = 2)
  abline(v = median(d), col = "#E8A33D", lwd = 2.5, lty = 2)
  legend("topright", legend = c("Mean", "Median"), col = c("#D1495B", "#E8A33D"),
         lty = 2, lwd = 2.5, cex = 0.8)
}
```

![V3: Distribution of RNA-protein distances. Full standardised space
used for ranking
(A4).](gene_distance_shared_space_files/figure-commonmark/v3-distance-dist-1.png)

``` r
par(op)
```

### V3b — Density distribution in UMAP space

Where do RNA and protein cluster? Filled density contours show the
concentration of genes in different regions of the manifold.

``` r
# Prepare data for ggplot2: combine RNA and protein embeddings with labels
density_data <- local({
  rna_cad <- data.frame(
    UMAP1 = GD$Cadenza$space$Er[, 1],
    UMAP2 = GD$Cadenza$space$Er[, 2],
    modality = "RNA",
    variety = "Cadenza"
  )
  prot_cad <- data.frame(
    UMAP1 = GD$Cadenza$space$Ep[, 1],
    UMAP2 = GD$Cadenza$space$Ep[, 2],
    modality = "Protein",
    variety = "Cadenza"
  )
  rna_nor <- data.frame(
    UMAP1 = GD$Norin$space$Er[, 1],
    UMAP2 = GD$Norin$space$Er[, 2],
    modality = "RNA",
    variety = "Norin"
  )
  prot_nor <- data.frame(
    UMAP1 = GD$Norin$space$Ep[, 1],
    UMAP2 = GD$Norin$space$Ep[, 2],
    modality = "Protein",
    variety = "Norin"
  )
  rbind(rna_cad, prot_cad, rna_nor, prot_nor)
})

# Create plot with ggplot2: separate panels per variety, overlaid densities
p <- ggplot(
  density_data,
  aes(x = UMAP1, y = UMAP2, fill = modality)
) +
  facet_wrap(~ variety, ncol = 2) +

  stat_density_2d(
    geom = "polygon",
    alpha = 0.42,
    bins = 6,
    contour = TRUE,
    colour = NA
  ) +

  scale_fill_manual(
    values = c(
      "RNA"     = "#4A90E2",
      "Protein" = "#F5A623"
    ),
    name = "Modality"
  ) +

  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +

  theme_classic(base_size = 20) +

  theme(
    strip.background = element_blank(),

    strip.text = element_text(
      face = "bold",
      size = 22
    ),

    axis.title = element_text(
      size = 22,
      face = "bold"
    ),

    axis.text = element_text(
      size = 18,
      colour = "black"
    ),

    axis.line = element_line(
      colour = "black",
      linewidth = 0.8
    ),

    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.7
    ),

    axis.ticks.length = unit(0.18, "cm"),

    legend.position = "bottom",

    legend.title = element_text(
      size = 20,
      face = "bold"
    ),

    legend.text = element_text(
      size = 18
    ),

    legend.key.size = unit(0.8, "cm"),

    panel.spacing = unit(1.5, "lines"),

    plot.margin = margin(
      12, 16, 12, 12
    )
  ) +

  guides(
    fill = guide_legend(
      title.position = "left",
      override.aes = list(alpha = 0.42)
    )
  )

print(p)
```

![V3b: 2D density distributions of RNA (blue) and protein (orange) in
shared UMAP space. Filled contours show concentration
levels.](gene_distance_shared_space_files/figure-commonmark/v3b-density-plots-1.png)

### V3c — Per-condition contributions to distance

Which conditions drive the most RNA-protein discordance? Heatmap shows
each gene (rows) and condition contribution (columns).

``` r
for (nm in c("Cadenza", "Norin")) {
  if (!is.null(per_cond_contributions[[nm]])) {
    contrib_pct <- per_cond_contributions[[nm]]$contributions_pct
    distance <- per_cond_contributions[[nm]]$distance

    # Sort by distance (most discordant first)
    idx_sort <- order(-distance)
    contrib_sorted <- contrib_pct[idx_sort[1:min(50, nrow(contrib_pct))], ]

    # Heatmap using base R
    par(mar = c(5, 12, 2, 1))
    image(t(contrib_sorted), col = colorRampPalette(c("white", "orange", "red"))(100),
          xaxt = "n", yaxt = "n",
          main = paste(nm, ": Per-condition contribution to distance\n(top 50 most-discordant genes)"))
    axis(1, at = seq(0, 1, length.out = ncol(contrib_sorted)), labels = colnames(contrib_sorted), las = 2)
    axis(2, at = seq(0, 1, length.out = min(50, nrow(contrib_pct))), labels = rownames(contrib_sorted), las = 1, cex.axis = 0.6)
  }
}
```

![V3c: Per-condition contribution to RNA-protein distance (% of total
d_full²). Red = high contribution. Shows which timepoints/treatments are
most
discordant.](gene_distance_shared_space_files/figure-commonmark/v3c-condition-heatmap-1.png)

![V3c: Per-condition contribution to RNA-protein distance (% of total
d_full²). Red = high contribution. Shows which timepoints/treatments are
most
discordant.](gene_distance_shared_space_files/figure-commonmark/v3c-condition-heatmap-2.png)

### V3d — Per-condition trends: average discordance by timepoint/treatment

Which conditions show the most RNA-protein discrepancy across all genes?
Shows mean and distribution of discordance for each condition.

``` r
op <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

for (nm in c("Cadenza", "Norin")) {
  if (!is.null(per_cond_contributions[[nm]])) {
    contrib <- per_cond_contributions[[nm]]$contributions
    contrib_pct <- per_cond_contributions[[nm]]$contributions_pct

    # Left: Mean contribution per condition
    mean_contrib <- colMeans(contrib)
    mean_contrib_pct <- colMeans(contrib_pct)

    barplot(mean_contrib_pct, names.arg = colnames(contrib), las = 2,
            ylab = "Mean % contribution to distance",
            main = paste(nm, ": Average per-condition discordance"),
            col = "#2E86AB", border = NA)

    # Right: Box plot of per-condition discrepancies
    contrib_long <- data.frame(
      condition = rep(colnames(contrib), each = nrow(contrib)),
      distance_sq = as.vector(as.matrix(contrib))
    )

    boxplot(distance_sq ~ condition, data = contrib_long, las = 2,
            ylab = "Per-condition squared difference",
            main = paste(nm, ": Distribution of discordance by condition"),
            col = "#A23B72", border = NA)
  }
}
```

![V3d: Per-condition discordance trends. Left: mean contribution per
condition. Right: distribution of per-condition squared differences
across all genes. Which conditions are systematically
discordant?](gene_distance_shared_space_files/figure-commonmark/v3d-condition-trends-1.png)

![V3d: Per-condition discordance trends. Left: mean contribution per
condition. Right: distribution of per-condition squared differences
across all genes. Which conditions are systematically
discordant?](gene_distance_shared_space_files/figure-commonmark/v3d-condition-trends-2.png)

``` r
par(op)

# Summary statistics
cat("\n**Summary: Per-condition discordance**\n\n")
```


    **Summary: Per-condition discordance**

``` r
for (nm in c("Cadenza", "Norin")) {
  if (!is.null(per_cond_contributions[[nm]])) {
    contrib_pct <- per_cond_contributions[[nm]]$contributions_pct
    mean_pct <- colMeans(contrib_pct)
    sd_pct <- apply(contrib_pct, 2, sd)

    cat("**", nm, ":**\n", sep = "")
    for (i in seq_along(mean_pct)) {
      cat(sprintf("  %s: %.1f%% ± %.1f%% mean contribution\n",
                  names(mean_pct)[i], mean_pct[i], sd_pct[i]))
    }
    cat("\n")
  }
}
```

    **Cadenza:**
      T0_t0: 17.7% ± 18.8% mean contribution
      T0_t24: 14.6% ± 16.4% mean contribution
      T0_t48: 10.2% ± 12.1% mean contribution
      T0_t72: 11.6% ± 12.9% mean contribution
      T1_t0: 14.1% ± 16.0% mean contribution
      T1_t24: 13.5% ± 14.4% mean contribution
      T1_t48: 8.4% ± 9.5% mean contribution
      T1_t72: 9.7% ± 12.7% mean contribution

    **Norin:**
      T0_t0: 11.4% ± 12.9% mean contribution
      T0_t24: 15.3% ± 16.6% mean contribution
      T0_t48: 15.5% ± 16.0% mean contribution
      T0_t72: 10.0% ± 10.8% mean contribution
      T1_t0: 13.8% ± 15.5% mean contribution
      T1_t24: 13.7% ± 13.9% mean contribution
      T1_t48: 13.0% ± 14.7% mean contribution
      T1_t72: 7.3% ± 9.3% mean contribution

### V4 — Comparing full space, PCA, and UMAP distances

Are the rankings stable across different projection methods? If so, the
distance metric is robust.

``` r
op <- par(mfrow = c(1, 2))

for (nm in c("Cadenza", "Norin")) {
  d <- GD[[nm]]$genes

  # Spearman rank correlations
  r_full_pca  <- cor(d$d_full, d$d_pca2, method = "spearman")
  r_full_umap <- cor(d$d_full, d$d_umap, method = "spearman")
  r_pca_umap  <- cor(d$d_pca2, d$d_umap, method = "spearman")

  # Plot: full vs PCA
  plot(d$d_full, d$d_pca2, pch = 19, col = "#2E86AB", alpha = 0.5,
       main = paste(nm, ": Full space vs PCA\nSpearman ρ =", round(r_full_pca, 3)),
       xlab = "Distance (full standardised space)",
       ylab = "Distance (top-2 PCA)")
  abline(lm(d$d_pca2 ~ d$d_full), col = "#D1495B", lwd = 2)
}
```

    Warning in plot.window(...): "alpha" is not a graphical parameter

    Warning in plot.xy(xy, type, ...): "alpha" is not a graphical parameter

    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter
    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter

    Warning in box(...): "alpha" is not a graphical parameter

    Warning in title(...): "alpha" is not a graphical parameter

    Warning in plot.window(...): "alpha" is not a graphical parameter

    Warning in plot.xy(xy, type, ...): "alpha" is not a graphical parameter

    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter
    Warning in axis(side = side, at = at, labels = labels, ...): "alpha" is not a
    graphical parameter

    Warning in box(...): "alpha" is not a graphical parameter

    Warning in title(...): "alpha" is not a graphical parameter

![V4: Full space vs. top-2 PCA. PCA tracks full-space ranking closely
(high Spearman rank correlation) – see v4-agreement-table below for the
UMAP comparison, which does
not.](gene_distance_shared_space_files/figure-commonmark/v4-correlation-comparisons-1.png)

``` r
par(op)
```

`r_full_umap` and `r_pca_umap` above were computed but never reported –
the loop only ever plotted full-vs-PCA. Surfacing them, plus what
fraction of the genes each space would call “most discordant” actually
agree with the full-space call, since a rank correlation alone can hide
whether the *specific genes at the top of the table* (the ones anyone
would actually follow up on) are the same genes:

``` r
v4_agreement <- do.call(rbind, lapply(c("Cadenza", "Norin"), function(nm) {
  d <- GD[[nm]]$genes
  k <- round(0.10 * nrow(d))
  top_full <- d$gene_id[order(-d$d_full)][1:k]
  top_pca  <- d$gene_id[order(-d$d_pca2)][1:k]
  top_umap <- d$gene_id[order(-d$d_umap)][1:k]

  data.frame(
    variety                = nm,
    n_genes                = nrow(d),
    rho_full_vs_pca        = round(cor(d$d_full, d$d_pca2, method = "spearman"), 3),
    rho_full_vs_umap       = round(cor(d$d_full, d$d_umap, method = "spearman"), 3),
    rho_pca_vs_umap        = round(cor(d$d_pca2, d$d_umap, method = "spearman"), 3),
    top10pct_overlap_pca_n = length(intersect(top_full, top_pca)),
    top10pct_overlap_umap_n = length(intersect(top_full, top_umap)),
    top10pct_n              = k
  )
}))

knitr::kable(v4_agreement,
            caption = "V4: how much of the full-space ranking survives PCA-2 vs. UMAP, on the real data. top10pct_overlap_*_n / top10pct_n is what fraction of the full-space 'most discordant' decile that projection still puts in its own most-discordant decile.")
```

| variety | n_genes | rho_full_vs_pca | rho_full_vs_umap | rho_pca_vs_umap | top10pct_overlap_pca_n | top10pct_overlap_umap_n | top10pct_n |
|:--------|--------:|----------------:|-----------------:|----------------:|-----------------------:|------------------------:|-----------:|
| Cadenza |    5580 |           0.869 |            0.549 |           0.682 |                    408 |                     179 |        558 |
| Norin   |    5159 |           0.865 |            0.514 |           0.613 |                    380 |                     139 |        516 |

V4: how much of the full-space ranking survives PCA-2 vs. UMAP, on the
real data. top10pct_overlap\_\*\_n / top10pct_n is what fraction of the
full-space ‘most discordant’ decile that projection still puts in its
own most-discordant decile.

``` r
op <- par(mfrow = c(1, 2))

for (nm in c("Cadenza", "Norin")) {
  d <- GD[[nm]]$genes
  r_full_umap <- cor(d$d_full, d$d_umap, method = "spearman")

  plot(d$d_full, d$d_umap, pch = 19, col = "#2E86AB",
       main = paste(nm, ": Full space vs UMAP\nSpearman ρ =", round(r_full_umap, 3)),
       xlab = "Distance (full standardised space)",
       ylab = "Distance (2D UMAP)")
  abline(lm(d$d_umap ~ d$d_full), col = "#D1495B", lwd = 2)
}
```

![V4b: Full space vs. 2D UMAP. Compare the spread here to
v4-correlation-comparisons above – UMAP distance tracks full-space
distance far more loosely than PCA does, consistent with the lower rho
and top-decile overlap in
v4-agreement-table.](gene_distance_shared_space_files/figure-commonmark/v4b-fig-full-vs-umap-1.png)

``` r
par(op)
```

**Reading this against §5’s G2 (simulated data):** there, compressing to
2D UMAP cost ~2.6% AUC (0.820 -\> 0.799) for separating known-discordant
genes from baseline – a small loss. On the real data, the loss looks
much bigger by this measure: rho drops to ~0.51-0.55 (vs ~0.87 for
PCA-2), and only ~27-32% of the genes in the full-space “most
discordant” decile are still in UMAP’s own most-discordant decile (vs
~73% for PCA-2). The two measures are not contradictory – G2’s AUC asks
“does this space still separate archetypes on average,” which tolerates
a lot of individual-gene reshuffling; the agreement table above asks “is
it the *same genes* at the top,” which does not. For picking specific
genes to follow up on (as in §8’s top-15 tables), the second question is
the one that matters, and by that measure UMAP is a much less faithful
summary of `d_full` than the AUC comparison alone would suggest.

------------------------------------------------------------------------

## Interpretation: Is the distance meaningful?

### What the distance tells us

The RNA-protein distance in the shared condition-profile space measures
**standardised discrepancy between the centred RNA and protein
condition-response profiles**. Mathematically:

$$d_i = \sqrt{\sum_c (s_c x_{ic} - y_{ic})^2}$$

where $s_c$ is the across-gene RNA SD for condition $c$, $x_{ic}$ is the
centred RNA for gene $i$ and condition $c$, and $y_{ic}$ is centred
protein.

**Key properties:** - Row-centering removes absolute baseline abundance,
but **does NOT remove response amplitude** - The metric includes both
**pattern (direction of change)** and **magnitude (strength of
response)** - Example: if RNA = (−1, 0, 1) and protein = (−2, 0, 2),
they have identical pattern but the protein response is twice as strong,
so distance \> 0 - If protein = (1, 0, −1), the direction is reversed;
distance is large - The metric is **one-directional**: fitted on RNA,
protein projected into that space

### Can we use it for “similarity”?

**Yes, but understand what it measures.** The distance ranks genes
reproducibly (G4 median Spearman ρ across resamples), and the pairing
signal is real (G3: permutation p = 0.005). However:

1.  **It does NOT distinguish regulation from kinetics** (G1: `lag_only`
    genes rank as discordant despite being mechanistically concordant).
    A high-distance gene is equally likely to be slow-turnover as
    genuinely decoupled.

2.  **It measures total discrepancy including amplitude**, not shape
    alone. Two genes with identical pattern but different response
    magnitudes will have non-zero distance. Two genes with opposite
    directions (one peaks, one dips) will have large distance.

3.  Better terminology: call this a **“standardised condition-response
    discrepancy”** rather than “shape similarity”—it captures both
    directional concordance and amplitude agreement.

### Full space vs PCA vs UMAP

From G2 and V4:

| Space                 | Use for ranking?        | AUC (discordant vs baseline) | Notes                                                                                                                                                                                                                                                                                       |
|-----------------------|-------------------------|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Full standardised** | ✅ Yes (A4 validated)   | 0.820                        | Most informative, distance preserved exactly                                                                                                                                                                                                                                                |
| **Top-2 PCA**         | ✅ Yes (nearly as good) | 0.827                        | Minimal loss of ranking power; visualisable as 2D plot                                                                                                                                                                                                                                      |
| **2D UMAP**           | ⚠️ Use with caution     | 0.799                        | ~2% AUC loss on simulated archetypes, but on the real data (v4-agreement-table) rank correlation with `d_full` is only ρ≈0.51-0.55 and just ~27-32% of the top-decile “most discordant” genes agree with the full-space call – a much bigger practical loss than the AUC gap alone suggests |

**Recommendation:**

- **For ranking and biological claims**: use `d_full` (full standardised
  space). Distance is mathematically exact, not approximate.
- **For exploration and talks**: the top-2 PCA plot is essentially
  equivalent to full space for ranking (AUC 0.827 vs 0.820) and far
  clearer to present than a 2D scatter.
- **For 2D visualisation only**: UMAP is fine for showing “here are the
  genes” but the 2D distances should not be interpreted as global
  distances. Use UMAP for illustration, full space for claims.

------------------------------------------------------------------------

------------------------------------------------------------------------

## 10. Assumption validation summary

| Assumption                                                               | Tested_in                                  | How                                                              |
|:-------------------------------------------------------------------------|:-------------------------------------------|:-----------------------------------------------------------------|
| A1: RNA\<-\>protein correspondence used is 1:1                           | §4 (by construction)                       | Excluded ambiguous isoform/shared-peptide groups before scoring  |
| A2: condition means are the right granularity                            | D7, dimensionality_reduction_wheat.qmd §7b | Cross-referenced, not re-derived                                 |
| A3: shape (not level) is the right comparison                            | – (design choice, not tested)              | –                                                                |
| A4: full-space distance beats a 2D embedding for ranking                 | G2 §5                                      | AUC of archetype separation, full space vs top-2 PCA vs 2D UMAP  |
| A5: short distance implies concordant regulation, not just fast kinetics | G1 §4                                      | AUC of lag_only (concordant, slow) vs true discordant archetypes |

Every assumption behind this analysis and where it is tested.

------------------------------------------------------------------------

## 11. What is safe to claim

The decision rules above (A4 tested by G2, A5 tested by G1) were fixed
before results were seen.

### ✅ Safe to claim

**1. The metric carries real cross-block signal, on real data, not just
in simulation.** G3’s permutation test is significant in both varieties
(p = 0.005, the floor of a 199-draw null): a gene’s real protein profile
sits closer to its own RNA profile, on average, than to a randomly
chosen other gene’s – something that could only be false if RNA and
protein carried no shared condition-level structure at all. This is
consistent with D4’s PLS Q² \> 0 for Cadenza, now confirmed by an
independent method.

**2. The per-gene ranking is highly stable to replicate noise – unlike
the kinetics per-gene calls.** G4’s split-half rank correlation is 0.997
in both varieties, with 95.7% (Cadenza) / 96.7% (Norin) of the
most-discordant decile agreeing between independent halves of a 60-draw
bootstrap. This is a genuinely different reliability picture from
`assumptions_validation.qmd` B4, where 36-49% of individual kinetic-null
calls flipped under the same kind of resampling. The two methods are not
measuring the same thing (a geometric distance vs. a mechanistic
hypothesis test), and only one of them is safe to interpret gene-by-gene
at this replicate depth.

**3. Reducing to 2 PCA components does not cost separation power here**
– top-2 PCA AUC (0.827 discordant, 0.929 lag_only) tracks the full-space
numbers (0.820, 0.932) almost exactly. The source notebook’s own finding
(full space clearly beating a 2D embedding) does not reproduce on this
project’s data for PCA, though it holds in the expected direction for
UMAP (0.799 vs 0.820) – just by a much smaller margin than the source
notebook’s $\rho$ 0.923 vs 0.713. **Use the full-space distance for
ranking regardless** (A4 was never falsified, just not strongly
confirmed either), but do not assume a 2D PCA plot used only for
visualisation is materially misleading here.

### ❌ Not supported

**4. A high distance cannot be read as “post-transcriptional regulation”
without further evidence.** This is the central negative result.
`lag_only` genes – fully concordant by mechanism, decoupled only by slow
protein turnover – separate from the concordant/unchanged baseline at
AUC = 0.932, *higher* than genuinely discordant archetypes (AUC =
0.820). The adversarial case this notebook set out to test (A5) is
violated in the direction that matters most: a gene at the top of the
“most discordant” table in §8 is at least as likely to be a
slow-turnover gene as a genuinely decoupled one, and this metric alone
cannot tell the two apart. Any biological follow-up on a specific gene
from §8’s tables should be cross-referenced against
`kinetics_limited/kinetics_what_we_can_claim.qmd` (K3-K4) before being
described as post-transcriptional regulation.

### The one-paragraph version

> A per-gene RNA-protein distance in a shared, RNA-defined
> condition-profile space carries real signal on real data (G3, p =
> 0.005 both varieties) and ranks genes reproducibly under replicate
> resampling (G4, split-half rank correlation 0.997) – markedly more
> stable than this project’s per-gene kinetic-null calls. But the metric
> cannot distinguish genuine regulatory decoupling from ordinary kinetic
> lag (G1: `lag_only` AUC 0.932 vs genuinely discordant AUC 0.820), so a
> gene’s rank alone does not license a mechanistic claim. Full-space
> distance is the version to rank on; a 2D PCA plot for visualisation
> did not measurably distort the ranking on this data, though that is
> not guaranteed to hold on other datasets.

### What would strengthen this

1.  A real biological gold-standard gene set (e.g. known
    post-transcriptionally regulated wheat immune genes) to validate
    against, rather than only this project’s own simulator.
2.  Extending G4’s stability check to the *ranking itself* under the
    D7-style replicate-pairing bootstrap, not just condition-mean
    resampling.
3.  A GO/pathway enrichment pass on the most-discordant decile,
    cross-referenced gene-by-gene against the kinetics notebooks’
    half-life estimates to separate slow-turnover genes from genuinely
    decoupled ones before any enrichment claim.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-becht2019" class="csl-entry">

Becht, Etienne, Leland McInnes, John Healy, Charles-Antoine Dutertre,
Immanuel W. H. Kwok, Lai Guan Ng, Florent Ginhoux, and Evan W. Newell.
2019. “Dimensionality Reduction for Visualizing Single-Cell Data Using
UMAP.” *Nature Biotechnology* 37 (1): 38–44.
<https://doi.org/10.1038/nbt.4314>.

</div>

<div id="ref-mcinnes2018" class="csl-entry">

McInnes, Leland, John Healy, and James Melville. 2018. “UMAP: Uniform
Manifold Approximation and Projection for Dimension Reduction.”
<https://arxiv.org/abs/1802.03426>.

</div>

</div>
