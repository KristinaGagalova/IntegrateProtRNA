## =============================================================================
## wheat_pipeline.R -- shared upstream steps for the wheat Case B analyses
##
## Loading, QC, missingness and temporal DE, factored out so that every
## analysis notebook under analysis/*/ runs the SAME upstream code rather than
## keeping its own copy. Behaviour is identical to the versions originally
## written inline in analysis/kinetics_limited/de_proteomics_wheat.qmd.
##
## Everything takes an explicit `design` argument (built by wheat_design())
## rather than reading a global, so two notebooks with different timepoint
## conventions cannot silently interfere with each other.
##
## Usage:
##   source(here::here("R", "wheat_pipeline.R"))
##   DESIGN <- wheat_design()
##   cad    <- load_variety("cadenza", DESIGN)
## =============================================================================

source(here::here("R", "utils.R"))

# ---------------------------------------------------------------------------
# Design constants
# ---------------------------------------------------------------------------

#' Build the design description shared by every wheat notebook.
#'
#' Defaults encode what the RNA-seq sample sheet confirms:
#'   T0 = NEG   = uninfected control (reference)
#'   T1 = PN143 = infected
#'   t0..t3     = 0 / 24 / 48 / 72 hours post inoculation
wheat_design <- function(treatment_ref    = "T0",
                         timepoint_values = c(0, 24, 48, 72),
                         n_reps           = 3) {

  treatments <- c(treatment_ref, setdiff(c("T0", "T1"), treatment_ref))

  cells <- as.vector(t(outer(
    treatments, timepoint_values,
    function(a, b) sprintf("%s_t%s", a, b)
  )))

  list(
    treatments = treatments,
    reference  = treatment_ref,
    timepoints = timepoint_values,
    n_reps     = n_reps,
    cells      = cells,
    time_unit  = "hours"
  )
}

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

#' Load one variety's RNA, protein, metadata and gene mapping.
#'
#' The LH-style sample IDs coincide between the RNA and protein matrices for a
#' given variety (both trace back to the same field-plot label), but per the
#' Case B design this is NOT treated as sample-level pairing anywhere
#' downstream: only treatment x timepoint CELL membership links the layers.
load_variety <- function(prefix, design, data_dir = here::here("data", "real")) {

  meta <- read.csv(
    file.path(data_dir, paste0(prefix, "_metadata.csv")),
    row.names = 1
  )

  meta$time_num <- design$timepoints[
    match(meta$timepoint, sort(unique(meta$timepoint)))
  ]

  meta$treatment <- factor(meta$treatment, levels = design$treatments)

  meta$cell <- factor(
    sprintf("%s_t%s", meta$treatment, meta$time_num),
    levels = design$cells
  )

  rna <- as.matrix(read.csv(
    file.path(data_dir, paste0(prefix, "-rnaseq.csv")),
    row.names = 1, check.names = FALSE
  ))

  prot <- as.matrix(read.csv(
    file.path(data_dir, paste0(prefix, "-prot.csv")),
    row.names = 1, check.names = FALSE
  ))

  gmap <- read.csv(
    file.path(data_dir, paste0(prefix, "_protein_gene_mapping.csv"))
  )

  stopifnot(
    all(colnames(rna)  %in% rownames(meta)),
    all(colnames(prot) %in% rownames(meta))
  )

  list(
    rna  = rna[,  rownames(meta), drop = FALSE],
    prot = prot[, rownames(meta), drop = FALSE],
    meta = meta,
    gmap = gmap
  )
}

# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

#' RNA QC: manual low-expression filter, then VST + voom-ready counts.
#'
#' The filter is the DESeq2-style rule (>= min_count in at least as many
#' samples as the smallest design group), not edgeR::filterByExpr's adaptive
#' heuristic -- it is more permissive and matches the criterion used in the
#' project's standalone limma DE workflow.
run_qc_rna <- function(counts, meta, design,
                       min_count = 10) {

  dge  <- edgeR::DGEList(counts = counts)
  keep <- rowSums(dge$counts >= min_count) >= design$n_reps

  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cf  <- dge$counts

  dge    <- edgeR::calcNormFactors(dge)
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 3)

  dds <- DESeq2::DESeqDataSetFromMatrix(cf, meta, ~cell)
  dds <- DESeq2::estimateSizeFactors(dds)

  vst <- SummarizedExperiment::assay(
    DESeq2::vst(dds, blind = TRUE, nsub = min(1000, sum(rowMeans(cf) > 5)))
  )

  list(
    counts = cf,
    logcpm = logcpm,
    vst    = vst,
    n_in   = nrow(counts),
    n_kept = nrow(cf)
  )
}

#' Protein QC.
#'
#' Zero is treated as "not quantified" (standard LFQ convention), not a
#' biological zero. The validity filter keeps a protein with >= n valid values
#' in AT LEAST ONE design cell -- a global ">= 50% valid" rule would delete
#' genuine on/off biology. No decoy/contaminant/unique-peptide filtering is
#' possible here: the supplied mapping files carry no such columns.
run_qc_prot <- function(intensity, meta, min_valid_in_cell = 2) {

  m <- intensity
  m[m == 0] <- NA
  m <- log2(m)

  valid_per_cell <- vapply(
    levels(meta$cell),
    function(cl) rowSums(!is.na(m[, meta$cell == cl, drop = FALSE])),
    numeric(nrow(m))
  )

  keep <- matrixStats::rowMaxs(valid_per_cell) >= min_valid_in_cell
  mf   <- m[keep, , drop = FALSE]

  # Median normalisation: corrects per-run loading, not biology
  shift <- apply(mf, 2, median, na.rm = TRUE) - median(mf, na.rm = TRUE)
  mn    <- sweep(mf, 2, shift, "-")

  list(
    norm      = mn,
    n_in      = nrow(intensity),
    n_kept    = nrow(mf),
    miss_rate = mean(is.na(mn))
  )
}

# ---------------------------------------------------------------------------
# Missingness
# ---------------------------------------------------------------------------

#' MNAR/MAR diagnosis and mechanism-specific imputation.
#'
#' A protein missing in EVERY replicate of a cell (but seen elsewhere) is MNAR
#' -- genuinely low/absent there, so it gets left-censored imputation. Missing
#' in only SOME replicates of its cell is MAR (stochastic) and gets kNN.
run_missingness <- function(prot_norm, meta,
                            downshift = 1.8,
                            width     = 0.3) {

  cells <- levels(meta$cell)
  is_na <- is.na(prot_norm)

  n_obs_cell <- vapply(
    cells,
    function(cl) rowSums(!is.na(prot_norm[, meta$cell == cl, drop = FALSE])),
    numeric(nrow(prot_norm))
  )

  lab <- matrix("obs", nrow(prot_norm), ncol(prot_norm),
                dimnames = dimnames(prot_norm))

  for (k in seq_along(cells)) {
    ii   <- which(meta$cell == cells[k])
    gone <- n_obs_cell[, k] == 0

    lab[gone, ii] <- "MNAR"
    lab[!gone, ii][is_na[!gone, ii]] <- "MAR"
  }

  lab[!is_na] <- "obs"

  mar_imp <- suppressWarnings(
    impute::impute.knn(prot_norm, k = 10, rowmax = 0.95, colmax = 0.95)$data
  )

  mnar_imp <- prot_norm

  for (j in seq_len(ncol(prot_norm))) {
    v   <- prot_norm[, j]
    nas <- is.na(v)
    if (!any(nas)) next

    mu  <- mean(v, na.rm = TRUE)
    sdv <- sd(v,   na.rm = TRUE)

    # Perseus-style down-shifted normal
    mnar_imp[nas, j] <- rnorm(sum(nas), mu - downshift * sdv, width * sdv)
  }

  mixed <- prot_norm
  mixed[lab == "MAR"]  <- mar_imp[lab == "MAR"]
  mixed[lab == "MNAR"] <- mnar_imp[lab == "MNAR"]

  stopifnot(!anyNA(mixed))

  list(
    mixed     = mixed,
    labels    = lab,
    mnar_frac = mean(lab == "MNAR"),
    mar_frac  = mean(lab == "MAR")
  )
}

# ---------------------------------------------------------------------------
# Temporal differential expression
# ---------------------------------------------------------------------------

#' Temporal DE for one layer, treatment-contrast parameterisation.
#'
#' Uses ~treatment*time with Control and t0 as reference levels, then
#' reconstructs the per-timepoint Infected-vs-Control contrasts (main effect,
#' plus the interaction term at each later timepoint). Equivalent to a
#' ~0+cell means model, but the bare `time` coefficients are also directly
#' interpretable as the CONTROL arm's own drift, which the A3 diagnostic uses.
run_de <- function(mat, meta, design, is_rna) {

  meta$treatment <- relevel(
    factor(meta$treatment, levels = design$treatments),
    ref = design$reference
  )

  meta$time_f <- relevel(
    factor(meta$time_num, levels = design$timepoints),
    ref = as.character(design$timepoints[1])
  )

  dm <- model.matrix(~ treatment * time_f, data = meta)
  colnames(dm) <- make.names(colnames(dm))

  if (is_rna) {
    dge <- edgeR::calcNormFactors(edgeR::DGEList(mat))
    v   <- limma::voom(dge, dm, plot = FALSE)
    fit <- limma::lmFit(v, dm)
  } else {
    fit <- limma::lmFit(mat, dm)
  }

  alt       <- setdiff(design$treatments, design$reference)
  main_coef <- paste0("treatment", alt)
  int_coefs <- paste0(main_coef, ".time_f", design$timepoints[-1])

  stopifnot(all(c(main_coef, int_coefs) %in% colnames(dm)))

  contrast_defs <- setNames(
    c(main_coef, paste(main_coef, int_coefs, sep = " + ")),
    paste0("tp_t", design$timepoints)
  )

  cm <- limma::makeContrasts(contrasts = contrast_defs, levels = dm)
  colnames(cm) <- names(contrast_defs)

  f2 <- limma::eBayes(limma::contrasts.fit(fit, cm),
                      trend = TRUE, robust = TRUE)

  lfc  <- f2$coefficients
  se   <- f2$stdev.unscaled * sqrt(f2$s2.post)
  padj <- apply(f2$p.value, 2, p.adjust, "BH")

  colnames(lfc) <- colnames(se) <- colnames(padj) <-
    paste0("t", design$timepoints)

  # as.numeric() strips rownames -- reattach explicitly, or every downstream
  # any_fdr[gene_subset] silently returns NA instead of erroring.
  ff    <- limma::classifyTestsF(f2, fstat.only = TRUE)
  p_any <- setNames(
    pf(as.numeric(ff), attr(ff, "df1"), attr(ff, "df2"), lower.tail = FALSE),
    rownames(f2)
  )

  # Control-arm drift: the bare time coefficients
  time_coefs <- paste0("time_f", design$timepoints[-1])
  ctrl_time  <- NULL

  if (all(time_coefs %in% colnames(fit$coefficients))) {
    ctrl_time <- fit$coefficients[, time_coefs, drop = FALSE]
    colnames(ctrl_time) <- paste0("t", design$timepoints[-1])
  }

  list(
    lfc       = lfc,
    se        = se,
    fdr       = padj,
    any_fdr   = p.adjust(p_any, "BH"),
    ctrl_time = ctrl_time
  )
}

#' Convenience wrapper: load -> QC -> impute -> DE for one variety.
#'
#' Returns everything the downstream notebooks need, including the t0-centred
#' response trajectories used for clustering and kinetic modelling.
prepare_variety <- function(prefix, design, seed = 1) {

  set.seed(seed)

  v  <- load_variety(prefix, design)
  qr <- run_qc_rna(v$rna,  v$meta, design)
  qp <- run_qc_prot(v$prot, v$meta)
  mi <- run_missingness(qp$norm, v$meta)

  de_rna <- run_de(qr$counts, v$meta, design, is_rna = TRUE)

  common  <- intersect(rownames(mi$mixed), rownames(qr$counts))
  de_prot <- run_de(mi$mixed[common, , drop = FALSE], v$meta, design,
                    is_rna = FALSE)

  # t0-centred response trajectories (Infected vs Control, relative to t0)
  Rl <- de_rna$lfc[common, , drop = FALSE]
  Pl <- de_prot$lfc
  Rs <- de_rna$se[common, , drop = FALSE]
  Ps <- de_prot$se

  Rl <- Rl - Rl[, 1]
  Pl <- Pl - Pl[, 1]

  list(
    label   = prefix,
    meta    = v$meta,
    qc_rna  = qr,
    qc_prot = qp,
    imputed = mi,
    de_rna  = de_rna,
    de_prot = de_prot,
    common  = common,
    Rl = Rl, Rs = Rs,
    Pl = Pl, Ps = Ps
  )
}
