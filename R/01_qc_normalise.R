#!/usr/bin/env Rscript
## =============================================================================
## 01_qc_normalise.R -- QC, filtering and normalisation for both omics layers
##
## RNA     : filterByExpr -> DESeq2 VST (for integration) + voom logCPM (for DE)
## Protein : decoy/contaminant removal -> unique-peptide filter ->
##           "at least n valid values in >= 1 design cell" -> log2 -> normalise
##
## Two deliberate choices worth defending in the methods section:
##
##  1. The protein validity filter is "at least n valid in AT LEAST ONE design
##     cell", not a global ">= 50% valid". A global rule deletes genuine on/off
##     biology (a protein present in all treated cells and absent from all
##     control cells has 50% missingness and is exactly what you are looking
##     for). This filter runs BEFORE any imputation.
##
##  2. VST is used for the integration matrix and voom for DE. Feeding raw
##     counts or plain logCPM to a Gaussian factor model gives factors that
##     track library size and mean-variance trend rather than biology.
##
## Inputs  : data/simulated/*  (or config paths$data_raw with the same names)
## Outputs : results/qc/*.tsv, figures/fig01_qc.pdf
## Usage   : Rscript R/01_qc_normalise.R [config/config.yaml] [--raw]
## =============================================================================

source("R/utils.R")
need_pkgs(c("limma", "edgeR", "DESeq2", "matrixStats"))

args    <- commandArgs(trailingOnly = TRUE)
cfgpath <- if (length(args) && !grepl("^--", args[1])) args[1] else "config/config.yaml"
cfg     <- load_config(cfgpath)
use_raw <- "--raw" %in% args
indir   <- if (use_raw) cfg$paths$data_raw else cfg$paths$data_sim
outdir  <- ensure_dir(cfg$paths$results, "qc")
figdir  <- ensure_dir(cfg$paths$figures)
set.seed(cfg$project$seed)
log_step("QC from ", indir)

## ------------------------------------------------------------------ load ----
counts  <- read_mat(file.path(indir, "rna_counts.tsv"))
cd_rna  <- read.delim(file.path(indir, "rna_coldata.tsv"), stringsAsFactors = FALSE)
prot    <- read_mat(file.path(indir, "protein_intensities.tsv"))
cd_prot <- read.delim(file.path(indir, "protein_coldata.tsv"), stringsAsFactors = FALSE)
pmeta   <- read.delim(file.path(indir, "protein_meta.tsv"), stringsAsFactors = FALSE)

for (d in c("cd_rna", "cd_prot")) {
  x <- get(d)
  x$cell      <- factor(x$cell, levels = cell_levels(cfg))
  x$treatment <- factor(x$treatment, levels = cfg$design$treatments)
  x$time_f    <- factor(x$time, levels = cfg$design$timepoints)
  rownames(x) <- x$sample_id
  assign(d, x)
}
counts <- counts[, cd_rna$sample_id, drop = FALSE]
prot   <- prot[,  cd_prot$sample_id, drop = FALSE]

## Case B assertion. If this ever fails you are NOT in Case B and should be
## using the sample-matched pipeline (MOFA/DIABLO on 24 paired rows).
shared <- intersect(colnames(counts), colnames(prot))
if (length(shared) > 0)
  warning("RNA and protein share ", length(shared), " sample IDs -- ",
          "if these are genuinely the same biological samples, use Case A.")
n_prot_in <- nrow(prot); n_rna_in <- nrow(counts)
log_step("RNA ", nrow(counts), "x", ncol(counts),
         " | protein ", nrow(prot), "x", ncol(prot))

## ============================================================== RNA-seq =====
keep <- edgeR::filterByExpr(counts, group = cd_rna$cell,
                            min.count = cfg$qc$rna_min_count,
                            min.total.count = 15)
log_step("RNA genes kept by filterByExpr: ", sum(keep), " / ", nrow(counts))
counts_f <- counts[keep, , drop = FALSE]

dge <- edgeR::calcNormFactors(edgeR::DGEList(counts_f, group = cd_rna$cell))
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 3)     # DE / plotting

dds <- DESeq2::DESeqDataSetFromMatrix(counts_f, cd_rna, ~ cell)
dds <- DESeq2::estimateSizeFactors(dds)
vst <- SummarizedExperiment::assay(
  DESeq2::vst(dds, blind = TRUE,
              nsub = min(1000, sum(rowMeans(counts_f) > 5))))            # integration
log_step("RNA VST done; size-factor range ",
         paste(round(range(DESeq2::sizeFactors(dds)), 2), collapse = "-"))

## ============================================================= proteomics ===
## 1. decoys, contaminants, site-only identifications
pmeta <- pmeta[match(rownames(prot), pmeta$group_id), , drop = FALSE]
bad <- (pmeta$reverse == "+") | (pmeta$potential_contaminant == "+") |
       (pmeta$only_identified_by_site == "+") |
       grepl("^(REV__|CON__)", pmeta$majority_protein_ids)
bad[is.na(bad)] <- TRUE
log_step("protein groups removed as decoy/contaminant/site-only: ", sum(bad))
prot  <- prot[!bad, , drop = FALSE]
pmeta <- pmeta[!bad, , drop = FALSE]

## 2. unique-peptide support
keep_pep <- pmeta$unique_peptides >= cfg$qc$prot_min_unique_peptides
log_step("protein groups removed for < ", cfg$qc$prot_min_unique_peptides,
         " unique peptides: ", sum(!keep_pep))
prot  <- prot[keep_pep, , drop = FALSE]
pmeta <- pmeta[keep_pep, , drop = FALSE]

## 3. log2 if needed
if (!isTRUE(cfg$qc$prot_input_is_log2)) {
  prot[prot <= 0] <- NA_real_
  prot <- log2(prot)
}

## 4. "at least n valid values in >= 1 design cell"  (see header note 1)
valid_per_cell <- vapply(levels(cd_prot$cell), function(cl)
  rowSums(!is.na(prot[, cd_prot$cell == cl, drop = FALSE])),
  numeric(nrow(prot)))
keep_valid <- matrixStats::rowMaxs(valid_per_cell) >= cfg$qc$prot_min_valid_in_a_group
log_step("protein groups removed by validity filter: ", sum(!keep_valid),
         " (kept ", sum(keep_valid), ")")
prot  <- prot[keep_valid, , drop = FALSE]
pmeta <- pmeta[keep_valid, , drop = FALSE]

## 5. normalisation -- corrects per-run loading, not biology
norm_method <- cfg$qc$prot_normalisation
prot_n <- switch(norm_method,
  median   = sweep(prot, 2, apply(prot, 2, median, na.rm = TRUE) -
                             median(prot, na.rm = TRUE), "-"),
  quantile = { m <- limma::normalizeQuantiles(prot); m },
  vsn      = { need_pkgs("vsn")
               m <- suppressMessages(vsn::justvsn(2^prot)); m },
  stop("unknown prot_normalisation: ", norm_method))
log_step("protein normalisation: ", norm_method)

## The observed/missing mask travels with the data. Downstream, Universe B
## analyses (concordance, kinetics) must use observed values ONLY -- imputed
## left-censored values manufacture anti-correlation and would produce false
## post-transcriptional regulation calls.
obs_mask <- !is.na(prot_n)

## ================================================================== QC fig ==
pca_scores <- function(m, n = 500) {
  m <- m[order(row_var(m), decreasing = TRUE, na.last = NA)[seq_len(min(n, nrow(m)))], ,
         drop = FALSE]
  m <- m[complete.cases(m), , drop = FALSE]
  p <- prcomp(t(m), scale. = FALSE)
  list(x = p$x, var = round(100 * p$sdev^2 / sum(p$sdev^2), 1))
}

open_pdf(file.path(figdir, "fig01_qc.pdf"), 11, 8.5)
op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))

for (nm in c("RNA (VST)", "Protein (norm)")) {
  m  <- if (nm == "RNA (VST)") vst else prot_n
  cd <- if (nm == "RNA (VST)") cd_rna else cd_prot
  p  <- pca_scores(m)
  cols <- c("#2C6FBB", "#D1495B")[as.integer(cd$treatment)]
  pch  <- c(15, 16, 17, 18)[as.integer(cd$time_f)]
  plot(p$x[, 1], p$x[, 2], col = cols, pch = pch, cex = 1.6,
       xlab = sprintf("PC1 (%.1f%%)", p$var[1]),
       ylab = sprintf("PC2 (%.1f%%)", p$var[2]),
       main = paste0("PCA: ", nm))
  legend("topright", bty = "n", cex = 0.75,
         legend = c(levels(cd$treatment), paste0("t", cfg$design$timepoints)),
         col = c("#2C6FBB", "#D1495B", rep("grey30", 4)),
         pch = c(16, 16, 15, 16, 17, 18))
}

## Missingness diagnosis: a strong negative slope of missingness on intensity
## is the signature of MNAR (left-censoring). This plot decides the imputation
## strategy in script 03 -- do not choose it by habit.
mean_int  <- rowMeans(prot_n, na.rm = TRUE)
miss_rate <- rowMeans(is.na(prot_n))
plot(mean_int, miss_rate, pch = 16, cex = 0.4, col = "#00000055",
     xlab = "mean log2 intensity", ylab = "missing rate",
     main = "Missingness vs abundance\n(negative slope => MNAR)")
ok <- is.finite(mean_int)
if (sum(ok) > 10) {
  lo <- lowess(mean_int[ok], miss_rate[ok], f = 0.4)
  lines(lo, col = "#D1495B", lwd = 2)
  cc <- cor(mean_int[ok], miss_rate[ok], method = "spearman")
  legend("topright", bty = "n", legend = sprintf("Spearman rho = %.2f", cc))
}

barplot(colSums(!is.na(prot_n)), las = 2, cex.names = 0.45,
        col = c("#2C6FBB", "#D1495B")[as.integer(cd_prot$treatment)],
        ylab = "proteins quantified", main = "Detection per MS run")

cr <- cor(prot_n, use = "pairwise.complete.obs")
image(seq_len(ncol(cr)), seq_len(ncol(cr)), cr, axes = FALSE,
      xlab = "", ylab = "", main = "Protein sample correlation",
      col = hcl.colors(32, "Blues", rev = TRUE))
box()

boxplot(as.data.frame(prot_n), las = 2, cex.axis = 0.45, outline = FALSE,
        col = c("#2C6FBB", "#D1495B")[as.integer(cd_prot$treatment)],
        main = paste0("Protein intensity after ", norm_method),
        ylab = "log2 intensity")
par(op); dev.off()
log_step("wrote ", file.path(figdir, "fig01_qc.pdf"))

## ================================================================== write ===
write_mat(vst,      file.path(outdir, "rna_vst.tsv"),        "gene_id")
write_mat(logcpm,   file.path(outdir, "rna_logcpm.tsv"),     "gene_id")
write_mat(counts_f, file.path(outdir, "rna_counts_filt.tsv"),"gene_id")
write_mat(prot_n,   file.path(outdir, "protein_norm.tsv"),   "group_id")
write_mat(obs_mask * 1, file.path(outdir, "protein_obs_mask.tsv"), "group_id")
write_tsv(pmeta,    file.path(outdir, "protein_meta_filt.tsv"))
write_tsv(cd_rna,   file.path(outdir, "rna_coldata.tsv"))
write_tsv(cd_prot,  file.path(outdir, "protein_coldata.tsv"))

qc_summary <- data.frame(
  metric = c("rna_genes_in", "rna_genes_kept", "prot_groups_in",
             "prot_groups_kept", "prot_missing_rate", "prot_norm"),
  value  = c(n_rna_in, nrow(counts_f), n_prot_in,
             nrow(prot_n), sprintf("%.3f", mean(is.na(prot_n))), norm_method))
write_tsv(qc_summary, file.path(outdir, "qc_summary.tsv"))
print(qc_summary)
log_step("01_qc_normalise done")
