#!/usr/bin/env Rscript
## =============================================================================
## 04_univariate_temporal_de.R -- per-layer temporal differential expression
##
## limma is used as the SINGLE engine for both layers (voom for counts,
## trend+robust eBayes for protein) so that effect sizes are directly
## comparable across omics. That comparability is what makes script 05
## possible; mixing DESeq2 shrunken LFCs with limma protein LFCs would put the
## two axes of the concordance plot on different scales.
##
## Design: ~ 0 + cell, cell = treatment_time (8 levels). From it:
##   tp_k     Trt vs Ctrl at timepoint k                      (4 contrasts)
##   ANY      F-test over the 4 timepoint contrasts           "any effect"
##   INT      F-test over (Trt-Ctrl)_k - (Trt-Ctrl)_0, k>0    "trajectory
##            differs by treatment" -- the true treatment x time interaction
##
## With 4 timepoints, an 8-level cell factor is preferred over a spline basis:
## it makes no smoothness assumption and every contrast is interpretable.
## A DESeq2 spline LRT is run alongside as an independent cross-check.
##
## KEY OUTPUT for script 05: logFC and its moderated standard error per gene
## per timepoint, for both layers. The SE is not optional -- with n=3 the
## concordance analysis has to be error-aware or it will read noise as
## post-transcriptional regulation.
##
## Usage: Rscript R/04_univariate_temporal_de.R [config/config.yaml]
## =============================================================================

source("R/utils.R")
need_pkgs(c("limma", "edgeR", "matrixStats"))

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- load_config(if (length(args)) args[1] else "config/config.yaml")
qcdir  <- file.path(cfg$paths$results, "qc")
mapdir <- file.path(cfg$paths$results, "mapping")
outdir <- ensure_dir(cfg$paths$results, "de")
figdir <- ensure_dir(cfg$paths$figures)
set.seed(cfg$project$seed)

tps  <- cfg$design$timepoints
trts <- cfg$design$treatments
ref  <- cfg$design$reference_treatment
alt  <- setdiff(trts, ref)
stopifnot(length(alt) == 1)
FDR  <- cfg$de$fdr

counts <- read_mat(file.path(qcdir, "rna_counts_filt.tsv"))
cd_rna <- read.delim(file.path(qcdir, "rna_coldata.tsv"), stringsAsFactors = FALSE)
prot   <- read_mat(file.path(qcdir, "protein_imputed_mixed.tsv"))
prot_obs <- read_mat(file.path(qcdir, "protein_norm.tsv"))
cd_prot  <- read.delim(file.path(qcdir, "protein_coldata.tsv"), stringsAsFactors = FALSE)
pairs    <- read.delim(file.path(mapdir, "universeB_pairs.tsv"), stringsAsFactors = FALSE)

for (d in c("cd_rna", "cd_prot")) {
  x <- get(d); x$cell <- factor(x$cell, levels = cell_levels(cfg))
  rownames(x) <- x$sample_id; assign(d, x)
}

## ------------------------------------------------- contrast construction ----
make_design <- function(cd) {
  d <- model.matrix(~ 0 + cell, cd)
  colnames(d) <- levels(cd$cell)
  d
}
cn <- function(tr, t) sprintf("%s_t%g", tr, t)

tp_contrasts <- setNames(
  sprintf("%s - %s", cn(alt, tps), cn(ref, tps)), paste0("tp_t", tps))
## Interaction: change in the treatment effect relative to baseline.
int_contrasts <- setNames(
  sprintf("(%s - %s) - (%s - %s)",
          cn(alt, tps[-1]), cn(ref, tps[-1]), cn(alt, tps[1]), cn(ref, tps[1])),
  paste0("int_t", tps[-1]))

fit_layer <- function(mat, cd, layer) {
  design <- make_design(cd)
  if (layer == "rna") {
    dge <- edgeR::calcNormFactors(edgeR::DGEList(mat))
    v   <- limma::voom(dge, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    trend <- FALSE
  } else {
    fit <- limma::lmFit(mat, design)
    trend <- TRUE
  }
  cm  <- limma::makeContrasts(contrasts = c(tp_contrasts, int_contrasts),
                              levels = design)
  colnames(cm) <- c(names(tp_contrasts), names(int_contrasts))
  f2  <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = trend, robust = TRUE)
  f2
}

## Moderated standard error of each contrast: stdev.unscaled * sqrt(s2.post).
## This is the posterior SD limma actually uses for its t-statistics.
contrast_se <- function(f2) f2$stdev.unscaled * sqrt(f2$s2.post)

## F-test across a set of contrast columns (any-effect / interaction tests).
f_over <- function(f2, cols) {
  ff <- limma::classifyTestsF(f2[, cols], fstat.only = TRUE)
  df1 <- attr(ff, "df1"); df2 <- attr(ff, "df2")
  p <- pf(as.numeric(ff), df1, df2, lower.tail = FALSE)
  data.frame(F = as.numeric(ff), P = p, FDR = p.adjust(p, "BH"),
             row.names = rownames(f2))
}

## ================================================================== fits =====
log_step("fitting RNA (voom + limma)")
fit_rna <- fit_layer(counts, cd_rna, "rna")

log_step("fitting protein (limma, trend + robust) on mixed-imputed matrix")
fit_prot_grp <- fit_layer(prot, cd_prot, "prot")

## Gene-level protein matrix for Universe B: collapse the imputed groups with
## the same linear-sum rule used in 02, so DE and mapping stay consistent.
collapse_to_gene <- function(m, pairs) {
  out <- t(vapply(seq_len(nrow(pairs)), function(i) {
    g <- strsplit(pairs$groups[i], ";", fixed = TRUE)[[1]]
    g <- g[g %in% rownames(m)]
    if (!length(g)) return(rep(NA_real_, ncol(m)))
    if (length(g) == 1) return(m[g, ])
    s <- colSums(2^m[g, , drop = FALSE], na.rm = TRUE)
    ## colSums(na.rm=TRUE) returns 0 when every group is missing, and log2(0)
    ## is -Inf, which complete.cases() does NOT catch. Mark it NA explicitly.
    s[colSums(!is.na(m[g, , drop = FALSE])) == 0] <- NA_real_
    log2(s)
  }, numeric(ncol(m))))
  rownames(out) <- pairs$gene_id; colnames(out) <- colnames(m)
  out[!is.finite(out)] <- NA_real_
  out[complete.cases(out), , drop = FALSE]
}
prot_gene <- collapse_to_gene(prot, pairs)
log_step("gene-level protein matrix: ", nrow(prot_gene), " genes")
fit_prot <- fit_layer(prot_gene, cd_prot, "prot")

## Sensitivity: observed-values-only (limma fits each row on its non-missing
## values). Compared against the imputed fit below.
prot_gene_obs <- collapse_to_gene(prot_obs, pairs)
## All 8 cell means must be estimable from observed data alone, otherwise the
## row yields partial-NA coefficients and breaks the eBayes variance trend.
obs_per_cell <- vapply(levels(cd_prot$cell), function(cl)
  rowSums(!is.na(prot_gene_obs[, cd_prot$cell == cl, drop = FALSE])),
  numeric(nrow(prot_gene_obs)))
keep_obs <- matrixStats::rowMins(obs_per_cell) >= 2
log_step("observed-only sensitivity fit uses ", sum(keep_obs), " / ",
         nrow(prot_gene_obs), " genes (>=2 observed values in every cell)")
fit_prot_obs <- fit_layer(prot_gene_obs[keep_obs, , drop = FALSE], cd_prot, "prot")

## ================================================== extract lfc / se / fdr ==
extract <- function(f2, tag) {
  lfc <- f2$coefficients[, names(tp_contrasts), drop = FALSE]
  se  <- contrast_se(f2)[, names(tp_contrasts), drop = FALSE]
  colnames(lfc) <- colnames(se) <- paste0("t", tps)
  padj <- apply(f2$p.value[, names(tp_contrasts), drop = FALSE], 2, p.adjust, "BH")
  colnames(padj) <- paste0("t", tps)
  list(lfc = lfc, se = se, fdr = padj,
       any = f_over(f2, names(tp_contrasts)),
       int = f_over(f2, names(int_contrasts)), tag = tag)
}
R <- extract(fit_rna,  "rna")
P <- extract(fit_prot, "prot")
Pg <- extract(fit_prot_grp, "prot_group")

for (X in list(R, P, Pg)) {
  write_mat(X$lfc, file.path(outdir, sprintf("%s_lfc.tsv", X$tag)),  "feature_id")
  write_mat(X$se,  file.path(outdir, sprintf("%s_se.tsv",  X$tag)),  "feature_id")
  write_mat(X$fdr, file.path(outdir, sprintf("%s_fdr.tsv", X$tag)),  "feature_id")
  write_tsv(data.frame(feature_id = rownames(X$any), X$any),
            file.path(outdir, sprintf("%s_Ftest_any.tsv", X$tag)))
  write_tsv(data.frame(feature_id = rownames(X$int), X$int),
            file.path(outdir, sprintf("%s_Ftest_interaction.tsv", X$tag)))
}

## --------------------------------------- DESeq2 spline LRT cross-check -----
if (has_pkg("DESeq2") && has_pkg("splines")) {
  log_step("DESeq2 spline LRT cross-check (~ treatment * ns(time, df=2))")
  cdd <- cd_rna; cdd$treatment <- factor(cdd$treatment, levels = trts)
  X   <- splines::ns(cdd$time, df = 2)
  colnames(X) <- c("s1", "s2"); cdd <- cbind(cdd, as.data.frame(X))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, cdd, ~ treatment * s1 + treatment * s2)
  dds <- DESeq2::DESeq(dds, test = "LRT",
                       reduced = ~ treatment + s1 + s2, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(dds))
  write_tsv(data.frame(gene_id = rownames(res), res),
            file.path(outdir, "rna_deseq2_spline_LRT.tsv"))
  agree <- table(limma = R$int$FDR < FDR,
                 deseq2 = ifelse(is.na(res$padj), FALSE, res$padj < FDR)[
                   match(rownames(R$int), rownames(res))])
  log_step("limma vs DESeq2 interaction agreement:")
  print(agree)
}

## ================================================================ summary ===
n_sig <- function(X) c(vapply(paste0("t", tps), function(k) sum(X$fdr[, k] < FDR), 0),
                       ANY = sum(X$any$FDR < FDR), INT = sum(X$int$FDR < FDR))
summ <- rbind(RNA = n_sig(R), Protein_gene = n_sig(P), Protein_group = n_sig(Pg))
cat("\n=== features at FDR < ", FDR, " ===\n", sep = ""); print(summ)
write_tsv(data.frame(layer = rownames(summ), summ), file.path(outdir, "de_summary.tsv"))

## Imputed vs observed-only sensitivity, on the shared genes.
common <- intersect(rownames(P$lfc), rownames(fit_prot_obs$coefficients))
lfc_obs <- fit_prot_obs$coefficients[common, names(tp_contrasts), drop = FALSE]
sens <- vapply(seq_along(tps), function(k)
  cor(P$lfc[common, k], lfc_obs[, k], use = "complete.obs"), 0)
names(sens) <- paste0("t", tps)
cat("\nlogFC correlation, mixed-imputed vs observed-only:\n"); print(round(sens, 3))
write_tsv(data.frame(timepoint = names(sens), cor = round(sens, 4)),
          file.path(outdir, "imputation_sensitivity_lfc.tsv"))

## ================================================================= figure ===
open_pdf(file.path(figdir, "fig02_univariate_de.pdf"), 12, 8)
op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3.5, 1))

bp <- barplot(t(summ[, paste0("t", tps)]), beside = TRUE, las = 1,
              col = hcl.colors(length(tps), "Blues3"),
              ylab = sprintf("features FDR<%.2f", FDR),
              main = "Treatment effect per timepoint")
legend("topleft", bty = "n", cex = 0.8, fill = hcl.colors(length(tps), "Blues3"),
       legend = paste0("t", tps))

for (k in seq_along(tps)[1:min(3, length(tps))]) {
  kk <- k + 1  # skip t0 where there is no treatment effect yet
  if (kk > length(tps)) break
  plot(R$lfc[, kk], -log10(pmax(R$fdr[, kk], 1e-300)), pch = 16, cex = 0.35,
       col = ifelse(R$fdr[, kk] < FDR, "#D1495B55", "#00000033"),
       xlab = "RNA log2FC", ylab = "-log10 FDR",
       main = sprintf("RNA volcano, t%g", tps[kk]))
  abline(h = -log10(FDR), lty = 2, col = "grey40")
}

## The overlap panel is deliberately the weakest-looking one: it sets up why
## per-layer DE is insufficient and integration is needed (fig 3).
sig_r <- rownames(R$any)[R$any$FDR < FDR]
sig_p <- rownames(P$any)[P$any$FDR < FDR]
univ  <- intersect(rownames(R$any), rownames(P$any))
ov <- c(`RNA only` = length(setdiff(intersect(sig_r, univ), sig_p)),
        Both       = length(intersect(intersect(sig_r, univ), sig_p)),
        `Prot only`= length(setdiff(intersect(sig_p, univ), sig_r)))
bb <- barplot(ov, col = c("#2C6FBB", "#8E7CC3", "#D1495B"), ylab = "genes",
              main = sprintf("Any-treatment-effect overlap\n(%d testable genes)", length(univ)))
text(bb, ov, ov, pos = 1, col = "white")

boxplot(list(RNA = as.vector(R$lfc[, -1]), Protein = as.vector(P$lfc[, -1])),
        col = c("#2C6FBB", "#D1495B"), outline = FALSE,
        ylab = "log2FC", main = "Effect-size scale by layer")
abline(h = 0, lty = 2)
par(op); dev.off()
log_step("04_univariate_temporal_de done -> ", outdir)
