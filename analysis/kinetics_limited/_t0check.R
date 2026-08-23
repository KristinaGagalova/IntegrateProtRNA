## Baseline (0 hpi) sanity check.
## At t0 the inoculation has only just happened, so Infected vs Control should
## be near-identical. Any large difference there is a confound, not biology.
setwd("C:/Claude/Projects/IntegrateProtRNA")
suppressMessages({library(limma); library(edgeR)})

meta_all <- read.csv("data/real/cadenza_metadata.csv", row.names = 1)

check <- function(prefix) {
  meta <- read.csv(sprintf("data/real/%s_metadata.csv", prefix), row.names = 1)
  meta$tp <- factor(meta$timepoint)
  meta$tr <- relevel(factor(meta$treatment), ref = "T0")

  rna <- as.matrix(read.csv(sprintf("data/real/%s-rnaseq.csv", prefix),
                            row.names = 1, check.names = FALSE))
  rna <- rna[, rownames(meta), drop = FALSE]

  keep <- rowSums(rna >= 10) >= 3
  dge  <- calcNormFactors(DGEList(rna[keep, ]))

  des <- model.matrix(~ tr * tp, data = meta)
  colnames(des) <- make.names(colnames(des))
  v   <- voom(dge, des, plot = FALSE)
  fit <- eBayes(lmFit(v, des), trend = TRUE, robust = TRUE)

  # the bare treatment coefficient == Infected vs Control AT t0
  tt <- topTable(fit, coef = "trT1", number = Inf, sort.by = "none")
  n_sig <- sum(tt$adj.P.Val < 0.01 & abs(tt$logFC) > 1)

  # for contrast: same test at the LAST timepoint
  cm  <- makeContrasts(last = trT1 + trT1.tpt3, levels = des)
  f2  <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
  tl  <- topTable(f2, coef = "last", number = Inf, sort.by = "none")
  n_last <- sum(tl$adj.P.Val < 0.01 & abs(tl$logFC) > 1)

  # how separated are the two arms at t0, in PCA space?
  lc <- cpm(dge, log = TRUE, prior.count = 3)
  t0 <- rownames(meta)[meta$timepoint == "t0"]
  p  <- prcomp(t(lc[order(apply(lc, 1, var), decreasing = TRUE)[1:2000], t0]))
  grp <- meta[t0, "treatment"]
  sep <- abs(diff(tapply(p$x[, 1], grp, mean))) / sd(p$x[, 1])

  data.frame(variety = prefix,
             genes_tested = sum(keep),
             DE_at_t0_0hpi = n_sig,
             DE_at_t3_72hpi = n_last,
             pc1_arm_separation_at_t0_SD = round(sep, 2))
}

res <- rbind(check("cadenza"), check("norin"))
print(res, row.names = FALSE)
