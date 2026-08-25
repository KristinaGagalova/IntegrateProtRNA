#!/usr/bin/env Rscript
## =============================================================================
## install_r_deps.R -- R dependencies for the IntegrateProtRNA pipeline.
##
##   Rscript env/install_r_deps.R          # everything
##   Rscript env/install_r_deps.R core     # CRAN + Bioc only (skip GitHub)
##
## Every pipeline script degrades gracefully when an optional package is
## absent, so a partial install still gets you a full run with fallbacks.
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
core_only <- length(args) && args[1] == "core"

cran <- c("yaml", "ggplot2", "pheatmap", "matrixStats", "reshape2", "cluster",
          "igraph", "RColorBrewer", "BiocManager", "remotes", "data.table",
          "OmicsPLS", "corpcor", "MASS", "patchwork", "ggrepel", "stringr",
          "tidyr", "dplyr", "tibble", "glmnet", "ranger", "readr")

bioc <- c("limma", "edgeR", "DESeq2", "vsn", "impute", "sva", "IHW",
          "MOFA2", "mixOmics", "ComplexHeatmap", "preprocessCore",
          "fgsea", "GO.db", "AnnotationDbi", "GOSemSim", "BiocParallel")

## imputeLCMD carries QRILC / MinProb; it was archived from CRAN at one point,
## so install it from source if the binary is unavailable.
cran_maybe <- c("imputeLCMD", "norm", "tmvtnorm")

github <- c("abodein/timeOmics")   # CRAN/Bioc build lags; GitHub is canonical

inst <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

repos <- "https://cloud.r-project.org"
need <- inst(cran)
if (length(need)) { message("CRAN: ", paste(need, collapse = ", "))
                    install.packages(need, repos = repos) }

need <- inst(cran_maybe)
if (length(need)) {
  message("CRAN (optional): ", paste(need, collapse = ", "))
  try(install.packages(need, repos = repos), silent = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = repos)
need <- inst(bioc)
if (length(need)) { message("Bioconductor: ", paste(need, collapse = ", "))
                    BiocManager::install(need, ask = FALSE, update = FALSE) }

if (!core_only) {
  for (g in github) {
    pkg <- sub(".*/", "", g)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("GitHub: ", g)
      try(remotes::install_github(g, upgrade = "never"), silent = TRUE)
    }
  }
}

## ------------------------------------------------------------------ report --
all_pkgs <- c(cran, cran_maybe, bioc, sub(".*/", "", github))
status <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
cat("\n==================== dependency status ====================\n")
for (i in order(!status, names(status)))
  cat(sprintf("  %-16s %s\n", names(status)[i], ifelse(status[i], "OK", "MISSING")))
missing <- names(status)[!status]
if (length(missing)) {
  cat("\nMissing (pipeline will use documented fallbacks):\n  ",
      paste(missing, collapse = ", "), "\n")
} else cat("\nAll dependencies present.\n")
