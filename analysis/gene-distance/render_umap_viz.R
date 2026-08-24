#!/usr/bin/env Rscript
"""
Render UMAP visualizations for gene-distance analysis.
Run this after the main gene_distance_shared_space.qmd analysis completes.
"""

if (!requireNamespace("here", quietly = TRUE)) stop("package 'here' required")

if (file.exists(here::here("renv", "activate.R"))) {
  Sys.setenv(RENV_PROJECT = here::here())
  source(here::here("renv", "activate.R"))
}

source(here::here("R", "utils.R"))
source(here::here("R", "wheat_pipeline.R"))

need_pkgs(c("limma", "edgeR", "DESeq2", "matrixStats", "impute", "uwot", "yaml"))

set.seed(20260823)

cat("\n=== UMAP Visualizations for Gene Distance Analysis ===\n")

# Load data
cat("Loading data...\n")
DESIGN <- wheat_design()
cad <- prepare_variety("cadenza", DESIGN)
nor <- prepare_variety("norin",   DESIGN)
REAL <- list(Cadenza = cad, Norin = nor)

# Source the functions from the main analysis
cat("Loading analysis functions...\n")

row_center <- function(M) M - rowMeans(M)

fit_col_scaler <- function(M) {
  mu  <- colMeans(M)
  sdv <- apply(M, 2, sd); sdv[sdv < 1e-8] <- 1
  function(X) sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
}

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

# Main analysis function
gene_distance_real <- function(v, n_perm = 199, seed = 1) {
  common <- intersect(rownames(v$qc_rna$vst), rownames(v$imputed$mixed))
  R_cond <- cell_means(v$qc_rna$vst[common, , drop = FALSE],    v$meta)
  P_cond <- cell_means(v$imputed$mixed[common, , drop = FALSE], v$meta)

  sp     <- build_shared_space(R_cond, P_cond)
  d_full <- euclid(sp$Rs, sp$Ps)
  d_pca2 <- euclid(sp$Zr, sp$Zp, k = 2)

  # UMAP fit on RNA, transform protein
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

  list(genes = data.frame(gene_id = common, d_full = d_full,
                          d_pca2 = d_pca2, d_umap = d_umap),
       space = list(Er = Er, Ep = Ep),
       obs = obs, null = null, p_perm = p_perm, n_genes = length(common))
}

# Run analysis
cat("Computing gene distances and UMAP embeddings...\n")
GD <- setNames(Map(function(v, s) gene_distance_real(v, seed = s), REAL, seq_along(REAL)), names(REAL))

# Plotting functions
plot_shared_umap <- function(space_list, dist_df, variety_name, output_file) {
  Er <- space_list$Er
  Ep <- space_list$Ep
  distances <- dist_df$d_umap

  png(output_file, width = 1400, height = 500, res = 100)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Left: colored by distance
  color_scale <- colorRampPalette(c("#2C6FBB", "#E8A33D", "#D1495B"))(100)
  color_idx <- pmin(99, pmax(1, round((distances - min(distances)) /
                                       (max(distances) - min(distances)) * 99) + 1))

  plot(Er[, 1], Er[, 2], type = "n", xlab = "UMAP1", ylab = "UMAP2",
       main = paste(variety_name, ": RNA-Protein Shared UMAP\n(colored by distance)"),
       cex.main = 1.1)

  points(Er[, 1], Er[, 2], col = color_scale[color_idx], pch = 19,
         cex = 2, alpha = 0.6)
  points(Ep[, 1], Ep[, 2], col = color_scale[color_idx], pch = 17,
         cex = 2, alpha = 0.6)

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
  dev.off()
}

plot_top_genes_umap <- function(space_list, dist_df, variety_name, output_file, n_top = 10) {
  Er <- space_list$Er
  Ep <- space_list$Ep
  distances <- dist_df$d_umap

  png(output_file, width = 800, height = 700, res = 100)
  par(mar = c(4, 4, 3, 1))

  plot(Er[, 1], Er[, 2], type = "n", xlab = "UMAP1", ylab = "UMAP2",
       main = paste(variety_name, ": Top", n_top, "Most/Least Discordant Genes"),
       cex.main = 1.1)

  points(Er[, 1], Er[, 2], col = rgb(0.5, 0.5, 0.5, 0.2), pch = 19, cex = 1.5)
  points(Ep[, 1], Ep[, 2], col = rgb(0.5, 0.5, 0.5, 0.2), pch = 17, cex = 1.5)

  # Top N most discordant
  top_idx <- order(-distances)[seq_len(n_top)]
  points(Er[top_idx, 1], Er[top_idx, 2], col = "#D1495B", pch = 19, cex = 3.5)
  points(Ep[top_idx, 1], Ep[top_idx, 2], col = "#D1495B", pch = 17, cex = 3.5)

  for (i in top_idx) {
    lines(c(Er[i, 1], Ep[i, 1]), c(Er[i, 2], Ep[i, 2]),
          col = "#D1495B", lwd = 1.5, alpha = 0.5)
  }

  # Bottom N most concordant
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

  dev.off()
}

plot_distance_dist <- function(GD, output_file) {
  png(output_file, width = 1200, height = 400, res = 100)
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

  par(op)
  dev.off()
}

# Create output directory
out_dir <- here::here("analysis", "gene-distance", "figs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Generate plots
cat("Generating visualizations...\n")

for (nm in c("Cadenza", "Norin")) {
  cat("  ", nm, "...\n")
  prefix <- tolower(nm)

  plot_shared_umap(
    GD[[nm]]$space, GD[[nm]]$genes, nm,
    file.path(out_dir, paste0(prefix, "_shared_umap.png"))
  )

  plot_top_genes_umap(
    GD[[nm]]$space, GD[[nm]]$genes, nm,
    file.path(out_dir, paste0(prefix, "_top_genes.png")),
    n_top = 10
  )
}

plot_distance_dist(GD, file.path(out_dir, "distance_distribution.png"))

cat("✓ All visualizations saved to", out_dir, "\n\n")

# Save distance tables
cat("Saving distance tables...\n")
for (nm in c("Cadenza", "Norin")) {
  prefix <- tolower(nm)
  sorted <- GD[[nm]]$genes[order(-GD[[nm]]$genes$d_full), ]
  write.csv(sorted, file.path(out_dir, paste0(prefix, "_gene_distances.csv")), row.names = FALSE)
  cat("  ", paste0(prefix, "_gene_distances.csv"), "\n")
}

cat("\nDone!\n")
