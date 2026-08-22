## utils.R -- shared helpers for the IntegrateProtRNA pipeline (Case B)
## Sourced by every R/NN_*.R script.

suppressWarnings(suppressMessages({
  library(yaml)
}))

## ---------------------------------------------------------------- config ----
load_config <- function(path = "config/config.yaml") {
  if (!file.exists(path)) stop("config not found: ", path,
                               " (run scripts from the project root)")
  cfg <- yaml::read_yaml(path)
  cfg$design$timepoints <- as.numeric(cfg$design$timepoints)
  ## YAML 1.1 parses "2.0e7" as a STRING (it requires an explicit exponent
  ## sign). Coerce anything that looks numeric so a config typo cannot
  ## silently propagate into a model.
  coerce <- function(x) {
    if (is.list(x)) return(lapply(x, coerce))
    if (is.character(x) && length(x) == 1L) {
      n <- suppressWarnings(as.numeric(x))
      if (!is.na(n)) return(n)
    }
    x
  }
  cfg <- coerce(cfg)
  cfg
}

ensure_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

## Fail loudly and early rather than half-way through a long run.
need_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("missing R packages: ", paste(missing, collapse = ", "),
         "\n  install with: Rscript env/install_r_deps.R", call. = FALSE)
  }
  invisible(TRUE)
}

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

## ---------------------------------------------------------------- design ----
## The design cell is the unit that links the two omics layers in Case B.
make_coldata <- function(cfg, layer = c("rna", "prot")) {
  layer <- match.arg(layer)
  d <- expand.grid(
    rep       = seq_len(cfg$design$n_reps),
    time      = cfg$design$timepoints,
    treatment = cfg$design$treatments,
    stringsAsFactors = FALSE
  )
  d <- d[order(d$treatment, d$time, d$rep), ]
  d$cell      <- sprintf("%s_t%g", d$treatment, d$time)
  d$sample_id <- sprintf("%s_%s_r%d", layer, d$cell, d$rep)
  d$treatment <- factor(d$treatment, levels = cfg$design$treatments)
  d$time_f    <- factor(d$time, levels = cfg$design$timepoints)
  d$cell      <- factor(d$cell, levels = cell_levels(cfg))
  d$layer     <- layer
  rownames(d) <- d$sample_id
  d[, c("sample_id", "layer", "treatment", "time", "time_f", "rep", "cell")]
}

cell_levels <- function(cfg) {
  as.vector(t(outer(cfg$design$treatments, cfg$design$timepoints,
                    function(a, b) sprintf("%s_t%g", a, b))))
}

## Design-cell means: the pseudo-samples used for Case B integration.
cell_means <- function(mat, coldata, na.rm = TRUE) {
  cells <- levels(coldata$cell)
  out <- vapply(cells, function(cl) {
    idx <- which(coldata$cell == cl)
    rowMeans(mat[, idx, drop = FALSE], na.rm = na.rm)
  }, numeric(nrow(mat)))
  out[is.nan(out)] <- NA_real_
  colnames(out) <- cells
  out
}

## Non-parametric replicate bootstrap of design-cell means. Every Case B
## integration statistic must be run through this to get an honest CI --
## the point estimate on its own hides the n=3 uncertainty completely.
boot_cell_means <- function(mat, coldata, B = 200, seed = 1L, FUN = identity) {
  set.seed(seed)
  cells <- levels(coldata$cell)
  idx_by_cell <- lapply(cells, function(cl) which(coldata$cell == cl))
  replicate(B, {
    bidx <- unlist(lapply(idx_by_cell, function(ii) sample(ii, length(ii), TRUE)))
    bcol <- coldata[bidx, , drop = FALSE]
    FUN(cell_means(mat[, bidx, drop = FALSE], bcol))
  }, simplify = FALSE)
}

## ------------------------------------------------------------- transforms ---
logit  <- function(p) log(p / (1 - p))
expit  <- function(x) 1 / (1 + exp(-x))

## Row-wise variance ignoring NAs, used for HVG/HVP selection.
row_var <- function(m) {
  mu <- rowMeans(m, na.rm = TRUE)
  n  <- rowSums(!is.na(m))
  v  <- rowSums((m - mu)^2, na.rm = TRUE) / pmax(n - 1, 1)
  v[n < 2] <- NA_real_
  v
}

top_variable <- function(m, n) {
  v <- row_var(m)
  head(names(sort(v, decreasing = TRUE, na.last = NA)), n)
}

## Exact solution of dP/dt = ks*R(t) - kd*P, with R linear on each grid step.
## Used by both the simulator and the R-side kinetic sanity checks, so the
## generative model and the fitted model share one implementation.
integrate_protein <- function(t_grid, R_grid, ks, kd, P0 = NULL) {
  stopifnot(length(t_grid) == length(R_grid), kd > 0)
  if (is.null(P0)) P0 <- ks * R_grid[1] / kd     # pre-treatment steady state
  P <- numeric(length(t_grid)); P[1] <- P0
  for (i in seq_len(length(t_grid) - 1L)) {
    dt <- t_grid[i + 1L] - t_grid[i]
    a  <- R_grid[i]
    b  <- if (dt > 0) (R_grid[i + 1L] - R_grid[i]) / dt else 0
    alpha <- ks * a / kd - ks * b / kd^2
    beta  <- ks * b / kd
    P[i + 1L] <- alpha + beta * dt + (P[i] - alpha) * exp(-kd * dt)
  }
  P
}

## -------------------------------------------------------------------- io ----
write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  invisible(path)
}

write_mat <- function(m, path, idcol = "feature_id") {
  df <- data.frame(id = rownames(m), m, check.names = FALSE)
  names(df)[1] <- idcol
  write_tsv(df, path)
}

read_mat <- function(path) {
  df <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  m  <- as.matrix(df[, -1, drop = FALSE])
  rownames(m) <- df[[1]]
  m
}

log_step <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
}

## ------------------------------------------------------------------ plot ----
open_pdf <- function(path, width = 9, height = 7) {
  ensure_dir(dirname(path)); grDevices::pdf(path, width = width, height = height)
}

## Vectorised over genes (rows). Same exact piecewise-linear solution as
## integrate_protein(), but recursing over time with whole columns at once.
integrate_protein_mat <- function(t_grid, R_mat, ks, kd, P0 = NULL) {
  stopifnot(ncol(R_mat) == length(t_grid), all(kd > 0))
  n <- nrow(R_mat)
  if (is.null(P0)) P0 <- ks * R_mat[, 1] / kd
  P <- matrix(NA_real_, n, length(t_grid),
              dimnames = list(rownames(R_mat), NULL))
  P[, 1] <- P0
  for (i in seq_len(length(t_grid) - 1L)) {
    dt    <- t_grid[i + 1L] - t_grid[i]
    a     <- R_mat[, i]
    b     <- (R_mat[, i + 1L] - a) / dt
    alpha <- ks * a / kd - ks * b / kd^2
    beta  <- ks * b / kd
    P[, i + 1L] <- alpha + beta * dt + (P[, i] - alpha) * exp(-kd * dt)
  }
  P
}
