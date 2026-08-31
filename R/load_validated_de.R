## =============================================================================
## load_validated_de.R -- RNA differential expression from the project's
##                        validated standalone limma workflow
##
## WHY THIS EXISTS
## ---------------
## `run_de()` in wheat_pipeline.R re-derives the RNA treatment contrasts inside
## this project. Doing so disagreed 2--6x with the project's standalone limma
## workflow (`norinXcadenza-shared`, `DE-varieties-limma`) on counts that match
## exactly -- documented as UNRESOLVED in assumptions_validation.qmd 2b. Rather
## than leave two disagreeing RNA DE results in one project, the standalone
## workflow's published per-timepoint contrasts are read directly. It is the
## result whose numbers appear elsewhere in the project, so adopting it removes
## the disagreement by construction.
##
## `run_de()` is still used for the PROTEIN layer, and for the control-arm time
## coefficients (`ctrl_time`) that the A3 diagnostic needs -- see the note in
## prepare_variety().
##
## TWO GAPS IN THE STANDALONE CSVs, AND HOW EACH IS HANDLED
## --------------------------------------------------------
##  * Standard errors are not exported. Recovered exactly as |logFC / t|: the
##    moderated t is defined as logFC/SE, so this inverts it losslessly.
##  * No omnibus F-statistic is exported. `any_fdr` is therefore the minimum
##    BH-adjusted p across the four timepoint contrasts, NOT an F-test. This is
##    anti-conservative relative to an F-test, so the RNA-responsive gate built
##    on it is more permissive than in earlier renders. Disclose in Methods.
## =============================================================================

#' Path to the standalone limma workflow's wheat DE results.
#'
#' Override with `options(integrateprotrna.standalone_de = "/path/to/dir")`, or
#' by setting the INTEGRATEPROTRNA_STANDALONE_DE environment variable, so this
#' is not pinned to one machine.
standalone_de_dir <- function() {
  p <- getOption("integrateprotrna.standalone_de",
                 Sys.getenv("INTEGRATEPROTRNA_STANDALONE_DE", unset = NA))
  if (is.na(p) || !nzchar(p)) {
    p <- file.path(Sys.getenv("USERPROFILE", Sys.getenv("HOME")),
                   "Documents", "GitHub_projects", "norinXcadenza-shared",
                   "results", "rnaseq", "DE", "wheat")
  }
  p
}

#' Load the standalone workflow's per-timepoint RNA DE contrasts.
#'
#' @param variety "Cadenza" or "Norin" (capitalised as in the filenames).
#' @param timepoints Sampled timepoints in HOURS. Filenames are in days post
#'   inoculation, so these are divided by 24 to build the path.
#' @param path Directory holding `limma_<Variety>_Infected_vs_Control_<N>dpi.csv`.
#'
#' @return A list shaped like `run_de()`'s return value: `lfc`, `se`, `fdr`
#'   (gene x timepoint matrices) and `any_fdr` (named vector). `ctrl_time` is
#'   NOT provided -- the standalone workflow does not export the bare `time`
#'   coefficients. See prepare_variety().
load_validated_de <- function(variety,
                              timepoints = c(0, 24, 48, 72),
                              path = standalone_de_dir()) {

  stopifnot(variety %in% c("Cadenza", "Norin"))
  if (!dir.exists(path)) {
    stop("Standalone DE directory not found: ", path, "\n",
         "Set options(integrateprotrna.standalone_de = ...) or the ",
         "INTEGRATEPROTRNA_STANDALONE_DE environment variable.")
  }

  genes <- NULL
  lfc_mat <- se_mat <- fdr_mat <- NULL

  for (i in seq_along(timepoints)) {
    dpi   <- timepoints[i] / 24
    fname <- file.path(path, sprintf("limma_%s_Infected_vs_Control_%ddpi.csv",
                                     variety, dpi))
    if (!file.exists(fname)) stop("Standalone DE file not found: ", fname)

    res <- read.csv(fname, stringsAsFactors = FALSE)
    colnames(res) <- c("gene_id", "logFC", "AveExpr", "t", "p_value", "fdr", "B")
    res <- res[order(res$gene_id), ]

    if (is.null(genes)) {
      genes   <- res$gene_id
      lfc_mat <- se_mat <- fdr_mat <-
        matrix(NA_real_, nrow = length(genes), ncol = length(timepoints))
    } else if (!identical(genes, res$gene_id)) {
      stop("Gene ID mismatch at timepoint ", timepoints[i], " (", variety, ")")
    }

    lfc_mat[, i] <- res$logFC
    fdr_mat[, i] <- res$fdr

    ## SE = |logFC / t|; exact, since t is defined as logFC/SE.
    se <- rep(NA_real_, nrow(res))
    ok <- res$t != 0 & is.finite(res$t)
    se[ok] <- abs(res$logFC[ok] / res$t[ok])
    se_mat[, i] <- se
  }

  dimnames(lfc_mat) <- dimnames(se_mat) <- dimnames(fdr_mat) <-
    list(genes, sprintf("t%d", timepoints))

  ## Defensive: the ODE fit interpolates each gene's RNA trajectory and needs
  ## >= 2 non-NA points. Complete rows are expected; this guards a malformed file.
  keep      <- rowSums(!is.na(lfc_mat)) >= 2
  n_removed <- sum(!keep)
  lfc_mat <- lfc_mat[keep, , drop = FALSE]
  se_mat  <- se_mat[keep, , drop = FALSE]
  fdr_mat <- fdr_mat[keep, , drop = FALSE]

  ## See the header note: this is min-of-adjusted-p, not an F-test.
  any_fdr <- apply(fdr_mat, 1, min, na.rm = TRUE)
  any_fdr[is.infinite(any_fdr)] <- NA_real_

  message(sprintf("Validated DE (%s): %d genes x %d timepoints%s",
                  variety, nrow(lfc_mat), ncol(lfc_mat),
                  if (n_removed > 0)
                    sprintf(" (%d dropped, <2 valid logFC)", n_removed) else ""))

  list(lfc = lfc_mat, se = se_mat, fdr = fdr_mat, any_fdr = any_fdr)
}
