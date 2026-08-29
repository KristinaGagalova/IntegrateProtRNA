#!/usr/bin/env Rscript
# Method Convergence Figure (Fig 5) - Integration of Kinetics + Distance + ML
# This script combines results from three complementary analyses
# Output: Publication-quality scatter plot showing method agreement

library(ggplot2)
library(dplyr)
library(tidyr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

OUTPUT_FILE <- "fig-method-convergence-cadenza.png"
VARIETY <- "Cadenza"
WIDTH <- 12
HEIGHT <- 8
DPI <- 300

# ==============================================================================
# STEP 1: Load Kinetics Results
# ==============================================================================

# Expected: cad_arch from de_proteomics_wheat.qmd
# If not in environment, load from saved RData or extract from markdown output

load_kinetics <- function() {
  # OPTION A: If de_proteomics_wheat.qmd was just rendered
  # (cad_arch already in environment from: cad_arch <- run_kinetic_archetypes(...))

  if (exists("cad_arch")) {
    message("✓ Kinetics results loaded from environment (cad_arch)")

    kinetics_df <- data.frame(
      gene_id = cad_arch$gene_id,
      amplitude = cad_arch$amplitude_a,
      regulated = cad_arch$regulated,
      rna_responds = cad_arch$rna_responds,
      lrt_fdr = cad_arch$lrt_fdr,
      archetype = cad_arch$archetype
    )

    return(kinetics_df)
  } else {
    message("⚠️  cad_arch not found in environment")
    message("   To use this script: run de_proteomics_wheat.qmd first")
    return(NULL)
  }
}

# ==============================================================================
# STEP 2: Load Distance Results
# ==============================================================================

load_distance <- function() {
  # Expected: results from gene-distance analysis
  # Need to extract d_full (distance score) from gene_distance_shared_space.qmd output

  # PLACEHOLDER: This would require parsing the markdown or loading saved data
  message("⚠️  Distance results extraction requires parsing gene_distance output")
  message("   For now, create mock distance data")

  # MOCK DATA (replace with actual distance extraction)
  distance_df <- data.frame(
    gene_id = character(0),
    distance = numeric(0)
  )

  if (nrow(distance_df) == 0) {
    message("   Distance data not available in this environment")
  }

  return(distance_df)
}

# ==============================================================================
# STEP 3: Load ML Predictions
# ==============================================================================

load_ml_predictions <- function() {
  # Expected: per-gene R² from predict_protein_from_rna.ipynb
  # Best model's cross-validated R² score on unseen genes

  # PLACEHOLDER: This requires parsing the Python notebook or loading saved results

  message("⚠️  ML predictions extraction requires parsing Python output")
  message("   For now, create mock ML confidence data")

  # MOCK DATA (replace with actual ML extraction)
  ml_df <- data.frame(
    gene_id = character(0),
    ml_r2 = numeric(0)
  )

  if (nrow(ml_df) == 0) {
    message("   ML predictions not available in this environment")
  }

  return(ml_df)
}

# ==============================================================================
# STEP 4: Combine and Compute Convergence
# ==============================================================================

compute_convergence <- function(kinetics_df, distance_df, ml_df) {

  if (is.null(kinetics_df) || nrow(kinetics_df) == 0) {
    message("ERROR: Kinetics data required")
    return(NULL)
  }

  # Start with kinetics results
  combined <- kinetics_df %>%
    filter(rna_responds == TRUE) %>%
    select(gene_id, amplitude, regulated, lrt_fdr, archetype)

  # Add distance scores (if available)
  if (!is.null(distance_df) && nrow(distance_df) > 0) {
    combined <- combined %>%
      left_join(distance_df, by = "gene_id")
    combined$distance_flag <- if_else(is.na(combined$distance), 0,
                                      if_else(combined$distance > quantile(combined$distance, 0.75, na.rm = TRUE), 1, 0))
  } else {
    combined$distance <- NA
    combined$distance_flag <- 0
  }

  # Add ML predictions (if available)
  if (!is.null(ml_df) && nrow(ml_df) > 0) {
    combined <- combined %>%
      left_join(ml_df, by = "gene_id")
    combined$ml_flag <- if_else(is.na(combined$ml_r2), 0,
                                if_else(combined$ml_r2 < quantile(combined$ml_r2, 0.25, na.rm = TRUE), 1, 0))
  } else {
    combined$ml_r2 <- NA
    combined$ml_flag <- 0
  }

  # Compute convergence: how many methods agree (0-3)
  combined <- combined %>%
    mutate(
      kinetics_flag = as.integer(regulated),
      convergence = kinetics_flag + distance_flag + ml_flag,
      convergence_label = case_when(
        convergence == 3 ~ "All 3 methods",
        convergence == 2 ~ "2 methods agree",
        convergence == 1 ~ "1 method flags",
        TRUE ~ "No agreement"
      ),
      convergence_label = factor(convergence_label,
                                 levels = c("All 3 methods", "2 methods agree", "1 method flags", "No agreement"))
    )

  return(combined)
}

# ==============================================================================
# STEP 5: Create Convergence Figure
# ==============================================================================

create_figure <- function(combined_df, output_file, variety) {

  if (is.null(combined_df) || nrow(combined_df) == 0) {
    message("ERROR: No data to plot")
    return(NULL)
  }

  # Filter to genes with both metrics
  plot_data <- combined_df %>%
    filter(!is.na(amplitude))

  message(sprintf("Plotting %d genes", nrow(plot_data)))

  # Create scatter plot
  p <- ggplot(plot_data, aes(x = abs(amplitude), y = distance, color = convergence_label, size = convergence)) +
    geom_point(alpha = 0.6, stroke = 0.5) +
    scale_color_manual(
      values = c(
        "All 3 methods" = "#e74c3c",
        "2 methods agree" = "#f39c12",
        "1 method flags" = "#95a5a6",
        "No agreement" = "#ecf0f1"
      ),
      name = "Method Agreement"
    ) +
    scale_size_continuous(name = "Convergence Score", range = c(2, 6)) +
    labs(
      title = sprintf("Method Convergence: Post-Transcriptional Modification Candidates (%s)", variety),
      x = "Kinetics Amplitude |a| (RNA-protein divergence)",
      y = "Distance Score (RNA-protein profile discordance)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "top",
      panel.grid.major = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  # Save figure
  ggsave(output_file, p, width = WIDTH, height = HEIGHT, dpi = DPI, bg = "white")
  message(sprintf("✓ Figure saved: %s", output_file))

  return(p)
}

# ==============================================================================
# STEP 6: Summary Statistics
# ==============================================================================

print_summary <- function(combined_df) {

  if (is.null(combined_df) || nrow(combined_df) == 0) return(NULL)

  message("\n=== METHOD CONVERGENCE SUMMARY ===\n")
  message(sprintf("Total genes analyzed: %d", nrow(combined_df)))

  conv_summary <- combined_df %>%
    group_by(convergence) %>%
    summarise(count = n(), pct = 100 * n() / nrow(combined_df), .groups = "drop") %>%
    arrange(desc(convergence))

  message("\nMethod agreement distribution:")
  for (i in 1:nrow(conv_summary)) {
    row <- conv_summary[i, ]
    message(sprintf("  %d methods agree: %5d genes (%5.1f%%)", row$convergence, row$count, row$pct))
  }

  # High-confidence candidates
  high_conf <- combined_df %>%
    filter(convergence == 3) %>%
    arrange(desc(abs(amplitude)))

  message(sprintf("\n🎯 HIGH-CONFIDENCE PTM CANDIDATES (All 3 methods agree): %d genes", nrow(high_conf)))
  if (nrow(high_conf) > 0) {
    message("   Top 10:")
    for (i in 1:min(10, nrow(high_conf))) {
      row <- high_conf[i, ]
      message(sprintf("   %2d. %-25s | Amplitude: %6.2f | Distance: %s",
                      i, row$gene_id, row$amplitude,
                      if_else(is.na(row$distance), "NA", sprintf("%6.2f", row$distance))))
    }
  }
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

message("Creating Method Convergence Figure (Fig 5)")
message("=" %+% strrep("=", 58))

# Load data
kinetics <- load_kinetics()
distance <- load_distance()
ml <- load_ml_predictions()

if (!is.null(kinetics)) {
  # Compute convergence
  combined <- compute_convergence(kinetics, distance, ml)

  if (!is.null(combined)) {
    # Create figure
    create_figure(combined, OUTPUT_FILE, VARIETY)

    # Print summary
    print_summary(combined)

    message("\n✓ Method Convergence Figure complete!")
  }
} else {
  message("\nTo generate the convergence figure, follow these steps:")
  message("  1. Run de_proteomics_wheat.qmd in R to generate cad_arch results")
  message("  2. Extract distance results from gene-distance analysis")
  message("  3. Extract ML R² scores from predict_protein_from_rna.ipynb")
  message("  4. Run this script in R environment where those data are available")
}

