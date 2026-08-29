# Extract Top Genes for Supplementary Table S1
# Run this in R/RStudio after loading kinetics results from de_proteomics_wheat.qmd
# Output: formatted CSV ready for publication

# ==============================================================================
# STEP 1: Load kinetics results from the main analysis
# ==============================================================================

# Option A: If you've already rendered de_proteomics_wheat.qmd in this session
#          (cad_arch, nor_arch, and cad_fit, nor_fit are in environment)

# Option B: Otherwise, source the kinetics functions and re-run the fitting
source(here::here("R", "07_bayesian_kinetics.R"))

# Assume you've already loaded:
# - cad_arch: Cadenza kinetics archetype results
# - cad_fit: Cadenza Bayesian fit summary
# - cad_rna_resp: indices of Cadenza genes with RNA response
# - cad_Rl, cad_Rs: Cadenza RNA logFC and SE matrices

# ==============================================================================
# STEP 2: Extract Cadenza top genes by regulation amplitude
# ==============================================================================

# Build results table for Cadenza
cad_results <- data.frame(
  Gene_ID = cad_arch$gene_id,
  Regulation_Amplitude = cad_arch$amplitude_a,
  Regulated_LRT = cad_arch$regulated,           # Profile-likelihood test
  RNA_Response = cad_arch$rna_responds,         # Did RNA change?
  Half_Life_h = cad_arch$t_half_M1_h,          # Protein half-life (hours)
  Identifiable = cad_arch$identifiable,         # Is half-life in resolvable range?
  LRT_p_value = cad_arch$lrt_p,
  LRT_FDR = cad_arch$lrt_fdr,
  Archetype = cad_arch$archetype,
  stringsAsFactors = FALSE
)

# Filter to genes with RNA response (more interpretable)
cad_top <- subset(cad_results, RNA_Response == TRUE)

# Sort by absolute amplitude (most regulated first)
cad_top <- cad_top[order(-abs(cad_top$Regulation_Amplitude)), ]

# Extract top 100 genes for supplementary table
cad_top_100 <- head(cad_top, 100)

# ==============================================================================
# STEP 3: Add RNA measurement uncertainty bootstrap estimates
# ==============================================================================

# If you ran the bootstrap uncertainty analysis (from assumptions_validation.qmd),
# load those results here. If not, note "Not yet available" in the table.

# Assuming bootstrap results are in 'rna_uncertainty' dataframe:
# cad_top_100 <- merge(cad_top_100,
#                      rna_uncertainty[, c("Gene_ID", "Bootstrap_SD", "Call_Flip_Rate")],
#                      by = "Gene_ID", all.x = TRUE)

# For now, add placeholder column with guidance
cad_top_100$RNA_Uncertainty_Call_Flip_Pct <- NA  # Placeholder
cad_top_100$Notes <- "See assumptions_validation.md §B4 for RNA uncertainty analysis"

# ==============================================================================
# STEP 4: Repeat for Norin (with caveat about baseline asymmetry)
# ==============================================================================

nor_results <- data.frame(
  Gene_ID = nor_arch$gene_id,
  Regulation_Amplitude = nor_arch$amplitude_a,
  Regulated_LRT = nor_arch$regulated,
  RNA_Response = nor_arch$rna_responds,
  Half_Life_h = nor_arch$t_half_M1_h,
  Identifiable = nor_arch$identifiable,
  LRT_p_value = nor_arch$lrt_p,
  LRT_FDR = nor_arch$lrt_fdr,
  Archetype = nor_arch$archetype,
  stringsAsFactors = FALSE
)

nor_top <- subset(nor_results, RNA_Response == TRUE)
nor_top <- nor_top[order(-abs(nor_top$Regulation_Amplitude)), ]
nor_top_100 <- head(nor_top, 100)

nor_top_100$RNA_Uncertainty_Call_Flip_Pct <- NA
nor_top_100$Notes <- "PROVISIONAL: See assumptions_validation.md §A0 for baseline asymmetry caveat"

# ==============================================================================
# STEP 5: Format and save for publication
# ==============================================================================

# Column order for publication
col_order <- c("Gene_ID", "Archetype", "Regulation_Amplitude", "Regulated_LRT",
               "RNA_Response", "Half_Life_h", "Identifiable", "LRT_FDR",
               "RNA_Uncertainty_Call_Flip_Pct", "Notes")

cad_top_100_pub <- cad_top_100[, col_order]
nor_top_100_pub <- nor_top_100[, col_order]

# Rename columns for human readability
colnames(cad_top_100_pub) <- c(
  "Gene ID",
  "Archetype",
  "Amplitude (a)",
  "Regulated (LRT)",
  "RNA Response",
  "Protein Half-Life (h)",
  "Half-Life Identifiable",
  "LRT FDR",
  "RNA Uncertainty: Call Flip (%)",
  "Notes"
)

colnames(nor_top_100_pub) <- colnames(cad_top_100_pub)

# ==============================================================================
# STEP 6: Save CSV files to results/kinetics_limited/
# ==============================================================================

results_dir <- here::here("results", "kinetics_limited")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(cad_top_100_pub,
          file = file.path(results_dir, "Supplementary_Table_S1_Cadenza_Top100_Genes.csv"),
          row.names = FALSE)

write.csv(nor_top_100_pub,
          file = file.path(results_dir, "Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv"),
          row.names = FALSE)

cat("✓ Supplementary tables exported to results/kinetics_limited/:\n")
cat("  - Supplementary_Table_S1_Cadenza_Top100_Genes.csv\n")
cat("  - Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv\n")

# ==============================================================================
# STEP 7: Print summary statistics for manuscript text
# ==============================================================================

cat("\n=== SUMMARY STATISTICS FOR MANUSCRIPT ===\n")
cat(sprintf("Cadenza: %d genes with RNA response\n", nrow(cad_top)))
cat(sprintf("  - Regulated (LRT): %d (%.1f%%)\n",
            sum(cad_top$Regulated_LRT), 100*mean(cad_top$Regulated_LRT)))
cat(sprintf("  - Amplitude mean: %.2f (SD: %.2f)\n",
            mean(abs(cad_top$Regulation_Amplitude)),
            sd(abs(cad_top$Regulation_Amplitude))))
cat(sprintf("  - Half-lives in resolvable range [12h-72h]: %d (%.1f%%)\n",
            sum(cad_top$Identifiable), 100*mean(cad_top$Identifiable)))

cat("\n")
cat(sprintf("Norin: %d genes with RNA response (PROVISIONAL)\n", nrow(nor_top)))
cat(sprintf("  - Regulated (LRT): %d (%.1f%%)\n",
            sum(nor_top$Regulated_LRT), 100*mean(nor_top$Regulated_LRT)))
cat(sprintf("  - Amplitude mean: %.2f (SD: %.2f)\n",
            mean(abs(nor_top$Regulation_Amplitude)),
            sd(abs(nor_top$Regulation_Amplitude))))
cat(sprintf("  - Half-lives in resolvable range [12h-72h]: %d (%.1f%%)\n",
            sum(nor_top$Identifiable), 100*mean(nor_top$Identifiable)))

# ==============================================================================
# STEP 8: Optional - Create summary table for main text
# ==============================================================================

summary_table <- data.frame(
  Variety = c("Cadenza", "Norin"),
  Total_Genes_Analyzed = c(nrow(cad_results), nrow(nor_results)),
  RNA_Responsive = c(nrow(cad_top), nrow(nor_top)),
  Post_Transcriptionally_Regulated = c(sum(cad_top$Regulated_LRT), sum(nor_top$Regulated_LRT)),
  Mean_Regulation_Amplitude = c(mean(abs(cad_top$Regulation_Amplitude)),
                                mean(abs(nor_top$Regulation_Amplitude))),
  Identifiable_Half_Lives = c(sum(cad_top$Identifiable), sum(nor_top$Identifiable))
)

print(summary_table)

write.csv(summary_table, file = "Table_Summary_Statistics.csv", row.names = FALSE)
cat("\n✓ Summary table saved: Table_Summary_Statistics.csv\n")

