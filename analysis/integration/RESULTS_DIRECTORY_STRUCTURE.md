# Results Directory Structure

## 📂 Organization Principle

**Results mirror Analysis:** Each `analysis/*/` directory has a corresponding `results/*/` directory containing CSV and tabular outputs.

```
analysis/                          results/
├── kinetics_limited/              ├── kinetics_limited/
│   ├── *.qmd (code)              │   └── *.csv (outputs)
│   └── *.md (markdown)           │
│                                 │
├── integration/                   ├── integration/
│   ├── *.qmd (code)              │   └── *.csv (outputs)
│   └── *.R (scripts)             │
│                                 │
├── gene-distance/                 ├── gene-distance/
│   └── *.qmd (code)              │   └── *.csv (outputs)
│                                 │
├── predictions/                   ├── predictions/
│   └── *.ipynb (code)            │   └── *.csv (outputs)
│                                 │
└── dimensionality-reduction/      └── dimensionality-reduction/
    └── *.qmd (code)                  └── *.csv (outputs)
```

---

## 📊 Current Results Structure

### `results/kinetics_limited/`
**Outputs from:** `analysis/kinetics_limited/de_proteomics_wheat.qmd`, `EXTRACT_TOP_GENES_R_SCRIPT.R`

Expected files when scripts run:
- `Supplementary_Table_S1_Cadenza_Top100_Genes.csv` — Top 100 regulated genes (Cadenza)
- `Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv` — Top 100 genes (Norin, provisional)
- `Table_Summary_Statistics.csv` — Summary statistics for methods section

### `results/integration/`
**Outputs from:** `analysis/integration/RUN_CONVERGENCE_ANALYSIS.R`, `method_convergence_figure.qmd`

Expected files when scripts run:
- `high_confidence_ptm_candidates_top50.csv` — Top 50 convergent candidates
- `method_convergence_all_genes.csv` — Full convergence rankings for all genes
- `convergence_summary_statistics.csv` — Method agreement summary statistics

### `results/gene-distance/`
**Outputs from:** `analysis/gene-distance/gene_distance_shared_space.qmd`

Already present:
- `cadenza_gene_distances_8condition.csv` — RNA-protein distance scores (Cadenza)
- `norin_gene_distances_8condition.csv` — RNA-protein distance scores (Norin)
- `cadenza_ml_predictions.csv` — ML proxy predictions (Cadenza)
- `norin_ml_predictions.csv` — ML proxy predictions (Norin)

### `results/predictions/`
**Outputs from:** `analysis/predictions/predict_protein_from_rna.ipynb`

To be populated when:
- Cross-validated R² per gene extracted from ML notebook
- Per-gene model performance metrics exported

### `results/dimensionality-reduction/`
**Outputs from:** `analysis/dimensionality-reduction/dimensionality_reduction_wheat.qmd`

To be populated when dimensionality reduction outputs are exported as CSV

---

## 🚀 Workflow: Scripts → Results

### Step 1: Run Kinetics Analysis
```r
quarto::quarto_render("analysis/kinetics_limited/de_proteomics_wheat.qmd")
```
✓ Generates figures and loads `cad_arch`, `nor_arch` into environment

### Step 2: Extract Gene Tables
```r
source("analysis/integration/EXTRACT_TOP_GENES_R_SCRIPT.R")
```
Outputs to `results/kinetics_limited/`:
- `Supplementary_Table_S1_Cadenza_Top100_Genes.csv`
- `Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv`

### Step 3: Generate Convergence Analysis
```r
source("analysis/integration/RUN_CONVERGENCE_ANALYSIS.R")
```
Outputs to `results/integration/`:
- `high_confidence_ptm_candidates_top50.csv`
- `method_convergence_all_genes.csv`
- `convergence_summary_statistics.csv`

### Step 4: Render Quarto (Optional)
```bash
quarto render analysis/integration/method_convergence_figure.qmd
```
Also outputs to `results/integration/`:
- Same CSV files as Step 3

---

## 📋 Files Created After Running Scripts

```
results/kinetics_limited/
├── Supplementary_Table_S1_Cadenza_Top100_Genes.csv
├── Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv
└── Table_Summary_Statistics.csv

results/integration/
├── high_confidence_ptm_candidates_top50.csv
├── method_convergence_all_genes.csv
└── convergence_summary_statistics.csv

results/gene-distance/
├── cadenza_gene_distances_8condition.csv (already present)
├── norin_gene_distances_8condition.csv (already present)
├── cadenza_ml_predictions.csv (already present)
└── norin_ml_predictions.csv (already present)
```

---

## 🔄 Adding Results from New Analyses

When you run new analyses and want to export results:

1. **Create the results directory** (if not exists):
   ```r
   output_dir <- here::here("results", "[analysis_name]")
   dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
   ```

2. **Export CSV files** to that directory:
   ```r
   write.csv(results_df, 
             file.path(output_dir, "output_filename.csv"),
             row.names = FALSE)
   ```

3. **Update this documentation** with the new expected outputs

---

## 📖 Manuscript Integration

### Pull results for manuscript sections:

**Methods:**
- Copy summary statistics from `results/kinetics_limited/Table_Summary_Statistics.csv`
- Include convergence breakdown from `results/integration/convergence_summary_statistics.csv`

**Supplementary Material:**
- Table S1: Use `Supplementary_Table_S1_Cadenza_Top100_Genes.csv`
- Table S2: Use `high_confidence_ptm_candidates_top50.csv`
- Optionally: Full rankings from `method_convergence_all_genes.csv`

**Discussion:**
- Reference high-confidence candidates from `high_confidence_ptm_candidates_top50.csv`
- Cross-reference with gene-distance scores from `results/gene-distance/cadenza_gene_distances_8condition.csv`

---

## ✅ Checklist: Results Organization

- [ ] `results/kinetics_limited/` exists
- [ ] `results/integration/` exists
- [ ] `results/predictions/` exists
- [ ] `results/dimensionality-reduction/` exists
- [ ] Scripts updated to output to correct results/ directories
- [ ] All scripts pushed to nectar
- [ ] Ran Step 1-3 and verified CSV files in results/
- [ ] Results ready for manuscript integration

---

**Status:** ✅ Results directory structure established and scripts configured

**Last Updated:** 2026-08-29
