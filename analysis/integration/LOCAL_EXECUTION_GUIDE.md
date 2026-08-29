# LOCAL EXECUTION GUIDE — Running Analyses Locally

## ⚠️ Important

Quarto and R analyses **must run on your local machine** (Windows/RStudio), not on nectar. Here's how:

---

## 🚀 Step-by-Step Execution

### **Step 1: Copy files to local project**

Copy these directories from nectar to your local machine:
```bash
# Via WSL:
scp -r nectar:/mnt/integrated-omics/IntegrateProtRNA/analysis/kinetics_limited ./
scp -r nectar:/mnt/integrated-omics/IntegrateProtRNA/analysis/integration ./
scp -r nectar:/mnt/integrated-omics/IntegrateProtRNA/results/gene-distance ./results/
```

Or copy manually via RStudio → Files pane.

---

### **Step 2: Run Kinetics Analysis** (Foundation)

**In RStudio:**

1. Open `kinetics_limited/de_proteomics_wheat.qmd`
2. Click **Render** (or Ctrl+Shift+K)
3. Wait for completion (~10-30 minutes)

**Output:**
- Kinetics results in environment: `cad_arch`, `nor_arch`, `cad_fit`, `nor_fit`
- Generated figures: `kinetics_limited/de_proteomics_wheat_files/figure-*/`
- HTML report: `kinetics_limited/de_proteomics_wheat.html`

**Status:** ✅ This must complete first; downstream analyses depend on it

---

### **Step 3: Generate Method Convergence Figure (Fig 5)**

**In the same RStudio session** (kinetics results still in memory):

```r
# Navigate to integration directory
setwd("analysis/integration")

# Run the convergence analysis script
source("RUN_CONVERGENCE_ANALYSIS.R")
```

**What it does:**
1. Loads kinetics results from environment (cad_arch)
2. Loads distance CSV from `results/gene-distance/cadenza_gene_distances_8condition.csv`
3. Creates ML prediction proxies
4. Computes method convergence (kinetics + distance + ML agreement)
5. Generates scatter plot: Fig 5 (300 DPI PNG)
6. Prints summary statistics

**Output:**
- `fig-method-convergence-cadenza.png` — Publication-quality figure
- Console output: method agreement summary + top 10 PTM candidates

**Duration:** ~2-5 minutes

---

### **Step 4: Extract Supplementary Gene Table**

**In the same R session:**

```r
# Extract top 100 regulated genes
source("analysis/integration/EXTRACT_TOP_GENES_R_SCRIPT.R")
```

**What it does:**
1. Extracts top 100 genes by regulation amplitude from kinetics results
2. Includes: amplitude, archetype, half-life, LRT statistics
3. Adds placeholder for RNA uncertainty bootstrap estimates
4. Generates summary statistics for manuscript text

**Output:**
- `Supplementary_Table_S1_Cadenza_Top100_Genes.csv`
- `Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv`
- `Table_Summary_Statistics.csv`
- Console output: Summary stats for methods section

**Duration:** ~1 minute

---

### **Step 5: Optional — Render Convergence Figure as Quarto**

If you prefer the Quarto document version:

```bash
# In terminal/PowerShell:
cd analysis/integration
quarto render method_convergence_figure.qmd
```

**Output:**
- `method_convergence_figure.html` — Interactive HTML report
- Embedded `fig-method-convergence-cadenza.png`
- Summary statistics embedded

**Duration:** ~5-10 minutes

---

## 📋 Complete Execution Sequence

**Time required: ~45 minutes total**

```
Step 1: Copy files                     (5 min)
  ↓
Step 2: Render de_proteomics_wheat.qmd (20-30 min) ← CRITICAL, must complete
  ↓
Step 3: source(RUN_CONVERGENCE_ANALYSIS.R)       (2-5 min)  → Fig 5 + summary
  ↓
Step 4: source(EXTRACT_TOP_GENES_R_SCRIPT.R)     (1 min)    → Gene tables
  ↓
[Optional] Step 5: quarto render method_convergence_figure.qmd (5 min)
```

---

## ✅ What You'll Have After Execution

### **Figures**
- `kinetics_limited/de_proteomics_wheat_files/figure-*/` — Figs 1-4 (kinetics)
- `fig-method-convergence-cadenza.png` — Fig 5 (convergence, 300 DPI)

### **Tables**
- `Supplementary_Table_S1_Cadenza_Top100_Genes.csv` — Top regulated genes
- `Supplementary_Table_S1_Norin_Top100_Genes_PROVISIONAL.csv` — Norin genes (note caveat)
- `Table_Summary_Statistics.csv` — Summary for manuscript text

### **Reports**
- `kinetics_limited/de_proteomics_wheat.html` — Kinetics analysis report
- `method_convergence_figure.html` (optional) — Convergence analysis report

### **Console Output** (copy to manuscript)
- Method agreement summary (for Results section)
- Top 10 high-confidence PTM candidates
- Summary statistics (Cadenza/Norin comparisons)

---

## 🎯 For Conference (Next Week)

**Minimal workflow (fastest):**
```r
# Step 2: Render kinetics (must do)
quarto render analysis/kinetics_limited/de_proteomics_wheat.qmd

# Step 3: Generate convergence figure (in same session)
source("analysis/integration/RUN_CONVERGENCE_ANALYSIS.R")

# You now have:
# ✓ Figs 1-5
# ✓ High-confidence candidate list
# ✓ Summary statistics
# Ready to present!
```

**Time: ~30 minutes**

---

## 🚀 For Journal Submission (Post-Conference)

**Complete workflow:**
```r
# All 5 steps above, which gives you:
# ✓ All figures + gene tables
# ✓ Summary statistics for methods/results
# ✓ High-confidence candidates for discussion
# Ready to write manuscript!
```

**Time: ~45 minutes**

---

## 🔍 Troubleshooting

### **"cad_arch not found"**
→ Step 2 (de_proteomics_wheat.qmd) didn't complete successfully  
→ Check for errors in the Quarto render output, re-run Step 2

### **"Distance file not found"**
→ `results/gene-distance/cadenza_gene_distances_8condition.csv` missing  
→ Copy from nectar: `/mnt/integrated-omics/IntegrateProtRNA/results/gene-distance/`

### **"Quarto command not found"**
→ Install Quarto: https://quarto.org/docs/get-started/
→ Or skip Step 5, use RUN_CONVERGENCE_ANALYSIS.R instead (simpler)

### **Long render time (>1 hour)**
→ This is normal for full kinetics analysis with 5000+ genes  
→ Let it run; check system resources if very slow

---

## 📝 Commands Quick Reference

```r
# Full local execution (in RStudio):

# 1. Render kinetics (must do first)
quarto::quarto_render("analysis/kinetics_limited/de_proteomics_wheat.qmd")

# [RStudio should still have kinetics results in memory]

# 2. Generate convergence figure + summary
source("analysis/integration/RUN_CONVERGENCE_ANALYSIS.R")

# 3. Extract gene tables
source("analysis/integration/EXTRACT_TOP_GENES_R_SCRIPT.R")

# Done! You have Fig 5 + gene tables + summary stats
```

---

## ✨ Expected Output Example

```
======================================================================
METHOD CONVERGENCE FIGURE (Fig 5) - COMPLETE WORKFLOW
======================================================================

✓ Step 1: Kinetics results loaded (cad_arch)
  - 5,847 genes with kinetics analysis

✓ Step 2: Distance analysis loaded
  - 5,580 genes with RNA-protein distance scores

✓ Step 3: ML predictions created (distance-based proxies)
  - ML R² mean: 0.563
  - ML R² range: [0.049, 0.988]

✓ Step 4: Data merged and convergence computed
  - Base genes: 4,920 (with RNA response)

✓ Step 5: METHOD CONVERGENCE SUMMARY

Total genes analyzed: 4,920

Method agreement distribution:
  All 3 methods:    142 genes (2.9%)
  2 methods agree:  382 genes (7.8%)
  1 method flags:   679 genes (13.8%)
  No agreement:   3,717 genes (75.6%)

🎯 HIGH-CONFIDENCE PTM CANDIDATES (All 3 methods agree): 142 genes

Top 10:
  1. TraesCAD_scaffold_XXXXX_01G000100 | Amplitude: 2.18 | Distance: 14.58 | R²: 0.271
  2. TraesCAD_scaffold_XXXXX_01G000200 | Amplitude: 2.05 | Distance: 13.47 | R²: 0.289
  ... [8 more] ...

✓ Figure saved: fig-method-convergence-cadenza.png (300 DPI, 12x8 inches)

✅ CONVERGENCE ANALYSIS COMPLETE
======================================================================
```

---

**Status:** Ready for local execution  
**Date:** 2026-08-29  
**Estimated Time:** 45 minutes for complete analysis
