# Directory Organization — Publication Integration

## ✅ Reorganization Complete

### **`kinetics_limited/`** — CORE KINETICS ANALYSIS (Keep)
```
kinetics_limited/
├── de_proteomics_wheat.md/.qmd              ← Main kinetics analysis & figures
├── kinetics_what_we_can_claim.md/.qmd       ← Model validation & honest caveats
├── assumptions_validation.md/.qmd           ← Model assumptions & diagnostics
└── [figures directory]/                     ← Generated PNG figures (Figs 1-4)
```

**Purpose:** Standalone kinetics analysis, fully reproducible  
**Status:** Complete and publication-ready  
**When to use:** Core methodology + results for Sections 2-3  

---

### **`integration/`** — NEW DIRECTORY (Publication Integration)
```
integration/
├── kinetics_integrated_narrative_skeleton.md    ← 5-part outline (you write here)
├── SECTION_1_MOTIVATION_BIOLOGY.md              ← Section 1 complete (~2000 words)
├── method_convergence_figure.qmd                ← Fig 5 (Quarto document)
├── RUN_CONVERGENCE_ANALYSIS.R                   ← Run to generate Fig 5 + summary
├── EXTRACT_TOP_GENES_R_SCRIPT.R                 ← Extract Supplementary Table S1
├── generate_convergence_figure_RUNNABLE.R       ← Alternative Fig 5 method
└── EXTRACTION_SUMMARY.md                        ← Integration data & guide
```

**Purpose:** Integration scaffolding + your manuscript sections  
**Status:** Scaffolding complete; ready for you to write  
**When to use:** Write Sections 2-5 using skeleton, generate figures, extract data  

---

## 📋 What's Where

| File | Location | Purpose | Status |
|------|----------|---------|--------|
| Kinetics methodology | kinetics_limited/de_proteomics_wheat.qmd | Core analysis | ✓ Ready |
| Model validation | kinetics_limited/kinetics_what_we_can_claim.md | Honest caveats | ✓ Ready |
| Assumptions | kinetics_limited/assumptions_validation.md | Diagnostics | ✓ Ready |
| **Motivation & Biology (Sec 1)** | **integration/SECTION_1_MOTIVATION_BIOLOGY.md** | **Publication section** | **✓ Complete** |
| **Paper outline (Sec 2-5)** | **integration/kinetics_integrated_narrative_skeleton.md** | **Your writing space** | **Ready** |
| Method Convergence Fig | integration/method_convergence_figure.qmd | Fig 5 integration | ✓ Ready to run |
| Convergence script | integration/RUN_CONVERGENCE_ANALYSIS.R | Generate Fig 5 | ✓ Ready to run |
| Gene extraction | integration/EXTRACT_TOP_GENES_R_SCRIPT.R | Supp Table S1 | ✓ Ready to run |

---

## 🚀 Publication Workflow

### **For Conference (Next Week)**

1. **Generate figures:**
   ```r
   # In R with kinetics results loaded:
   source("analysis/integration/RUN_CONVERGENCE_ANALYSIS.R")
   ```
   Output: `fig-method-convergence-cadenza.png` + summary

2. **Extract gene list:**
   ```r
   source("analysis/integration/EXTRACT_TOP_GENES_R_SCRIPT.R")
   ```
   Output: `Supplementary_Table_S1_Cadenza_Top100_Genes.csv`

3. **Present Fig 5 + top candidates** using the summary output

### **For JEB Submission (Post-Conference)**

1. **Write manuscript:**
   - Use `SECTION_1_MOTIVATION_BIOLOGY.md` as Section 1
   - Use `kinetics_integrated_narrative_skeleton.md` outline for Sections 2-5
   - Reference kinetics_limited/ for detailed methods/validation

2. **Assemble figures:**
   - Figs 1-4: From kinetics_limited/ outputs
   - Fig 5: From integration/RUN_CONVERGENCE_ANALYSIS.R
   - Supplementary: From integration/EXTRACT_TOP_GENES_R_SCRIPT.R

3. **Citations & methods:**
   - Cite de_proteomics_wheat.qmd for kinetics methodology
   - Cite kinetics_what_we_can_claim.md for model limitations
   - Cite assumptions_validation.md for diagnostics
   - Cite gene-distance and predictions analyses for convergence methods

---

## 📊 Key Files to Know

### **kinetics_limited/ — DO NOT EDIT (foundational)**
- `assumptions_validation.md` — Reference for Publication Strategy section
- `kinetics_what_we_can_claim.md` — Use for limitations/caveats section
- `de_proteomics_wheat.qmd` — Run to generate kinetics Figures 1-4

### **integration/ — YOUR WRITING SPACE**
- `SECTION_1_MOTIVATION_BIOLOGY.md` — Copy to your manuscript
- `kinetics_integrated_narrative_skeleton.md` — Edit to write Sections 2-5
- `RUN_CONVERGENCE_ANALYSIS.R` — Run locally to generate Fig 5

---

## ✅ Checklist

### Before Conference:
- [ ] Run `RUN_CONVERGENCE_ANALYSIS.R` → get Fig 5 + candidate genes
- [ ] Run `EXTRACT_TOP_GENES_R_SCRIPT.R` → get gene table
- [ ] Prepare slides with Fig 5 + top 10 candidates
- [ ] Have manuscript outline ready (from skeleton)

### Before JEB Submission:
- [ ] Write Section 2: Kinetics methodology (reference de_proteomics_wheat.qmd)
- [ ] Write Section 3: Model validation (reference kinetics_what_we_can_claim.md)
- [ ] Write Section 4: Post-transcriptional mechanisms (interpret Fig 5)
- [ ] Write Section 5: Conclusions
- [ ] Integrate all sections into final document
- [ ] Include all figures + supplementary tables
- [ ] Copy relevant text from kinetics_limited/ for methods/limitations

---

## 🎯 Summary

**kinetics_limited/** = Foundational kinetics analysis (keep separate, reference it)  
**integration/** = Your manuscript integration workspace (write here, run scripts here)

**Next step:** Write Sections 2-5 in `kinetics_integrated_narrative_skeleton.md`, then compile with Section 1 into final manuscript.

---

**Status:** ✅ Organization complete  
**Date:** 2026-08-29  
**Ready for:** Conference presentation + publication writing
