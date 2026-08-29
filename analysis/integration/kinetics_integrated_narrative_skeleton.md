# Integrated Kinetics Analysis: RNA-Protein Divergence in Wheat Infection

## Structure for Publication (Biology-First Narrative)

### **PART 1: MOTIVATION & BIOLOGY**

#### 1.1 Background
- Wheat (*Triticum aestivum*) varieties Cadenza (susceptible) and Norin (resistant) infected with PN143
- Kinetics allow us to ask: does protein response lag RNA, and by how much?
- Post-transcriptional regulation (PTR) may explain defense strategy differences

**Figure/Table**: Table 1 — Experimental design (varieties, timepoints 0-72h, replicates)

---

#### 1.2 RNA and Protein Respond Differently
- **Finding:** RNA shows robust infection response (424–5745 DE genes by timepoint)
- **Finding:** Protein response is more muted (6–532 DE proteins)
- **Implication:** Post-transcriptional mechanisms are at play

**Figure**: Fig 1A — Heatmap: RNA logFC and protein logFC side-by-side, 0-72h  
**Figure**: Fig 1B — Scatter: RNA amplitude (abs logFC) vs protein amplitude by condition  
**Rationale**: Shows the biological divergence clearly to plant scientists

---

### **PART 2: KINETICS METHODOLOGY & FINDINGS**

#### 2.1 Fitting Kinetics: Why and How
- First-order model: dP/dt = ks·R(t) − kd·P(t)
- Estimates post-transcriptional lag, protein half-life per gene
- Models regulation as protein response departure from RNA trajectory (amplitude *a*: if a ≠ 1, post-transcriptional regulation)

**Key message for audience:** Kinetics gives us a testable null hypothesis (protein predicts from RNA) and tells us which genes break it.

---

#### 2.2 Main Findings: Population-Level Post-Transcriptional Regulation
- **Robust claim:** ~66–75% of genes show post-transcriptional regulation signal (excess over no-regulation null)
- **Robust claim:** Cadenza shows stronger regulation signature than Norin (but see caveats)
- **Safe interpretation:** Infection-response proteome is actively regulated post-transcriptionally, not just passively following RNA

**Figure**: Fig 2 — Amplitude distribution (a parameter) across Cadenza vs Norin  
**Figure**: Fig 3 — Kinetic archetypes pie chart: % regulated, buffered, lag, unchanged  
**Rationale**: Communicates aggregate findings that are robust despite model limitations

---

### **PART 3: WHAT CAN & CANNOT BE CLAIMED**

#### 3.1 Honest Assessment of Model & Data Limitations
**Key Limitation 1 (SHOWSTOPPER - B4):** RNA measurement error flips 36–49% of individual calls  
- Per-gene "regulated/not" lists are unstable  
- **Consequence:** Individual targets require orthogonal validation

**Key Limitation 2 (SHOWSTOPPER - A1):** Model shape mismatch  
- First-order kinetics don't capture true trajectory shape  
- Opposite bias direction between Cadenza and Norin  
- **Consequence:** Cannot compare kinetic *mechanisms* across varieties (only regulation extent)

**Key Limitation 3 (MANAGEABLE - Identifiability):** Only 23–27% of genes have resolvable half-lives  
- 50% above detection window, 24% pinned at floor  
- **Consequence:** Wheat proteome half-life *distribution* cannot be estimated

**Figure**: Fig 4 — Identifiability plot: which genes are resolvable vs unidentifiable  
**Rationale**: Transparency builds trust with reviewers

---

#### 3.2 Hypothesis-Generating Targets
- Individual gene lists framed as **candidates for validation**, not results
- "These 500 genes show kinetic divergence signals; further work needed"
- **Validation strategy:** qPCR, Northern, or independent proteomics cohort

**Table**: Supplementary Table S1 — Top 100 genes by regulation amplitude (with caveats)

---

### **PART 4: CANDIDATE POST-TRANSCRIPTIONAL MECHANISMS**

#### 4.1 Post-Transcriptional Modifications (PTM) Hypothesis
- Genes diverging in kinetics may show:
  - Translational control (selective protein synthesis)
  - Proteolysis regulation (differential half-life)
  - Protein localization (sequestration, export)
  
- **Approach:** Combine kinetics (here) + distance (RNA-protein trajectory concordance) + ML (predictability) to rank confidence

**Figure**: Fig 5 — Scatter: kinetics amplitude vs distance score, colored by ML confidence  
**Rationale**: Shows data integration without overselling

---

### **PART 5: CONCLUSIONS & NEXT STEPS**

- Post-transcriptional regulation is a major layer in wheat defense
- Individual genes require validation; population signal is robust
- Distance + ML + kinetics converge on candidate PTM regulators
- Functional studies needed on validated targets (mutants, proteolysis assays)

---

## Figure Selection & Generation Strategy

### **Main Figures (Publication Quality)**

| Figure | Content | Source | Status |
|---|---|---|---|
| **Fig 1A** | RNA/Protein heatmaps (Cadenza focus, Norin supplementary) | de_proteomics_wheat.md | ✓ Exists? |
| **Fig 1B** | Scatter: RNA vs Protein amplitude | de_proteomics_wheat.md | ✓ Exists? |
| **Fig 2** | Amplitude (*a*) distribution: Cadenza vs Norin | de_proteomics_wheat.md | ✓ Exists? |
| **Fig 3** | Kinetic archetype proportions (pie/bar chart) | de_proteomics_wheat.md | ✓ Exists? |
| **Fig 4** | Identifiability: resolvable vs unresolvable genes | assumptions_validation.md | ✓ fig-identifiability-1.png |
| **Fig 5** | Method convergence: kinetics + distance + ML | kinetics_what_we_can_claim.md | ? Needs creation |

### **Supplementary Figures**

| Figure | Content | Source | Status |
|---|---|---|---|
| **Supp Fig 1** | A1 residuals by timepoint (model misfit) | assumptions_validation.md | ✓ fig-a1-residuals-1.png |
| **Supp Fig 2** | A0 baseline asymmetry (Cadenza vs Norin) | assumptions_validation.md | ✓ fig-a0-baseline-1.png |
| **Supp Fig 3** | B4 RNA uncertainty (bootstrap flip rate) | assumptions_validation.md | ✓ fig-b4-1.png |
| **Supp Fig 4** | MNAR dominance test | assumptions_validation.md | ✓ fig-mnar-1.png |

---

## Writing Workflow

### **Phase 1: Scaffold (this document)**
- [x] Narrative structure
- [x] Figure strategy  
- [x] Caveats placement

### **Phase 2: Refine Figures**
- [ ] Check existing figures for publication quality
- [ ] Generate missing figures (Fig 5)
- [ ] Add captions
- [ ] Verify color-blind accessible palettes

### **Phase 3: Write Narrative**
- [ ] Section 1 (Motivation): ~500 words
- [ ] Section 2 (Methodology & Findings): ~800 words
- [ ] Section 3 (Limitations & Claims): ~600 words
- [ ] Section 4 (Mechanisms): ~300 words
- [ ] Section 5 (Conclusions): ~200 words

### **Phase 4: Integrate with Assumptions**
- [ ] Cross-link to Publication Strategy section
- [ ] Embed assumption summaries as text boxes
- [ ] Add caution statements throughout

---

## Key Messages for Plant Scientists (Cadenza-Focused, Norin Provisional)

1. **"Infection-triggered proteome is actively post-transcriptionally regulated"** ← aggregate, robust
2. **"500+ genes show kinetic divergence signals; candidates for functional follow-up"** ← hypothesis-generating
3. **"Protein half-lives vary 12–72h window; slower movers evade detection"** ← honest about identifiability
4. **"Model assumptions tested; per-gene calls require orthogonal validation"** ← methodological rigor

