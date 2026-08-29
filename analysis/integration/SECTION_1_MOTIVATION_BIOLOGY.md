# SECTION 1: Motivation & Biology
## Understanding Post-Transcriptional Regulation in Wheat Infection Response

### Introduction

Fungal pathogens trigger rapid reprogramming of plant gene expression, but the relationship between **transcript abundance and protein abundance** during infection remains poorly understood. While transcriptomics has revealed which genes respond to pathogen challenge, proteomics shows that protein levels often lag behind transcript changes—sometimes dramatically. This mismatch between RNA and protein trajectories suggests that **post-transcriptional mechanisms** (translational control, protein degradation, sequestration) actively shape the infection response.

To understand these mechanisms, we need tools that quantify the *magnitude and kinetics* of RNA-protein divergence. Here, we apply Bayesian kinetic modeling to paired time-course transcriptome and proteome data from two wheat (*Triticum aestivum*) varieties infected with the fungal pathogen *Magnaporthe oryzae* strain PN143. By fitting first-order kinetic models to protein trajectories predicted from RNA, we identify genes whose proteins depart significantly from RNA-expected trajectories, revealing which genes are likely subject to active post-transcriptional regulation.

---

### Biological Context: Wheat Varieties and Infection

**Wheat varieties:**
- **Cadenza** (susceptible variety): Develops strong disease symptoms; used here as primary analytical focus
- **Norin** (resistant variety): Shows enhanced defense responses; included for comparison

**Infection model:**
- Pathogen: *Magnaporthe oryzae* strain PN143 (causal agent of rice blast; causes systemic necrosis in wheat)
- Inoculation: Point infection at 3-week-old flag leaf
- Sampling: Four timepoints spanning 72 hours post-infection (hpi)
  - 0 hpi (uninfected control)
  - 24 hpi (early response)
  - 48 hpi (peak response)
  - 72 hpi (late response)
- Replication: Three biological replicates per treatment × timepoint combination

This design captures the full trajectory of an active infection response in a susceptible host, providing the temporal resolution needed to distinguish transcriptional lag from post-transcriptional regulation.

---

### The RNA-Protein Divergence Puzzle

**What transcriptomics shows:**
Infection triggers a robust transcriptional response in both varieties:
- Cadenza: 305–17,846 differentially expressed (DE) genes across timepoints (using stringent thresholds: FDR < 0.01, |logFC| > 1)
- Norin: 5,745–15,102 DE genes (approximately 2–4× higher response amplitude than Cadenza)

This suggests Norin mounts a more vigorous transcriptional defense, consistent with its resistant phenotype.

**What proteomics shows:**
The infection-triggered proteome response is remarkably *muted* relative to transcriptomics:
- Cadenza: 6–532 DE proteins across timepoints
- Norin: 6–86 DE proteins
- **Key observation:** Protein DE is **10–100× smaller than RNA DE**, despite both layers measuring the same biological samples

**The interpretive challenge:**
This RNA-protein mismatch could reflect three non-exclusive mechanisms:

1. **Translational buffering:** Protein synthesis rates are kept low despite high RNA levels (selective translation of defense genes, while bulk translation is suppressed)

2. **Protein stability differences:** Some gene products are rapidly degraded (high turnover), while others accumulate (slow degradation), independent of transcript levels

3. **Kinetic lag:** Protein synthesis and degradation are slower processes than transcription, so protein responses lag RNA responses by hours

**Why kinetics matter:**
By measuring protein trajectories over 72 hours, we can distinguish these scenarios. A gene whose protein trajectory **lags** its RNA trajectory (delayed response but similar magnitude) points to kinetic lag. A gene whose protein trajectory **diverges in amplitude** (protein much lower than RNA predicts) points to translational buffering or rapid degradation. Genes whose proteins remain **stable despite RNA fluctuations** suggest selective protein stabilization.

---

### Hypothesis: Post-Transcriptional Regulation is a Major Defense Layer

We hypothesize that **post-transcriptional mechanisms are actively employed** to fine-tune the infection response in wheat. Specifically:

- **Defense genes** (kinases, transcription factors, secondary metabolism enzymes) may be selectively translated and stabilized to maximize impact per transcript
- **Metabolic genes** (biosynthesis, ribosomal proteins) may be buffered (kept at moderate levels despite high transcripts) to conserve resources during stress
- **Regulatory proteins** (degradation and localization signals) may be rapidly turned over to enable dynamic response to changing pathogen signals

This hypothesis predicts that many genes will show kinetic divergence (RNA-protein decoupling), and these genes will be enriched for functions suggesting active regulation (kinases, TFs, protein turnover machinery).

---

### Kinetics as a Discovery Tool

**First-order kinetic model:**

For each gene, we fit:
$$\frac{dP}{dt} = k_s \cdot R(t) - k_d \cdot P(t)$$

Where:
- $P(t)$ = protein abundance at time $t$
- $R(t)$ = RNA abundance at time $t$ (measured; treated as known)
- $k_s$ = protein synthesis rate (inferred)
- $k_d$ = protein degradation rate (inferred)

**Key insight:** If this model fits well, then protein trajectories are *explained* by RNA trajectories alone. If the model **fails** (residuals are systematic, or the fit requires amplitude adjustment $a \neq 1$), then **unmeasured post-transcriptional regulation** is present.

**Regulation detection:**
A gene is flagged as "regulated" when the data significantly favor a model where the protein response *amplitude* differs from the RNA-predicted response:
- Amplitude $a < 1$: Protein is buffered below RNA-expected level (translational inhibition or rapid degradation)
- Amplitude $a > 1$: Protein is amplified above RNA-expected level (enhanced translation or stabilization)
- Amplitude $a \approx 1$: Protein follows RNA (kinetic lag only; no post-transcriptional regulation)

**Population-level interpretation:**
By aggregating across thousands of genes, we can estimate what fraction of the proteome undergoes active regulation, independent of per-gene mechanism.

---

### Roadmap: Structure of This Analysis

This document presents:

1. **Data loading and quality control** — Confirm RNA and protein measurements are reliable and comparable across samples

2. **Kinetics fitting and archetype classification** — Fit the kinetic model to each gene, classify genes as "regulated," "buffered," "lagged," or "unchanged"

3. **Population-level findings** — Quantify how widespread post-transcriptional regulation is, and whether it differs between Cadenza and Norin

4. **Honest assessment of model limitations** — Validate assumptions, quantify uncertainties (especially RNA measurement error), and establish which claims are robust vs. conditional

5. **Hypothesis-generating candidate lists** — Identify individual genes showing strong regulation signals for functional follow-up

6. **Integration with complementary methods** — (Later) combine kinetics with distance metrics and machine learning to cross-validate regulation candidates

---

### Key Claims This Section Supports

By the end of this analysis, we will be able to claim:

✅ **Robust (population-level):**
- Post-transcriptional regulation affects a substantial minority (66–75% excess genes over no-regulation null) of the infection-response proteome
- Regulation is detectable at the aggregate level with high confidence, independent of individual gene uncertainties

⚠️ **Conditional (individual gene level):**
- Individual gene targets are hypothesis-generating and require orthogonal validation (qPCR, Northern, functional assays)
- Per-gene regulation calls are sensitive to RNA measurement noise and should not be published as definitive without uncertainty propagation

❌ **Not claimable:**
- Specific kinetic rate constants (*k_s*, *k_d*) for the proteome — identifiability constraints prevent reliable estimation for slow-turnover proteins
- Mechanistic comparisons between varieties — model misfit (systematic residual structure) differs between Cadenza and Norin in ways that confound biological interpretation

---

### Figure 1: Experimental Design & Overview

**[Figure 1A]** Heatmap layout:
- Rows: genes (sorted by regulation amplitude)
- Columns: treatment × timepoint cells (Cadenza left, Norin right for reference)
- Left block: RNA log₂FC (centered at T0 control baseline)
- Right block: Protein log₂FC (same scale for direct comparison)
- Color scale: blue (low) to red (high), white = no change

**Interpretation:** Visual scanning should show that protein heatmaps are "smoother" and show less dynamic range than RNA heatmaps, supporting the RNA-protein mismatch observation.

**[Figure 1B]** Scatter plot:
- X-axis: RNA amplitude (abs log₂FC, treatment effect)
- Y-axis: Protein amplitude (abs log₂FC, treatment effect)
- Points: individual genes
- Color: by regulation archetype (regulated = red, lagged = blue, buffered = orange, unchanged = gray)

**Interpretation:** Genes above the diagonal show protein amplification (rare); below show buffering (common). Spread of points indicates widespread decoupling.

---

### Next Steps

With this biological context established, we proceed to:
- Section 2: Describe kinetics methodology and report aggregate findings
- Section 3: Present model validation and honest caveats
- Section 4: Discuss candidate post-transcriptional mechanisms
- Section 5: Conclusions and implications for functional studies

