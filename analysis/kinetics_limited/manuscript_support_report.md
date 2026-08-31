# Manuscript support report — kinetics analysis

**Scope.** What the kinetics analysis in `de_proteomics_wheat.qmd` can and cannot
support in a submitted manuscript: the wording to use, the revisions required
before submission, and the literature to ground it in.

**Provenance.** All numbers below are read from the render of
`de_proteomics_wheat.qmd` dated 2026-08-29 18:03, which is the **first render
that uses the validated standalone limma DE** (`norinXcadenza-shared`,
`DE-varieties-limma`) for the RNA treatment contrasts. Numbers in the current
notebook prose predate that change and are superseded — see §0.

---

## 0. What changed in this render, and what it invalidates

The RNA differential expression is now read from the project's standalone limma
workflow rather than re-derived inside the kinetics pipeline. This resolves the
§2b reconciliation failure by construction: there is one RNA DE result in the
project, not two that disagree 2–6×.

Three consequences that change what the manuscript can say:

**Universe B grew.** Cadenza 5,580 and Norin 5,159 matched genes (RNA + protein,
both QC-passing), up from ~4,000. The earlier figure was depressed by a gene-set
misalignment: `cad_common` was built from the local QC'd gene list only and then
used to index the validated DE matrices, injecting all-NA rows. Fixed with a
three-way intersect.

**The Norin baseline problem is far worse than previously reported.** At 0 hpi,
when infection cannot yet have acted:

| Variety | DE at 0 hpi (FDR < 0.05) | Rate | Previously reported |
|:--------|-------------------------:|-----:|--------------------:|
| Cadenza | 432 / 59,931 | **0.72 %** | 0.4 % |
| Norin   | 10,499 / 59,502 | **17.64 %** | 5.8 % |

A 24× asymmetry. Every statement in the notebook prose citing "Cadenza 0.4 %,
Norin 5.8 %" must be updated. This strengthens the existing "weight Cadenza"
recommendation, but it also means **Norin's 0 hpi arm cannot be treated as a
clean baseline at all** — and the t0-centring step (`Rl <- Rl - Rl[,1]`) is
absorbing a large real signal in Norin, not a batch artefact. State this rather
than glossing it.

**The RNA-responsive gate changed definition.** The standalone workflow does not
export an omnibus F-statistic, so the gate is now the minimum BH-adjusted p
across the four timepoint contrasts, paired as before with |log2FC| > 0.5. A
min-of-adjusted-p is anti-conservative relative to an F-test, so the responsive
set is more permissive than in earlier renders. This is independent of the DE
change itself and must be disclosed in Methods. It is reversible: the F-test
gate can be taken from the local fit the same way `ctrl_time` now is.

---

## 1. Overall judgment

The analysis is coherent **as a design-cell-level RNA–protein integration
study**. It is not yet defensible **as a proteome-wide measurement of protein
turnover, or as a mechanistic account of post-transcriptional regulation**.

That distinction should drive the framing of the whole paper. What the
experiment measures well is whether transcript trajectories carry information
about protein trajectories across an infection time course in two wheat
varieties with unpaired samples. What it cannot measure is *which*
post-transcriptional process is responsible, or a reliable per-gene half-life.

---

## 2. The four defensible conclusions

Each is supported by the analysis's own diagnostics.

**C1 — RNA trajectories predict held-out protein trajectories.**
Leave-one-design-cell-out Q², against a permutation null:

| Variety | Q² (comp 1) | perm q95 | p |
|:--------|------------:|---------:|---:|
| Cadenza | **0.269** | 0.095 | 0.015 |
| Norin   | **0.084** | 0.077 | 0.045 |

This is the most trustworthy number in the analysis. It is computed on design
cells held out of the fit, so it cannot be inflated by the in-sample artefact
that makes raw cross-block correlation meaningless here (a PLS on 8
pseudo-samples reaches r ≈ 0.72 on *permuted* data). The claim is about
**between-cell** coupling; within-cell covariance is not estimable in a Case B
design by any method.

**C2 — A fixed first-order RNA-only model does not adequately describe many
protein trajectories.** Roughly 23 % of RNA-responsive genes reject the kinetic
null in both varieties (Cadenza 1,050 / 4,559; Norin 896 / 3,796). See §4 for
why this must **not** be worded as evidence for post-transcriptional regulation.

**C3 — The experiment cannot distinguish altered translation from altered
degradation.** Absolute `k_s` cancels algebraically in fold-change space and is
absorbed by the per-protein MS response factor. The amplitude parameter `a`
conflates the two. This is a structural property of the design, not a power
problem; no amount of extra depth or replication fixes it.

**C4 — Per-gene calls are exploratory.** RNA measurement error is not propagated
into the kinetic model; propagating it by parametric bootstrap flips 36 % of
Cadenza and 49 % of Norin regulated/not-regulated calls. Any named gene is a
hypothesis requiring orthogonal validation (PRM, western, qPCR).

---

## 3. Claim wording — what to write instead

| Current framing | Recommended manuscript wording |
|:----------------|:-------------------------------|
| RNA and protein responses are coupled | "RNA trajectories predicted held-out protein trajectories at the treatment-by-timepoint design-cell level." |
| Post-transcriptional regulation is present | "Protein responses showed departures from an RNA-predicted reference trajectory." |
| Protein half-lives can be estimated | "Half-life-like parameters were estimable for a subset of genes whose trajectories contained sufficient kinetic information." |
| Translation or degradation changed | **Do not claim this from the current design.** |
| Gene X is post-transcriptionally regulated | "Gene X is a candidate for regulation beyond transcript abundance." |
| The first-order model fits | "The first-order model was used as a reference model but showed systematic residual misspecification." |
| Cadenza and Norin externally validate one another | "The two varieties provide independent cross-variety consistency evidence." |

The recurring move is from *mechanism* to *departure from a reference model*.
That is the honest description of what an amplitude parameter fitted to
fold-change trajectories measures.

---

## 4. The central issue: C2 is over-claimed as written

The excess of regulated calls over the simulated null is **not yet specific
evidence for post-transcriptional regulation**, because the null omits features
that are demonstrably present in the real data:

- control-arm drift (measured — see below);
- structured trajectory-shape errors (A1 residuals are systematically non-zero
  per timepoint, in opposite directions per variety);
- RNA measurement uncertainty (accepted by the fitter, never used);
- correlated residuals across genes (residual PC1 absorbs 62–69 % of variance
  against ~25 % expected under independence);
- time-varying degradation;
- omitted transcription→translation delay and translational saturation.

Any of these produces departures from a fixed-rate, independent-error null
**without any post-transcriptional regulation being present**. The null is
therefore not the right comparator for the claim being made.

### The control-drift evidence, stated precisely

A3 (control arm at quasi-steady state) is measurably violated:

| Variety | n | median abs. drift | > 1 log2FC | > 2 log2FC |
|:--------|--:|------------------:|-----------:|-----------:|
| Cadenza | 4,300 | 0.67 | 32.0 % | 8.6 % |
| Norin   | 3,533 | 0.69 | 34.1 % | 11.7 % |

Whether that drift *drives* the regulated calls differs by variety. Adjusted
logistic regression of `regulated` on control drift, controlling for the gene's
own RNA response amplitude (odds ratio per 1 SD):

| Variety | term | OR | 95 % CI | p |
|:--------|:-----|---:|:--------|---:|
| Cadenza | control drift | 1.060 | 0.979–1.146 | 0.145 |
| Cadenza | RNA amplitude | 1.056 | 0.975–1.141 | 0.176 |
| Norin   | control drift | **1.316** | 1.204–1.439 | < 1e-15 |
| Norin   | RNA amplitude | 0.777 | 0.701–0.859 | 1.2e-06 |

**Cadenza is clean on this test; Norin is not.** In Norin there is a real ~30 %
increase in the odds of being called regulated per SD of control drift, in
exactly the direction A3 violation predicts. Two honest caveats: drift and
treatment amplitude are collinear by construction (Spearman 0.45–0.52), so the
two coefficients cannot be cleanly separated — note that the RNA-amplitude
coefficient itself flips direction between varieties, a symptom of precisely
that; and association is not causation, since a gene with an unstable control
arm may also be genuinely more regulated.

### Wording to use for C2

> The data show more departures from a fixed-rate, RNA-only reference trajectory
> than expected under an independent-error null. However, because the control arm
> drifts and the first-order model is misspecified, these departures are
> interpreted as **transcript–protein discordance** rather than as direct
> evidence for altered translation or degradation.

---

## 5. Required revisions, in priority order

**R1 — Resolve the DE discrepancy.** *Done in this render*, by adopting the
standalone workflow's contrasts (§0). Remaining: update all notebook prose
citing the old counts, and disclose the gate change in Methods. Until §0 is
propagated through the text, **do not report DE counts**.

**R2 — Make the null control-aware.** The ratio model `dp/dt = k_d(r − p)`
assumes a stable control trajectory; the measured drift violates it. Fit the
infected and control arms jointly, or add an explicit time-varying baseline
term. This is the single change that would most strengthen C2.
Supports: `Lee11`, `Mar23c`.

**R3 — Replace the analytic LRT with a parametric bootstrap.** The χ²
approximation is weak here because the half-life is selected on a grid (a
non-regular model-selection step), RNA is estimated rather than known, the four
timepoint contrasts are correlated, and protein variances are estimated. The
bootstrap must repeat the *full* fitting procedure — grid selection included —
and propagate RNA uncertainty, protein uncertainty, timepoint covariance, and
imputation. Supports: `Tch14`, `Teo17`.

**R4 — Use the full contrast covariance.** The four treatment contrasts are
built from shared fitted coefficients, so their errors are correlated. The
current weighting `W = 1/max(SE, 0.05)²` is diagonal and uses marginal SEs only.
This understates uncertainty in a way that biases the LRT toward rejection.

**R5 — Make cross-validation genuinely nested.** Feature selection, imputation,
scaling, and dimension reduction must all occur inside each held-out fold.
Relatedly, the pipeline documentation is inconsistent about whether imputed
values enter the kinetic analysis — resolve this explicitly in Methods, because
it changes how C1 should be read.

**R6 — Demote the Bayesian model to sensitivity analysis.** Fixing the half-life
at its profile-likelihood optimum yields amplitude intervals *conditional on
that half-life*, understating joint uncertainty. The hierarchical model flags
roughly twice as many genes as the LRT (Cadenza 1,913 vs 951 of 4,300), and that
gap is largely the removed uncertainty, not extra evidence. Report it as a
sensitivity analysis, and treat genes flagged by **both** methods as the
high-confidence set. Supports: `Kuc18`.

**R7 — Soften the MNAR language.** The negative abundance–missingness
correlation (Spearman −0.74 / −0.75) supports *abundance-dependent* missingness.
It does not establish that each completely-missing cell is MNAR.
Supports: `Ahl19`, `Lazar16`, `Val17`, `Tya16`.

---

## 6. Two results whose numbers must be corrected in the text

### Half-lives are less resolvable than the current text says

The notebook reports 23–27 % of fitted genes in the estimable band.
Cross-tabulating the `identifiable` flag against where each fit landed on the
half-life grid gives a harsher picture:

| | Cadenza (n = 5,580) | Norin (n = 5,159) |
|:--|--:|--:|
| Not identifiable | 2,732 (49.0 %) | 2,651 (51.4 %) |
| Identifiable, but pinned at the 12 h grid floor | 1,304 (23.4 %) | 908 (17.6 %) |
| **Identifiable and genuinely inside the resolvable band** | **392 (7.0 %)** | **643 (12.5 %)** |
| Fit landed above the 72 h identifiability bound | 2,376 (42.6 %) | 1,952 (37.8 %) |

A fit at the 12 h floor means "≤ 12 h", not "12 h" — it is censored, and
counting it as identifiable overstates resolution. **Only 7 % of Cadenza and
13 % of Norin genes have a half-life this design actually resolves.** State that,
not the 23–27 %. The fitted `t_half` distribution is a property of the sampling
window (0/24/48/72 h), not of the wheat proteome.

### Cross-variety agreement is real, modest, and stronger than the last render

Reciprocal-best-hit orthologs, both varieties analysed independently with no
shared samples:

- 1:1 RBH pairs: 73,971; testable in both Universe Bs: 2,505; RNA-responsive in
  both: 1,819.
- Continuous amplitude correlation: Spearman **0.04** — near zero.
- Categorical archetype agreement: χ² = **727.7**, df = 42, p < 2e-16,
  Cramér's V = **0.220** (up from 0.14 in the previous render).

Explain the near-zero continuous correlation rather than burying it: each
variety's per-gene amplitude is itself a noisy n = 3 point estimate, and two
noisy proxies for the same underlying quantity correlate far more weakly than
the true quantities do. The categorical call is more robust to that noise, and
it is the one carrying the signal. This is genuine external validation — the
varieties share no samples and were processed separately — but V = 0.22 is a
modest effect and should be described as such.

---

## 7. References to add

> **Note on these keys.** The keys below are reproduced exactly as supplied,
> with the supplied statement of what each supports. They have **not** been
> expanded into full citations here: doing so from memory risks attaching the
> wrong author, year, or journal to a claim in a submitted manuscript. Resolve
> each key against the reference manager and add the entries to
> `analysis/kinetics_limited/references.bib` before citing. That file already
> backs the `@schwanhausser2011` citation used in the A1 section.

### Dynamic kinetic modelling

| Key | Supports |
|:----|:---------|
| `Tch14` | Simple first-order ODEs; profile shape, noise, identifiability. Closest simple-ODE precedent for this model. |
| `Teo13` | PECA's decomposition of dynamic RNA- and protein-level regulation. |
| `Teo17` | PECAplus: change-point analysis, smoothing, uncertainty propagation. |
| `Kuc18` | External half-life information and time-varying half-life models. Precedent for R6 and for `concordance.half_life_prior_file`. |
| `Jov15` | Direct pathogen-response measurements using pulsed SILAC. **Highest priority** — shows directly why total RNA and protein abundance cannot separate synthesis from degradation without pulse labelling. Cite in support of C3. |
| `Lee11` | Growth arrest, dilution, and nonstationary control physiology altering protein trajectories. Cite in support of R2. |
| `Mar23c` | Flux analysis outside equilibrium assumptions; moving beyond fixed-rate models. |
| `Liu17` | Why direct synthesis measurements add information beyond RNA or ribosome occupancy. |
| `Vel14` | Biological versus technical RNA–protein discordance. |
| `Haj10` | Statistical assessment of transcript–protein concordance in plants. |

### Plant and crop evidence

| Key | Supports |
|:----|:---------|
| `Abu23b` | Direct plant-immune turnover measurements separating synthesis and degradation. **Highest priority** — also the source for the limitation that early regulatory changes occur before 24 h, which this design's first post-baseline timepoint cannot see. |
| `Li17b` | Plant protein half-lives from hours to months; dependence on tissue state. |
| `Ish15` | Stable-isotope measurement of plant protein synthesis and degradation. |
| `Cao21b` | Direct protein-turnover measurements in developing wheat grain. |
| `Bai21b` | Delayed protein responses during Arabidopsis seed germination. |
| `Gal13` | Selective translation and protein turnover in plants. |
| `Xu17` | Translational reprogramming during plant immune activation. |
| `Yoo19` | Translational efficiency changes during effector-triggered immunity. |
| `Cai21` | Ribosome-associated translation changes during fungal infection in cotton. |
| `Wan24d` | RNA-seq and ribosome profiling during rice virus infection. |
| `Kag21` | Translational control in *Fusarium graminearum* during wheat infection. |
| `Sta17` | Technical and biological sources of transcript–protein discordance in potato. |
| `Jaw17` | Transcript–protein comparison during barley fungal infection. |
| `Zha19` | Wheat stripe-rust transcriptome–proteome integration. |
| `Lem24` | Recent barley–*Ramularia* transcriptome–proteome integration. |
| `Moe12` | Pathogen-induced changes in polysome-associated plant mRNAs. |

> **Use the turnover papers cautiously.** Half-lives vary strongly with tissue,
> development stage, growth conditions, and treatment. `Li17b`, `Ish15` and
> `Cao21b` support **hierarchical or weakly informative priors**; they do not
> license transferring fixed gene-specific rates from Arabidopsis, or from wheat
> grain, to infected wheat leaves.

### Missingness and proteomics statistics

| Key | Supports |
|:----|:---------|
| `Ahl19` | Probabilistic dropout modelling without imputation. |
| `Wol22` | Peptide aggregation, variance moderation, complex proteomics designs. |
| `Val17` | Benchmarking of filtering and imputation methods. |
| `Chu21` | Normalisation and uncertainty issues in multi-species RNA-seq. |
| `Lazar16` | Mixed MAR/MNAR mechanisms in label-free proteomics. Cite in support of R7. |
| `Tya16` | Perseus-style downshifted-normal imputation — cite if that scheme is retained. |

---

## 8. Software and data resources

| Resource | Role in this manuscript |
|:---------|:------------------------|
| **PECAplus** | Dynamic-method comparator. Provides source, tutorials, example data, and a Perseus plugin. **Designed for paired dynamic omics**, so it is *not* a drop-in solution for this Case B design — cite as a comparator and say why it does not apply. |
| **proDA** | Recommended no-imputation sensitivity analysis for the protein DE: models intensity-dependent dropout directly. Addresses R7. |
| **prolfqua** | Peptide-to-protein aggregation, formula-based designs, variance moderation, protein- vs peptide-level model comparison. |
| **MSstats** | Worth using if peptide- or precursor-level data are available; supports label-free, DIA, DDA, targeted. |
| **ProteinTurnover** | For a future isotope-labelling arm; useful in Discussion for explaining what direct turnover measurement would add. |
| **PRIDE / ProteomeXchange** | Source of cereal, plant–pathogen, and isotope-labelling datasets for external priors or benchmarking. |
| **IWGSC / Ensembl Plants** | **Record the exact wheat assembly, annotation, and protein release used for RNA–protein mapping.** Currently unstated anywhere in the pipeline. IWGSC RefSeq v2.1 and its annotation are available; Ensembl Plants provides wheat gene/transcript/protein annotations. A reproducibility requirement, not optional. |

---

## 9. Suggested closing paragraph

> We integrated unpaired transcriptome and proteome measurements at the
> treatment-by-timepoint design-cell level. RNA trajectories predicted held-out
> protein trajectories better than a permutation null, indicating reproducible
> cross-omic association at the population level. However, protein trajectories
> frequently departed from a fixed-rate first-order RNA-predicted model, and
> these departures were sensitive to control-arm drift, missing-value handling,
> and uncertainty in the RNA trajectories. We therefore interpret the results as
> evidence for transcript–protein discordance and regulation beyond transcript
> abundance, rather than as direct estimates of protein turnover, translation,
> or degradation. Individual candidates should be treated as hypotheses
> requiring orthogonal validation.

---

## 10. Complete reference bundle

```
Tch14, Teo13, Teo17, Kuc18, Jov15, Lee11, Mar23c, Liu17, Vel14, Haj10,
Abu23b, Li17b, Ish15, Cao21b, Bai21b, Gal13, Xu17, Yoo19, Cai21, Wan24d,
Kag21, Sta17, Jaw17, Zha19, Lem24, Moe12,
Ahl19, Wol22, Val17, Chu21, Lazar16, Tya16
```

32 keys. No unrelated general multi-omics papers.

---

## 11. Pre-submission checklist

- [ ] **R1** — propagate §0 through all notebook prose; remove every
      "0.4 % / 5.8 %" and "2–6× discrepancy" statement; disclose the
      RNA-responsive gate change in Methods.
- [ ] **R2** — control-aware null (joint arms, or time-varying baseline term).
- [ ] **R3** — parametric bootstrap replacing the analytic LRT.
- [ ] **R4** — full contrast covariance in the weighting matrix.
- [ ] **R5** — nested CV; resolve the imputation-in-kinetics ambiguity.
- [ ] **R6** — Bayesian model demoted to sensitivity analysis.
- [ ] **R7** — MNAR language softened.
- [ ] Half-life resolution restated as 7 % / 13 %, not 23–27 % (§6).
- [ ] Cross-variety Cramér's V updated to 0.220 (§6).
- [ ] Norin 0 hpi rate updated to 17.6 %; consequences for t0-centring stated (§0).
- [ ] Wheat assembly / annotation / protein release recorded (§8).
- [ ] All 32 keys resolved into `references.bib` (§7).
