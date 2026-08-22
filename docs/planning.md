# Analysis planning — Pillars 2–5

Original design reasoning for the IntegrateProtRNA analysis, written before any
code existed. Kept as a reference for *why* the pipeline in `R/` and the
summary in [`PIPELINE.md`](PIPELINE.md) are built the way they are. Pillars
2–3 describe Case A (sample-matched) tooling as background; once the design
was confirmed as **Case B** (unmatched RNA/protein samples), the constraints
in `PIPELINE.md` supersede the sample-level integration tools described here
(MOFA+/MEFISTO/DIABLO/O2PLS at the sample level are not valid under Case B —
see `PIPELINE.md` for what replaces them). Pillar 4's custom-model proposals
are a superset of what has been implemented; only 4A (the kinetic model) has
a built equivalent so far, in `R/05_concordance_archetypes.R`.

---

## Pillar 2 — Mapping and pre-processing

### 2.1 The key structural insight

For all latent-variable integration (MOFA+, MEFISTO, DIABLO, O2PLS), the
shared dimension is **samples, not features**. Views may have entirely
disjoint feature sets. So:

> Do **not** subset to the RNA∩protein intersection before running
> MOFA/DIABLO. That throws away ~70% of your transcriptome for no reason.

Feature-level matching is required only for the *pairwise* analyses
(archetypes, kinetics, lag). Keep two parallel feature universes:

- **Universe A (full):** all genes + all protein groups → factor models,
  network models, pathway integration.
- **Universe B (matched core):** 1:1 gene↔protein pairs → archetypes, ODE
  kinetics, delay estimation, concordance scatterplots.

### 2.2 ID mapping — do it in the protein→gene direction

Transcript-level → gene-level first (`tximport` with
`countsFromAbundance="lengthScaledTPM"` for DE, or scaled TPM). Then anchor
everything on **Ensembl gene ID** as the join key, and map proteins *onto*
it, not the reverse. Rationale: the ambiguity is concentrated on the MS side
(protein groups, shared peptides, isoforms), so resolve it there.

Pipeline:

1. Fix a genome/annotation build and a UniProt release. **Record both.**
   Version skew between your FASTA search database and your RNA annotation
   is the #1 silent source of mapping loss.
2. From MaxQuant `proteinGroups.txt` / DIA-NN `report.pg_matrix.tsv`: drop
   `Reverse`, `Potential contaminant`, `Only identified by site`. Require ≥1
   (preferably ≥2) unique peptides.
3. Expand the semicolon-delimited `Majority protein IDs` into a full
   accession list per group. Do **not** silently take the leading razor
   protein — record it, but keep the list.
4. Map UniProt accessions → Ensembl gene via `UniProt.ws` / the UniProt
   ID-mapping REST API, plus `biomaRt` as a cross-check. Prefer UniProt's
   own cross-references over biomaRt where they disagree (biomaRt's UniProt
   links are lossier). Handle deprecated/merged accessions via UniProt's
   secondary-accession table.

### 2.3 1-to-many, many-to-1, many-to-many — the bipartite-graph solution

Build a bipartite graph: gene nodes ↔ protein-group nodes, edges = mapping
evidence. Take **connected components**. Then:

- Component with 1 gene, 1 group → **clean pair**, goes into Universe B.
  (Expect 60–80%.)
- Multiple transcripts/isoforms → 1 group → collapse to gene level (you're
  gene-level on the RNA side anyway). Clean pair.
- 1 gene → multiple groups (isoform-resolved, differential PTM forms) →
  either sum/median the groups, or keep separate and analyse as
  isoform-level. Decide by whether unique peptides distinguish them.
- Multiple genes ↔ multiple groups (paralogue families sharing peptides:
  tubulins, histones, HLA, keratins) → **do not force a pairing.** Collapse
  the whole component into a single "feature cluster" (sum RNA, sum
  protein) and flag it, or exclude from Universe B and analyse only in
  Universe A. Report the count of such components in the methods.

Also handle **homologous / paralogous gene assignment** explicitly: if
you're crossing species or working in a non-model organism, run
OrthoFinder/eggNOG-mapper once, freeze the orthogroup table, and treat
orthogroups as the join key rather than gene symbols. Gene symbols are the
worst possible join key — non-unique, non-stable, and Excel-corruptible
(`SEPT7`, `MARCH1`, `1-Mar`). Never join on symbols.

### 2.4 Non-overlapping features — a 3-tier strategy, not deletion

1. **Tier 1 (matched):** full cross-omics analysis.
2. **Tier 2 (RNA-only):** genes with no protein evidence. Analyse
   univariately; then ask *why* no protein — low abundance (check RNA level
   distribution vs matched set), secreted/membrane, small/lysine-poor (few
   tryptic peptides). Distinguish "not expressed" from "not detectable" —
   this is a technical censoring problem, and conflating the two produces
   false "post-transcriptional repression" claims. Model detectability with
   a logistic regression on peptide-count/hydrophobicity/length covariates
   and report it.
3. **Tier 3 (protein-only):** rarer; often mapping failures. Re-inspect
   before interpreting.

Reintegrate all three tiers at the **pathway level** — a pathway can be
supported by RNA-only, protein-only, and matched members. `ReactomeGSA`
supports genuine multi-omics comparative pathway analysis; `PathIntegrate`
(Ebbels lab) does pathway-level multi-omics with PLS and is a good fit here.

### 2.5 Missing values — the decision that will get scrutinised

Proteomics missingness is a mixture: **MNAR** (left-censored, below
detection) and **MAR** (stochastic MS/MS sampling, ID transfer failure). DIA
data is more MAR than DDA. Diagnose before choosing:

- Plot mean intensity vs. missingness rate per protein. A strong negative
  slope ⇒ MNAR-dominant.
- Plot missingness pattern heatmap by run order/batch ⇒ technical MAR.

Then, in order of preference:

1. **Filter first, hard.** Keep a protein if it has ≥2 (I'd use ≥2 of 3)
   valid values in **at least one** treatment×time group. This "at least n
   in one group" rule preserves on/off biology that a global "≥50% valid"
   rule deletes. Do this *before* any imputation.
2. **Avoid imputation where the model can handle NAs.** MOFA+/MEFISTO
   marginalise missing values in the ELBO natively — a genuine advantage
   and a reason to prefer them here. `MSstats` and `limma` handle unbalanced
   missingness within the linear model. `msImpute` explicitly distinguishes
   MAR/MNAR modes.
3. **If you must impute:** mixed strategy — MNAR proteins (missing in all
   reps of a group, present in others) get left-censored imputation
   (`QRILC`, `MinProb`, or Perseus-style down-shifted normal: width 0.3,
   downshift 1.8 SD); MAR proteins get kNN or `missForest`. `imputeLCMD` and
   the `DEP` package implement this split. Benchmark alternatives with
   `NAguideR`.
4. **Sensitivity analysis is non-negotiable.** Re-run your headline results
   under ≥3 imputation schemes (including complete-case) and show in a
   supplementary figure that conclusions are stable. Reviewers of
   proteomics manuscripts ask for this reliably.

Critical trap: imputing MNAR values with a low constant then computing
RNA–protein correlations **manufactures anti-correlation** and false
"post-transcriptional repression." Exclude imputed values from the
correlation/kinetics analyses (Universe B) even if you keep them for factor
models.

### 2.6 Normalisation and scaling for integration

- **RNA:** `DESeq2::vst()` (or `limma::voom` if you stay in the limma
  world). You need a homoscedastic continuous matrix for factor models — do
  *not* feed raw counts or plain logCPM to MOFA's Gaussian likelihood.
  (MOFA+ does offer a Poisson likelihood, but VST + Gaussian is better
  behaved at n=24.)
- **Protein:** log₂, then `vsn` or median/quantile normalisation. Check for
  MS batch/run-order drift; correct with `limma::removeBatchEffect` passing
  your full design matrix (`~treatment*time`) as `design=` so you don't
  regress out biology. Do not run ComBat without the covariate-preserving
  design.
- **Per-view scaling:** center features, then scale each *view* to unit
  total variance (`MOFA2::prepare_mofa(scale_views=TRUE)`). Without this,
  the higher-variance block (usually RNA, ~15k features) dominates every
  factor.
- **Never center or scale within timepoint.** That erases the temporal
  signal you're trying to model. Surprisingly common error.
- **Feature pre-selection:** for RNA, top 3,000–5,000 by variance (after
  VST) for factor models; keep the full set for univariate DE and pathway
  analysis. Report the selection and show factor robustness to the cutoff.

---

## Pillar 3 — Published tools, and why each fits

### 3.1 MEFISTO — the primary integration engine (Case A)

MEFISTO (Velten et al., *Nat Methods* 2022) extends MOFA+ by placing
**Gaussian-process priors on the latent factors over a continuous
covariate** — here, time.

- **Model:** $Y^{(m)} \approx Z W_m^\top + \epsilon^{(m)}$, with views
  $m \in \{\text{RNA}, \text{protein}\}$, $Z \in \mathbb{R}^{24 \times K}$
  latent factors. MOFA+ adds a two-level sparsity prior (ARD across factors
  per view + spike-and-slab across features) so each factor is
  automatically labelled as view-shared or view-specific. MEFISTO adds
  $z_k(t) \sim \mathcal{GP}(0, K_k(t,t'))$ with a squared-exponential
  kernel, and learns a per-factor **smoothness parameter** $s_k \in [0,1]$
  — factor $k$ is temporally structured if $s_k \to 1$, and
  unstructured/noise-like if $s_k \to 0$.
- **Why it fits:** (i) it handles missing values natively — no imputation
  needed for the factor model; (ii) it explicitly decomposes variance into
  *shared RNA+protein* vs *RNA-only* vs *protein-only* factors, which is a
  direct quantitative answer to Pillar 1's decoupling question; (iii)
  treating the two treatments as MEFISTO **groups** with a shared GP lets it
  learn factor trajectories per treatment on a common latent axis, and its
  **alignment/warping** feature can detect that treatment shifts the
  response in *time*, not just amplitude — i.e. delayed response detection
  at the factor level; (iv) unsupervised, so it can't overfit to the
  treatment label.
- **Inference:** variational Bayes maximising the ELBO, with sparse GP
  (inducing points) — though with T=4 you'll use exact GP.
- **Caveat to state:** 4 timepoints is thin for GP hyperparameter
  estimation. Fix or tightly prior the lengthscale, and check factor
  stability across random seeds (MOFA is non-convex; run ≥10 seeds, keep the
  best ELBO, and report factor correlation across seeds).
- **Case B note:** requires sample-matched rows. Superseded in this project
  by the design-cell PLS + leave-one-cell-out Q² approach in
  `R/06_integration_caseB.R` — see `PIPELINE.md`.

### 3.2 DIABLO (mixOmics) — supervised signature (Case A)

- **Model:** sparse generalised canonical correlation analysis
  (sGCCA/RGCCA). Maximises $\sum_{i,j} c_{ij}\,\text{cov}(X_i a_i, X_j a_j)$
  subject to $\ell_1$ penalties on the loadings $a_i$, where $C$ is a
  user-specified **design matrix** encoding how strongly you want each pair
  of blocks correlated vs. how strongly you want them to discriminate the
  outcome.
- **Why it fits:** you have a designed experiment with a known label. DIABLO
  gives a compact, sparse, cross-validated multi-omics **signature**
  discriminating treatment — exactly what a manuscript needs as a "here is
  the integrated biomarker panel" figure, plus the circos/network plots
  that make cross-omics correlation legible.
- **Handling time:** DIABLO has no temporal model. Two valid uses: (a)
  outcome = 8-level `treatment:time` factor; (b) fit DIABLO separately per
  timepoint and compare signature drift across time — often more
  interpretable.
- **Caveat:** supervised + n=24 ⇒ overfitting risk is real. Use `perf()`
  with repeated stratified M-fold CV (M=3 or 4, ≥50 repeats — LOOCV is
  high-variance here), tune `keepX` by CV, and **always report a
  permutation null** (permute treatment labels, re-run the entire tuning +
  CV, compare error rate distributions). A DIABLO signature without a
  permutation null will not survive review.
- **Case B note:** also requires sample-matched rows.

### 3.3 timeOmics — trajectory clustering across blocks

Designed precisely for longitudinal multi-omics. Three stages: (1)
per-feature **LMM-spline** modelling (`lmms`) to get smooth profiles and
filter out features whose "trajectory" is a straight line/noise
(Breusch–Godfrey test); (2) profile scaling; (3) **multiblock sPLS**
clustering of RNA and protein features into shared trajectory clusters, with
cluster membership by loading sign.

- **Why it fits:** it directly produces "these genes and these proteins move
  together over time" — the clusters are the integration.
- **Caveat:** with T=4 you are at the floor for spline modelling. Restrict
  to linear or `df=2` natural splines; do not let `lmms` fit cubic splines
  through 4 points. In Case B (unmatched), this is the strongest published
  fit, since it works on fitted profiles rather than matched samples — not
  yet implemented in `R/`.

### 3.4 O2PLS — the right choice for exactly two blocks

`OmicsPLS` (Bouhaddani et al.). Decomposes each block into **joint**,
**orthogonal** (block-specific), and **residual** variation:

$$X = T_{joint}W^\top + T^\perp_X P_X^\top + E, \qquad Y = U_{joint}C^\top + U^\perp_Y P_Y^\top + F$$

Symmetric (unlike PLS, no X→Y asymmetry), SVD-based, fast, and gives a clean
scalar: *what fraction of transcriptome variance is shared with the
proteome, and vice versa?* That number, computed per timepoint, is a
compelling quantification of temporal decoupling. Two blocks is its exact
design point. (Case A; the design-cell PLS in `R/06_integration_caseB.R` is
the Case B analogue of the "shared vs specific variance" question this
answers.)

### 3.5 Univariate temporal baselines (reviewers expect these)

- **RNA:** `DESeq2` LRT, `~treatment*ns(time, df=2)` vs
  `~treatment + ns(time, df=2)`, testing the interaction — i.e. "does the
  temporal trajectory differ by treatment?" Alternatively treat time as an
  ordered factor with 4 levels and test the interaction; with T=4 and n=3
  this is often more robust and more interpretable than splines.
  `maSigPro` is the classic two-step polynomial regression alternative
  purpose-built for multi-series time-course RNA-seq.
- **Both layers:** `limma` with the same design (`voom` for RNA,
  `eBayes(trend=TRUE, robust=TRUE)` for protein) gives a single consistent
  framework across omics — a real advantage for comparability. **This is
  the approach implemented** in `R/04_univariate_temporal_de.R`.
- **`ImpulseDE2`:** models impulse-shaped (transient) responses in
  case-control time series. Fits this design well if the treatment response
  is transient; it explicitly classifies genes as transient vs. monotonic,
  which feeds directly into the archetype table.
- **`MSstats`** for the proteomics side if you have peptide/feature-level
  data — its linear mixed model handles missingness and technical
  replication properly, better than protein-matrix-level limma.

### 3.6 What to skip

`iClusterPlus`, `intNMF`, `SNF` — sample-clustering methods for discovering
subtypes. This is a designed experiment with known groups; there are no
subtypes to discover. `MINT` — cross-study integration, not this problem.

---

## Pillar 4 — Custom models

Design constraint restated: **n=24.** Any architecture whose parameters
scale with sample count will memorise. The proposals below each dodge this
differently, ordered by expected return on effort. **Status:** only 4A has a
built equivalent so far (`R/05_concordance_archetypes.R`'s kinetic null +
amplitude LRT). 4B–4E are proposals, not yet implemented — and a
substantially overlapping condition-aligned autoencoder + GNN-adjacent
codebase already exists on `nectar` (`rnaprot/`, see `PIPELINE.md`/project
memory), so 4B/4C should be evaluated against that before building new code.

### 4A — Hierarchical Bayesian kinetic state-space model *(highest value, lowest risk — implemented as a point-estimate profile-likelihood version; full Bayesian version not yet built)*

Because it *is* the biology, and it's the analysis MOFA cannot do.

**Model.** For gene $g$, treatment $c$, discretising the ODE:

$$P_{g,c}(t_{i+1}) = P_{g,c}(t_i) + \Delta t\left[k^{syn}_{g,c} R_{g,c}(t_i) - k^{deg}_{g,c} P_{g,c}(t_i)\right]$$

Observation model, working in log space for positivity:

$$\log P^{obs}_{g,c,r}(t_i) \sim \mathcal{N}\!\left(\log P_{g,c}(t_i),\ \sigma^2_g\right), \qquad \log R^{obs} \sim \text{NB or } \mathcal{N}(\log R_{g,c}(t_i), \tau^2_g)$$

with $R_{g,c}(t)$ itself a latent smooth function (GP or monotone spline)
fitted to the RNA observations — this handles RNA measurement noise properly
rather than plugging in point estimates.

**Partial pooling — the part that makes it identifiable:**

$$\log k^{syn}_{g,c} = \mu^{syn} + \alpha^{syn}_g + \delta^{syn}_g \cdot \mathbb{1}[c = \text{treated}], \quad \alpha^{syn}_g \sim \mathcal{N}(0, \sigma^2_{syn}), \quad \delta^{syn}_g \sim \text{Horseshoe}(\lambda)$$

The **horseshoe (or regularised horseshoe) prior** on the treatment effects
$\delta$ encodes "most genes' rate constants don't change; a few change a
lot." This is what buys identifiability from 4 timepoints — genes borrow
strength from each other for the baseline rates, and only
genuinely-perturbed genes escape shrinkage.

**Loss / objective:** negative log posterior; fit by NUTS in Stan or
NumPyro (NumPyro + GPU is much faster at this scale). Realistically: run on
1,000–3,000 well-quantified matched genes, not 15,000. Consider ADVI/SVI for
a first pass, then NUTS on the final gene set.

**Output → figures:** posterior distributions of $\delta^{syn}_g$ and
$\delta^{deg}_g$; a 2D "kinetic volcano" of Δsynthesis vs Δdegradation with
genes coloured by archetype; GO enrichment of the degradation-shifted set.
Priors on $k^{deg}$ can be informed by published half-life datasets to
further stabilise fits — implemented as the optional
`concordance.half_life_prior_file` in the current pipeline (validated: lifts
regulated-gene recall 27% → 72% on simulated data — see `PIPELINE.md`).

**Validation:** posterior predictive checks; simulation-based calibration;
leave-one-timepoint-out predictive log-likelihood vs. a null model with
$\delta \equiv 0$.

### 4B — Graph neural network with genes as observations *(not yet built; overlaps existing nectar code)*

The key idea: stop treating samples as observations. Treat **genes** as
observations. You then have ~3,000–10,000 training examples instead of 24,
and deep learning becomes legitimate.

**Setup.** Node = gene. Node feature vector = concatenation of its RNA
trajectory and protein trajectory across both conditions:
$x_g = [\Delta R_g(t_1..t_4),\ \Delta P_g(t_1..t_4),\ \text{covariates}] \in \mathbb{R}^{8+p}$,
where covariates include mean abundance, CDS/UTR length, codon adaptation
index, predicted miRNA target sites, peptide-detectability score. Edges from
**STRING** (high-confidence, ≥0.7) + a TF–target regulatory network + CORUM
complex co-membership, as a multi-relational graph (use R-GCN or a
heterogeneous GAT to keep edge types distinct).

**Task — self-supervised, mask-and-predict:** mask the protein trajectory of
a random 20% of nodes; predict it from RNA trajectory + graph neighbourhood.

$$\mathcal{L} = \underbrace{\frac{1}{|M|}\sum_{g \in M}\left\|\hat{P}_g - P_g\right\|^2}_{\text{masked reconstruction}} + \lambda_1 \underbrace{\sum_{(g,h)\in E} w_{gh}\|h_g - h_h\|^2}_{\text{graph Laplacian smoothness}} + \lambda_2 \|\Theta\|^2$$

**Why this is scientifically interesting, not just an exercise:**

1. The **residual** $P_g - \hat{P}_g^{\text{RNA-only baseline}}$, after
   conditioning on the network, *is* an estimate of post-transcriptional
   regulation not explained by mRNA or by pathway context. Genes with large
   residuals are post-transcriptional regulation candidates, now with a
   proper null.
2. Attention weights on edges (with GAT) reveal *which neighbours* carry
   predictive information — interpretable as regulatory or
   complex-membership dependence.
3. It **imputes protein trajectories for Tier-2 RNA-only genes**, extending
   proteome coverage in a principled, validatable way.

**Architecture:** 2–3 GAT layers, hidden 64–128, 4 heads, dropout 0.3–0.5,
layer norm, residual connections. Deeper than 3 layers over-smooths.
PyTorch Geometric or DGL.

**Validation (essential):** hold out genes, not samples. Compare against
three baselines — (i) mean protein trajectory, (ii) ridge regression on RNA
trajectory alone, (iii) same GNN on a **degree-preserving rewired graph**.
If the real graph doesn't beat the rewired graph, the network prior is
contributing nothing and that should be reported as such.

### 4C — Multi-view VAE with shared temporal latent dynamics *(not yet built; a condition-aligned autoencoder already exists in `rnaprot/unpaired.py` on nectar — evaluate that first)*

Two encoders → a shared latent $z$ with temporal structure → two decoders.

$$q(z \mid x^R, x^P) \propto p(z)\, q_R(z \mid x^R)\, q_P(z \mid x^P) \quad \text{(Product-of-Experts)}$$

PoE is the right fusion choice specifically because it handles **missing
views and missing features gracefully** — an absent protein measurement
just drops that expert's contribution, mirroring how totalVI/MultiVI handle
partial modalities.

**Loss:**

$$\mathcal{L} = \underbrace{\mathbb{E}_q[\log p_{NB}(x^R \mid z)]}_{\text{NB likelihood, RNA counts}} + \underbrace{\mathbb{E}_q[m \odot \log p_{\mathcal{N}}(x^P \mid z)]}_{\text{masked Gaussian, protein}} - \beta\, D_{KL}(q\|p) - \gamma \underbrace{\sum_{t}\|z_{t+1} - z_t\|^2}_{\text{temporal smoothness}} - \eta \underbrace{\mathcal{L}_{\text{align}}(z_R, z_P)}_{\text{cross-view alignment}}$$

- The **mask $m$** on the protein term means you never impute — missing
  entries simply contribute nothing to the loss. This is the cleanest
  possible answer to Pillar 2.5.
- $\mathcal{L}_{\text{align}}$: InfoNCE contrastive loss treating the RNA
  and protein views of the *same sample* as positives and other samples as
  negatives. With n=24 the negative pool is tiny; a soft-CCA or MMD
  alignment penalty is the safer alternative. (Note: under Case B there is
  no "same sample" across views — this term must instead be a
  condition-mean alignment, as `rnaprot/unpaired.py`'s
  `L_condition_latent_alignment` already implements.)
- Temporal prior: replace the i.i.d. $\mathcal{N}(0,I)$ prior with a **GP
  prior over $z(t)$** (i.e. a neural MEFISTO) or a latent neural ODE
  $\frac{dz}{dt} = f_\theta(z)$. The neural ODE version is elegant — it
  learns a continuous latent dynamical system and lets you
  extrapolate/interpolate between the 4 timepoints — but is the most
  data-hungry option.

**Regularisation for n=24 — mandatory, all of it:** latent dim 4–8; reduce
input to ~2,000 HVG + all quantified proteins; hidden layers ≤128 with
strong dropout; early stopping on leave-one-sample-out; heavy weight decay;
and **ensemble over ≥20 seeds**, reporting only latent structure stable
across seeds.

**Honest framing for the manuscript:** report the VAE as a complementary
nonlinear view, validated by (i) correlating its latent factors with
MEFISTO/PLS factors — if they agree, the nonlinear model found the same
biology and adds confidence; if they disagree, investigate whether the
extra structure is real nonlinearity or overfitting — and (ii) permutation
nulls with shuffled time labels.

### 4D — Bolt-on: learned per-gene translation delay *(not yet built)*

A small module worth adding to 4B or 4C. Predict protein at time $t_i$ as a
soft-attention-weighted combination of RNA at all earlier timepoints:

$$\hat{P}_g(t_i) = \sum_{j \le i} \alpha_{g,ij}\, W R_g(t_j), \qquad \alpha_{g,i\cdot} = \text{softmax}\big(f_\theta(R_g, \text{covariates})\big)$$

The learned attention distribution's centre of mass is an interpretable
**per-gene delay** $\tau_g$. Output: a genome-wide "translation delay atlas"
with functional enrichment of fast vs. slow genes. Small model, few
parameters, directly interpretable, novel figure.

### 4E — The pragmatic win: gene-level supervised learning with interpretable features *(not yet built)*

Don't overlook this. Observations = genes (thousands). Target = protein
log₂FC trajectory (or the archetype class). Features = RNA trajectory +
sequence/structural covariates (5′/3′UTR length, codon usage bias, GC,
miRNA site counts, RBP motifs, predicted half-life, complex membership,
abundance). Fit **XGBoost / random forest**, interpret with **SHAP**.

This answers "*what determines whether a gene is buffered?*" with a
rigorous, small-n-safe, fully interpretable model, and it often produces
the single most citable figure in a paper like this. Group-aware CV (split
by chromosome or gene family) to avoid leakage from paralogues.

---

## Pillar 5 — Publishable workflow and figures

### Pipeline

```
S0  Design audit          → confirm sample matching (Case A/B); confirm no repeated measures;
                            fix annotation + UniProt versions; power/limitations statement
S1  QC & normalisation    → RNA: tximport → filter → VST | Protein: contaminant/decoy filter,
                            ≥2 unique peptides, "≥2 valid in ≥1 group", log2 + vsn
                            → batch/run-order diagnostics → design-preserving correction
S2  ID mapping            → bipartite graph → connected components → Universe A (all) +
                            Universe B (clean 1:1 pairs); orthogroups if cross-species
S3  Missingness           → MAR/MNAR diagnosis → mixed imputation (Universe A only)
                            → 3-scheme sensitivity analysis
S4  Univariate baseline   → limma (both layers), ~treatment*time interaction; ImpulseDE2;
                            FDR via BH + IHW
S5  Concordance/archetype → per-timepoint RNA-protein logFC concordance; archetype
                            classification; lag-vs-regulation separation
S6  Latent integration    → MEFISTO (primary) + O2PLS (joint/specific variance) +
                            DIABLO (supervised signature, with permutation null)
S7  Trajectory clustering → timeOmics multiblock sPLS clusters
S8  Kinetic modelling     → hierarchical Bayesian ODE → Δk_syn, Δk_deg posteriors
S9  Custom ML             → GNN residual analysis / delay atlas / XGBoost+SHAP
S10 Functional layer      → ReactomeGSA or PathIntegrate multi-omics pathway analysis;
                            CORUM stoichiometry test; DGCA network rewiring
S11 Validation            → permutation nulls, imputation sensitivity, seed stability,
                            orthogonal validation (PRM/western/qPCR on 5-10 targets)
S12 Reproducibility       → Snakemake/Nextflow + renv/conda; GEO + PRIDE deposition
```

Under Case B, S6 is replaced by design-cell PLS + leave-one-cell-out Q² (see
`R/06_integration_caseB.R`, `PIPELINE.md`); S0–S5 map directly onto
`R/00`–`R/05`. S7–S12 are not yet implemented.

### Manuscript figures

**Fig 1 — Design & QC.** (a) design schematic; (b) PCA per omics, coloured
by treatment, shaped by time, with variance-explained; (c) sample–sample
correlation heatmaps; (d) UpSet plot of RNA/protein feature overlap after
mapping; (e) protein missingness: intensity-vs-missingness plot establishing
MNAR, plus the pattern heatmap. *(implemented: `fig01_qc.pdf`,
`fig01b_mapping.pdf`, `fig01c_missingness.pdf`)*

**Fig 2 — Temporal DE per layer.** (a) counts of DE features per
timepoint/contrast, both layers; (b) clustered heatmaps of temporal profiles
per omics; (c) representative trajectories; (d) overlap of RNA-DE vs
protein-DE genes — deliberately underwhelming, and that's the setup for Fig
3. *(implemented: `fig02_univariate_de.pdf`)*

**Fig 3 — Concordance & decoupling (the money figure).** (a) 2×4 grid of
ΔRNA vs ΔProtein scatterplots, one per timepoint per contrast, with ρ and a
fitted slope; (b) global coupling coefficient trajectory per treatment, with
CI — showing decoupling emerging over time; (c) stacked bar of archetype
proportions across timepoints; (d) exemplar trajectories for each
archetype; (e) GO enrichment per archetype. *(implemented, minus GO
enrichment: `fig03_concordance.pdf`)*

**Fig 4 — MEFISTO integration.** (a) variance decomposition matrix (factor
× view); (b) factor trajectories over time by treatment, GP posterior mean
± 95% CI, with smoothness score $s_k$ annotated; (c) top feature weights per
factor per view; (d) factor–pathway enrichment; (e) if warping is used, the
learned time alignment between treatments — a literal, quantitative "the
response is delayed by X" result. Add O2PLS joint-vs-specific variance as a
panel or inset. *(Case B equivalent implemented as `fig04_integration_caseB.pdf`
— design-cell PLS trajectories with bootstrap CIs and permutation-null Q²,
not MEFISTO factors)*

**Fig 5 — Trajectory clusters.** timeOmics multiblock clusters, RNA and
protein members overlaid per cluster; cluster silhouette; per-cluster
functional annotation; a cross-omics correlation circos for the top
cluster. *(not yet implemented)*

**Fig 6 — Kinetics.** (a) example ODE fits with posterior bands; (b)
kinetic volcano: Δk_syn vs Δk_deg, points sized by posterior certainty; (c)
genes with credible degradation-rate shifts, enriched pathways; (d) inferred
half-life distributions vs. published reference half-lives (external sanity
check). *(partially implemented within `fig03_concordance.pdf`'s
half-life/amplitude panels; full Bayesian posterior version not yet built)*

**Fig 7 — Delay atlas.** (a) distribution of learned $\tau_g$; (b)
functional enrichment of fast- vs slow-responding genes; (c) exemplar
RNA/protein trajectory pairs with the fitted delay overlaid. *(not yet
implemented)*

**Fig 8 — Network / ML.** (a) GNN residual modules on the STRING
subnetwork, nodes coloured by unexplained protein change; (b) CORUM
complex-level variance-ratio test showing stoichiometric buffering; (c)
SHAP summary from the XGBoost buffering model — what sequence/structural
features predict buffering. *(not yet implemented)*

**Fig 9 — Validation.** (a) orthogonal validation of 5–10 targets; (b)
permutation-null distributions vs. observed for DIABLO error rate, GNN
performance, factor variance; (c) imputation-scheme sensitivity; (d)
seed-stability of MEFISTO factors. *(imputation sensitivity implemented via
`results/de/imputation_sensitivity_lfc.tsv`; permutation nulls implemented
for the Case B integration Q²; remainder not yet built)*

### Statistical reporting that will pre-empt reviewer objections

- State n=3/cell and the resulting power limits explicitly, once, clearly.
  Don't bury it.
- Report FDR everywhere; consider IHW with mean abundance as the covariate
  for the proteomics side.
- Every ML/latent result gets a permutation null. No exceptions.
- Distinguish "not detected" from "not expressed" every time you interpret
  a Tier-2 gene.
- Frame all causal language as temporal precedence, never causation.
- Report annotation build, UniProt release, tool versions, and the exact
  filtering thresholds in a methods table.

### Recommended execution order

Each stage's output de-risks the next, and you can stop at any point with a
publishable paper:

1. **S0–S3** (mapping + missingness) — unglamorous, but every downstream
   result depends on it, and errors here are invisible until they've
   contaminated everything.
2. **S4–S5** (univariate + concordance) — gives you Fig 2 and Fig 3. **This
   alone is a paper.**
3. **S6–S7** (integration + trajectory clustering) — Figs 4–5. Now it's a
   *good* paper.
4. **S8** (Bayesian kinetics) — Fig 6. Now it's a mechanistic paper, and
   this is the highest-value custom work per unit effort.
5. **S9** (GNN / delay atlas / SHAP) — Figs 7–8. Methodological novelty.
   Attempt after 1–4 are solid, and be prepared to report it as
   negative/inconclusive if the nulls don't separate.
