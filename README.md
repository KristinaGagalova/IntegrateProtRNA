# IntegrateProtRNA

RNA-seq + quantitative proteomics integration pipeline for a time-course
treatment study: 2 treatments × 4 timepoints × 3 biological replicates,
**24 independent samples per omics layer**.

**RNA and protein are measured on different biological samples** (unpaired /
"Case B" design), not split aliquots of the same 24 samples. That single fact
rules out sample-level integration (MOFA+/MEFISTO, DIABLO, O2PLS all require
matched rows) and shapes every method choice in this pipeline — see
[`docs/PIPELINE.md`](docs/PIPELINE.md) for the full reasoning and the
validated results on simulated ground-truth data.

## Quick start

```bash
Rscript env/install_r_deps.R      # CRAN + Bioconductor + timeOmics
./run_all.sh                      # simulate ground truth, run all 6 steps
```

Outputs land in `results/{qc,mapping,de,concordance,integration}/` and
`figures/`.

To run on real data instead of the simulator, populate `data/raw/` with the
same file layout as `data/simulated/` (see `R/00_simulate_data.R` for the
exact schema) and run `./run_all.sh --raw`.

## Pipeline

| Step | Script | Purpose |
|---|---|---|
| 00 | `R/00_simulate_data.R` | Ground-truth simulator: two independent 24-sample datasets from one shared latent biology |
| 00b | `R/00b_check_simulation.R` | Asserts the simulator produces the biology it claims |
| 01 | `R/01_qc_normalise.R` | QC, filtering, VST/voom (RNA), decoy/contaminant/validity filtering + normalisation (protein) |
| 02 | `R/02_id_mapping.R` | Bipartite gene↔protein-group graph → Universe A (all features) / Universe B (clean 1:1 pairs) |
| 03 | `R/03_missingness.R` | MNAR/MAR diagnosis, mechanism-specific imputation, sensitivity schemes |
| 04 | `R/04_univariate_temporal_de.R` | limma temporal DE for both layers, with moderated SEs |
| 05 | `R/05_concordance_archetypes.R` | Kinetic null model + amplitude LRT → regulatory archetypes |
| 06 | `R/06_integration_caseB.R` | Design-cell PLS with replicate bootstrap and leave-one-cell-out Q² |

Full write-up, including the three key validated findings (kinetic
buffering-vs-lag separation, the value of an external half-life prior, and
why in-sample cross-block correlation is worthless at 8 pseudo-samples):
[`docs/PIPELINE.md`](docs/PIPELINE.md).

## Claude Code setup (optional)

This repo was developed with [Claude Code](https://claude.com/product/claude-code)
configured so that local file edits stay inside this folder and remote
compute is restricted to one directory on a specific SSH host. If you want
that setup:

- `.claude/settings.json` / `.claude/hooks/validate_ssh.py` — client-side
  permission rules and an SSH command validator (not included in this repo;
  recreate locally if wanted, since they reference machine-specific paths)
- `docs/server-side-ssh-setup.md` — the part that actually matters: a
  forced-command SSH key on the server so the restriction holds even if the
  client-side hook is wrong or bypassed

Fill in your own host/user/path — this repo's copies have those redacted.
