#!/usr/bin/env Rscript
## Build a DEMO external half-life prior from the simulator truth, with noise
## added so it behaves like a published turnover dataset rather than an oracle.
## In a real study this file comes from the literature or a pulse-SILAC arm.
source("R/utils.R")
cfg <- load_config(); set.seed(7)
tr  <- read.delim(file.path(cfg$paths$data_sim, "truth.tsv"), stringsAsFactors = FALSE)
hp  <- data.frame(gene_id  = tr$gene_id,
                  t_half_h = exp(log(tr$t_half_h) + rnorm(nrow(tr), 0, 0.5)),
                  sd_log   = 0.6)
write_tsv(hp, file.path(cfg$paths$data_sim, "half_life_prior.tsv"))
log_step("demo half-life prior written (log-scale noise sd = 0.5)")
