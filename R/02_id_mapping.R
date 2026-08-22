#!/usr/bin/env Rscript
## =============================================================================
## 02_id_mapping.R -- transcript/protein ID reconciliation via a bipartite graph
##
## The mapping between genes and protein groups is many-to-many, and the
## ambiguity is concentrated on the MS side (shared peptides, isoforms, razor
## proteins). So we map protein -> gene, never the reverse, and we resolve the
## ambiguity structurally rather than by picking the leading razor protein and
## hoping.
##
## Method: build a bipartite graph  gene <--> protein_group  and take its
## CONNECTED COMPONENTS. Each component is then one of four types:
##
##   simple    1 gene  : 1 group    -> clean pair, Universe B
##   split     1 gene  : N groups   -> collapse groups (sum linear intensity),
##                                     then a clean pair, Universe B
##   shared    N genes : 1 group    -> peptides cannot distinguish the genes;
##                                     NOT a pair. Feature-cluster, Universe A
##   complex   N genes : M groups   -> paralogue family; feature-cluster only
##
## Two universes come out of this, and they are used for different things:
##   Universe A  all features, both layers, NO pairing required
##               -> factor models, pathway analysis, network models
##   Universe B  clean 1:1 gene<->protein pairs
##               -> concordance, archetypes, kinetics, delay estimation
##
## Do NOT subset to the intersection before running a factor model: the shared
## dimension there is the design cell, not the feature.
##
## Never join on gene symbols. They are non-unique, non-stable, and
## Excel-corruptible (SEPT7, MARCH1 -> dates). Everything here joins on
## stable accessions and Ensembl gene IDs.
##
## Usage: Rscript R/02_id_mapping.R [config/config.yaml]
## =============================================================================

source("R/utils.R")
need_pkgs("igraph")

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- load_config(if (length(args)) args[1] else "config/config.yaml")
qcdir  <- file.path(cfg$paths$results, "qc")
outdir <- ensure_dir(cfg$paths$results, "mapping")
figdir <- ensure_dir(cfg$paths$figures)

prot   <- read_mat(file.path(qcdir, "protein_norm.tsv"))
pmeta  <- read.delim(file.path(qcdir, "protein_meta_filt.tsv"), stringsAsFactors = FALSE)
rna    <- read_mat(file.path(qcdir, "rna_vst.tsv"))
idmap  <- read.delim(file.path(cfg$paths$data_sim, "id_map_raw.tsv"),
                     stringsAsFactors = FALSE)
genes  <- rownames(rna)
log_step("mapping ", nrow(prot), " protein groups onto ", length(genes), " genes")

## ------------------------------------------------- 1. explode accessions ----
## MaxQuant / DIA-NN pack a whole protein group into one semicolon-delimited
## field. Taking only the leading entry silently discards the ambiguity we are
## trying to characterise, so expand the full list.
acc_long <- do.call(rbind, lapply(seq_len(nrow(pmeta)), function(i) {
  a <- trimws(strsplit(pmeta$majority_protein_ids[i], ";", fixed = TRUE)[[1]])
  a <- a[nzchar(a)]
  if (!length(a)) return(NULL)
  data.frame(group_id = pmeta$group_id[i], accession = a, stringsAsFactors = FALSE)
}))
## Strip UniProt isoform suffixes (P12345-2 -> P12345) for the gene-level join,
## but keep the original so isoform-resolved analysis stays possible later.
acc_long$accession_base <- sub("-[0-9]+$", "", acc_long$accession)
log_step("accessions after expansion: ", nrow(acc_long),
         " (", length(unique(acc_long$accession_base)), " unique)")

## ------------------------------------------------------ 2. accession->gene --
idmap <- idmap[!is.na(idmap$gene_id) & nzchar(idmap$gene_id), , drop = FALSE]
idmap <- idmap[!duplicated(idmap$accession), , drop = FALSE]
acc_long$gene_id <- idmap$gene_id[match(acc_long$accession_base, idmap$accession)]

n_unmapped <- sum(is.na(acc_long$gene_id))
log_step("accessions with no gene mapping: ", n_unmapped,
         sprintf(" (%.1f%%)", 100 * n_unmapped / nrow(acc_long)))

## Keep only genes that survived RNA filtering; the rest cannot be paired.
edges <- acc_long[!is.na(acc_long$gene_id) & acc_long$gene_id %in% genes,
                  c("gene_id", "group_id")]
edges <- unique(edges)

## ------------------------------------------------- 3. bipartite components --
g <- igraph::graph_from_data_frame(
  data.frame(from = paste0("G:", edges$gene_id),
             to   = paste0("P:", edges$group_id)), directed = FALSE)
comp <- igraph::components(g)
memb <- comp$membership

nodes <- data.frame(node = names(memb), comp = as.integer(memb),
                    stringsAsFactors = FALSE)
nodes$kind <- ifelse(startsWith(nodes$node, "G:"), "gene", "group")
nodes$id   <- sub("^[GP]:", "", nodes$node)

cs <- aggregate(kind ~ comp, nodes, function(k) c(sum(k == "gene"), sum(k == "group")))
comp_tab <- data.frame(comp = cs$comp, n_gene = cs$kind[, 1], n_group = cs$kind[, 2])
comp_tab$type <- with(comp_tab, ifelse(
  n_gene == 1 & n_group == 1, "simple",
  ifelse(n_gene == 1 & n_group > 1, "split",
  ifelse(n_gene > 1 & n_group == 1, "shared", "complex"))))
log_step("component types:")
print(table(comp_tab$type))

## ------------------------------------------- 4. gene-level protein matrix ---
## For simple and split components the gene is unambiguous, so collapse its
## protein groups. Summing LINEAR intensities (not averaging log2) is the
## right operation: isoform groups partition the same protein pool, and the
## quantity of interest is total abundance.
lin_sum <- function(m) {
  if (nrow(m) == 1) return(m[1, ])
  s <- colSums(2^m, na.rm = TRUE)
  s[colSums(!is.na(m)) == 0] <- NA_real_
  log2(s)
}

pairable <- comp_tab$comp[comp_tab$type %in% c("simple", "split")]
pair_rows <- lapply(pairable, function(cc) {
  nd  <- nodes[nodes$comp == cc, ]
  gid <- nd$id[nd$kind == "gene"]
  grp <- nd$id[nd$kind == "group"]
  grp <- grp[grp %in% rownames(prot)]
  if (!length(grp)) return(NULL)
  list(gene_id = gid, groups = paste(grp, collapse = ";"),
       n_groups = length(grp), values = lin_sum(prot[grp, , drop = FALSE]))
})
pair_rows <- Filter(Negate(is.null), pair_rows)

prot_gene <- do.call(rbind, lapply(pair_rows, `[[`, "values"))
rownames(prot_gene) <- vapply(pair_rows, `[[`, character(1), "gene_id")
universeB <- data.frame(
  gene_id  = rownames(prot_gene),
  groups   = vapply(pair_rows, `[[`, character(1), "groups"),
  n_groups = vapply(pair_rows, `[[`, numeric(1),   "n_groups"),
  stringsAsFactors = FALSE)

## ------------------------------------------------- 5. ambiguous clusters ----
## shared / complex components get one FEATURE-CLUSTER id each. They stay in
## Universe A (factor models, pathway analysis) but are excluded from every
## pairwise analysis, and the exclusion count is reported in the methods.
ambig <- comp_tab$comp[comp_tab$type %in% c("shared", "complex")]
clusters <- do.call(rbind, lapply(ambig, function(cc) {
  nd <- nodes[nodes$comp == cc, ]
  data.frame(cluster_id = sprintf("FC%04d", cc),
             type       = comp_tab$type[comp_tab$comp == cc],
             genes      = paste(sort(nd$id[nd$kind == "gene"]),  collapse = ";"),
             groups     = paste(sort(nd$id[nd$kind == "group"]), collapse = ";"),
             stringsAsFactors = FALSE)
}))

## --------------------------------------------------------- 6. tiering ------
tier1 <- universeB$gene_id                                    # matched
tier2 <- setdiff(genes, unique(edges$gene_id))                # RNA-only
tier3 <- setdiff(rownames(prot), unique(edges$group_id))      # protein-only
log_step("Tier 1 matched pairs : ", length(tier1))
log_step("Tier 2 RNA-only genes: ", length(tier2))
log_step("Tier 3 protein-only  : ", length(tier3))

## Tier 2 is NOT evidence of post-transcriptional repression. Most of it is
## detectability, so model that explicitly rather than interpreting silence.
## With only abundance available here we test the abundance component; add
## peptide count, hydrophobicity and sequence length when you have them.
rna_mean <- rowMeans(read_mat(file.path(qcdir, "rna_logcpm.tsv")))
det <- data.frame(gene_id = genes,
                  detected = as.integer(genes %in% unique(edges$gene_id)),
                  rna_mean = rna_mean[genes])
det_fit <- glm(detected ~ rna_mean, binomial, det)
det$p_detect <- predict(det_fit, type = "response")
log_step("detectability model: odds ratio per log2CPM = ",
         round(exp(coef(det_fit))[2], 3),
         "  (Tier 2 is dominated by abundance if this is >> 1)")

## ------------------------------------------------------------- 7. figure ----
open_pdf(file.path(figdir, "fig01b_mapping.pdf"), 11, 4.2)
op <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1))

bp <- barplot(table(comp_tab$type), col = "#2C6FBB", ylab = "components",
              main = "Bipartite component types")
text(bp, table(comp_tab$type), table(comp_tab$type), pos = 1, col = "white")

barplot(c(`Tier 1\nmatched` = length(tier1), `Tier 2\nRNA-only` = length(tier2),
          `Tier 3\nprot-only` = length(tier3)),
        col = c("#2C6FBB", "#7FB2E5", "#D1495B"), ylab = "features",
        main = "Feature coverage tiers")

o <- order(det$rna_mean)
plot(det$rna_mean[o], det$detected[o], pch = 16, cex = 0.3, col = "#00000033",
     xlab = "mean RNA log2 CPM", ylab = "protein detected (0/1)",
     main = "Detectability vs abundance\n(Tier 2 is mostly a detection limit)")
lines(det$rna_mean[o], det$p_detect[o], col = "#D1495B", lwd = 2)
par(op); dev.off()

## -------------------------------------------------------------- 8. write ----
write_mat(prot_gene, file.path(outdir, "protein_gene_level.tsv"), "gene_id")
write_tsv(universeB,  file.path(outdir, "universeB_pairs.tsv"))
write_tsv(comp_tab,   file.path(outdir, "components.tsv"))
if (!is.null(clusters)) write_tsv(clusters, file.path(outdir, "feature_clusters.tsv"))
write_tsv(data.frame(gene_id = tier2), file.path(outdir, "tier2_rna_only.tsv"))
write_tsv(data.frame(group_id = tier3), file.path(outdir, "tier3_protein_only.tsv"))
write_tsv(det, file.path(outdir, "detectability.tsv"))
write_tsv(data.frame(
  metric = c("accessions_expanded", "accessions_unmapped", "components",
             "simple", "split", "shared", "complex",
             "universeB_pairs", "tier2_rna_only", "tier3_protein_only"),
  value = c(nrow(acc_long), n_unmapped, nrow(comp_tab),
            sum(comp_tab$type == "simple"), sum(comp_tab$type == "split"),
            sum(comp_tab$type == "shared"), sum(comp_tab$type == "complex"),
            nrow(universeB), length(tier2), length(tier3))),
  file.path(outdir, "mapping_summary.tsv"))

log_step("02_id_mapping done -> ", outdir)
