# -*- coding: utf-8 -*-
# P1 pseudobulk: DESCRIPTIVE aggregation only (no formal statistical testing)
# Note: This dataset lacks independent biological replicates.
# Results are presented as directional support, not as statistical evidence.

set.seed(123)
suppressMessages({library(Seurat)})

# ---- Load data (relative path; user places rds files in data/) ----
base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
dnv <- readRDS(file.path(base, "denv_cc.rds"))
zik <- readRDS(file.path(base, "zika_clustered.rds"))

counts_dnv <- LayerData(dnv, assay="RNA", layer="counts")
counts_zik <- LayerData(zik, assay="RNA", layer="counts")

meta_dnv <- dnv@meta.data
meta_zik <- zik@meta.data

# ---- Build group labels ----
# Infection threshold consistent with P2: viral_log > 0.5
groups_dnv <- data.frame(
  cell = rownames(meta_dnv),
  virus = "DENV",
  time = meta_dnv$time,
  condition = ifelse(meta_dnv$viral_log > 0.5, "infected", "uninfected")
)
groups_zik <- data.frame(
  cell = rownames(meta_zik),
  virus = "ZIKV",
  time = meta_zik$time,
  condition = ifelse(meta_zik$viral_log > 0.5, "infected", "uninfected")
)
groups <- rbind(groups_dnv, groups_zik)
groups$group <- paste0(groups$virus, "_t", groups$time, "_", groups$condition)

common_genes <- intersect(rownames(counts_dnv), rownames(counts_zik))
unique_groups <- unique(groups$group)

# ---- Pseudobulk aggregation (per virus, no cross-virus mixing) ----
pb_counts <- matrix(0, nrow=length(common_genes), ncol=length(unique_groups),
                    dimnames=list(common_genes, unique_groups))
for (g in unique_groups) {
  cells <- groups$cell[groups$group == g]
  virus <- groups$virus[groups$group == g][1]
  use_counts <- if (virus == "DENV") counts_dnv else counts_zik
  pb_counts[, g] <- rowSums(use_counts[common_genes, cells, drop=FALSE])
}

# ---- Key genes for descriptive comparison ----
key_genes <- c("CXCL10","VCAM1","NFKB2","RELB","SOCS2","HSPA6","ATF3","DDIT3","HERPUD1","TNF","IL6","CXCL8","CCL5")

# ---- Report: counts per cell at each time point ----
cat("\n=== Pseudobulk descriptive analysis (no statistical testing) ===\n")
cat("Note: No biological replicates available. Results are directional only.\n\n")

for (g in key_genes) {
  if (!(g %in% rownames(pb_counts))) next
  cat(sprintf("\n--- %s ---\n", g))
  for (grp in sort(unique_groups)) {
    n_cells <- sum(groups$group == grp)
    cpc <- pb_counts[g, grp] / n_cells
    cat(sprintf("  %-25s n=%4d  counts/cell=%.3f\n", grp, n_cells, cpc))
  }
}

cat("\n=== DONE ===\n")
