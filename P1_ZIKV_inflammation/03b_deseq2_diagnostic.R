# -*- coding: utf-8 -*-
# 03b_deseq2_diagnostic.R — DESeq2 dispersion diagnostics for pseudobulk
# NOTE: This is an EXPLORATORY diagnostic script.
# The dataset lacks independent biological replicates, so DESeq2 formal testing
# is NOT performed. This script only reports dispersion estimates to document
# why formal testing is unreliable.

set.seed(123)
suppressMessages({library(Seurat); library(DESeq2)})

base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README.")
dnv <- readRDS(file.path(base, "denv_cc.rds"))
zik <- readRDS(file.path(base, "zika_clustered.rds"))

# ---- Build pseudobulk count matrix (same groups as 03_pseudobulk.R) ----
counts_dnv <- LayerData(dnv, assay="RNA", layer="counts")
counts_zik <- LayerData(zik, assay="RNA", layer="counts")

meta_dnv <- dnv@meta.data
meta_zik <- zik@meta.data

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

pb_counts <- matrix(0, nrow=length(common_genes), ncol=length(unique_groups),
                    dimnames=list(common_genes, unique_groups))
for (g in unique_groups) {
  cells <- groups$cell[groups$group == g]
  virus <- groups$virus[groups$group == g][1]
  use_counts <- if (virus == "DENV") counts_dnv else counts_zik
  pb_counts[, g] <- rowSums(use_counts[common_genes, cells, drop=FALSE])
}

# ---- DESeq2 dispersion diagnostic (NO formal testing) ----
cat("\n=== DESeq2 dispersion diagnostic (exploratory only) ===\n")
cat("Design: ~ virus + time + condition (4 time points x 2 viruses x 2 conditions)\n")
cat("Note: No biological replicates. dispersion estimates are descriptive.\n\n")

coldata <- data.frame(
  group = unique_groups,
  virus = sapply(strsplit(unique_groups, "_"), `[`, 1),
  time = as.numeric(gsub("t", "", sapply(strsplit(unique_groups, "_"), `[`, 2))),
  condition = sapply(strsplit(unique_groups, "_"), `[`, 3)
)
rownames(coldata) <- coldata$group

dds <- DESeqDataSetFromMatrix(countData=round(pb_counts[, unique_groups]),
                              colData=coldata,
                              design = ~ virus + time + condition)
# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat(sprintf("Genes after filtering: %d\n", nrow(dds)))

# Estimate dispersions (no formal test)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

disp <- dispersions(dds)
cat(sprintf("Dispersion estimates:\n"))
cat(sprintf("  Median dispersion: %.3f\n", median(disp)))
cat(sprintf("  Mean dispersion:   %.3f\n", mean(disp)))
cat(sprintf("  %% genes with inestimable dispersion (disp > 10): %.1f%%\n",
            mean(disp > 10) * 100))
cat(sprintf("  Residual degrees of freedom: %d\n",
            ncol(dds) - length(coefficients(dds))))

cat("\n=== Conclusion: dispersion estimates are unreliable due to lack of biological ===\n")
cat("=== replicates. Results are reported descriptively in the manuscript.        ===\n")

# Save diagnostic output
sink("P1_ZIKV_inflammation/deseq2_diagnostic.txt")
cat("DESeq2 dispersion diagnostic (exploratory, no formal testing)\n")
cat("=============================================================\n\n")
cat(sprintf("Genes after filtering: %d\n", nrow(dds)))
cat(sprintf("Median dispersion: %.3f\n", median(disp)))
cat(sprintf("Mean dispersion: %.3f\n", mean(disp)))
cat(sprintf("%% genes with inestimable dispersion (disp > 10): %.1f%%\n",
            mean(disp > 10) * 100))
cat(sprintf("Residual degrees of freedom: %d\n",
            ncol(dds) - length(coefficients(dds))))
cat("\nNo biological replicates available; DESeq2 results are not reported as statistical evidence.\n")
sink()
cat("\nSaved deseq2_diagnostic.txt\n")
