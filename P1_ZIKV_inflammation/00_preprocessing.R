# -*- coding: utf-8 -*-
# 00_preprocessing.R — QC, normalization, dimensionality reduction, clustering
# Run this FIRST to produce the processed Seurat objects (denv_cc.rds, zika_clustered.rds)
# from raw count matrices.
#
# Input:  data/DENV_counts.csv, data/DENV_metadata.csv (or .rds)
#         data/Zika_counts.csv, data/Zika_metadata.csv
# Output: data/denv_cc.rds, data/zika_clustered.rds

set.seed(123)
suppressMessages({library(Seurat); library(ggplot2)})

base <- "data"

# ---- Load raw data ----
# The exact file format depends on how you downloaded GSE110496.
# This script assumes you have placed raw count matrices in data/.
# If starting from the Zanini et al. supplementary matrices, convert to:
#   rows = genes, columns = cells, values = raw counts

load_raw <- function(virus_name) {
  if (virus_name == "DENV") {
    counts_file <- file.path(base, "DENV_counts.csv")
    meta_file <- file.path(base, "DENV_metadata.csv")
  } else {
    counts_file <- file.path(base, "Zika_counts.csv")
    meta_file <- file.path(base, "Zika_metadata.csv")
  }
  if (!file.exists(counts_file))
    stop("Error: ", counts_file, " not found. Place raw count matrix in data/.")
  if (!file.exists(meta_file))
    stop("Error: ", meta_file, " not found. Place metadata in data/.")

  counts <- read.csv(counts_file, row.names=1, check.names=FALSE)
  meta <- read.csv(meta_file, row.names=1, check.names=FALSE)
  counts <- as(as.matrix(counts), "dgCMatrix")
  list(counts=counts, meta=meta)
}

# ---- QC filter ----
qc_filter <- function(seu) {
  n_before <- ncol(seu)
  # Remove low-quality cells: 200-7500 detected genes, <10% mitochondrial
  seu <- subset(seu, subset = nFeature_RNA >= 200 & nFeature_RNA <= 7500 &
                           percent.mt < 10)
  n_after <- ncol(seu)
  cat(sprintf("  QC: %d -> %d cells (removed %d)\n", n_before, n_after, n_before - n_after))
  return(seu)
}

# ---- Process one virus ----
process_virus <- function(virus_name) {
  cat(sprintf("=== Processing %s ===\n", virus_name))
  raw <- load_raw(virus_name)

  seu <- CreateSeuratObject(counts=raw$counts, meta.data=raw$meta)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")

  # QC
  seu <- qc_filter(seu)

  # Normalize
  seu <- NormalizeData(seu, normalization.method="LogNormalize", scale.factor=10000)
  seu <- FindVariableFeatures(seu, nfeatures=2000)
  seu <- ScaleData(seu, vars.to.regress=NULL)

  # PCA
  seu <- RunPCA(seu, npcs=50, verbose=FALSE)

  # UMAP
  seu <- RunUMAP(seu, dims=1:20, n.neighbors=30, min.dist=0.3)

  # Clustering (Louvain, resolution 0.5)
  seu <- FindNeighbors(seu, dims=1:20)
  seu <- FindClusters(seu, resolution=0.5)

  # Cell cycle scoring (human genes)
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  seu <- CellCycleScoring(seu, s.features=s_genes, g2m.features=g2m_genes)

  # Save
  outfile <- if (virus_name == "DENV") file.path(base, "denv_cc.rds")
             else file.path(base, "zika_clustered.rds")
  saveRDS(seu, outfile)
  cat(sprintf("  Saved: %s (%d cells)\n", outfile, ncol(seu)))
  return(seu)
}

# ---- Run ----
denv <- process_virus("DENV")
zika <- process_virus("Zika")

cat("\n=== Preprocessing complete. Proceed to run_all.R for analysis. ===\n")
