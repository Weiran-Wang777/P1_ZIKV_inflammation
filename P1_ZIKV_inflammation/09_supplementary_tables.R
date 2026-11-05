# ============================================================
# P1 Supplementary Tables export
# Run from P1_ZIKV_inflammation_code/ root: Rscript P1_ZIKV_inflammation/09_supplementary_tables.R
# ============================================================
set.seed(123)
suppressMessages({
  library(Seurat); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)
})

outdir <- "P1_ZIKV_inflammation/supplementary"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

denv <- readRDS("data/denv_cc.rds")
zika <- readRDS("data/zika_clustered.rds")

# ---- Auto-detect high-viral cluster (same as 01_qc_cluster.R) ----
find_high_viral <- function(seu) {
  viral_by_cluster <- aggregate(seu$viral_log, by = list(cluster = seu$seurat_clusters),
                                FUN = median, na.rm = TRUE)
  hc <- viral_by_cluster$cluster[which.max(viral_by_cluster$x)]
  return(list(cluster=hc, median=max(viral_by_cluster$x)))
}
hc_d <- find_high_viral(denv)
hc_z <- find_high_viral(zika)
cat(sprintf("DENV high-viral: cluster %s, median viral_log=%.2f\n", hc_d$cluster, hc_d$median))
cat(sprintf("ZIKV high-viral: cluster %s, median viral_log=%.2f\n", hc_z$cluster, hc_z$median))

# ================= S1: Full DEG tables =================
Idents(denv) <- "seurat_clusters"
Idents(zika) <- "seurat_clusters"

deg_d <- FindMarkers(denv, ident.1=hc_d$cluster, ident.2=setdiff(levels(denv), hc_d$cluster),
                     logfc.threshold=0, min.pct=0.01, only.pos=FALSE)
deg_z <- FindMarkers(zika, ident.1=hc_z$cluster, ident.2=setdiff(levels(zika), hc_z$cluster),
                     logfc.threshold=0, min.pct=0.01, only.pos=FALSE)

deg_d$gene <- rownames(deg_d)
deg_z$gene <- rownames(deg_z)
deg_d$virus <- "DENV"
deg_z$virus <- "ZIKV"

# Upregulated as defined in manuscript (log2FC>0.5, padj<0.05)
deg_d$upregulated <- deg_d$avg_log2FC > 0.5 & deg_d$p_val_adj < 0.05
deg_z$upregulated <- deg_z$avg_log2FC > 0.5 & deg_z$p_val_adj < 0.05

write.csv(deg_d, file.path(outdir, "TableS1_DENV_DEG.csv"), row.names=FALSE)
write.csv(deg_z, file.path(outdir, "TableS1_ZIKV_DEG.csv"), row.names=FALSE)
cat(sprintf("S1 exported: DENV=%d genes, ZIKV=%d genes\n", nrow(deg_d), nrow(deg_z)))

# ================= S2: GO BP + KEGG full tables =================
deg_d_up <- rownames(deg_d)[deg_d$upregulated]
deg_z_up <- rownames(deg_z)[deg_z$upregulated]
only_z <- setdiff(deg_z_up, deg_d_up)
only_d <- setdiff(deg_d_up, deg_z_up)
shared <- intersect(deg_d_up, deg_z_up)
cat(sprintf("Up: DENV=%d ZIKV=%d shared=%d onlyD=%d onlyZ=%d\n",
            length(deg_d_up), length(deg_z_up), length(shared), length(only_d), length(only_z)))

run_enrich <- function(genes, type="go") {
  eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if (nrow(eg) < 5) return(NULL)
  if (type=="go") enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                           pAdjustMethod="BH", readable=TRUE)
  else enrichKEGG(gene=eg$ENTREZID, organism="hsa", pAdjustMethod="BH")
}

go_d <- run_enrich(deg_d_up)
go_z <- run_enrich(deg_z_up)
kegg_z <- run_enrich(only_z, type="kegg")
kegg_d <- run_enrich(only_d, type="kegg")

go_dd <- as.data.frame(go_d); go_dd$Virus <- "DENV"
go_zd <- as.data.frame(go_z); go_zd$Virus <- "ZIKV"
go_all <- rbind(go_dd, go_zd)
write.csv(go_all, file.path(outdir, "TableS2_GO_BP_all.csv"), row.names=FALSE)
cat(sprintf("S2 GO BP: DENV=%d terms, ZIKV=%d terms\n", nrow(go_dd), nrow(go_zd)))

kz <- as.data.frame(kegg_z); kz$Virus <- "ZIKV-specific"
if (!is.null(kegg_d)) {
  kd <- as.data.frame(kegg_d); kd$Virus <- "DENV-specific"
  kegg_all <- rbind(kz, kd)
} else {
  kegg_all <- kz
}
write.csv(kegg_all, file.path(outdir, "TableS2_KEGG_all.csv"), row.names=FALSE)
cat(sprintf("S2 KEGG: ZIKV-specific=%d terms, DENV-specific=%d terms\n",
            nrow(kz), ifelse(is.null(kegg_d), 0, nrow(kd))))

# ================= S3: Cluster summary + QC =================
cluster_summary <- function(seu, virus, hc_cluster) {
  meta <- seu@meta.data
  df <- data.frame(
    Virus = virus,
    Cluster = as.character(meta$seurat_clusters),
    n_cells = 1
  )
  agg <- aggregate(n_cells ~ Virus + Cluster, data=df, FUN=sum)
  viral_by_cl <- aggregate(viral_log ~ seurat_clusters, data=meta, FUN=median, na.rm=TRUE)
  nGene_by_cl <- aggregate(nFeature_RNA ~ seurat_clusters, data=meta, FUN=median, na.rm=TRUE)
  nUMI_by_cl <- aggregate(nCount_RNA ~ seurat_clusters, data=meta, FUN=median, na.rm=TRUE)
  mt_by_cl <- aggregate(percent.mt ~ seurat_clusters, data=meta, FUN=median, na.rm=TRUE)

  agg$median_viral_log <- viral_by_cl$viral_log[match(agg$Cluster, viral_by_cl$seurat_clusters)]
  agg$median_nGene <- nGene_by_cl$nFeature_RNA[match(agg$Cluster, nGene_by_cl$seurat_clusters)]
  agg$median_nUMI <- nUMI_by_cl$nCount_RNA[match(agg$Cluster, nUMI_by_cl$seurat_clusters)]
  agg$median_percent_mt <- mt_by_cl$percent.mt[match(agg$Cluster, mt_by_cl$seurat_clusters)]
  agg$is_high_viral <- agg$Cluster == hc_cluster
  return(agg)
}

s3_d <- cluster_summary(denv, "DENV", hc_d$cluster)
s3_z <- cluster_summary(zika, "ZIKV", hc_z$cluster)
s3_all <- rbind(s3_d, s3_z)
write.csv(s3_all, file.path(outdir, "TableS3_cluster_summary.csv"), row.names=FALSE)
cat(sprintf("S3 cluster summary: %d rows\n", nrow(s3_all)))

# Time-point x infection state breakdown
tp_breakdown <- function(seu, virus) {
  meta <- seu@meta.data
  # Try to find time point and MOI columns
  cat(sprintf("  Meta columns: %s\n", paste(colnames(meta), collapse=", ")))
  # Build from available columns
  result <- data.frame(Virus=virus)
  if ("time" %in% colnames(meta)) {
    result <- aggregate(viral_log ~ time + (viral_log > 0.5), data=meta, FUN=length)
    colnames(result) <- c("time", "infected", "n_cells")
    result$Virus <- virus
  }
  return(result)
}

tp_d <- tryCatch(tp_breakdown(denv, "DENV"), error=function(e) NULL)
tp_z <- tryCatch(tp_breakdown(zika, "ZIKV"), error=function(e) NULL)
if (!is.null(tp_d) && !is.null(tp_z)) {
  write.csv(rbind(tp_d, tp_z), file.path(outdir, "TableS3_timepoint_breakdown.csv"), row.names=FALSE)
}

# Viral-load matched comparison numbers
matched <- function(seu, virus) {
  meta <- seu@meta.data
  inf <- meta[meta$viral_log > 0.5 & meta$viral_log >= 7 & meta$viral_log <= 9, ]
  data.frame(Virus=virus, n_matched=nrow(inf),
             median_viral_log=median(inf$viral_log))
}
mat_d <- matched(denv, "DENV")
mat_z <- matched(zika, "ZIKV")
write.csv(rbind(mat_d, mat_z), file.path(outdir, "TableS3_viral_load_matched.csv"), row.names=FALSE)

cat("\nALL SUPPLEMENTARY TABLES EXPORTED to P1_ZIKV_inflammation/supplementary/\n")
