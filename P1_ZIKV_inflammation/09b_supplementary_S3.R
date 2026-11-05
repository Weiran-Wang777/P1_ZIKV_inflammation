# ============================================================
# P1 Supplementary Table S3 only (S1/S2 already exported)
# ============================================================
set.seed(123)
suppressMessages(library(Seurat))

outdir <- "P1_ZIKV_inflammation/supplementary"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

denv <- readRDS("data/denv_cc.rds")
zika <- readRDS("data/zika_clustered.rds")

find_high_viral <- function(seu) {
  vbc <- aggregate(seu$viral_log, by = list(cluster = seu$seurat_clusters), FUN = median, na.rm = TRUE)
  return(list(cluster=vbc$cluster[which.max(vbc$x)], median=max(vbc$x)))
}
hc_d <- find_high_viral(denv)
hc_z <- find_high_viral(zika)

# ---- S3: Cluster summary ----
cluster_summary <- function(seu, virus, hc_cluster) {
  meta <- seu@meta.data
  clusters <- sort(unique(meta$seurat_clusters))
  rows <- lapply(clusters, function(cl) {
    sub <- meta[meta$seurat_clusters == cl, ]
    data.frame(
      Virus = virus,
      Cluster = as.character(cl),
      n_cells = nrow(sub),
      median_viral_log = round(median(sub$viral_log, na.rm=TRUE), 2),
      median_nGene = median(sub$nFeature_RNA, na.rm=TRUE),
      median_nUMI = median(sub$nCount_RNA, na.rm=TRUE),
      is_high_viral = (as.character(cl) == hc_cluster)
    )
  })
  do.call(rbind, rows)
}

s3_d <- cluster_summary(denv, "DENV", hc_d$cluster)
s3_z <- cluster_summary(zika, "ZIKV", hc_z$cluster)
s3_all <- rbind(s3_d, s3_z)
write.csv(s3_all, file.path(outdir, "TableS3_cluster_summary.csv"), row.names=FALSE)
cat("S3 cluster summary exported:\n")
print(s3_all)

# ---- S3b: Time point x infection state breakdown ----
tp_breakdown <- function(seu, virus) {
  meta <- seu@meta.data
  meta$infected <- ifelse(meta$viral_log > 0.5, "infected", "uninfected")
  agg <- aggregate(row.names(meta) ~ time + MOI + infected, data=meta, FUN=length)
  colnames(agg) <- c("time", "MOI", "infected", "n_cells")
  agg$Virus <- virus
  agg
}

tp_d <- tp_breakdown(denv, "DENV")
tp_z <- tp_breakdown(zika, "ZIKV")
tp_all <- rbind(tp_d, tp_z)
write.csv(tp_all, file.path(outdir, "TableS3_timepoint_breakdown.csv"), row.names=FALSE)
cat("\nS3 time-point breakdown exported:", nrow(tp_all), "rows\n")

# ---- S3c: Viral-load matched ----
matched <- function(seu, virus) {
  meta <- seu@meta.data
  inf <- meta[meta$viral_log > 0.5 & meta$viral_log >= 7 & meta$viral_log <= 9, ]
  data.frame(Virus=virus, n_matched=nrow(inf),
             median_viral_log=round(median(inf$viral_log), 2))
}
mat_d <- matched(denv, "DENV")
mat_z <- matched(zika, "ZIKV")
write.csv(rbind(mat_d, mat_z), file.path(outdir, "TableS3_viral_load_matched.csv"), row.names=FALSE)
cat("S3 viral-load matched exported\n")

cat("\nS3 ALL DONE\n")
