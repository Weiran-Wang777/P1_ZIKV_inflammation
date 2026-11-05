# -*- coding: utf-8 -*-
# P1 figures v2 — Emerging Microbes & Infections style.
# Regenerate Fig1-4 for the ZIKV-vs-DENV inflammation paper.

set.seed(123)
suppressMessages({library(Seurat); library(ggplot2); library(patchwork)
                  library(dplyr); library(clusterProfiler); library(org.Hs.eg.db)})

# ---- Load data (relative path; user places rds in data/) ----
base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
denv <- readRDS(file.path(base, "denv_cc.rds"))
zika <- readRDS(file.path(base, "zika_clustered.rds"))
outdir <- "P1_ZIKV_inflammation/figures/"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

# ---- Auto-detect high-viral cluster ----
find_high_viral <- function(seu) {
  viral_by_cluster <- aggregate(seu$viral_log, by = list(cluster = seu$seurat_clusters),
                                FUN = median, na.rm = TRUE)
  hc <- viral_by_cluster$cluster[which.max(viral_by_cluster$x)]
  message("High-viral cluster: ", hc, " (median viral_log = ", round(max(viral_by_cluster$x), 2), ")")
  return(hc)
}
hc_d <- find_high_viral(denv)
hc_z <- find_high_viral(zika)

theme_pub <- function(base=12) {
  theme_classic(base_size=base) +
  theme(text=element_text(family="sans"),
        axis.text=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        legend.key=element_blank(),
        legend.title=element_text(size=base-2),
        plot.title=element_text(face="bold", size=base+1, hjust=0.5),
        plot.margin=margin(6,6,6,6))
}
save_fig <- function(p, name, w, h) {
  # Cap width at 7.5 inches (190mm = T&F double-column max)
  if (w > 7.5) { h <- h * (7.5 / w); w <- 7.5 }
  ggsave(paste0(outdir, name, ".png"), p, width=w, height=h, dpi=600, bg="white")
  ggsave(paste0(outdir, name, ".pdf"), p, width=w, height=h, device=cairo_pdf, bg="white")
  cat("Saved", name, " (", round(w*25.4), "x", round(h*25.4), "mm)\n", sep="")
}
cols_virus <- c("DENV"="#C0392B", "ZIKV"="#1A6FA0")

# ================= Fig1: UMAP clusters + viral load (2x2) =================
Idents(denv) <- "seurat_clusters"; Idents(zika) <- "seurat_clusters"
umap_d <- as.data.frame(Embeddings(denv, "umap")); umap_d$cluster <- denv$seurat_clusters; umap_d$viral <- denv$viral_log
umap_z <- as.data.frame(Embeddings(zika, "umap")); umap_z$cluster <- zika$seurat_clusters; umap_z$viral <- zika$viral_log

pal_cl <- c("#4C72B0","#DD8452","#55A868","#C44E52","#8172B3","#937860")
f1a <- ggplot(umap_d, aes(umap_1, umap_2, colour=cluster)) +
  geom_point(size=0.8, alpha=0.85) +
  scale_colour_manual(values=pal_cl, name="Cluster") +
  labs(title="DENV", x="UMAP 1", y="UMAP 2") + theme_pub(12) +
  theme(legend.position="right", legend.key.size=unit(0.4,"cm"))
f1b <- ggplot(umap_d, aes(umap_1, umap_2, colour=viral)) +
  geom_point(size=0.8) +
  scale_color_gradientn(colors=c("grey92","#F5B7B1","#E74C3C","#641E16"),
                        values=c(0,0.08,0.25,1), name="viral\nload") +
  labs(title="DENV", x="UMAP 1", y="UMAP 2") + theme_pub(12) +
  theme(legend.position="right", legend.key.height=unit(0.6,"cm"))
f1c <- ggplot(umap_z, aes(umap_1, umap_2, colour=cluster)) +
  geom_point(size=0.8, alpha=0.85) +
  scale_colour_manual(values=pal_cl[1:5], name="Cluster") +
  labs(title="ZIKV", x="UMAP 1", y="UMAP 2") + theme_pub(12) +
  theme(legend.position="right", legend.key.size=unit(0.4,"cm"))
f1d <- ggplot(umap_z, aes(umap_1, umap_2, colour=viral)) +
  geom_point(size=0.8) +
  scale_color_gradientn(colors=c("grey92","#AED6F1","#2E86C1","#1B4F72"),
                        values=c(0,0.08,0.25,1), name="viral\nload") +
  labs(title="ZIKV", x="UMAP 1", y="UMAP 2") + theme_pub(12) +
  theme(legend.position="right", legend.key.height=unit(0.6,"cm"))
save_fig((f1a | f1b) / (f1c | f1d) + plot_annotation(tag_levels='A'), "Fig1_UMAP_viral_load", 10.5, 8.5)

# ================= DEG: high-viral cluster vs rest =================
deg_d <- FindMarkers(denv, ident.1=hc_d, ident.2=setdiff(levels(denv), hc_d), logfc.threshold=0, min.pct=0.01, only.pos=FALSE)
deg_z <- FindMarkers(zika, ident.1=hc_z, ident.2=setdiff(levels(zika), hc_z), logfc.threshold=0, min.pct=0.01, only.pos=FALSE)
deg_d_up <- rownames(deg_d)[deg_d$avg_log2FC > 0.5 & deg_d$p_val_adj < 0.05]
deg_z_up <- rownames(deg_z)[deg_z$avg_log2FC > 0.5 & deg_z$p_val_adj < 0.05]
shared <- intersect(deg_d_up, deg_z_up)
only_d <- setdiff(deg_d_up, deg_z_up)
only_z <- setdiff(deg_z_up, deg_d_up)
cat(sprintf("Up: DENV=%d ZIKV=%d shared=%d onlyD=%d onlyZ=%d\n",
    length(deg_d_up), length(deg_z_up), length(shared), length(only_d), length(only_z)))

# ================= Fig2: DEG overlap + top DEG bar =================
overlap_df <- data.frame(Category=c("ZIKV-specific","Shared","DENV-specific"),
                         Count=c(length(only_z), length(shared), length(only_d)))
overlap_df$Category <- factor(overlap_df$Category, levels=c("ZIKV-specific","Shared","DENV-specific"))
f2a <- ggplot(overlap_df, aes(Category, Count, fill=Category)) +
  geom_col(width=0.6) +
  geom_text(aes(label=Count), vjust=-0.4, size=4.5) +
  scale_fill_manual(values=c("ZIKV-specific"="#1A6FA0","Shared"="#55A868","DENV-specific"="#C0392B"), guide="none") +
  labs(x=NULL, y="Number of upregulated genes", title="High-viral cluster DEG overlap") +
  theme_pub(12) + theme(axis.text.x=element_text(size=10))

er_shared <- c("HSPA6","ATF3","DDIT3","TRIB3","HERPUD1","CTH",
               "DNAJB9","SQSTM1","ERO1B","HSPA5","DNAJC3","PDIA4")
top <- er_shared[er_shared %in% shared]
if (length(top) < 2) {
  shared_df <- data.frame(gene=shared, DENV=deg_d[shared,"avg_log2FC"], ZIKV=deg_z[shared,"avg_log2FC"])
  shared_df$mean <- rowMeans(shared_df[,2:3]); shared_df <- shared_df[order(-shared_df$mean),]
  top <- head(shared_df$gene, 10)
}
top_df <- rbind(data.frame(gene=top, Virus="DENV", fc=deg_d[top,"avg_log2FC"]),
                data.frame(gene=top, Virus="ZIKV", fc=deg_z[top,"avg_log2FC"]))
top_df$gene <- factor(top_df$gene, levels=rev(top))
f2b <- ggplot(top_df, aes(gene, fc, fill=Virus)) +
  geom_col(position="dodge", width=0.7) +
  scale_fill_manual(values=cols_virus) +
  coord_flip() +
  labs(x=NULL, y="avg log2FC", title="Shared ER-stress genes in high-viral clusters") +
  theme_pub(11) + theme(axis.text.y=element_text(size=9))
save_fig((f2a / f2b) + plot_annotation(tag_levels='A'), "Fig2_DEG_overlap", 9, 8)

# ================= GO / KEGG =================
run_enrich <- function(genes, ont="BP", type="go") {
  eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if (nrow(eg) < 5) return(NULL)
  if (type=="go") enrichGO(gene=eg$ENTREZID, OrgDb=org.Hs.eg.db, ont=ont, pAdjustMethod="BH", readable=TRUE)
  else enrichKEGG(gene=eg$ENTREZID, organism="hsa", pAdjustMethod="BH")
}
go_d <- run_enrich(deg_d_up); go_z <- run_enrich(deg_z_up)
kegg_z <- run_enrich(only_z, type="kegg")
kegg_d <- run_enrich(only_d, type="kegg")
cat(sprintf("DENV-specific KEGG terms: %d\n", ifelse(is.null(kegg_d), 0, nrow(as.data.frame(kegg_d)))))

go_dd <- as.data.frame(go_d); go_zd <- as.data.frame(go_z)
if (nrow(go_dd)>10) go_dd <- go_dd[1:10,]
if (nrow(go_zd)>10) go_zd <- go_zd[1:10,]
go_dd$Virus <- "DENV"; go_zd$Virus <- "ZIKV"
go_all <- rbind(go_dd[,c("Description","p.adjust","Virus")], go_zd[,c("Description","p.adjust","Virus")])
top_go <- union(head(go_dd$Description,6), head(go_zd$Description,6))
go_plot <- go_all[go_all$Description %in% top_go,]
go_plot$Description <- factor(go_plot$Description, levels=rev(unique(go_plot$Description)))
f3a <- ggplot(go_plot, aes(-log10(p.adjust), Description, fill=Virus)) +
  geom_col(position="dodge", width=0.7) +
  scale_fill_manual(values=cols_virus) +
  labs(x="-log10(adjusted P)", y=NULL, title="Top GO BP terms in high-viral clusters") +
  theme_pub(11) + theme(axis.text.y=element_text(size=8.5))

kz <- as.data.frame(kegg_z)
key_kegg <- c("TNF signaling pathway",
              "NF-kappa B signaling pathway",
              "PD-L1 expression and PD-1 checkpoint pathway in cancer",
              "C-type lectin receptor signaling pathway",
              "IL-17 signaling pathway")
kz <- kz[kz$Description %in% key_kegg, ]
kz <- kz[order(kz$Count, decreasing=TRUE), ]
kz$Description <- factor(kz$Description, levels=rev(kz$Description))
f3b <- ggplot(kz, aes(Count, Description)) +
  geom_col(fill="#7A5195", width=0.65) +
  geom_text(aes(label=Count), hjust=-0.2, size=4) +
  labs(x="Gene count", y=NULL, title="ZIKV-specific inflammatory KEGG pathways") +
  theme_pub(11) + theme(axis.text.y=element_text(size=9))
save_fig((f3a / f3b) + plot_annotation(tag_levels='A'), "Fig3_GO_KEGG", 9, 9)

# ================= Fig4: immune genes DotPlot + ZIKV UMAP =================
imm <- c("CXCL10","VCAM1","NFKB1","NFKB2","RELB","ICAM1","SOCS2","TNF","IL6",
         "ISG15","MX1","HSPA6","ATF3","DDIT3","TRIB3")
imm <- intersect(imm, rownames(denv))
f4a <- DotPlot(denv, features=imm, idents=hc_d, cols=c("grey90","#F1948A","#C0392B")) +
  RotatedAxis() + ggtitle("DENV high-viral") + theme_pub(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9))
f4b <- DotPlot(zika, features=imm, idents=hc_z, cols=c("grey90","#85C1E9","#1A6FA0")) +
  RotatedAxis() + ggtitle("ZIKV high-viral") + theme_pub(11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9))

feat <- intersect(c("CXCL10","VCAM1","NFKB2","SOCS2"), rownames(zika))
fp <- list()
for (g in feat) {
  fp[[g]] <- FeaturePlot(zika, features=g, cols=c("grey92","#85C1E9","#1A6FA0")) +
    labs(title=g, x=NULL, y=NULL) +
    theme_pub(11) +
    theme(axis.text=element_blank(), axis.ticks=element_blank())
}
f4c <- Reduce(`|`, fp)
save_fig((f4a / f4b / f4c) + plot_annotation(tag_levels='A'), "Fig4_immune_genes", 9, 11)

cat("ALL P1 FIGURES REGENERATED (EMI style)\n")
