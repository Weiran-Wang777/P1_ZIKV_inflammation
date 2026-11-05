# Fig 5: Robustness analysis - time-course + viral-load-matched comparison
set.seed(123)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
dnv <- readRDS(file.path(base, "denv_cc.rds"))
zik <- readRDS(file.path(base, "zika_clustered.rds"))

inflammatory_genes <- c("CXCL10","VCAM1","NFKB2","RELB","SOCS2")

mean_expr <- function(obj, gene, cells) {
  if (gene %in% rownames(obj)) {
    mean(GetAssayData(obj, layer="data")[gene, cells], na.rm=TRUE)
  } else { NA }
}

# ===== Panel A: ZIKV inflammatory genes over time =====
meta <- zik@meta.data
uninf <- rownames(meta[meta$MOI == 0, , drop=FALSE])

time_data <- data.frame()
for (tp in c(4,12,24,48)) {
  tp_inf <- rownames(meta[meta$viral_log > 0.5 & meta$time==tp & meta$MOI >= 1, , drop=FALSE])
  for (g in inflammatory_genes) {
    e_inf <- mean_expr(zik, g, tp_inf)
    e_uninf <- mean_expr(zik, g, uninf)
    time_data <- rbind(time_data, data.frame(
      time=tp, gene=g, delta=e_inf - e_uninf, n=length(tp_inf)
    ))
  }
}

pA <- ggplot(time_data, aes(x=factor(time), y=delta, color=gene, group=gene)) +
  geom_line(linewidth=0.8) +
  geom_point(size=2.5) +
  labs(x="Time post-infection (h)", y="Mean expression (delta vs uninfected)",
       title="ZIKV: inflammatory gene induction over time") +
  theme_bw(base_size=11) +
  theme(legend.title=element_blank(),
        plot.title=element_text(face="bold", size=12),
        panel.grid.minor=element_blank())

# ===== Panel B: Viral-load-matched DENV vs ZIKV =====
vl_data <- data.frame()
for (virus_name in c("DENV","ZIKV")) {
  obj <- if (virus_name=="DENV") dnv else zik
  m <- obj@meta.data
  matched <- rownames(m[m$viral_log > 0.5 & m$viral_log >= 7 & m$viral_log <= 9 & m$MOI >= 1, , drop=FALSE])
  for (g in inflammatory_genes) {
    e <- mean_expr(obj, g, matched)
    vl_data <- rbind(vl_data, data.frame(
      virus=virus_name, gene=g, mean_expr=e, n=length(matched)
    ))
  }
}

pB <- ggplot(vl_data, aes(x=gene, y=mean_expr, fill=virus)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.6) +
  labs(x="", y="Mean normalized expression",
       title="Viral-load-matched comparison (viral_log 7-9)") +
  scale_fill_manual(values=c("DENV"="#C0392B", "ZIKV"="#1A6FA0")) +
  theme_bw(base_size=11) +
  theme(legend.title=element_blank(),
        plot.title=element_text(face="bold", size=12),
        axis.text.x=element_text(angle=30, hjust=1),
        panel.grid.minor=element_blank())

combined <- pA + pB + plot_annotation(tag_levels="A")
out_pdf <- "P1_ZIKV_inflammation/figures/Fig5_robustness.pdf"
out_png <- "P1_ZIKV_inflammation/figures/Fig5_robustness.png"
ggsave(out_pdf, combined, width=7.5, height=3.375, device=cairo_pdf)
ggsave(out_png, combined, width=7.5, height=3.375, dpi=600)
cat("Saved:", out_pdf, "\n")
