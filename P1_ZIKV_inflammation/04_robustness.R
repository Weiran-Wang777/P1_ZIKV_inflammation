# P1 robustness analysis: time-point stratification + viral-load matching + uninfected controls
set.seed(123)
suppressPackageStartupMessages(library(Seurat))

base <- "data"
if (!file.exists(file.path(base, "denv_cc.rds")))
  stop("Error: data/denv_cc.rds not found. See README for data download instructions.")
if (!file.exists(file.path(base, "zika_clustered.rds")))
  stop("Error: data/zika_clustered.rds not found. See README for data download instructions.")
dnv <- readRDS(file.path(base, "denv_cc.rds"))
zik <- readRDS(file.path(base, "zika_clustered.rds"))

inflammatory_genes <- c("CXCL10","VCAM1","NFKB2","RELB","SOCS2")
er_stress_genes <- c("HSPA6","ATF3","DDIT3","HERPUD1")

mean_expr <- function(obj, gene, cells) {
  if (gene %in% rownames(obj)) {
    mean(GetAssayData(obj, layer="data")[gene, cells], na.rm=TRUE)
  } else { NA }
}

cat("========================================\n")
cat("ROBUSTNESS ANALYSIS FOR P1\n")
cat("========================================\n\n")

# ---- 1. Time-point stratification ----
# Note: MOI >= 1 combines both MOI=1 and MOI=10 infected cells for DENV.
cat("=== 1. Time-point stratification (infected cells only, MOI>=1) ===\n\n")
for (virus_name in c("DENV","ZIKV")) {
  obj <- if (virus_name=="DENV") dnv else zik
  meta <- obj@meta.data
  infected <- rownames(meta[meta$viral_log > 0.5 & meta$MOI >= 1, , drop=FALSE])
  uninfected <- rownames(meta[meta$MOI == 0, , drop=FALSE])
  cat(sprintf("--- %s: %d infected cells, %d uninfected ---\n", virus_name, length(infected), length(uninfected)))
  for (tp in c(4,12,24,48)) {
    tp_cells <- intersect(infected, rownames(meta[meta$time==tp, , drop=FALSE]))
    cat(sprintf("  Time %dh: %d cells\n", tp, length(tp_cells)))
    for (g in c(inflammatory_genes, er_stress_genes)) {
      if (g %in% rownames(obj)) {
        e_inf <- mean_expr(obj, g, tp_cells)
        e_uninf <- mean_expr(obj, g, uninfected)
        cat(sprintf("    %-10s infected=%.3f  uninfected=%.3f  delta=%.3f\n", g, e_inf, e_uninf, e_inf-e_uninf))
      }
    }
  }
  cat("\n")
}

# ---- 2. Viral-load-matched comparison ----
cat("=== 2. Viral-load-matched DENV vs ZIKV comparison ===\n\n")
for (virus_name in c("DENV","ZIKV")) {
  obj <- if (virus_name=="DENV") dnv else zik
  meta <- obj@meta.data
  matched <- rownames(meta[meta$viral_log > 0.5 & meta$viral_log >= 7 & meta$viral_log <= 9 & meta$MOI >= 1, , drop=FALSE])
  cat(sprintf("--- %s: %d cells with viral_log in [7,9] ---\n", virus_name, length(matched)))
  for (g in c(inflammatory_genes, er_stress_genes)) {
    if (g %in% rownames(obj)) {
      e <- mean_expr(obj, g, matched)
      cat(sprintf("    %-10s mean=%.3f\n", g, e))
    }
  }
  cat("\n")
}

# ---- 3. Uninfected control ----
cat("=== 3. Uninfected (MOI=0) baseline ===\n\n")
for (virus_name in c("DENV","ZIKV")) {
  obj <- if (virus_name=="DENV") dnv else zik
  meta <- obj@meta.data
  uninf <- rownames(meta[meta$MOI == 0, , drop=FALSE])
  cat(sprintf("--- %s: %d uninfected cells ---\n", virus_name, length(uninf)))
  for (g in c(inflammatory_genes, er_stress_genes)) {
    if (g %in% rownames(obj)) {
      e <- mean_expr(obj, g, uninf)
      cat(sprintf("    %-10s mean=%.3f\n", g, e))
    }
  }
  cat("\n")
}

# ---- 4. ZIKV inflammatory genes by time point ----
cat("=== 4. ZIKV inflammatory genes by time point (vs ZIKV uninfected) ===\n\n")
meta <- zik@meta.data
uninf <- rownames(meta[meta$MOI == 0, , drop=FALSE])
for (tp in c(4,12,24,48)) {
  tp_inf <- rownames(meta[meta$viral_log > 0.5 & meta$time==tp & meta$MOI >= 1, , drop=FALSE])
  cat(sprintf("ZIKV %dh (%d cells):\n", tp, length(tp_inf)))
  for (g in inflammatory_genes) {
    if (g %in% rownames(zik)) {
      e_inf <- mean_expr(zik, g, tp_inf)
      e_uninf <- mean_expr(zik, g, uninf)
      cat(sprintf("  %-10s inf=%.3f  uninf=%.3f  delta=%.3f\n", g, e_inf, e_uninf, e_inf-e_uninf))
    }
  }
}

# ---- 5. DENV MOI comparison ----
cat("\n=== 5. DENV MOI=1 vs MOI=10 (high-viral cells) ===\n\n")
meta <- dnv@meta.data
for (moi in c(1,10)) {
  moi_cells <- rownames(meta[meta$viral_log > 0.5 & meta$MOI == moi, , drop=FALSE])
  cat(sprintf("DENV MOI=%d (%d infected cells):\n", moi, length(moi_cells)))
  for (g in inflammatory_genes) {
    if (g %in% rownames(dnv)) {
      e <- mean_expr(dnv, g, moi_cells)
      cat(sprintf("  %-10s mean=%.3f\n", g, e))
    }
  }
}

cat("\n=== DONE ===\n")
