# ============================================================
# run_all.R — One-click reproduction of P1 (ZIKV inflammation)
# Run this from the repository root:  Rscript run_all.R
# ============================================================
set.seed(123)
cat("========================================\n")
cat("P1: ZIKV-specific NF-kB/TNF inflammatory signature\n")
cat("========================================\n\n")

cat(">>> Step 0: Preprocessing (QC, clustering, cell-cycle scoring) ...\n")
cat("    [Only run if you do NOT have data/denv_cc.rds and data/zika_clustered.rds]\n")
cat("    [Reconstructed from GEO raw matrices (GSE110496) via 00_preprocessing.R]\n")
source("P1_ZIKV_inflammation/00_preprocessing.R")

cat("\n>>> Step 1: QC, clustering, DEG, GO/KEGG, Fig1-4 ...\n")
source("P1_ZIKV_inflammation/01_qc_cluster.R")

cat("\n>>> Step 2: Pseudobulk descriptive analysis ...\n")
source("P1_ZIKV_inflammation/03_pseudobulk.R")

cat("\n>>> Step 3: DESeq2 dispersion diagnostic (exploratory) ...\n")
source("P1_ZIKV_inflammation/03b_deseq2_diagnostic.R")

cat("\n>>> Step 4: Robustness analysis (text output) ...\n")
source("P1_ZIKV_inflammation/04_robustness.R")

cat("\n>>> Step 5: Fig5 robustness figure ...\n")
source("P1_ZIKV_inflammation/04b_fig5_robustness.R")

cat("\n>>> Step 6: Recording session info (after all packages loaded) ...\n")
source("P1_ZIKV_inflammation/00_sessionInfo.R")

cat("\n========================================\n")
cat("ALL P1 ANALYSES COMPLETE.\n")
cat("Figures saved to P1_ZIKV_inflammation/figures/\n")
cat("========================================\n")
