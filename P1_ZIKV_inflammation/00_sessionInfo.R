# 00_sessionInfo.R — Run this to record package versions for reproducibility
# Run at the end of the analysis and save output to session_info.txt

set.seed(123)
sink("P1_ZIKV_inflammation/session_info.txt")
cat("Session info for ZIKV-DENV scRNA reanalysis\n")
cat("===========================================\n\n")
print(sessionInfo())
cat("\n\nKey package versions:\n")
for (pkg in c("Seurat", "dplyr", "ggplot2", "clusterProfiler", "org.Hs.eg.db")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %s: %s\n", pkg, as.character(packageVersion(pkg))))
  }
}
sink()
cat("Saved session_info.txt\n")
