# ============================================================
# Common header for all analysis scripts
# ============================================================
# NOTE: This file documents the shared helper functions and settings.
# Individual scripts are self-contained and define these functions inline
# so they can be run independently. This file is provided for reference only.
# ============================================================
# Set random seed for reproducibility
set.seed(123)

# Load libraries
suppressMessages({
  library(Seurat)
  library(ggplot2)
})

# ---- Helper: auto-detect high-viral cluster (no hardcoded numbers) ----
# Returns the cluster with highest median viral_log
find_high_viral_cluster <- function(seu_obj) {
  viral_by_cluster <- aggregate(
    seu_obj$viral_log,
    by = list(cluster = seu_obj$seurat_clusters),
    FUN = median, na.rm = TRUE
  )
  high_cluster <- viral_by_cluster$cluster[which.max(viral_by_cluster$x)]
  message("Auto-detected high-viral cluster: ", high_cluster,
          " (median viral_log = ", round(max(viral_by_cluster$x), 2), ")")
  return(high_cluster)
}

# ---- Helper: cell-cycle scoring with HUMAN gene set ----
# Uses Seurat's built-in human cc.genes (correct for Huh7 human cells)
score_cell_cycle <- function(seu_obj) {
  # Use human cell cycle genes (cc.genes is human by default in Seurat)
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  seu_obj <- CellCycleScoring(seu_obj,
                              s.features = s_genes,
                              g2m.features = g2m_genes)
  return(seu_obj)
}

# ---- Helper: infection classification (consistent threshold) ----
# All scripts use viral_log > 0.5 = infected
is_infected <- function(seu_obj) {
  seu_obj$viral_log > 0.5
}
