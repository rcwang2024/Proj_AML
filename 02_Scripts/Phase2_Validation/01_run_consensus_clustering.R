#!/usr/bin/env Rscript
# ==============================================================================
# Phase 2.1: Consensus Clustering on Diagnostic Cohort
# ==============================================================================
# Objective:
#   1. Load the diagnostic-only expression matrix.
#   2. Filter for the top 5,000 most variable genes by MAD.
#   3. Perform consensus clustering using ConsensusClusterPlus.
#   4. Extract and save cluster assignments for k=2 through k=10.
# ==============================================================================

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
  library(matrixStats)
  library(tidyverse)
})

set.seed(42)
setwd("d:/Proj_AML")

# Create directories
dir.create("03_Results/06_Molecular_Subtypes", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/03_Consensus_Clustering/ConsensusPlots", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 2.1: CONSENSUS CLUSTERING ON DIAGNOSTIC COHORT\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load and Filter Expression Data
# ------------------------------------------------------------------------------
cat("STEP 1: Loading expression data...\n")
expr_data <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")
cat(sprintf("Loaded expression matrix: %d genes × %d samples\n\n", nrow(expr_data), ncol(expr_data)))

# Calculate MAD for each gene
cat("Calculating MAD per gene...\n")
gene_mads <- rowMads(as.matrix(expr_data))
names(gene_mads) <- rownames(expr_data)

# Get top 5000 variable genes
top_genes <- names(sort(gene_mads, decreasing = TRUE)[1:5000])
writeLines(top_genes, "03_Results/06_Molecular_Subtypes/top5000_variable_genes.txt")
cat("✓ Saved top 5000 variable genes list.\n\n")

# Subset and scale
expr_subset <- expr_data[top_genes, ]
expr_scaled <- t(scale(t(expr_subset)))

# ------------------------------------------------------------------------------
# STEP 2: Run Consensus Clustering
# ------------------------------------------------------------------------------
cat("STEP 2: Running ConsensusClusterPlus (k=2-10, 1000 reps)...\n")
cat("This may take a moment...\n")

# Run ConsensusClusterPlus
consensus_results <- ConsensusClusterPlus(
  d = as.matrix(expr_scaled),
  maxK = 10,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1.0,
  title = "04_Figures/03_Consensus_Clustering/ConsensusPlots",
  clusterAlg = "hc",
  distance = "euclidean",
  writeTable = FALSE,
  plot = "png",
  seed = 42,
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2"
)

# Save raw consensus results
saveRDS(consensus_results, "04_Figures/03_Consensus_Clustering/ConsensusPlots/consensus_results.rds")
cat("✓ Saved raw consensus results (consensus_results.rds).\n\n")

# ------------------------------------------------------------------------------
# STEP 3: Extract and Save Assignments
# ------------------------------------------------------------------------------
cat("STEP 3: Extracting cluster assignments...\n")

all_k_assignments <- list()

for (k in 2:10) {
  classes <- consensus_results[[k]]$consensusClass
  df <- data.frame(
    sample_id = names(classes),
    cluster = as.integer(classes),
    row.names = names(classes),
    stringsAsFactors = FALSE
  )
  all_k_assignments[[paste0("k", k)]] <- df
}

# Save all assignments
saveRDS(all_k_assignments, "03_Results/06_Molecular_Subtypes/all_k_cluster_assignments.rds")
cat("✓ Saved all k cluster assignments (all_k_cluster_assignments.rds).\n")

# Save k=2 assignments as the default classification
k2_df <- all_k_assignments[["k2"]]
write.csv(k2_df, "03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv", row.names = FALSE)
cat("✓ Saved default k=2 cluster assignments (sample_cluster_assignments.csv).\n\n")

cat("==============================================================================\n")
cat("CONSENSUS CLUSTERING COMPLETE\n")
cat("==============================================================================\n")
