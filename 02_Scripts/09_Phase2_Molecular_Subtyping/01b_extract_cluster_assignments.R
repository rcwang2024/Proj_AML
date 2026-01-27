#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 2.1b: Extract Cluster Assignments
# ==============================================================================
# This script extracts cluster assignments from completed consensus clustering
# ==============================================================================

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("EXTRACTING CLUSTER ASSIGNMENTS FROM CONSENSUS CLUSTERING\n")
cat("==============================================================================\n\n")

# Load the consensus clustering results
results_dir <- "04_Figures/03_Consensus_Clustering/ConsensusPlots"

cat("Loading consensus clustering results...\n\n")

# Read the results object saved by ConsensusClusterPlus
results_file <- file.path(results_dir, "consensus_results.rds")

if (!file.exists(results_file)) {
  cat("ERROR: consensus_results.rds not found. Checking for alternative files...\n")
  print(list.files(results_dir, full.names = TRUE))
  stop("Cannot find consensus clustering results")
}

results <- readRDS(results_file)

# Determine optimal k by examining consensus scores
max_k <- length(results)
cat(sprintf("Found results for k=2 to k=%d\n\n", max_k))

# Calculate mean consensus score for each k
consensus_scores <- numeric(max_k - 1)

for (k in 2:max_k) {
  consensus_matrix <- results[[k]]$consensusMatrix
  consensus_scores[k-1] <- mean(consensus_matrix[lower.tri(consensus_matrix)])
}

# Also calculate delta area (change in consensus)
delta_consensus <- c(NA, diff(consensus_scores))

cat("Consensus scores by k:\n")
quality_table <- data.frame(
  k = 2:max_k,
  mean_consensus = round(consensus_scores, 3),
  delta_consensus = round(delta_consensus, 3)
)
print(quality_table)
cat("\n")

# Select optimal k (largest drop in delta or k=4-6 range typical for AML)
# Rule: Choose k where delta is largest negative in k=3-7 range
# k=4 typically works well for AML molecular subtypes
delta_search <- abs(delta_consensus[2:6])  # k=3 to k=7
optimal_k_idx <- which.max(delta_search)
optimal_k <- optimal_k_idx + 2  # Add 2 because we started at k=3

# Override if clusters are very imbalanced - try k=4 or k=5
# Check if any cluster has <5% of samples
test_k <- optimal_k
for (test_val in c(4, 5, 6)) {
  test_assignments <- results[[test_val]]$consensusClass
  test_sizes <- table(test_assignments)
  min_pct <- min(test_sizes) / sum(test_sizes) * 100
  if (min_pct >= 5) {
    optimal_k <- test_val
    break
  }
}

cat(sprintf("Optimal k selected: %d\n", optimal_k))
cat(sprintf("Mean consensus at k=%d: %.3f\n\n", optimal_k, consensus_scores[optimal_k-1]))

# Extract cluster assignments for optimal k
cluster_assignments <- results[[optimal_k]]$consensusClass
sample_ids <- names(cluster_assignments)

# Calculate per-sample consensus scores
consensus_matrix <- results[[optimal_k]]$consensusMatrix
sample_consensus <- numeric(length(sample_ids))

for (i in 1:length(sample_ids)) {
  # Average consensus with samples in same cluster
  same_cluster <- cluster_assignments == cluster_assignments[i]
  if (sum(same_cluster) > 1) {
    sample_consensus[i] <- mean(consensus_matrix[i, same_cluster])
  } else {
    sample_consensus[i] <- 1.0
  }
}

# Create assignment dataframe
assignments_df <- data.frame(
  sample_id = sample_ids,
  cluster = cluster_assignments,
  consensus_score = round(sample_consensus, 3),
  stringsAsFactors = FALSE
)

# Sort by cluster then sample ID
assignments_df <- assignments_df[order(assignments_df$cluster, assignments_df$sample_id), ]

# Save cluster assignments
dir.create("03_Results/06_Molecular_Subtypes", recursive = TRUE, showWarnings = FALSE)

write.csv(assignments_df,
          "03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv",
          row.names = FALSE)

cat(sprintf("✓ Saved: sample_cluster_assignments.csv (%d samples)\n\n", nrow(assignments_df)))

# Summary by cluster
cat("Cluster sizes:\n")
cluster_summary <- table(assignments_df$cluster)
print(cluster_summary)
cat("\n")

# Save quality metrics summary
quality_summary <- data.frame(
  optimal_k = optimal_k,
  n_samples = nrow(assignments_df),
  mean_consensus = round(consensus_scores[optimal_k-1], 3),
  min_cluster_size = min(cluster_summary),
  max_cluster_size = max(cluster_summary)
)

write.csv(quality_summary,
          "03_Results/06_Molecular_Subtypes/clustering_quality_summary.csv",
          row.names = FALSE)

cat("✓ Saved: clustering_quality_summary.csv\n\n")

# Save all k results for reference
all_k_assignments <- list()
for (k in 2:max_k) {
  all_k_assignments[[paste0("k", k)]] <- data.frame(
    sample_id = names(results[[k]]$consensusClass),
    cluster = results[[k]]$consensusClass,
    stringsAsFactors = FALSE
  )
}

saveRDS(all_k_assignments,
        "03_Results/06_Molecular_Subtypes/all_k_cluster_assignments.rds")

cat("✓ Saved: all_k_cluster_assignments.rds\n\n")

cat("==============================================================================\n")
cat("CLUSTER ASSIGNMENT EXTRACTION COMPLETE\n")
cat("==============================================================================\n\n")

cat(sprintf("Optimal k: %d clusters\n", optimal_k))
cat(sprintf("Total samples: %d\n", nrow(assignments_df)))
cat(sprintf("Mean consensus: %.3f\n\n", consensus_scores[optimal_k-1]))

cat("NEXT STEP: Proceed with Phase 2.2 (Pathway Enrichment)\n\n")
