#!/usr/bin/env Rscript
# ==============================================================================
# Final cluster assignment with outlier handling
# ==============================================================================

results <- readRDS("04_Figures/03_Consensus_Clustering/ConsensusPlots/consensus_results.rds")

cat("==============================================================================\n")
cat("FINALIZING CLUSTER ASSIGNMENTS\n")
cat("==============================================================================\n\n")

# Based on analysis, k=5 gives 2 major clusters plus small groups
# k=7 adds one more decent sized cluster (5%)
# Strategy: Use k=5 but merge tiny clusters (<2%) into nearest major cluster

optimal_k <- 5

assignments <- results[[optimal_k]]$consensusClass
sizes <- table(assignments)

cat(sprintf("Starting with k=%d:\n", optimal_k))
for (i in 1:optimal_k) {
  cat(sprintf("  Cluster %d: %d samples (%.1f%%)\n",
              i, sizes[i], sizes[i]/sum(sizes)*100))
}
cat("\n")

# Identify tiny clusters (<2% or <15 samples)
tiny_threshold <- max(15, sum(sizes) * 0.02)
tiny_clusters <- which(sizes < tiny_threshold)

cat(sprintf("Tiny clusters to reassign (<%d samples): %s\n\n",
            tiny_threshold, paste(tiny_clusters, collapse=", ")))

# For each tiny cluster, reassign to nearest major cluster
# based on consensus with other clusters
consensus_matrix <- results[[optimal_k]]$consensusMatrix
final_assignments <- assignments

if (length(tiny_clusters) > 0) {
  for (tiny_clust in tiny_clusters) {
    tiny_samples <- which(assignments == tiny_clust)

    # Calculate mean consensus with each other cluster
    consensus_with_others <- numeric(optimal_k)
    for (other_clust in 1:optimal_k) {
      if (other_clust != tiny_clust) {
        other_samples <- which(assignments == other_clust)
        # Mean consensus between tiny cluster and this cluster
        consensus_with_others[other_clust] <- mean(
          consensus_matrix[tiny_samples, other_samples, drop=FALSE]
        )
      }
    }

    # Assign to cluster with highest consensus
    best_match <- which.max(consensus_with_others)
    final_assignments[tiny_samples] <- best_match

    cat(sprintf("Reassigning cluster %d (%d samples) to cluster %d (consensus=%.3f)\n",
                tiny_clust, length(tiny_samples), best_match,
                consensus_with_others[best_match]))
  }
  cat("\n")
}

# Renumber clusters 1, 2, 3, ... (removing gaps)
unique_clusters <- sort(unique(final_assignments))
cluster_map <- setNames(1:length(unique_clusters), unique_clusters)
final_assignments <- cluster_map[as.character(final_assignments)]

final_sizes <- table(final_assignments)
final_k <- length(final_sizes)

cat(sprintf("Final clustering: k=%d\n", final_k))
for (i in 1:final_k) {
  cat(sprintf("  Cluster %d: %d samples (%.1f%%)\n",
              i, final_sizes[i], final_sizes[i]/sum(final_sizes)*100))
}
cat("\n")

# Calculate per-sample consensus scores with final clusters
sample_ids <- names(results[[optimal_k]]$consensusClass)
if (is.null(sample_ids)) {
  sample_ids <- colnames(results[[optimal_k]]$consensusMatrix)
}
sample_consensus <- numeric(length(sample_ids))

for (i in 1:length(sample_ids)) {
  same_cluster <- final_assignments == final_assignments[i]
  if (sum(same_cluster) > 1) {
    sample_consensus[i] <- mean(consensus_matrix[i, same_cluster])
  } else {
    sample_consensus[i] <- 1.0
  }
}

# Create final dataframe
assignments_df <- data.frame(
  sample_id = sample_ids,
  cluster = final_assignments,
  consensus_score = round(sample_consensus, 3),
  stringsAsFactors = FALSE
)

# Sort by cluster then sample ID
assignments_df <- assignments_df[order(assignments_df$cluster, assignments_df$sample_id), ]

# Save
dir.create("03_Results/06_Molecular_Subtypes", recursive = TRUE, showWarnings = FALSE)

write.csv(assignments_df,
          "03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv",
          row.names = FALSE)

cat(sprintf("✓ Saved: sample_cluster_assignments.csv (%d samples in %d clusters)\n\n",
            nrow(assignments_df), final_k))

# Summary metrics
quality_summary <- data.frame(
  optimal_k = final_k,
  n_samples = nrow(assignments_df),
  mean_consensus = round(mean(sample_consensus), 3),
  min_cluster_size = min(final_sizes),
  max_cluster_size = max(final_sizes),
  min_cluster_pct = round(min(final_sizes)/sum(final_sizes)*100, 1),
  max_cluster_pct = round(max(final_sizes)/sum(final_sizes)*100, 1)
)

write.csv(quality_summary,
          "03_Results/06_Molecular_Subtypes/clustering_quality_summary.csv",
          row.names = FALSE)

cat("✓ Saved: clustering_quality_summary.csv\n\n")

cat("==============================================================================\n")
cat("READY TO PROCEED WITH DOWNSTREAM ANALYSES\n")
cat("==============================================================================\n")
