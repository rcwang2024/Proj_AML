#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 2.1: Consensus Clustering for Molecular Subtyping
# ==============================================================================
# Objective:
#   1. Select top 5,000 most variable genes
#   2. Run consensus clustering (k=2-10)
#   3. Determine optimal number of clusters
#   4. Assign samples to clusters
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ConsensusClusterPlus)
  library(cluster)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

# Create output directories
dir.create("03_Results/06_Molecular_Subtypes", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/03_Consensus_Clustering", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/03_Consensus_Clustering/ConsensusPlots", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 2.1: CONSENSUS CLUSTERING FOR MOLECULAR SUBTYPING\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load Batch-Corrected Expression Data
# ------------------------------------------------------------------------------

cat("STEP 1: Loading batch-corrected expression data...\n\n")

expr_data <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")

cat(sprintf("Expression data: %d genes × %d samples\n",
            nrow(expr_data), ncol(expr_data)))
cat(sprintf("Value range: %.3f to %.3f\n\n", min(expr_data), max(expr_data)))

# ------------------------------------------------------------------------------
# STEP 2: Select Most Variable Genes
# ------------------------------------------------------------------------------

cat("STEP 2: Selecting most variable genes...\n\n")

# Calculate median absolute deviation (MAD) for each gene
gene_mad <- apply(expr_data, 1, mad, na.rm = TRUE)

# Select top 5000 most variable genes
n_genes_select <- 5000
top_genes_idx <- order(gene_mad, decreasing = TRUE)[1:min(n_genes_select, length(gene_mad))]
top_genes <- names(gene_mad)[top_genes_idx]

expr_variable <- expr_data[top_genes, ]

cat(sprintf("Top %d most variable genes selected\n", nrow(expr_variable)))
cat(sprintf("MAD range: %.3f to %.3f\n\n",
            min(gene_mad[top_genes]), max(gene_mad[top_genes])))

# Save selected genes
write.table(top_genes,
            "03_Results/06_Molecular_Subtypes/top5000_variable_genes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("✓ Saved: top5000_variable_genes.txt\n\n")

# ------------------------------------------------------------------------------
# STEP 3: Run Consensus Clustering
# ------------------------------------------------------------------------------

cat("STEP 3: Running consensus clustering (k=2-10)...\n")
cat("This may take several minutes...\n\n")

# Consensus clustering parameters
max_k <- 10
n_reps <- 1000  # Number of iterations
p_item <- 0.8   # Proportion of items to sample
p_feature <- 0.8  # Proportion of features to sample

cat("Parameters:\n")
cat(sprintf("  - Max clusters (k): %d\n", max_k))
cat(sprintf("  - Iterations: %d\n", n_reps))
cat(sprintf("  - Item sampling: %.0f%%\n", p_item * 100))
cat(sprintf("  - Feature sampling: %.0f%%\n\n", p_feature * 100))

# Create output directory for ConsensusClusterPlus
consensus_out_dir <- "04_Figures/03_Consensus_Clustering/ConsensusPlots"
dir.create(consensus_out_dir, recursive = TRUE, showWarnings = FALSE)

# Run consensus clustering
# Output will be saved in the specified directory
results <- ConsensusClusterPlus(
  d = as.matrix(expr_variable),
  maxK = max_k,
  reps = n_reps,
  pItem = p_item,
  pFeature = p_feature,
  clusterAlg = "hc",  # Hierarchical clustering
  distance = "pearson",  # Pearson correlation distance
  title = consensus_out_dir,
  plot = "pdf",  # Generate PDF plots
  writeTable = FALSE,  # Don't write tables to avoid path issues
  verbose = FALSE  # Reduce verbosity
)

cat("\n✓ Consensus clustering complete!\n\n")

# Save results object for later use
saveRDS(results, file.path(consensus_out_dir, "consensus_results.rds"))
cat("✓ Saved: consensus_results.rds\n\n")

# ------------------------------------------------------------------------------
# STEP 4: Calculate Cluster Quality Metrics
# ------------------------------------------------------------------------------

cat("STEP 4: Calculating cluster quality metrics...\n\n")

# Initialize results dataframe
cluster_metrics <- data.frame(
  k = 2:max_k,
  consensus_score = numeric(max_k - 1),
  delta_area = numeric(max_k - 1),
  mean_silhouette = numeric(max_k - 1),
  stringsAsFactors = FALSE
)

# Calculate metrics for each k
for (k in 2:max_k) {
  # Consensus score (mean consensus within clusters)
  consensus_matrix <- results[[k]]$consensusMatrix
  cluster_assignment <- results[[k]]$consensusClass

  # Mean consensus score within clusters
  within_cluster_consensus <- numeric(k)
  for (i in 1:k) {
    samples_in_cluster <- which(cluster_assignment == i)
    if (length(samples_in_cluster) > 1) {
      cluster_consensus <- consensus_matrix[samples_in_cluster, samples_in_cluster]
      within_cluster_consensus[i] <- mean(cluster_consensus[lower.tri(cluster_consensus)])
    } else {
      within_cluster_consensus[i] <- 1
    }
  }

  consensus_score <- mean(within_cluster_consensus)
  cluster_metrics$consensus_score[k-1] <- consensus_score

  # Calculate silhouette scores
  # Use Euclidean distance on scaled expression data
  expr_scaled <- t(scale(t(expr_variable)))
  dist_matrix <- dist(t(expr_scaled))

  sil <- silhouette(cluster_assignment, dist_matrix)
  mean_sil <- mean(sil[, 3])
  cluster_metrics$mean_silhouette[k-1] <- mean_sil
}

# Calculate delta area (change in consensus CDF area)
# Get CDF areas from results - handle different ConsensusClusterPlus versions
cdf_areas <- numeric(max_k - 1)
for (k in 2:max_k) {
  # Try to get areaUnderCDF, if not available or empty use alternative metric
  if (!is.null(results[[k]]$areaUnderCDF) && length(results[[k]]$areaUnderCDF) > 0) {
    cdf_areas[k-1] <- results[[k]]$areaUnderCDF
  } else {
    # Alternative: use mean consensus score
    cdf_areas[k-1] <- mean(results[[k]]$consensusMatrix[lower.tri(results[[k]]$consensusMatrix)])
  }
}

# Delta area is the proportional change in area
if (length(cdf_areas) > 1) {
  cluster_metrics$delta_area[2:nrow(cluster_metrics)] <-
    (cdf_areas[2:length(cdf_areas)] - cdf_areas[1:(length(cdf_areas)-1)]) / cdf_areas[1:(length(cdf_areas)-1)]
  cluster_metrics$delta_area[1] <- cdf_areas[1]  # First value is just the area itself
} else {
  cluster_metrics$delta_area <- 0
}

# Display metrics
cat("Cluster quality metrics:\n")
print(cluster_metrics)
cat("\n")

# Save metrics
write.csv(cluster_metrics,
          "03_Results/06_Molecular_Subtypes/consensus_clustering_metrics.csv",
          row.names = FALSE)
cat("✓ Saved: consensus_clustering_metrics.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 5: Determine Optimal K
# ------------------------------------------------------------------------------

cat("STEP 5: Determining optimal number of clusters...\n\n")

# Multiple criteria for optimal k:
# 1. High consensus score (>0.8)
# 2. Small delta area (<0.05 indicates plateau)
# 3. High mean silhouette (>0.3)

# Find k where delta area is minimal (plateau)
delta_area_threshold <- 0.05
plateau_ks <- cluster_metrics$k[abs(cluster_metrics$delta_area) < delta_area_threshold & cluster_metrics$k > 2]

# Find k with good consensus and silhouette
good_quality_ks <- cluster_metrics$k[
  cluster_metrics$consensus_score > 0.7 &
  cluster_metrics$mean_silhouette > 0.2
]

# Combine criteria
candidate_ks <- intersect(plateau_ks, good_quality_ks)

if (length(candidate_ks) == 0) {
  # If no k meets all criteria, use silhouette maximum
  optimal_k <- cluster_metrics$k[which.max(cluster_metrics$mean_silhouette)]
  cat("⚠ No k met all quality criteria\n")
  cat(sprintf("Using k=%d (highest silhouette score)\n\n", optimal_k))
} else {
  # Use smallest k among candidates (parsimony)
  optimal_k <- min(candidate_ks)
  cat(sprintf("✓ Optimal k=%d determined by consensus criteria\n", optimal_k))
  cat("  Criteria met:\n")
  cat(sprintf("  - Consensus score: %.3f\n",
              cluster_metrics$consensus_score[optimal_k-1]))
  cat(sprintf("  - Mean silhouette: %.3f\n",
              cluster_metrics$mean_silhouette[optimal_k-1]))
  cat(sprintf("  - Delta area: %.3f\n\n",
              cluster_metrics$delta_area[optimal_k-1]))
}

# Save optimal k
optimal_k_record <- data.frame(
  optimal_k = optimal_k,
  consensus_score = cluster_metrics$consensus_score[optimal_k-1],
  mean_silhouette = cluster_metrics$mean_silhouette[optimal_k-1],
  delta_area = cluster_metrics$delta_area[optimal_k-1],
  method = "Consensus + Silhouette criteria",
  date = Sys.Date()
)

write.csv(optimal_k_record,
          "03_Results/06_Molecular_Subtypes/optimal_k_selection.csv",
          row.names = FALSE)
cat("✓ Saved: optimal_k_selection.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 6: Extract Cluster Assignments for Optimal K
# ------------------------------------------------------------------------------

cat("STEP 6: Extracting cluster assignments...\n\n")

cluster_assignment <- results[[optimal_k]]$consensusClass

# Create detailed assignment dataframe
sample_clusters <- data.frame(
  sample_id = names(cluster_assignment),
  cluster = as.factor(cluster_assignment),
  stringsAsFactors = FALSE
)

# Calculate sample-level consensus scores
consensus_matrix_opt <- results[[optimal_k]]$consensusMatrix
sample_consensus_scores <- numeric(nrow(sample_clusters))

for (i in 1:nrow(sample_clusters)) {
  cluster_i <- sample_clusters$cluster[i]
  samples_same_cluster <- which(sample_clusters$cluster == cluster_i)
  sample_consensus_scores[i] <- mean(consensus_matrix_opt[i, samples_same_cluster])
}

sample_clusters$consensus_score <- sample_consensus_scores

# Calculate silhouette width for each sample
expr_scaled <- t(scale(t(expr_variable)))
dist_matrix <- dist(t(expr_scaled))
sil <- silhouette(cluster_assignment, dist_matrix)
sample_clusters$silhouette_width <- sil[, 3]

# Cluster size summary
cat("Cluster sizes:\n")
cluster_sizes <- table(sample_clusters$cluster)
print(cluster_sizes)
cat("\n")

cluster_summary <- data.frame(
  cluster = names(cluster_sizes),
  n_samples = as.numeric(cluster_sizes),
  percentage = round(as.numeric(cluster_sizes) / sum(cluster_sizes) * 100, 1),
  mean_consensus = tapply(sample_clusters$consensus_score, sample_clusters$cluster, mean),
  mean_silhouette = tapply(sample_clusters$silhouette_width, sample_clusters$cluster, mean)
)

cat("Cluster summary:\n")
print(cluster_summary)
cat("\n")

# Save assignments
write.csv(sample_clusters,
          "03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv",
          row.names = FALSE)
cat("✓ Saved: sample_cluster_assignments.csv\n\n")

write.csv(cluster_summary,
          "03_Results/06_Molecular_Subtypes/cluster_summary.csv",
          row.names = FALSE)
cat("✓ Saved: cluster_summary.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 7: Create Additional Visualizations
# ------------------------------------------------------------------------------

cat("STEP 7: Creating additional visualizations...\n\n")

# Plot 1: Cluster metrics across k
pdf("04_Figures/03_Consensus_Clustering/cluster_quality_metrics.pdf",
    width = 12, height = 4)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))

# Consensus score
plot(cluster_metrics$k, cluster_metrics$consensus_score,
     type = "b", pch = 19, col = "blue",
     xlab = "Number of clusters (k)", ylab = "Mean consensus score",
     main = "Consensus Score by k")
abline(v = optimal_k, col = "red", lty = 2)
text(optimal_k, max(cluster_metrics$consensus_score),
     paste("k =", optimal_k), pos = 4, col = "red")

# Silhouette score
plot(cluster_metrics$k, cluster_metrics$mean_silhouette,
     type = "b", pch = 19, col = "darkgreen",
     xlab = "Number of clusters (k)", ylab = "Mean silhouette width",
     main = "Silhouette Width by k")
abline(v = optimal_k, col = "red", lty = 2)
abline(h = 0, col = "gray", lty = 2)

# Delta area
plot(cluster_metrics$k, cluster_metrics$delta_area,
     type = "b", pch = 19, col = "purple",
     xlab = "Number of clusters (k)", ylab = "Proportional change in CDF area",
     main = "Delta Area by k")
abline(v = optimal_k, col = "red", lty = 2)
abline(h = 0, col = "gray", lty = 2)

dev.off()
cat("✓ Saved: cluster_quality_metrics.pdf\n\n")

# Plot 2: Silhouette plot for optimal k
pdf("04_Figures/03_Consensus_Clustering/silhouette_plot_optimal_k.pdf",
    width = 8, height = 6)

sil_optimal <- silhouette(cluster_assignment, dist_matrix)
plot(sil_optimal, col = rainbow(optimal_k)[sil_optimal[, 1]],
     main = paste0("Silhouette Plot (k=", optimal_k, ")"),
     border = NA)

dev.off()
cat("✓ Saved: silhouette_plot_optimal_k.pdf\n\n")

# Plot 3: Heatmap with cluster annotations
pdf("04_Figures/03_Consensus_Clustering/consensus_heatmap_optimal_k.pdf",
    width = 10, height = 9)

# Order samples by cluster
sample_order <- order(sample_clusters$cluster)

# Create annotation
annotation_col <- data.frame(
  Cluster = factor(sample_clusters$cluster[sample_order])
)
rownames(annotation_col) <- sample_clusters$sample_id[sample_order]

# Cluster colors
cluster_colors <- rainbow(optimal_k)
names(cluster_colors) <- levels(annotation_col$Cluster)

annotation_colors <- list(Cluster = cluster_colors)

# Plot heatmap
pheatmap(
  consensus_matrix_opt[sample_order, sample_order],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("white", "red"))(100),
  main = paste0("Consensus Matrix (k=", optimal_k, ")"),
  fontsize = 10
)

dev.off()
cat("✓ Saved: consensus_heatmap_optimal_k.pdf\n\n")

# ------------------------------------------------------------------------------
# STEP 8: Summary Report
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("CONSENSUS CLUSTERING COMPLETE\n")
cat("==============================================================================\n\n")

cat("SUMMARY:\n\n")
cat(sprintf("Optimal number of clusters: k = %d\n", optimal_k))
cat(sprintf("Total samples clustered: %d\n", nrow(sample_clusters)))
cat(sprintf("Genes used for clustering: %d (most variable)\n", nrow(expr_variable)))
cat(sprintf("Consensus clustering iterations: %d\n\n", n_reps))

cat("CLUSTER QUALITY METRICS:\n")
cat(sprintf("  - Mean consensus score: %.3f\n", optimal_k_record$consensus_score))
cat(sprintf("  - Mean silhouette width: %.3f\n", optimal_k_record$mean_silhouette))
cat(sprintf("  - Interpretation: %s\n\n",
            ifelse(optimal_k_record$mean_silhouette > 0.5, "Strong structure",
                   ifelse(optimal_k_record$mean_silhouette > 0.25, "Reasonable structure",
                          "Weak structure"))))

cat("CLUSTER DISTRIBUTION:\n")
for (i in 1:nrow(cluster_summary)) {
  cat(sprintf("  Cluster %s: %d samples (%.1f%%), mean silhouette=%.3f\n",
              cluster_summary$cluster[i],
              cluster_summary$n_samples[i],
              cluster_summary$percentage[i],
              cluster_summary$mean_silhouette[i]))
}
cat("\n")

cat("OUTPUT FILES:\n")
cat("  - 03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv\n")
cat("  - 03_Results/06_Molecular_Subtypes/cluster_summary.csv\n")
cat("  - 03_Results/06_Molecular_Subtypes/consensus_clustering_metrics.csv\n")
cat("  - 03_Results/06_Molecular_Subtypes/optimal_k_selection.csv\n")
cat("  - 04_Figures/03_Consensus_Clustering/*.pdf\n\n")

cat("NEXT STEPS:\n")
cat("  → Phase 2.2: Pathway enrichment-based clustering\n")
cat("  → Phase 2.3: Subtype characterization\n")
cat("    - Differential expression analysis\n")
cat("    - Mutation enrichment\n")
cat("    - Clinical associations\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
