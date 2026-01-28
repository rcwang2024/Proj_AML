# TASK 2: Re-evaluate k=3, k=4, k=5 Clustering Solutions
# Compare different numbers of clusters to validate k=2 is optimal

library(tidyverse)
library(survival)
library(cluster)

cat("=== TASK 2: EVALUATING CLUSTERING SOLUTIONS k=2 to k=5 ===\n\n")

# Load all clustering solutions
cat("Loading clustering solutions...\n")
all_clusters <- readRDS("03_Results/06_Molecular_Subtypes/all_k_cluster_assignments.rds")
cat("Available k values:", paste(names(all_clusters), collapse=", "), "\n\n")

# Load survival data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
cat("Loaded survival data:", nrow(survival_data), "samples\n\n")

# Load expression data for silhouette calculation
expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
cat("Loaded expression data:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n\n")

# Use top 5000 variable genes (same as used for clustering)
library(matrixStats)
gene_mad <- rowMads(as.matrix(expr_data))
top_genes <- order(gene_mad, decreasing = TRUE)[1:5000]
expr_subset <- expr_data[top_genes, ]

cat("Using top 5000 most variable genes for silhouette calculation\n\n")

# Initialize comparison results
comparison_results <- data.frame()

cat("=== COMPARING k=2 THROUGH k=5 ===\n\n")

for (k in 2:5) {
  cat("--- k =", k, "---\n")

  # Get cluster assignments
  clusters_k <- all_clusters[[k]]

  if (is.null(clusters_k)) {
    cat("⚠️ No data for k =", k, "\n\n")
    next
  }

  # 1. Cluster sizes
  cluster_sizes <- table(clusters_k$cluster)
  cat("Cluster sizes:\n")
  print(cluster_sizes)

  min_size <- min(cluster_sizes)
  pct_min <- min_size / sum(cluster_sizes) * 100

  if (pct_min < 5) {
    cat("⚠️ WARNING: Smallest cluster is only", round(pct_min, 1), "% of samples\n")
  }
  cat("\n")

  # 2. Survival differences
  cat("Testing survival differences...\n")
  surv_data_k <- survival_data %>%
    inner_join(clusters_k, by = c("sample_id" = "sample_id")) %>%
    rename(cluster_k = cluster.y) %>%
    select(-cluster.x)

  cat("Matched samples with survival:", nrow(surv_data_k), "\n")

  if (nrow(surv_data_k) > 0) {
    surv_obj <- Surv(time = surv_data_k$OS_months, event = surv_data_k$OS_event)

    # Log-rank test
    survdiff_result <- survdiff(surv_obj ~ cluster_k, data = surv_data_k)
    pvalue_survival <- 1 - pchisq(survdiff_result$chisq,
                                    df = length(unique(surv_data_k$cluster_k)) - 1)

    cat("Survival p-value:", format(pvalue_survival, scientific = TRUE), "\n")

    # Calculate median survival per cluster
    fit_k <- survfit(surv_obj ~ cluster_k, data = surv_data_k)
    medians <- summary(fit_k)$table[, "median"]
    cat("Median survival by cluster (months):\n")
    print(medians)

    # Range of median survival
    survival_range <- max(medians, na.rm = TRUE) - min(medians, na.rm = TRUE)
    cat("Survival range:", round(survival_range, 1), "months\n")
  } else {
    pvalue_survival <- NA
    survival_range <- NA
  }

  cat("\n")

  # 3. Silhouette score (clustering quality)
  cat("Calculating silhouette score...\n")

  # Match samples between expression and clusters
  common_samples <- intersect(colnames(expr_subset), clusters_k$sample_id)
  expr_matched <- expr_subset[, common_samples]
  cluster_vec <- clusters_k$cluster[match(common_samples, clusters_k$sample_id)]

  # Calculate distance matrix
  dist_matrix <- dist(t(expr_matched))

  # Calculate silhouette
  sil <- silhouette(cluster_vec, dist_matrix)
  mean_sil <- mean(sil[, "sil_width"])
  cat("Mean silhouette score:", round(mean_sil, 3), "\n")

  # Check for negative silhouettes (poorly assigned samples)
  neg_sil <- sum(sil[, "sil_width"] < 0)
  pct_neg <- neg_sil / length(cluster_vec) * 100
  cat("Negative silhouettes:", neg_sil, "(", round(pct_neg, 1), "%)\n")

  cat("\n")

  # Store results
  comparison_results <- rbind(comparison_results, data.frame(
    k = k,
    n_clusters = k,
    min_cluster_size = min_size,
    pct_min_cluster = pct_min,
    survival_pvalue = pvalue_survival,
    survival_range_months = survival_range,
    mean_silhouette = mean_sil,
    pct_negative_silhouette = pct_neg,
    stringsAsFactors = FALSE
  ))
}

# Summary comparison
cat("\n=== SUMMARY COMPARISON ===\n")
print(comparison_results, row.names = FALSE)
cat("\n")

# Determine recommendation
cat("=== RECOMMENDATION ===\n\n")

# Criteria:
# 1. No tiny clusters (<5%)
# 2. Significant survival differences (p < 0.05)
# 3. Good silhouette score (>0.1)
# 4. Maximize survival discrimination

valid_k <- comparison_results %>%
  filter(pct_min_cluster >= 5) %>%
  filter(survival_pvalue < 0.05 | is.na(survival_pvalue)) %>%
  filter(mean_silhouette > 0.1)

if (nrow(valid_k) > 0) {
  # Among valid solutions, choose one with best survival discrimination
  best_k <- valid_k %>%
    arrange(desc(survival_range_months)) %>%
    slice(1) %>%
    pull(k)

  cat("Recommended k =", best_k, "\n\n")

  cat("Rationale:\n")
  row_best <- comparison_results %>% filter(k == best_k)
  cat("  - All clusters ≥ 5% of samples: ✓\n")
  cat("  - Survival p-value:", format(row_best$survival_pvalue, scientific = TRUE), "\n")
  cat("  - Survival range:", round(row_best$survival_range_months, 1), "months\n")
  cat("  - Mean silhouette:", round(row_best$mean_silhouette, 3), "\n")

  if (best_k == 2) {
    cat("\n✓ k=2 is OPTIMAL (simplest robust solution)\n")
  } else {
    cat("\n→ k=", best_k, "shows better performance than k=2\n", sep="")
    cat("  Consider using k=", best_k, " for more detailed subtyping\n", sep="")
  }
} else {
  cat("No clear winner based on strict criteria\n")
  cat("k=2 remains most robust choice (balance of simplicity and biology)\n")
}

cat("\n")

# Additional insights
cat("=== ADDITIONAL INSIGHTS ===\n\n")

cat("Silhouette interpretation:\n")
cat("  >0.7 = Strong structure\n")
cat("  >0.5 = Reasonable structure\n")
cat("  >0.25 = Weak structure\n")
cat("  <0.25 = No substantial structure\n\n")

for (i in 1:nrow(comparison_results)) {
  k_val <- comparison_results$k[i]
  sil_val <- comparison_results$mean_silhouette[i]

  if (sil_val > 0.5) {
    cat("k=", k_val, ": Strong to reasonable structure (silhouette=", round(sil_val, 3), ")\n", sep="")
  } else if (sil_val > 0.25) {
    cat("k=", k_val, ": Weak structure (silhouette=", round(sil_val, 3), ")\n", sep="")
  } else {
    cat("k=", k_val, ": Very weak structure (silhouette=", round(sil_val, 3), ")\n", sep="")
  }
}

# Save comparison
write.csv(comparison_results,
          "03_Results/11_Extended_Analysis/clustering_k_comparison.csv",
          row.names = FALSE)

cat("\n✓ Saved: 03_Results/11_Extended_Analysis/clustering_k_comparison.csv\n")

cat("\n### Task 2 COMPLETE ###\n")
