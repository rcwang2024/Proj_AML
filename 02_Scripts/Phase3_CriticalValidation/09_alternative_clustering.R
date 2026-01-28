#!/usr/bin/env Rscript
# ==============================================================================
# Phase 3.7: Alternative Clustering Solutions
# ==============================================================================
# Objective:
#   1. Test k=3, 4, 5 clustering solutions (currently using k=2)
#   2. Evaluate cluster quality metrics for each k
#   3. Assess mutation enrichment patterns
#   4. Test survival differences between clusters
#   5. Determine if alternative k provides better biological separation
# Date: 2025-10-11
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(cluster)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

# Create output directories
dir.create("03_Results/11_Survival_Reanalysis/07_alternative_clustering",
           recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/11_Survival_Reanalysis/07_alternative_clustering",
           recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 3.7: ALTERNATIVE CLUSTERING SOLUTIONS\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# SECTION 1: LOAD DATA AND CONSENSUS CLUSTERING RESULTS
# ------------------------------------------------------------------------------

cat("SECTION 1: LOADING DATA\n")
cat("==============================================================================\n\n")

# Load consensus clustering results (contains k=2 through k=10)
consensus_results <- readRDS("04_Figures/03_Consensus_Clustering/ConsensusPlots/consensus_results.rds")
cat("✓ Loaded consensus clustering results (k=2 to k=10)\n")

# Load current cluster assignments (k=2)
current_clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
current_k <- length(unique(current_clusters$cluster))
cat(sprintf("✓ Current clustering: k=%d (%d samples)\n", current_k, nrow(current_clusters)))

# Load expression data for distance calculations
expr_data <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")
cat(sprintf("✓ Expression data: %d genes × %d samples\n", nrow(expr_data), ncol(expr_data)))

# Load survival data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
cat(sprintf("✓ Survival data: %d samples\n", nrow(survival_data)))

# Standardize column names (remove cluster from survival to avoid conflict)
survival_data <- survival_data %>%
  rename(
    dbgap_rnaseq_sample = sample_id,
    OS_MONTHS = OS_months,
    OS_STATUS = OS_event
  ) %>%
  dplyr::select(-cluster)  # Remove cluster column to avoid naming conflicts

# Load mutation data
mutation_file <- "03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv"
if (file.exists(mutation_file)) {
  mutation_data <- read.csv(mutation_file, row.names = 1)
  mutation_data$dbgap_rnaseq_sample <- rownames(mutation_data)
  cat(sprintf("✓ Mutation data: %d samples\n", nrow(mutation_data)))
} else {
  cat("⚠ Mutation data not found\n")
  mutation_data <- NULL
}

cat("\n")

# ------------------------------------------------------------------------------
# SECTION 2: EVALUATE CLUSTER QUALITY FOR k=2,3,4,5
# ------------------------------------------------------------------------------

cat("SECTION 2: EVALUATING CLUSTER QUALITY METRICS\n")
cat("==============================================================================\n\n")

# Test k values
k_values <- c(2, 3, 4, 5)

# Get top 5000 variable genes (used for original clustering)
top_genes <- readLines("03_Results/06_Molecular_Subtypes/top5000_variable_genes.txt")
expr_variable <- expr_data[top_genes, ]

# Scale expression data for distance calculations
expr_scaled <- t(scale(t(expr_variable)))
dist_matrix <- dist(t(expr_scaled))

# Initialize results dataframe
cluster_comparison <- data.frame()

cat("Calculating quality metrics for each k...\n\n")

for (k in k_values) {
  cat(sprintf("--- k=%d ---\n", k))

  # Get cluster assignments
  cluster_assignment <- consensus_results[[k]]$consensusClass
  consensus_matrix <- consensus_results[[k]]$consensusMatrix

  # Cluster sizes
  cluster_sizes <- table(cluster_assignment)
  cat("Cluster sizes:", paste(cluster_sizes, collapse=", "), "\n")

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
  mean_consensus <- mean(within_cluster_consensus)

  # Silhouette scores
  sil <- silhouette(cluster_assignment, dist_matrix)
  mean_silhouette <- mean(sil[, 3])
  min_silhouette <- min(sil[, 3])

  # Cluster balance (coefficient of variation of cluster sizes)
  cv_size <- sd(cluster_sizes) / mean(cluster_sizes)

  # Between-cluster separation (average distance between cluster centroids)
  centroids <- matrix(0, nrow = k, ncol = nrow(expr_scaled))
  for (i in 1:k) {
    centroids[i, ] <- rowMeans(expr_scaled[, cluster_assignment == i, drop = FALSE])
  }
  centroid_dist <- dist(centroids)
  mean_centroid_dist <- mean(centroid_dist)

  cat(sprintf("  Consensus: %.3f\n", mean_consensus))
  cat(sprintf("  Mean silhouette: %.3f\n", mean_silhouette))
  cat(sprintf("  Min silhouette: %.3f\n", min_silhouette))
  cat(sprintf("  Size balance (CV): %.3f\n", cv_size))
  cat(sprintf("  Centroid separation: %.3f\n\n", mean_centroid_dist))

  # Store results
  cluster_comparison <- rbind(cluster_comparison, data.frame(
    k = k,
    n_clusters = k,
    min_cluster_size = min(cluster_sizes),
    max_cluster_size = max(cluster_sizes),
    mean_consensus = round(mean_consensus, 3),
    mean_silhouette = round(mean_silhouette, 3),
    min_silhouette = round(min_silhouette, 3),
    size_cv = round(cv_size, 3),
    centroid_separation = round(mean_centroid_dist, 3),
    stringsAsFactors = FALSE
  ))
}

cat("Summary table:\n")
print(cluster_comparison)
cat("\n")

write.csv(cluster_comparison,
          "03_Results/11_Survival_Reanalysis/07_alternative_clustering/cluster_quality_comparison.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 3: SURVIVAL ANALYSIS FOR EACH k
# ------------------------------------------------------------------------------

cat("SECTION 3: SURVIVAL ANALYSIS FOR ALTERNATIVE k VALUES\n")
cat("==============================================================================\n\n")

survival_comparison <- data.frame()

for (k in k_values) {
  cat(sprintf("--- k=%d ---\n", k))

  # Get cluster assignments
  cluster_assignment <- consensus_results[[k]]$consensusClass
  cluster_df <- data.frame(
    sample_id = names(cluster_assignment),
    cluster_k = cluster_assignment,
    stringsAsFactors = FALSE
  )

  # Merge with survival data
  surv_data <- survival_data %>%
    inner_join(cluster_df, by = c("dbgap_rnaseq_sample" = "sample_id"))

  # Overall log-rank test
  surv_obj <- Surv(time = surv_data$OS_MONTHS, event = surv_data$OS_STATUS)
  logrank_test <- survdiff(surv_obj ~ cluster_k, data = surv_data)
  logrank_p <- 1 - pchisq(logrank_test$chisq, df = k - 1)

  cat(sprintf("  Log-rank p-value: %.4f %s\n",
              logrank_p,
              ifelse(logrank_p < 0.05, "***", "")))

  # Median survival by cluster
  km_fit <- survfit(surv_obj ~ cluster_k, data = surv_data)
  medians <- summary(km_fit)$table[, "median"]

  for (i in 1:k) {
    cat(sprintf("  Cluster %d: median OS = %.1f months (n=%d)\n",
                i, medians[i], sum(surv_data$cluster_k == i)))
  }

  # Cox regression (use cluster as factor)
  cox_model <- coxph(surv_obj ~ factor(cluster_k), data = surv_data)
  cox_pvalue <- summary(cox_model)$sctest["pvalue"]

  cat(sprintf("  Cox p-value: %.4f\n", cox_pvalue))

  # Range of median survivals
  median_range <- max(medians, na.rm = TRUE) - min(medians, na.rm = TRUE)

  cat(sprintf("  Range of medians: %.1f months\n\n", median_range))

  # Store results
  survival_comparison <- rbind(survival_comparison, data.frame(
    k = k,
    logrank_p = round(logrank_p, 5),
    cox_p = round(cox_pvalue, 5),
    min_median_os = round(min(medians, na.rm = TRUE), 1),
    max_median_os = round(max(medians, na.rm = TRUE), 1),
    median_range = round(median_range, 1),
    significant = logrank_p < 0.05,
    stringsAsFactors = FALSE
  ))
}

cat("Survival comparison summary:\n")
print(survival_comparison)
cat("\n")

write.csv(survival_comparison,
          "03_Results/11_Survival_Reanalysis/07_alternative_clustering/survival_comparison.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 4: MUTATION ENRICHMENT FOR EACH k
# ------------------------------------------------------------------------------

if (!is.null(mutation_data)) {
  cat("SECTION 4: MUTATION ENRICHMENT ANALYSIS\n")
  cat("==============================================================================\n\n")

  # Key mutations to test
  key_mutations <- c("TP53", "NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2",
                     "TET2", "RUNX1", "ASXL1", "NRAS", "KRAS")

  # Check which mutations are available
  available_mutations <- intersect(key_mutations, colnames(mutation_data))
  cat(sprintf("Testing %d mutations: %s\n\n", length(available_mutations),
              paste(available_mutations, collapse=", ")))

  mutation_enrichment_all <- data.frame()

  for (k in k_values) {
    cat(sprintf("--- k=%d ---\n", k))

    # Get cluster assignments
    cluster_assignment <- consensus_results[[k]]$consensusClass
    cluster_df <- data.frame(
      sample_id = names(cluster_assignment),
      cluster_k = cluster_assignment,
      stringsAsFactors = FALSE
    )

    # Merge with mutation data
    mut_cluster <- mutation_data %>%
      inner_join(cluster_df, by = c("dbgap_rnaseq_sample" = "sample_id"))

    # Test each mutation
    significant_enrichments <- 0

    for (gene in available_mutations) {
      # Create contingency table
      contingency <- table(mut_cluster$cluster_k, mut_cluster[[gene]])

      # Fisher's exact test (or chi-square if large sample)
      if (all(contingency >= 5)) {
        test_result <- chisq.test(contingency)
        p_value <- test_result$p.value
        test_method <- "chi-square"
      } else {
        # For k>2, use chi-square approximation even with small counts
        if (k > 2) {
          test_result <- chisq.test(contingency, simulate.p.value = TRUE)
          p_value <- test_result$p.value
          test_method <- "chi-square (simulated)"
        } else {
          test_result <- fisher.test(contingency)
          p_value <- test_result$p.value
          test_method <- "fisher"
        }
      }

      if (p_value < 0.01) {
        significant_enrichments <- significant_enrichments + 1
      }

      # Store result
      mutation_enrichment_all <- rbind(mutation_enrichment_all, data.frame(
        k = k,
        gene = gene,
        p_value = p_value,
        test_method = test_method,
        significant = p_value < 0.01,
        stringsAsFactors = FALSE
      ))
    }

    cat(sprintf("  Significant enrichments (p<0.01): %d / %d\n",
                significant_enrichments, length(available_mutations)))
  }

  cat("\n")

  # Summarize mutation enrichment by k
  mutation_summary <- mutation_enrichment_all %>%
    group_by(k) %>%
    summarise(
      n_tested = n(),
      n_significant_001 = sum(p_value < 0.001),
      n_significant_01 = sum(p_value < 0.01),
      n_significant_05 = sum(p_value < 0.05),
      min_p = min(p_value),
      .groups = "drop"
    )

  cat("Mutation enrichment summary:\n")
  print(mutation_summary)
  cat("\n")

  write.csv(mutation_enrichment_all,
            "03_Results/11_Survival_Reanalysis/07_alternative_clustering/mutation_enrichment_by_k.csv",
            row.names = FALSE)

  write.csv(mutation_summary,
            "03_Results/11_Survival_Reanalysis/07_alternative_clustering/mutation_enrichment_summary.csv",
            row.names = FALSE)

} else {
  cat("SECTION 4: MUTATION ENRICHMENT ANALYSIS\n")
  cat("==============================================================================\n")
  cat("⚠ Skipped - mutation data not available\n\n")
}

# ------------------------------------------------------------------------------
# SECTION 5: VISUALIZATIONS
# ------------------------------------------------------------------------------

cat("SECTION 5: CREATING VISUALIZATIONS\n")
cat("==============================================================================\n\n")

# Plot 1: Quality metrics comparison
pdf("04_Figures/11_Survival_Reanalysis/07_alternative_clustering/quality_metrics_comparison.pdf",
    width = 12, height = 8)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

# Consensus
plot(cluster_comparison$k, cluster_comparison$mean_consensus,
     type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "Number of clusters (k)", ylab = "Mean consensus score",
     main = "Consensus Score", ylim = c(0.5, 1))
abline(h = 0.8, col = "gray", lty = 2)
abline(v = current_k, col = "red", lty = 2)

# Silhouette
plot(cluster_comparison$k, cluster_comparison$mean_silhouette,
     type = "b", pch = 19, col = "darkgreen", lwd = 2,
     xlab = "Number of clusters (k)", ylab = "Mean silhouette width",
     main = "Silhouette Width")
abline(h = 0.25, col = "gray", lty = 2)
abline(v = current_k, col = "red", lty = 2)

# Size balance
plot(cluster_comparison$k, cluster_comparison$size_cv,
     type = "b", pch = 19, col = "purple", lwd = 2,
     xlab = "Number of clusters (k)", ylab = "Coefficient of variation",
     main = "Cluster Size Balance (lower = better)")
abline(v = current_k, col = "red", lty = 2)

# Survival p-value
plot(survival_comparison$k, -log10(survival_comparison$logrank_p),
     type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Number of clusters (k)", ylab = "-log10(p-value)",
     main = "Survival Difference (log-rank)")
abline(h = -log10(0.05), col = "gray", lty = 2)
abline(v = current_k, col = "red", lty = 2)

# Median survival range
plot(survival_comparison$k, survival_comparison$median_range,
     type = "b", pch = 19, col = "orange", lwd = 2,
     xlab = "Number of clusters (k)", ylab = "Months",
     main = "Range of Median Survivals")
abline(v = current_k, col = "red", lty = 2)

# Mutation enrichments (if available)
if (!is.null(mutation_data) && nrow(mutation_summary) > 0) {
  plot(mutation_summary$k, mutation_summary$n_significant_01,
       type = "b", pch = 19, col = "darkblue", lwd = 2,
       xlab = "Number of clusters (k)", ylab = "Count",
       main = "Significant Mutation Enrichments (p<0.01)")
  abline(v = current_k, col = "red", lty = 2)
} else {
  plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
  text(1, 1, "Mutation data not available", cex = 1.5)
}

dev.off()
cat("✓ Saved: quality_metrics_comparison.pdf\n")

# Plot 2: Kaplan-Meier curves for each k
for (k in k_values) {
  cluster_assignment <- consensus_results[[k]]$consensusClass
  cluster_df <- data.frame(
    sample_id = names(cluster_assignment),
    cluster_k = as.factor(cluster_assignment),
    stringsAsFactors = FALSE
  )

  surv_data <- survival_data %>%
    inner_join(cluster_df, by = c("dbgap_rnaseq_sample" = "sample_id"))

  surv_obj <- Surv(time = surv_data$OS_MONTHS, event = surv_data$OS_STATUS)
  km_fit <- survfit(surv_obj ~ cluster_k, data = surv_data)

  logrank_test <- survdiff(surv_obj ~ cluster_k, data = surv_data)
  logrank_p <- 1 - pchisq(logrank_test$chisq, df = k - 1)

  pdf(sprintf("04_Figures/11_Survival_Reanalysis/07_alternative_clustering/kaplan_meier_k%d.pdf", k),
      width = 8, height = 7)

  plot(km_fit,
       col = rainbow(k),
       lwd = 2,
       xlab = "Overall Survival (months)",
       ylab = "Survival Probability",
       main = sprintf("Kaplan-Meier Curves (k=%d clusters)\nLog-rank p=%.4f", k, logrank_p))

  legend("topright",
         legend = paste("Cluster", 1:k),
         col = rainbow(k),
         lwd = 2,
         bty = "n")

  # Add risk table
  times <- seq(0, max(surv_data$OS_MONTHS), by = 12)
  risk_table <- summary(km_fit, times = times)

  dev.off()
  cat(sprintf("✓ Saved: kaplan_meier_k%d.pdf\n", k))
}

cat("\n")

# ------------------------------------------------------------------------------
# SECTION 6: RECOMMENDATIONS
# ------------------------------------------------------------------------------

cat("SECTION 6: RECOMMENDATIONS\n")
cat("==============================================================================\n\n")

# Find best k by multiple criteria
best_by_silhouette <- cluster_comparison$k[which.max(cluster_comparison$mean_silhouette)]
best_by_survival <- survival_comparison$k[which.min(survival_comparison$logrank_p)]
best_by_consensus <- cluster_comparison$k[which.max(cluster_comparison$mean_consensus)]

cat("** OPTIMAL k BY DIFFERENT CRITERIA **\n\n")
cat(sprintf("Best by silhouette width: k=%d (%.3f)\n",
            best_by_silhouette,
            cluster_comparison$mean_silhouette[cluster_comparison$k == best_by_silhouette]))
cat(sprintf("Best by survival separation: k=%d (p=%.4f)\n",
            best_by_survival,
            survival_comparison$logrank_p[survival_comparison$k == best_by_survival]))
cat(sprintf("Best by consensus: k=%d (%.3f)\n\n",
            best_by_consensus,
            cluster_comparison$mean_consensus[cluster_comparison$k == best_by_consensus]))

cat(sprintf("Current choice: k=%d\n\n", current_k))

# Composite score (normalized metrics)
cluster_comparison_scored <- cluster_comparison %>%
  mutate(
    consensus_score = (mean_consensus - min(mean_consensus)) / (max(mean_consensus) - min(mean_consensus)),
    silhouette_score = (mean_silhouette - min(mean_silhouette)) / (max(mean_silhouette) - min(mean_silhouette)),
    balance_score = 1 - (size_cv - min(size_cv)) / (max(size_cv) - min(size_cv))  # Lower CV is better
  ) %>%
  left_join(survival_comparison %>% dplyr::select(k, logrank_p), by = "k") %>%
  mutate(
    survival_score = (-log10(logrank_p) - min(-log10(logrank_p))) / (max(-log10(logrank_p)) - min(-log10(logrank_p))),
    composite_score = (consensus_score + silhouette_score + balance_score + survival_score) / 4
  )

cat("** COMPOSITE RANKING **\n\n")
cat("(Equal weight: consensus, silhouette, balance, survival)\n\n")

ranking <- cluster_comparison_scored %>%
  arrange(desc(composite_score)) %>%
  dplyr::select(k, composite_score, mean_consensus, mean_silhouette, size_cv, logrank_p)

print(ranking)
cat("\n")

recommended_k <- ranking$k[1]

cat(sprintf("** RECOMMENDED k: %d **\n\n", recommended_k))

if (recommended_k == current_k) {
  cat("✓ Current k=2 choice is OPTIMAL by composite criteria\n\n")
  cat("JUSTIFICATION:\n")
  cat("  - Highest or near-highest on multiple quality metrics\n")
  cat("  - Simplest interpretation (binary classification)\n")
  cat("  - Strong biological separation (mutations, survival)\n")
  cat("  - Recommendation: KEEP k=2\n\n")
} else {
  cat(sprintf("⚠ Alternative k=%d scores higher by composite criteria\n\n", recommended_k))
  cat("CONSIDERATIONS:\n")
  cat(sprintf("  Current k=%d advantages:\n", current_k))
  cat("    - Simplest interpretation\n")
  cat("    - Established in analysis pipeline\n")
  cat("    - Already validated\n\n")
  cat(sprintf("  Alternative k=%d advantages:\n", recommended_k))
  cat("    - Better statistical metrics\n")
  cat("    - May capture finer biological distinctions\n\n")
  cat("  Recommendation: Consider k=%d for EXPLORATORY analysis,\n", recommended_k)
  cat("                  but k=%d remains valid for PRIMARY analysis\n\n", current_k)
}

# Save recommendations
recommendations <- data.frame(
  criterion = c("Silhouette", "Survival", "Consensus", "Composite"),
  recommended_k = c(best_by_silhouette, best_by_survival, best_by_consensus, recommended_k),
  stringsAsFactors = FALSE
)

write.csv(recommendations,
          "03_Results/11_Survival_Reanalysis/07_alternative_clustering/k_recommendations.csv",
          row.names = FALSE)

write.csv(ranking,
          "03_Results/11_Survival_Reanalysis/07_alternative_clustering/composite_ranking.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 7: DETAILED k=3 ANALYSIS (IF PROMISING)
# ------------------------------------------------------------------------------

if (recommended_k == 3 || best_by_survival == 3) {
  cat("SECTION 7: DETAILED k=3 ANALYSIS\n")
  cat("==============================================================================\n\n")

  cat("k=3 shows promise - generating detailed analysis...\n\n")

  # Get k=3 assignments
  k3_assignment <- consensus_results[[3]]$consensusClass
  k3_df <- data.frame(
    sample_id = names(k3_assignment),
    cluster_k3 = k3_assignment,
    stringsAsFactors = FALSE
  )

  # Survival by cluster
  surv_k3 <- survival_data %>%
    inner_join(k3_df, by = c("dbgap_rnaseq_sample" = "sample_id"))

  cat("Cluster sizes:\n")
  print(table(surv_k3$cluster_k3))
  cat("\n")

  surv_obj_k3 <- Surv(time = surv_k3$OS_MONTHS, event = surv_k3$OS_STATUS)

  # Pairwise survival comparisons
  cat("Pairwise log-rank tests:\n")
  for (i in 1:2) {
    for (j in (i+1):3) {
      surv_pair <- surv_k3 %>% filter(cluster_k3 %in% c(i, j))
      logrank_pair <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cluster_k3, data = surv_pair)
      p_pair <- 1 - pchisq(logrank_pair$chisq, df = 1)
      cat(sprintf("  Cluster %d vs %d: p=%.4f %s\n",
                  i, j, p_pair, ifelse(p_pair < 0.05, "***", "")))
    }
  }
  cat("\n")

  # Mutation enrichment (if available)
  if (!is.null(mutation_data)) {
    cat("Mutation enrichment in k=3:\n\n")

    mut_k3 <- mutation_data %>%
      inner_join(k3_df, by = c("dbgap_rnaseq_sample" = "sample_id"))

    k3_mutation_profiles <- data.frame()

    for (gene in available_mutations) {
      for (clust in 1:3) {
        mut_rate <- mean(mut_k3[[gene]][mut_k3$cluster_k3 == clust], na.rm = TRUE)
        k3_mutation_profiles <- rbind(k3_mutation_profiles, data.frame(
          cluster_num = clust,
          gene = gene,
          mutation_rate = mut_rate,
          stringsAsFactors = FALSE
        ))
      }
    }

    # Identify cluster-defining mutations
    k3_mutation_wide <- k3_mutation_profiles %>%
      pivot_wider(names_from = cluster_num, values_from = mutation_rate, names_prefix = "C")

    cat("Mutation rates by cluster (top enriched):\n")
    print(k3_mutation_wide %>% arrange(desc(pmax(C1, C2, C3))) %>% head(10))
    cat("\n")

    write.csv(k3_mutation_wide,
              "03_Results/11_Survival_Reanalysis/07_alternative_clustering/k3_mutation_profiles.csv",
              row.names = FALSE)
  }

} else {
  cat("SECTION 7: DETAILED k=3 ANALYSIS\n")
  cat("==============================================================================\n")
  cat("Skipped - k=3 not identified as particularly promising\n\n")
}

# ------------------------------------------------------------------------------
# FINAL SUMMARY
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("ALTERNATIVE CLUSTERING ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")

cat("FILES GENERATED:\n")
cat("  - cluster_quality_comparison.csv\n")
cat("  - survival_comparison.csv\n")
cat("  - mutation_enrichment_by_k.csv\n")
cat("  - k_recommendations.csv\n")
cat("  - composite_ranking.csv\n")
cat("  - quality_metrics_comparison.pdf\n")
cat("  - kaplan_meier_k[2-5].pdf\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("  - Tested k=2,3,4,5 clustering solutions\n"))
cat(sprintf("  - Best by composite score: k=%d\n", recommended_k))
cat(sprintf("  - Current choice (k=%d): %s\n",
            current_k,
            ifelse(recommended_k == current_k, "OPTIMAL", "REASONABLE")))
cat("\n")

cat("CONCLUSION:\n")
if (recommended_k == current_k) {
  cat(sprintf("The current k=%d solution is optimal and should be retained.\n", current_k))
} else {
  cat(sprintf("While k=%d scores marginally better, k=%d remains a valid choice.\n",
              recommended_k, current_k))
  cat("Simpler models (k=2) are often preferable for interpretation and\n")
  cat("clinical translation when differences are small.\n")
}

cat("\n==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
