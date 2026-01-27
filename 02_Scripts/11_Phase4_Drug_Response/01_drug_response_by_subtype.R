#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 4: Drug Response Integration
# ==============================================================================
# Objective:
#   1. Prepare drug response data
#   2. Analyze subtype-drug associations
#   3. Identify subtype-specific drug sensitivities
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(pheatmap)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

dir.create("03_Results/09_Drug_Response", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/06_Drug_Response", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 4: DRUG RESPONSE INTEGRATION\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load Drug Response Data
# ------------------------------------------------------------------------------

cat("STEP 1: Loading drug response data...\n\n")

drug_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt",
                  data.table = FALSE)

cat(sprintf("Raw drug data: %d rows\n", nrow(drug_raw)))
cat("Columns:", paste(colnames(drug_raw)[1:10], collapse = ", "), "...\n\n")

# Load gold standard cohort
gold_data <- readRDS("03_Results/05_Analysis_Ready_Data/gold_standard_cohort.rds")

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Load sample mapping
sample_mapping <- read.csv("03_Results/01_Processed_Data/master_sample_id_mapping.csv")

# ------------------------------------------------------------------------------
# STEP 2: Process Drug Response Matrix
# ------------------------------------------------------------------------------

cat("STEP 2: Creating drug response matrix...\n\n")

# Pivot drug data to wide format (samples × drugs)
# Each row in raw data is sample-drug-AUC

drug_matrix <- drug_raw %>%
  dplyr::select(dbgap_rnaseq_sample, inhibitor, auc) %>%
  filter(!is.na(auc)) %>%
  pivot_wider(names_from = inhibitor, values_from = auc,
              values_fn = mean)  # Take mean if duplicate

drug_matrix_df <- as.data.frame(drug_matrix)
rownames(drug_matrix_df) <- drug_matrix_df$dbgap_rnaseq_sample
drug_matrix_df$dbgap_rnaseq_sample <- NULL

cat(sprintf("Drug response matrix: %d samples × %d drugs\n\n",
            nrow(drug_matrix_df), ncol(drug_matrix_df)))

# ------------------------------------------------------------------------------
# STEP 3: Match with Cluster Assignments
# ------------------------------------------------------------------------------

cat("STEP 3: Matching drug data with clusters...\n\n")

# Find samples with both drug data and clusters
common_samples <- intersect(rownames(drug_matrix_df), clusters$sample_id)

drug_clustered <- drug_matrix_df[common_samples, ]
cluster_assignments <- clusters$cluster[match(common_samples, clusters$sample_id)]

cat(sprintf("Samples with both drug and cluster data: %d\n\n", length(common_samples)))

# ------------------------------------------------------------------------------
# STEP 4: Analyze Drug Response by Cluster
# ------------------------------------------------------------------------------

cat("STEP 4: Analyzing drug response by cluster...\n\n")

n_clusters <- length(unique(cluster_assignments))

# For each drug, test if there's a difference across clusters
drug_cluster_results <- data.frame(
  drug = character(),
  kruskal_pvalue = numeric(),
  mean_auc_cluster = character(),
  stringsAsFactors = FALSE
)

drugs_to_test <- colnames(drug_clustered)

for (drug in drugs_to_test) {
  drug_values <- drug_clustered[[drug]]

  # Remove NAs
  valid_idx <- !is.na(drug_values)

  if (sum(valid_idx) >= 30) {  # At least 30 samples
    # Kruskal-Wallis test
    kw_test <- kruskal.test(drug_values[valid_idx] ~ cluster_assignments[valid_idx])

    # Mean AUC per cluster
    mean_aucs <- tapply(drug_values[valid_idx], cluster_assignments[valid_idx], mean, na.rm = TRUE)

    drug_cluster_results <- rbind(drug_cluster_results, data.frame(
      drug = drug,
      n_samples = sum(valid_idx),
      kruskal_pvalue = kw_test$p.value,
      mean_auc_overall = mean(drug_values[valid_idx], na.rm = TRUE),
      sd_auc_overall = sd(drug_values[valid_idx], na.rm = TRUE),
      mean_auc_cluster = paste(round(mean_aucs, 2), collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values
drug_cluster_results$adj_pvalue <- p.adjust(drug_cluster_results$kruskal_pvalue, method = "BH")

# Significant drugs
sig_drugs <- drug_cluster_results %>%
  filter(adj_pvalue < 0.10) %>%
  arrange(adj_pvalue)

cat(sprintf("Drugs with differential response across clusters (FDR<0.10): %d\n\n",
            nrow(sig_drugs)))

if (nrow(sig_drugs) > 0) {
  cat("Top significant drugs:\n")
  print(head(sig_drugs[, c("drug", "n_samples", "kruskal_pvalue", "adj_pvalue")], 10))
  cat("\n")
}

write.csv(drug_cluster_results,
          "03_Results/09_Drug_Response/drug_cluster_associations.csv",
          row.names = FALSE)
cat("✓ Saved: drug_cluster_associations.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 5: Heatmap of Drug Sensitivities
# ------------------------------------------------------------------------------

cat("STEP 5: Creating drug sensitivity heatmap...\n\n")

# Select top variable drugs or significant drugs
if (nrow(sig_drugs) >= 10) {
  drugs_to_plot <- sig_drugs$drug[1:min(30, nrow(sig_drugs))]
} else {
  # Use most variable drugs
  drug_var <- apply(drug_clustered, 2, var, na.rm = TRUE)
  drugs_to_plot <- names(sort(drug_var, decreasing = TRUE))[1:min(30, length(drug_var))]
}

# Create mean AUC matrix (clusters × drugs)
auc_by_cluster <- matrix(NA, nrow = n_clusters, ncol = length(drugs_to_plot))
rownames(auc_by_cluster) <- paste0("Cluster_", 1:n_clusters)
colnames(auc_by_cluster) <- drugs_to_plot

for (k in 1:n_clusters) {
  cluster_samples <- common_samples[cluster_assignments == k]
  for (drug in drugs_to_plot) {
    auc_by_cluster[k, drug] <- mean(drug_clustered[cluster_samples, drug], na.rm = TRUE)
  }
}

# Z-score normalize for visualization
auc_by_cluster_z <- scale(auc_by_cluster)

# Plot heatmap
pdf("04_Figures/06_Drug_Response/drug_sensitivity_heatmap.pdf",
    width = 10, height = 8)

pheatmap(
  t(auc_by_cluster_z),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "Drug Sensitivity by Molecular Subtype\n(Z-scored AUC, lower=more sensitive)",
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 0
)

dev.off()

cat("✓ Saved: drug_sensitivity_heatmap.pdf\n\n")

# ------------------------------------------------------------------------------
# STEP 6: Identify Subtype-Specific Sensitivities
# ------------------------------------------------------------------------------

cat("STEP 6: Identifying subtype-specific drug sensitivities...\n\n")

subtype_drug_recommendations <- list()

for (k in 1:n_clusters) {
  cluster_drugs <- auc_by_cluster[k, ]

  # Lower AUC = more sensitive
  # Find drugs where this cluster is most sensitive
  most_sensitive <- names(sort(cluster_drugs))[1:min(10, length(cluster_drugs))]

  subtype_drug_recommendations[[paste0("Cluster_", k)]] <- list(
    most_sensitive = most_sensitive,
    mean_auc_sensitive = cluster_drugs[most_sensitive]
  )

  cat(sprintf("Cluster %d - Most sensitive drugs (lowest AUC):\n", k))
  n_show <- min(5, length(most_sensitive))
  for (i in 1:n_show) {
    cat(sprintf("  %d. %s (AUC=%.2f)\n",
                i, most_sensitive[i], cluster_drugs[most_sensitive[i]]))
  }
  cat("\n")
}

saveRDS(subtype_drug_recommendations,
        "03_Results/09_Drug_Response/subtype_drug_recommendations.rds")
cat("✓ Saved: subtype_drug_recommendations.rds\n\n")

# ------------------------------------------------------------------------------
# Summary Report
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("DRUG RESPONSE ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")

cat("SUMMARY:\n")
cat(sprintf("  - Samples with drug data: %d\n", nrow(drug_clustered)))
cat(sprintf("  - Drugs analyzed: %d\n", ncol(drug_clustered)))
cat(sprintf("  - Molecular subtypes: %d\n", n_clusters))
cat(sprintf("  - Drugs with differential response (FDR<0.10): %d\n\n", nrow(sig_drugs)))

cat("KEY FINDINGS:\n")
if (nrow(sig_drugs) > 0) {
  cat(sprintf("  - Identified %d drugs with subtype-specific responses\n", nrow(sig_drugs)))
  cat("  - Each subtype has distinct drug sensitivity profile\n")
  cat("  - Can guide precision medicine treatment selection\n\n")
} else {
  cat("  - No drugs reached FDR<0.10 significance\n")
  cat("  - May need larger sample size or different drugs\n\n")
}

cat("OUTPUT FILES:\n")
cat("  - 03_Results/09_Drug_Response/drug_cluster_associations.csv\n")
cat("  - 03_Results/09_Drug_Response/subtype_drug_recommendations.rds\n")
cat("  - 04_Figures/06_Drug_Response/drug_sensitivity_heatmap.pdf\n\n")

cat("NEXT STEPS:\n")
cat("  → Build predictive models for drug response\n")
cat("  → Validate findings in independent cohorts\n")
cat("  → Develop clinical decision support tool\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
