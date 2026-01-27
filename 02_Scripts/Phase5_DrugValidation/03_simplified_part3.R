#!/usr/bin/env Rscript
# ==============================================================================
# Phase 5 Part 3: SIMPLIFIED THREE-WAY INTERACTIONS
# ==============================================================================
# Simplified approach using step-wise R² comparison
# Avoids complex formula syntax issues
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("PART 3: SIMPLIFIED THREE-WAY INTERACTIONS\n")
cat("Question: Do clusters add R² beyond mutations?\n")
cat("==============================================================================\n\n")

# Load data
drug_results <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")
drug_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt", data.table = FALSE)
drug_auc <- drug_raw %>%
  dplyr::select(dbgap_rnaseq_sample, inhibitor, auc) %>%
  filter(!is.na(auc)) %>%
  pivot_wider(names_from = inhibitor, values_from = auc, values_fn = mean)
drug_auc_df <- as.data.frame(drug_auc)
rownames(drug_auc_df) <- drug_auc_df$dbgap_rnaseq_sample
drug_auc_df$dbgap_rnaseq_sample <- NULL
drug_auc <- drug_auc_df

clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

common_samples <- intersect(rownames(drug_auc), clusters$sample_id)
cluster_assign <- clusters$cluster[match(common_samples, clusters$sample_id)]
drug_auc_filt <- drug_auc[common_samples, ]

# Prepare mutation data (select key mutations)
key_mutations <- c("NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2", "TET2", "TP53",
                   "RUNX1", "ASXL1", "NRAS", "KRAS")

# Filter to mutations that actually exist
available_mutations <- intersect(key_mutations, colnames(mutations))
cat(sprintf("Available mutations: %s\n\n", paste(available_mutations, collapse = ", ")))

mutation_data <- mutations[common_samples, available_mutations, drop = FALSE]

# Create analysis dataset
analysis_data <- data.frame(
  sample_id = common_samples,
  cluster = as.numeric(cluster_assign),
  mutation_data,
  check.names = FALSE
)

# Get top 20 drugs
sig_drugs <- drug_results %>% filter(fdr < 0.05) %>% arrange(wilcoxon_pvalue)
top_drugs <- sig_drugs$drug[1:min(20, nrow(sig_drugs))]

cat(sprintf("Testing %d drugs...\n\n", length(top_drugs)))

interaction_results <- data.frame()

for (drug_name in top_drugs) {

  # Get AUC values
  auc_values <- drug_auc_filt[[drug_name]]
  valid_idx <- !is.na(auc_values)

  if (sum(valid_idx) < 50) next

  # Create test dataset
  test_data <- analysis_data[valid_idx, ]
  test_data$auc <- auc_values[valid_idx]

  # Remove mutations with <10 occurrences
  mut_counts <- colSums(test_data[, available_mutations, drop = FALSE], na.rm = TRUE)

  if (drug_name == "Venetoclax") {
    cat(sprintf("    Debug - mut_counts: %s\n", paste(names(mut_counts), "=", mut_counts, collapse=", ")))
  }

  valid_muts <- names(mut_counts)[mut_counts >= 10]

  if (drug_name == "Venetoclax") {
    cat(sprintf("    Debug - valid_muts after filter: %s\n", paste(valid_muts, collapse=", ")))
  }

  if (length(valid_muts) == 0) next

  tryCatch({

    # Model 1: Mutations only
    # Make sure valid_muts are in test_data
    valid_muts_in_data <- intersect(valid_muts, colnames(test_data))
    if (length(valid_muts_in_data) == 0) {
      stop("No valid mutations found in test data")
    }

    X_mutations <- as.matrix(test_data[, valid_muts_in_data, drop = FALSE])
    model1 <- lm(test_data$auc ~ X_mutations)
    r2_mutations <- summary(model1)$r.squared

    # Model 2: Cluster only
    model2 <- lm(test_data$auc ~ test_data$cluster)
    r2_cluster <- summary(model2)$r.squared

    # Model 3: Mutations + Cluster
    X_mutations_cluster <- cbind(X_mutations, cluster = test_data$cluster)
    model3 <- lm(test_data$auc ~ X_mutations_cluster)
    r2_additive <- summary(model3)$r.squared

    # Test if adding cluster improves R²
    anova_result <- anova(model1, model3)
    p_cluster_adds_value <- anova_result$`Pr(>F)`[2]

    # Calculate improvement
    r2_improvement <- r2_additive - r2_mutations
    pct_improvement <- 100 * r2_improvement / r2_mutations

    interaction_results <- rbind(interaction_results, data.frame(
      drug = drug_name,
      n_samples = sum(valid_idx),
      n_mutations = length(valid_muts_in_data),
      mutations_tested = paste(valid_muts_in_data, collapse = ";"),
      r2_mutations_only = r2_mutations,
      r2_cluster_only = r2_cluster,
      r2_mutations_plus_cluster = r2_additive,
      r2_improvement = r2_improvement,
      pct_improvement = pct_improvement,
      p_cluster_adds_value = p_cluster_adds_value,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  ✓ %-30s: R²_mut=%.3f, R²_mut+clus=%.3f, ΔR²=%.4f (p=%.3f)\n",
                substr(drug_name, 1, 30),
                r2_mutations, r2_additive, r2_improvement, p_cluster_adds_value))

  }, error = function(e) {
    cat(sprintf("  ✗ %-30s: %s\n", substr(drug_name, 1, 30), e$message))
  })
}

cat("\n")

if (nrow(interaction_results) > 0) {
  # FDR correction
  interaction_results$fdr_cluster <- p.adjust(interaction_results$p_cluster_adds_value, method = "BH")
  interaction_results <- interaction_results %>% arrange(p_cluster_adds_value)

  # Identify independent drugs
  independent_drugs <- interaction_results %>% filter(p_cluster_adds_value < 0.05)

  cat(sprintf("✓ Successfully tested: %d drugs\n", nrow(interaction_results)))
  cat(sprintf("✓ Cluster adds value (p<0.05): %d drugs\n", nrow(independent_drugs)))
  cat(sprintf("✓ Cluster adds value (FDR<0.05): %d drugs\n\n",
              sum(interaction_results$fdr_cluster < 0.05, na.rm = TRUE)))

  if (nrow(independent_drugs) > 0) {
    cat("⭐⭐⭐ CLUSTERS HAVE INDEPENDENT PREDICTIVE VALUE ⭐⭐⭐\n\n")
    cat("Drugs where cluster adds significant R² beyond mutations:\n\n")

    print_cols <- c("drug", "r2_mutations_only", "r2_mutations_plus_cluster",
                    "r2_improvement", "pct_improvement", "p_cluster_adds_value", "fdr_cluster")
    print(independent_drugs[, print_cols], row.names = FALSE)
    cat("\n")

    cat("Interpretation:\n")
    cat(sprintf("  • %d drugs show independent cluster effect after controlling for mutations\n",
                nrow(independent_drugs)))
    cat("  • Cluster membership provides additional predictive information\n")
    cat("  • These drugs may be particularly useful for cluster-based stratification\n\n")

  } else {
    cat("❌ CLUSTERS DO NOT ADD SIGNIFICANT VALUE BEYOND MUTATIONS\n\n")
    cat("Interpretation:\n")
    cat("  • Cluster effects on drug response are explained by mutation patterns\n")
    cat("  • Clusters are proxies for key mutations (NPM1, TP53, etc.)\n")
    cat("  • No independent predictive value for drug response\n\n")
  }

  # Show summary statistics
  cat("Overall R² Statistics:\n")
  cat(sprintf("  • Mean R² (mutations only):      %.3f (SD=%.3f)\n",
              mean(interaction_results$r2_mutations_only),
              sd(interaction_results$r2_mutations_only)))
  cat(sprintf("  • Mean R² (mutations + cluster): %.3f (SD=%.3f)\n",
              mean(interaction_results$r2_mutations_plus_cluster),
              sd(interaction_results$r2_mutations_plus_cluster)))
  cat(sprintf("  • Mean R² improvement:           %.4f (%.1f%%)\n\n",
              mean(interaction_results$r2_improvement),
              mean(interaction_results$pct_improvement)))

} else {
  cat("ERROR: No drugs could be tested\n\n")
}

write.csv(interaction_results,
          "03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv",
          row.names = FALSE)
cat("✓ Saved: drug_cluster_independence_SIMPLIFIED.csv\n\n")

cat("==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================================\n")
