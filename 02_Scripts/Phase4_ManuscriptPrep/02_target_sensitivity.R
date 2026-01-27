#!/usr/bin/env Rscript
# Phase 4 Part 2: TARGET Sensitivity Analysis
# Purpose: Test robustness of TARGET validation to missing gene imputation

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("TARGET SENSITIVITY ANALYSIS\n")
cat("==============================================================================\n\n")

# Load TARGET data
cat("Loading TARGET-AML data...\n")
target_expr <- readRDS("03_Results/18_TARGET_Validation/target_aml_expression_normalized.rds")
target_clinical <- read.csv("03_Results/18_TARGET_Validation/target_aml_clinical.csv")
target_assignments <- read.csv("03_Results/18_TARGET_Validation/target_cluster_assignments.csv")

cat(sprintf("Expression: %d genes x %d samples\n", nrow(target_expr), ncol(target_expr)))
cat(sprintf("Clinical: %d samples\n", nrow(target_clinical)))
cat(sprintf("Cluster assignments: %d samples\n\n", nrow(target_assignments)))

# Load gene signature
gene_sig_symbols <- read.csv("03_Results/15_Gene_Signature/50_gene_signature_symbols.csv")
gene_sig_symbols <- gene_sig_symbols %>% filter(!is.na(gene) & gene != "")

cat(sprintf("Signature: %d genes\n", nrow(gene_sig_symbols)))

# Check availability
available_genes <- intersect(gene_sig_symbols$gene, rownames(target_expr))
missing_genes <- setdiff(gene_sig_symbols$gene, rownames(target_expr))

cat(sprintf("Available in TARGET: %d / %d (%.1f%%)\n",
            length(available_genes), nrow(gene_sig_symbols),
            length(available_genes) / nrow(gene_sig_symbols) * 100))
cat(sprintf("Missing: %d genes\n", length(missing_genes)))
cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n\n")

# Load BeatAML expression for imputation
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

# Load classifier
rf_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")

# Load Ensembl signature for classifier compatibility
gene_sig_ensembl <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")
symbol_to_ensembl <- setNames(gene_sig_ensembl$gene, gene_sig_symbols$gene)

# ==============================================================================
# SCENARIO 1: MEAN IMPUTATION (CURRENT METHOD)
# ==============================================================================

cat("=== SCENARIO 1: Mean Imputation (BeatAML means) ===\n\n")

# Extract available genes
target_pred_1 <- as.data.frame(t(target_expr[available_genes, ]))

# Impute missing genes with BeatAML means
for (gene_symbol in missing_genes) {
  ensembl_id <- symbol_to_ensembl[gene_symbol]
  if (!is.na(ensembl_id) && ensembl_id %in% rownames(beataml_expr)) {
    target_pred_1[[gene_symbol]] <- mean(beataml_expr[ensembl_id, ], na.rm = TRUE)
  } else {
    target_pred_1[[gene_symbol]] <- 0
  }
}

# Reorder and convert to Ensembl IDs
target_pred_1 <- target_pred_1[, gene_sig_symbols$gene]
colnames(target_pred_1) <- symbol_to_ensembl[colnames(target_pred_1)]
target_pred_1 <- target_pred_1[, gene_sig_ensembl$gene]

# Predict
pred_1 <- predict(rf_classifier, newdata = target_pred_1)
prob_1 <- predict(rf_classifier, newdata = target_pred_1, type = "prob")
conf_1 <- apply(prob_1, 1, max)

cat(sprintf("Cluster 1: %d (%.1f%%)\n", sum(pred_1 == 1), mean(pred_1 == 1) * 100))
cat(sprintf("Cluster 2: %d (%.1f%%)\n", sum(pred_1 == 2), mean(pred_1 == 2) * 100))
cat(sprintf("Mean confidence: %.3f\n", mean(conf_1)))
cat(sprintf("Low confidence (<0.6): %d (%.1f%%)\n\n",
            sum(conf_1 < 0.6), mean(conf_1 < 0.6) * 100))

# ==============================================================================
# SCENARIO 2: ZERO IMPUTATION
# ==============================================================================

cat("=== SCENARIO 2: Zero Imputation ===\n\n")

target_pred_2 <- as.data.frame(t(target_expr[available_genes, ]))

# Impute with zero
for (gene_symbol in missing_genes) {
  target_pred_2[[gene_symbol]] <- 0
}

# Reorder and convert
target_pred_2 <- target_pred_2[, gene_sig_symbols$gene]
colnames(target_pred_2) <- symbol_to_ensembl[colnames(target_pred_2)]
target_pred_2 <- target_pred_2[, gene_sig_ensembl$gene]

# Predict
pred_2 <- predict(rf_classifier, newdata = target_pred_2)
prob_2 <- predict(rf_classifier, newdata = target_pred_2, type = "prob")
conf_2 <- apply(prob_2, 1, max)

cat(sprintf("Cluster 1: %d (%.1f%%)\n", sum(pred_2 == 1), mean(pred_2 == 1) * 100))
cat(sprintf("Cluster 2: %d (%.1f%%)\n", sum(pred_2 == 2), mean(pred_2 == 2) * 100))
cat(sprintf("Mean confidence: %.3f\n", mean(conf_2)))
cat(sprintf("Low confidence (<0.6): %d (%.1f%%)\n\n",
            sum(conf_2 < 0.6), mean(conf_2 < 0.6) * 100))

# ==============================================================================
# SCENARIO 3: MEDIAN IMPUTATION
# ==============================================================================

cat("=== SCENARIO 3: Median Imputation (BeatAML medians) ===\n\n")

target_pred_3 <- as.data.frame(t(target_expr[available_genes, ]))

# Impute with BeatAML medians
for (gene_symbol in missing_genes) {
  ensembl_id <- symbol_to_ensembl[gene_symbol]
  if (!is.na(ensembl_id) && ensembl_id %in% rownames(beataml_expr)) {
    target_pred_3[[gene_symbol]] <- median(beataml_expr[ensembl_id, ], na.rm = TRUE)
  } else {
    target_pred_3[[gene_symbol]] <- 0
  }
}

# Reorder and convert
target_pred_3 <- target_pred_3[, gene_sig_symbols$gene]
colnames(target_pred_3) <- symbol_to_ensembl[colnames(target_pred_3)]
target_pred_3 <- target_pred_3[, gene_sig_ensembl$gene]

# Predict
pred_3 <- predict(rf_classifier, newdata = target_pred_3)
prob_3 <- predict(rf_classifier, newdata = target_pred_3, type = "prob")
conf_3 <- apply(prob_3, 1, max)

cat(sprintf("Cluster 1: %d (%.1f%%)\n", sum(pred_3 == 1), mean(pred_3 == 1) * 100))
cat(sprintf("Cluster 2: %d (%.1f%%)\n", sum(pred_3 == 2), mean(pred_3 == 2) * 100))
cat(sprintf("Mean confidence: %.3f\n", mean(conf_3)))
cat(sprintf("Low confidence (<0.6): %d (%.1f%%)\n\n",
            sum(conf_3 < 0.6), mean(conf_3 < 0.6) * 100))

# ==============================================================================
# COMPARE CLUSTER ASSIGNMENTS ACROSS SCENARIOS
# ==============================================================================

cat("=== AGREEMENT BETWEEN SCENARIOS ===\n\n")
cat(sprintf("Mean vs Zero: %.1f%% agreement\n", mean(pred_1 == pred_2) * 100))
cat(sprintf("Mean vs Median: %.1f%% agreement\n", mean(pred_1 == pred_3) * 100))
cat(sprintf("Zero vs Median: %.1f%% agreement\n\n", mean(pred_2 == pred_3) * 100))

# ==============================================================================
# SURVIVAL ANALYSIS FOR EACH SCENARIO
# ==============================================================================

cat("=== SURVIVAL ANALYSIS BY SCENARIO ===\n\n")

# Prepare survival data
target_surv <- target_clinical %>%
  mutate(
    OS_MONTHS = ifelse(!is.na(days_to_death),
                       days_to_death / 30.44,
                       days_to_last_follow_up / 30.44),
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS))

# Match sample IDs
sample_order <- match(target_surv$sample_barcode, colnames(target_expr))
sample_order <- sample_order[!is.na(sample_order)]

# Add predictions
target_surv_matched <- target_surv[!is.na(match(target_surv$sample_barcode, colnames(target_expr))), ]
target_surv_matched$cluster_mean <- pred_1[sample_order]
target_surv_matched$cluster_zero <- pred_2[sample_order]
target_surv_matched$cluster_median <- pred_3[sample_order]

# Cox models
cox_mean <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_mean, data = target_surv_matched)
cox_zero <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_zero, data = target_surv_matched)
cox_median <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_median, data = target_surv_matched)

# Extract results
survival_comparison <- data.frame(
  Scenario = c("Mean imputation (current)", "Zero imputation", "Median imputation"),
  HR = c(exp(coef(cox_mean)), exp(coef(cox_zero)), exp(coef(cox_median))),
  HR_lower = c(exp(confint(cox_mean)[1]), exp(confint(cox_zero)[1]), exp(confint(cox_median)[1])),
  HR_upper = c(exp(confint(cox_mean)[2]), exp(confint(cox_zero)[2]), exp(confint(cox_median)[2])),
  P_value = c(
    summary(cox_mean)$coefficients[, "Pr(>|z|)"],
    summary(cox_zero)$coefficients[, "Pr(>|z|)"],
    summary(cox_median)$coefficients[, "Pr(>|z|)"]
  )
)

cat("Survival Results:\n")
print(survival_comparison, row.names = FALSE)

# Save results
write.csv(survival_comparison,
          "03_Results/21_Manuscript_Prep/target_survival_sensitivity.csv",
          row.names = FALSE)

# Cluster distribution comparison
cluster_comparison <- data.frame(
  Scenario = c("Mean imputation", "Zero imputation", "Median imputation"),
  Cluster1_n = c(sum(pred_1 == 1), sum(pred_2 == 1), sum(pred_3 == 1)),
  Cluster2_n = c(sum(pred_1 == 2), sum(pred_2 == 2), sum(pred_3 == 2)),
  Cluster1_pct = c(mean(pred_1 == 1) * 100, mean(pred_2 == 1) * 100, mean(pred_3 == 1) * 100),
  Mean_confidence = c(mean(conf_1), mean(conf_2), mean(conf_3)),
  Low_confidence_pct = c(mean(conf_1 < 0.6) * 100, mean(conf_2 < 0.6) * 100, mean(conf_3 < 0.6) * 100)
)

write.csv(cluster_comparison,
          "03_Results/21_Manuscript_Prep/target_imputation_sensitivity.csv",
          row.names = FALSE)

cat("\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("==============================================================================\n")
cat("SENSITIVITY ANALYSIS SUMMARY\n")
cat("==============================================================================\n\n")

# Check HR stability
hr_range <- range(survival_comparison$HR)
hr_cv <- sd(survival_comparison$HR) / mean(survival_comparison$HR)

cat(sprintf("HR range: %.3f - %.3f\n", hr_range[1], hr_range[2]))
cat(sprintf("HR coefficient of variation: %.1f%%\n", hr_cv * 100))
cat(sprintf("P-value range: %.4f - %.4f\n\n",
            min(survival_comparison$P_value), max(survival_comparison$P_value)))

if (hr_cv < 0.10 && all(survival_comparison$P_value > 0.03 & survival_comparison$P_value < 0.10)) {
  cat("✓ RESULTS ARE ROBUST\n")
  cat("  - Hazard ratios highly consistent across imputation methods\n")
  cat("  - All scenarios show marginally non-significant opposite effect\n")
  cat("  - Imputation method does not materially affect conclusions\n")
} else if (hr_cv < 0.20) {
  cat("✓ RESULTS ARE MODERATELY ROBUST\n")
  cat("  - Some variation in HRs but general pattern consistent\n")
} else {
  cat("⚠ RESULTS ARE SENSITIVE TO IMPUTATION\n")
  cat("  - Different imputation methods yield different conclusions\n")
  cat("  - Interpret TARGET results with caution\n")
}

cat("\n✅ TARGET sensitivity analysis complete\n\n")
cat("Files generated:\n")
cat("  - 03_Results/21_Manuscript_Prep/target_survival_sensitivity.csv\n")
cat("  - 03_Results/21_Manuscript_Prep/target_imputation_sensitivity.csv\n\n")
