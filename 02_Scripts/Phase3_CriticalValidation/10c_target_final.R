#!/usr/bin/env Rscript
# TARGET-AML Validation - FINAL with proper mapping

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(randomForest)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("TARGET-AML VALIDATION - FINAL\n")
cat("==============================================================================\n\n")

# Load expression data
cat("Loading expression data...\n")
expr <- readRDS("03_Results/18_TARGET_Validation/target_aml_expression_normalized.rds")
cat(sprintf("Expression: %d genes x %d samples\n\n", nrow(expr), ncol(expr)))

# Load clinical data
cat("Loading clinical data...\n")
clinical <- read.csv("03_Results/18_TARGET_Validation/target_aml_clinical.csv", stringsAsFactors = FALSE)
cat(sprintf("Clinical: %d samples\n\n", nrow(clinical)))

# Load UUID mapping
cat("Loading UUID-barcode mapping...\n")
mapping <- read.csv("03_Results/18_TARGET_Validation/uuid_barcode_mapping.csv", stringsAsFactors = FALSE)
cat(sprintf("Mapping: %d files\n\n", nrow(mapping)))

# Match expression samples (UUIDs) to case barcodes
expr_uuids <- colnames(expr)

# Create UUID to case mapping
uuid_to_case <- setNames(mapping$case_submitter_id, mapping$file_id)

# Map expression columns to case IDs
expr_cases <- uuid_to_case[expr_uuids]

# Find samples with both expression and clinical data
common_cases <- intersect(expr_cases[!is.na(expr_cases)], clinical$sample_barcode)

cat(sprintf("Matched samples: %d\n\n", length(common_cases)))

if (length(common_cases) < 50) {
  cat("⚠ WARNING: Too few matched samples for robust analysis\n")
  quit(save = "no", status = 1)
}

# Filter and match data
expr_idx <- which(expr_cases %in% common_cases)
expr_matched <- expr[, expr_idx]

# Rename columns to case IDs
colnames(expr_matched) <- expr_cases[expr_idx]

# For cases with multiple samples, take the first one
unique_cases <- unique(colnames(expr_matched))
expr_matched <- expr_matched[, match(unique_cases, colnames(expr_matched))]

# Match clinical data
clinical_matched <- clinical %>%
  filter(sample_barcode %in% unique_cases) %>%
  arrange(match(sample_barcode, colnames(expr_matched)))

cat(sprintf("Final matched dataset: %d samples\n", ncol(expr_matched)))
cat(sprintf("Events: %d (%.1f%%)\n", sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))
cat(sprintf("Age range: %.1f - %.1f years\n",
            min(clinical_matched$age_years), max(clinical_matched$age_years)))
cat(sprintf("Median follow-up: %.1f months\n\n", median(clinical_matched$OS_MONTHS)))

# ==============================================================================
# APPLY CLASSIFIER
# ==============================================================================

cat("APPLYING BEATAML CLASSIFIER\n")
cat("==============================================================================\n\n")

beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
signature_genes <- read.csv("03_Results/15_Gene_Signature/50_gene_signature_symbols.csv")
signature_genes <- signature_genes[!is.na(signature_genes$gene), , drop = FALSE]  # Remove NAs
classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")

# Load original Ensembl signature for classifier compatibility
signature_ensembl <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")

# Create mapping from symbols to Ensembl IDs
symbol_to_ensembl <- setNames(signature_ensembl$gene, signature_genes$gene)

# Check gene availability
available_genes <- intersect(signature_genes$gene, rownames(expr_matched))
cat(sprintf("Available genes: %d / 50 (%.1f%%)\n\n",
            length(available_genes),
            length(available_genes) / 50 * 100))

# Prepare prediction data with gene symbols first
pred_data_symbols <- data.frame(t(expr_matched[available_genes, ]))

# Convert column names from symbols to Ensembl IDs for classifier
pred_data <- pred_data_symbols
colnames(pred_data) <- symbol_to_ensembl[colnames(pred_data)]

# Impute missing genes (using Ensembl IDs)
missing_ensembl <- setdiff(signature_ensembl$gene, colnames(pred_data))
if (length(missing_ensembl) > 0) {
  cat(sprintf("Imputing %d missing genes with BeatAML means...\n", length(missing_ensembl)))
  for (ensembl_id in missing_ensembl) {
    if (ensembl_id %in% rownames(beataml_expr)) {
      pred_data[[ensembl_id]] <- mean(beataml_expr[ensembl_id, ])
    }
  }
}

# Reorder columns to match training (Ensembl IDs)
pred_data <- pred_data[, signature_ensembl$gene]

# Predict
cat("Predicting cluster assignments...\n\n")
cluster_pred <- predict(classifier, pred_data)
cluster_prob <- predict(classifier, pred_data, type = "prob")
cluster_conf <- apply(cluster_prob, 1, max)

cat("Cluster distribution:\n")
cat(sprintf("  Cluster 1: %d (%.1f%%)\n",
            sum(cluster_pred == 1), mean(cluster_pred == 1) * 100))
cat(sprintf("  Cluster 2: %d (%.1f%%)\n",
            sum(cluster_pred == 2), mean(cluster_pred == 2) * 100))
cat(sprintf("  Mean confidence: %.3f\n", mean(cluster_conf)))
cat(sprintf("  Low confidence (<0.6): %d (%.1f%%)\n\n",
            sum(cluster_conf < 0.6), mean(cluster_conf < 0.6) * 100))

# Add to clinical
clinical_matched$cluster <- as.numeric(cluster_pred)
clinical_matched$confidence <- cluster_conf

write.csv(clinical_matched,
          "03_Results/18_TARGET_Validation/target_cluster_assignments.csv",
          row.names = FALSE)
cat("✓ Saved cluster assignments\n\n")

# ==============================================================================
# SURVIVAL ANALYSIS
# ==============================================================================

cat("SURVIVAL ANALYSIS\n")
cat("==============================================================================\n\n")

surv_obj <- Surv(time = clinical_matched$OS_MONTHS,
                event = clinical_matched$OS_STATUS)

# Log-rank test
logrank_test <- survdiff(surv_obj ~ cluster, data = clinical_matched)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = 1)

# Cox regression
cox_model <- coxph(surv_obj ~ cluster, data = clinical_matched)
cox_summary <- summary(cox_model)
hr <- exp(coef(cox_model))
hr_ci <- exp(confint(cox_model))
cox_p <- cox_summary$coefficients[, "Pr(>|z|)"]

# Kaplan-Meier
km_fit <- survfit(surv_obj ~ cluster, data = clinical_matched)
km_table <- summary(km_fit)$table
medians <- km_table[, "median"]

cat("** SURVIVAL RESULTS **\n\n")
cat(sprintf("Sample size: %d (Cluster 1: %d, Cluster 2: %d)\n",
            nrow(clinical_matched),
            sum(clinical_matched$cluster == 1),
            sum(clinical_matched$cluster == 2)))
cat(sprintf("Events: %d (%.1f%%)\n\n", sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))

cat(sprintf("Log-rank test: p = %.4f %s\n",
            logrank_p,
            ifelse(logrank_p < 0.05, "*** SIGNIFICANT", "")))
cat(sprintf("Cox HR: %.3f (95%% CI: %.3f - %.3f)\n",
            hr, hr_ci[1], hr_ci[2]))
cat(sprintf("Cox p-value: %.4f\n\n", cox_p))

cat("Median survival by cluster:\n")
cat(sprintf("  Cluster 1: %.1f months\n", medians[1]))
cat(sprintf("  Cluster 2: %.1f months\n", medians[2]))
cat(sprintf("  Difference: %.1f months ", abs(medians[1] - medians[2])))
cat(ifelse(medians[1] > medians[2], "(C1 better)\n\n", "(C2 better)\n\n"))

# Save survival summary
survival_summary <- data.frame(
  cohort = "TARGET-AML",
  age_group = "Pediatric (0-30y)",
  n_samples = nrow(clinical_matched),
  n_events = sum(clinical_matched$OS_STATUS),
  cluster1_n = sum(clinical_matched$cluster == 1),
  cluster2_n = sum(clinical_matched$cluster == 2),
  cluster1_median_os = medians[1],
  cluster2_median_os = medians[2],
  hr = hr,
  hr_lower = hr_ci[1],
  hr_upper = hr_ci[2],
  cox_p = cox_p,
  logrank_p = logrank_p,
  significant = logrank_p < 0.05,
  stringsAsFactors = FALSE
)

write.csv(survival_summary,
          "03_Results/18_TARGET_Validation/target_survival_validation.csv",
          row.names = FALSE)
cat("✓ Saved survival summary\n\n")

# ==============================================================================
# CROSS-COHORT COMPARISON
# ==============================================================================

cat("CROSS-COHORT COMPARISON\n")
cat("==============================================================================\n\n")

# Load adult cohorts
adult_cohorts <- read.csv("03_Results/11_Survival_Reanalysis/07_meta_analysis_cohort_results.csv")

# Add age group if not present
if (!"age_group" %in% names(adult_cohorts)) {
  adult_cohorts$age_group <- "Adult"
}

# Combine all cohorts
all_cohorts <- rbind(
  adult_cohorts[, c("cohort", "age_group", "n_samples", "n_events", "hr", "hr_lower", "hr_upper")],
  survival_summary[, c("cohort", "age_group", "n_samples", "n_events", "hr", "hr_lower", "hr_upper")]
)

cat("All cohorts:\n")
print(all_cohorts)

# Test heterogeneity
weights <- 1 / ((log(all_cohorts$hr_upper) - log(all_cohorts$hr_lower)) / (2 * 1.96))^2
pooled_log_hr <- sum(log(all_cohorts$hr) * weights) / sum(weights)
pooled_hr <- exp(pooled_log_hr)
Q <- sum(weights * (log(all_cohorts$hr) - pooled_log_hr)^2)
Q_df <- nrow(all_cohorts) - 1
Q_p <- 1 - pchisq(Q, df = Q_df)
I2 <- max(0, (Q - Q_df) / Q * 100)

cat("\n\nHeterogeneity test:\n")
cat(sprintf("  Pooled HR: %.3f\n", pooled_hr))
cat(sprintf("  Cochran's Q: %.3f (df=%d), p = %.4f\n", Q, Q_df, Q_p))
cat(sprintf("  I² statistic: %.1f%%\n", I2))

if (Q_p < 0.10) {
  cat("  ⚠ Significant heterogeneity detected\n\n")
} else {
  cat("  ✓ No significant heterogeneity\n\n")
}

# Save combined results
write.csv(all_cohorts,
          "03_Results/18_TARGET_Validation/all_cohorts_comparison.csv",
          row.names = FALSE)

# ==============================================================================
# FIGURES
# ==============================================================================

cat("CREATING FIGURES\n")
cat("==============================================================================\n\n")

dir.create("04_Figures/18_TARGET_Validation", recursive = TRUE, showWarnings = FALSE)

# Kaplan-Meier plot
pdf("04_Figures/18_TARGET_Validation/target_kaplan_meier.pdf", width = 10, height = 8)
par(mar = c(5, 5, 4, 2))

plot(km_fit,
     col = c("blue", "red"),
     lwd = 3,
     lty = 1,
     xlab = "Overall Survival (months)",
     ylab = "Survival Probability",
     main = sprintf("TARGET-AML (Pediatric, Age 0-30)\nn=%d, Log-rank p=%.4f",
                    nrow(clinical_matched), logrank_p),
     cex.lab = 1.3,
     cex.axis = 1.2,
     cex.main = 1.4)

legend("topright",
       legend = c(sprintf("Cluster 1 (n=%d, median=%.1fm)",
                          sum(clinical_matched$cluster == 1), medians[1]),
                  sprintf("Cluster 2 (n=%d, median=%.1fm)",
                          sum(clinical_matched$cluster == 2), medians[2])),
       col = c("blue", "red"),
       lwd = 3,
       cex = 1.2,
       bty = "n")

# Add HR text
hr_text <- sprintf("HR = %.2f (95%% CI: %.2f - %.2f), p = %.4f",
                   hr, hr_ci[1], hr_ci[2], cox_p)
mtext(hr_text, side = 1, line = 4, cex = 1.1)

dev.off()
cat("✓ Saved: target_kaplan_meier.pdf\n")

# Forest plot - All cohorts
pdf("04_Figures/18_TARGET_Validation/forest_plot_all_cohorts.pdf", width = 12, height = 7)
par(mar = c(5, 10, 4, 2))

y_pos <- seq(nrow(all_cohorts), 1, -1)
colors <- ifelse(all_cohorts$age_group == "Pediatric", "darkgreen", "darkblue")

plot(all_cohorts$hr, y_pos,
     xlim = c(0.5, 3.5),
     ylim = c(0.5, nrow(all_cohorts) + 0.5),
     pch = 19,
     cex = 2,
     col = colors,
     xlab = "Hazard Ratio (95% CI)",
     ylab = "",
     yaxt = "n",
     main = "Molecular Subtypes: Prognostic Effect Across Cohorts",
     cex.lab = 1.3,
     cex.main = 1.4)

# Add CI lines
segments(all_cohorts$hr_lower, y_pos, all_cohorts$hr_upper, y_pos,
         lwd = 2, col = colors)

# Add reference line
abline(v = 1, lty = 2, col = "gray40", lwd = 2)
abline(v = pooled_hr, lty = 3, col = "red", lwd = 2)

# Add cohort labels with age group
labels <- paste0(all_cohorts$cohort, " (", all_cohorts$age_group, ")")
axis(2, at = y_pos, labels = labels, las = 1, cex.axis = 1.1)

# Add sample sizes and HR values
text(3.3, y_pos, sprintf("n=%d", all_cohorts$n_samples), cex = 1.0)
text(0.6, y_pos, sprintf("%.2f", all_cohorts$hr), cex = 0.9)

# Add legend
legend("topright",
       legend = c("Adult cohort", "Pediatric cohort",
                  sprintf("Pooled HR: %.2f", pooled_hr),
                  sprintf("I²: %.1f%%, p=%.3f", I2, Q_p)),
       col = c("darkblue", "darkgreen", "red", "white"),
       pch = c(19, 19, NA, NA),
       lty = c(NA, NA, 3, NA),
       lwd = c(NA, NA, 2, NA),
       cex = 1.1,
       bty = "n")

dev.off()
cat("✓ Saved: forest_plot_all_cohorts.pdf\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("==============================================================================\n")
cat("TARGET-AML VALIDATION COMPLETE\n")
cat("==============================================================================\n\n")

cat("COHORT:\n")
cat(sprintf("  Dataset: TARGET-AML (pediatric, age 0-30)\n"))
cat(sprintf("  Samples: %d\n", nrow(clinical_matched)))
cat(sprintf("  Events: %d (%.1f%%)\n", sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))
cat(sprintf("  Median age: %.1f years\n", median(clinical_matched$age_years)))
cat(sprintf("  Median follow-up: %.1f months\n\n", median(clinical_matched$OS_MONTHS)))

cat("VALIDATION RESULTS:\n")
cat(sprintf("  Log-rank p: %.4f %s\n", logrank_p,
            ifelse(logrank_p < 0.05, "(SIGNIFICANT)", "")))
cat(sprintf("  Cox HR: %.2f (95%% CI: %.2f - %.2f)\n", hr, hr_ci[1], hr_ci[2]))
cat(sprintf("  Cox p: %.4f\n\n", cox_p))

cat("INTERPRETATION:\n")

if (survival_summary$significant) {
  if (abs(log(hr) - log(1.35)) < 0.5) {
    cat("✓✓✓ STRONG VALIDATION\n")
    cat("  - Significant survival difference in pediatric AML\n")
    cat("  - Effect size consistent with adult cohorts (HR~1.35)\n")
    cat("  - Molecular subtypes show AGE-INDEPENDENT prognostic value\n")
    cat("  - Validates biological mechanism across age groups\n")
  } else {
    cat("✓ PARTIAL VALIDATION\n")
    cat("  - Significant survival difference detected\n")
    cat("  - Effect size differs from adult cohorts\n")
    cat("  - May indicate age-specific effect modification\n")
  }
} else {
  cat("✗ NO VALIDATION\n")
  cat("  - No significant survival difference in pediatric cohort\n")
  cat("  Possible reasons:\n")
  cat("    1. Age-specific biology (pediatric vs adult AML)\n")
  cat("    2. Treatment differences (intensive pediatric protocols)\n")
  cat("    3. Different mutation landscape in children\n")
  cat("    4. Insufficient statistical power\n")
}

cat("\n==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
