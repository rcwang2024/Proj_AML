#!/usr/bin/env Rscript
# TARGET-AML Validation - Continuation from saved expression data
# This script continues from the successfully saved expression matrix

suppressPackageStartupMessages({
  library(tidyverse)
  library(TCGAbiolinks)
  library(survival)
  library(randomForest)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("TARGET-AML VALIDATION - CONTINUATION\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# LOAD SAVED EXPRESSION DATA
# ------------------------------------------------------------------------------

cat("Loading saved expression data...\n")
expr_log <- readRDS("03_Results/18_TARGET_Validation/target_aml_expression_normalized.rds")
cat(sprintf("Expression: %d genes x %d samples\n\n", nrow(expr_log), ncol(expr_log)))

# ------------------------------------------------------------------------------
# PROCESS CLINICAL DATA
# ------------------------------------------------------------------------------

cat("SECTION 1: PROCESSING CLINICAL DATA\n")
cat("==============================================================================\n\n")

# Get clinical from GDC
target_clinical_gdc <- GDCquery_clinic(project = "TARGET-AML", type = "clinical")

# Fix duplicate names
colnames(target_clinical_gdc) <- make.unique(colnames(target_clinical_gdc))

cat(sprintf("Clinical records: %d\n", nrow(target_clinical_gdc)))

# Process clinical data carefully
target_clinical <- target_clinical_gdc %>%
  as_tibble() %>%
  mutate(
    sample_barcode = submitter_id,
    age_days = suppressWarnings(as.numeric(gsub(" days", "", age_at_diagnosis))),
    age_years = age_days / 365.25,
    sex = as.character(gender),
    vital_status = as.character(vital_status),
    days_death = suppressWarnings(as.numeric(gsub(" days", "", days_to_death))),
    days_followup = suppressWarnings(as.numeric(gsub(" days", "", days_to_last_follow_up)))
  ) %>%
  mutate(
    OS_DAYS = ifelse(!is.na(days_death), days_death, days_followup),
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0),
    OS_MONTHS = OS_DAYS / 30.44
  ) %>%
  filter(!is.na(OS_DAYS) & OS_DAYS > 0 & !is.na(age_years)) %>%
  dplyr::select(sample_barcode, age_years, sex, vital_status, OS_DAYS, OS_STATUS, OS_MONTHS)

cat(sprintf("Valid clinical data: %d samples\n", nrow(target_clinical)))
cat(sprintf("Events: %d deaths (%.1f%%)\n",
            sum(target_clinical$OS_STATUS),
            mean(target_clinical$OS_STATUS) * 100))
cat(sprintf("Age range: %.1f - %.1f years\n",
            min(target_clinical$age_years), max(target_clinical$age_years)))
cat(sprintf("Median follow-up: %.1f months\n\n", median(target_clinical$OS_MONTHS)))

# Save
write.csv(target_clinical,
          "03_Results/18_TARGET_Validation/target_aml_clinical.csv",
          row.names = FALSE)
cat("✓ Saved clinical data\n\n")

# ------------------------------------------------------------------------------
# MATCH EXPRESSION TO CLINICAL
# ------------------------------------------------------------------------------

cat("SECTION 2: MATCHING EXPRESSION TO CLINICAL DATA\n")
cat("==============================================================================\n\n")

# Get query for UUID-barcode mapping
query_target <- GDCquery(
  project = "TARGET-AML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Create UUID to barcode mapping
query_results <- query_target$results[[1]]
uuid_to_barcode <- setNames(query_results$cases, query_results$id)

# Map expression samples to barcodes
expr_sample_uuids <- colnames(expr_log)
expr_sample_barcodes <- uuid_to_barcode[expr_sample_uuids]

# Find samples with both expression and clinical
common_barcodes <- intersect(expr_sample_barcodes[!is.na(expr_sample_barcodes)],
                             target_clinical$sample_barcode)

cat(sprintf("Samples with both expression and clinical: %d\n", length(common_barcodes)))

if (length(common_barcodes) < 50) {
  cat("\n⚠ WARNING: Very few matched samples - barcode mapping may have issues\n\n")
}

# Filter data
expr_idx <- which(expr_sample_barcodes %in% common_barcodes)
expr_matched <- expr_log[, expr_idx]
colnames(expr_matched) <- expr_sample_barcodes[expr_idx]

clinical_matched <- target_clinical %>%
  filter(sample_barcode %in% common_barcodes) %>%
  arrange(match(sample_barcode, colnames(expr_matched)))

cat(sprintf("\nFinal matched dataset: %d samples\n", ncol(expr_matched)))
cat(sprintf("Events: %d (%.1f%%)\n\n", sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))

# ------------------------------------------------------------------------------
# APPLY CLASSIFIER
# ------------------------------------------------------------------------------

cat("SECTION 3: APPLYING BEATAML CLASSIFIER\n")
cat("==============================================================================\n\n")

# Load BeatAML reference
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
signature_genes <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")
classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")

cat(sprintf("Loaded 50-gene signature\n"))

# Check availability
available_genes <- intersect(signature_genes$gene, rownames(expr_matched))
cat(sprintf("Available genes: %d / 50 (%.1f%%)\n",
            length(available_genes),
            length(available_genes) / 50 * 100))

# Prepare prediction data
pred_data <- data.frame(t(expr_matched[available_genes, ]))

# Impute missing genes
missing_genes <- setdiff(signature_genes$gene, available_genes)
if (length(missing_genes) > 0) {
  cat(sprintf("Imputing %d missing genes...\n", length(missing_genes)))
  for (gene in missing_genes) {
    if (gene %in% rownames(beataml_expr)) {
      pred_data[[gene]] <- mean(beataml_expr[gene, ])
    }
  }
}

# Reorder
pred_data <- pred_data[, signature_genes$gene]

# Predict
cat("Predicting clusters...\n")
cluster_pred <- predict(classifier, pred_data)
cluster_prob <- predict(classifier, pred_data, type = "prob")
cluster_conf <- apply(cluster_prob, 1, max)

cat("\nCluster distribution:\n")
cat(sprintf("  Cluster 1: %d (%.1f%%)\n",
            sum(cluster_pred == 1), mean(cluster_pred == 1) * 100))
cat(sprintf("  Cluster 2: %d (%.1f%%)\n",
            sum(cluster_pred == 2), mean(cluster_pred == 2) * 100))
cat(sprintf("  Mean confidence: %.3f\n\n", mean(cluster_conf)))

# Add to clinical
clinical_matched$cluster <- as.numeric(cluster_pred)
clinical_matched$confidence <- cluster_conf

write.csv(clinical_matched,
          "03_Results/18_TARGET_Validation/target_cluster_assignments.csv",
          row.names = FALSE)
cat("✓ Saved cluster assignments\n\n")

# ------------------------------------------------------------------------------
# SURVIVAL ANALYSIS
# ------------------------------------------------------------------------------

cat("SECTION 4: SURVIVAL ANALYSIS\n")
cat("==============================================================================\n\n")

surv_obj <- Surv(time = clinical_matched$OS_MONTHS,
                event = clinical_matched$OS_STATUS)

# Log-rank
logrank_test <- survdiff(surv_obj ~ cluster, data = clinical_matched)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = 1)

# Cox
cox_model <- coxph(surv_obj ~ cluster, data = clinical_matched)
cox_summary <- summary(cox_model)
hr <- exp(coef(cox_model))
hr_ci <- exp(confint(cox_model))
cox_p <- cox_summary$coefficients[, "Pr(>|z|)"]

# KM
km_fit <- survfit(surv_obj ~ cluster, data = clinical_matched)
medians <- summary(km_fit)$table[, "median"]

cat("** SURVIVAL RESULTS **\n\n")
cat(sprintf("Log-rank: p = %.4f %s\n",
            logrank_p,
            ifelse(logrank_p < 0.05, "*** SIGNIFICANT", "")))
cat(sprintf("Cox HR: %.3f (95%% CI: %.3f - %.3f), p = %.4f\n\n",
            hr, hr_ci[1], hr_ci[2], cox_p))

cat("Median survival:\n")
cat(sprintf("  Cluster 1: %.1f months (n=%d)\n",
            medians[1], sum(clinical_matched$cluster == 1)))
cat(sprintf("  Cluster 2: %.1f months (n=%d)\n",
            medians[2], sum(clinical_matched$cluster == 2)))
cat(sprintf("  Difference: %.1f months\n\n", abs(medians[1] - medians[2])))

# Save summary
survival_summary <- data.frame(
  cohort = "TARGET-AML",
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
  significant = logrank_p < 0.05
)

write.csv(survival_summary,
          "03_Results/18_TARGET_Validation/target_survival_validation.csv",
          row.names = FALSE)
cat("✓ Saved survival results\n\n")

# ------------------------------------------------------------------------------
# CROSS-COHORT COMPARISON
# ------------------------------------------------------------------------------

cat("SECTION 5: CROSS-COHORT COMPARISON\n")
cat("==============================================================================\n\n")

# Load adult cohorts
beataml_tcga <- read.csv("03_Results/11_Survival_Reanalysis/07_meta_analysis_cohort_results.csv")

# Combine
all_cohorts <- rbind(
  beataml_tcga %>% dplyr::select(cohort, n_samples, n_events, hr, hr_lower, hr_upper),
  survival_summary %>% dplyr::select(cohort, n_samples, n_events, hr, hr_lower, hr_upper)
)

print(all_cohorts)

# Heterogeneity
weights <- 1 / ((log(all_cohorts$hr_upper) - log(all_cohorts$hr_lower)) / (2 * 1.96))^2
pooled_log_hr <- sum(log(all_cohorts$hr) * weights) / sum(weights)
Q <- sum(weights * (log(all_cohorts$hr) - pooled_log_hr)^2)
Q_df <- nrow(all_cohorts) - 1
Q_p <- 1 - pchisq(Q, df = Q_df)
I2 <- max(0, (Q - Q_df) / Q * 100)

cat("\n\nHeterogeneity:\n")
cat(sprintf("  Q: %.3f (p = %.4f)\n", Q, Q_p))
cat(sprintf("  I²: %.1f%%\n\n", I2))

# ------------------------------------------------------------------------------
# FIGURES
# ------------------------------------------------------------------------------

cat("SECTION 6: CREATING FIGURES\n")
cat("==============================================================================\n\n")

dir.create("04_Figures/18_TARGET_Validation", recursive = TRUE, showWarnings = FALSE)

# KM plot
pdf("04_Figures/18_TARGET_Validation/target_kaplan_meier.pdf", width = 8, height = 7)
plot(km_fit, col = c("blue", "red"), lwd = 2,
     xlab = "Overall Survival (months)", ylab = "Survival Probability",
     main = sprintf("TARGET-AML (Pediatric)\nn=%d, p=%.4f",
                    nrow(clinical_matched), logrank_p))
legend("topright",
       legend = c(sprintf("Cluster 1 (n=%d)", sum(clinical_matched$cluster == 1)),
                  sprintf("Cluster 2 (n=%d)", sum(clinical_matched$cluster == 2))),
       col = c("blue", "red"), lwd = 2)
dev.off()
cat("✓ Saved KM plot\n")

# Forest plot
pdf("04_Figures/18_TARGET_Validation/forest_plot_all_cohorts.pdf", width = 10, height = 6)
par(mar = c(5, 8, 4, 2))
y_pos <- seq(nrow(all_cohorts), 1, -1)
plot(all_cohorts$hr, y_pos, xlim = c(0.5, 3), ylim = c(0.5, nrow(all_cohorts) + 0.5),
     pch = 15, cex = 1.5, xlab = "Hazard Ratio (95% CI)", ylab = "", yaxt = "n",
     main = "Prognostic Effect Across Cohorts (Adult + Pediatric)")
segments(all_cohorts$hr_lower, y_pos, all_cohorts$hr_upper, y_pos, lwd = 2)
abline(v = 1, lty = 2, col = "gray")
axis(2, at = y_pos, labels = all_cohorts$cohort, las = 1)
text(2.8, y_pos, sprintf("n=%d", all_cohorts$n_samples), cex = 0.8)
dev.off()
cat("✓ Saved forest plot\n\n")

# ------------------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("TARGET-AML VALIDATION COMPLETE\n")
cat("==============================================================================\n\n")

cat(sprintf("Samples: %d pediatric AML\n", nrow(clinical_matched)))
cat(sprintf("Events: %d (%.1f%%)\n", sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))
cat(sprintf("Log-rank p: %.4f %s\n", logrank_p,
            ifelse(logrank_p < 0.05, "(SIGNIFICANT)", "")))
cat(sprintf("HR: %.2f (95%% CI: %.2f-%.2f)\n\n", hr, hr_ci[1], hr_ci[2]))

if (survival_summary$significant && abs(log(hr) - log(1.35)) < 0.5) {
  cat("✓✓✓ STRONG VALIDATION\n")
  cat("  Molecular subtypes validated in pediatric AML\n")
  cat("  Effect size consistent with adults\n")
  cat("  Age-independent biology confirmed\n")
} else if (survival_summary$significant) {
  cat("✓ PARTIAL VALIDATION\n")
  cat("  Significant but different effect size\n")
} else {
  cat("✗ NO VALIDATION\n")
  cat("  No significant survival difference\n")
}

cat("\n==============================================================================\n")
cat(sprintf("Completed: %s\n", Sys.time()))
cat("==============================================================================\n")
