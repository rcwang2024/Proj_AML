#!/usr/bin/env Rscript
# ==============================================================================
# Phase 3: TARGET-AML Validation (Part 4)
# ==============================================================================
# Objective:
#   1. Download TARGET-AML data from GDC
#   2. Process and normalize expression data
#   3. Apply BeatAML-trained classifier
#   4. Validate survival predictions
#   5. Compare with adult AML cohorts
# Date: 2025-10-11
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(survival)
  library(randomForest)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

# Create output directories
dir.create("03_Results/18_TARGET_Validation", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/18_TARGET_Validation", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 3.4: TARGET-AML VALIDATION\n")
cat("==============================================================================\n\n")

cat("TARGET-AML: Pediatric AML cohort (ages 0-20)\n")
cat("Data source: NCI Genomic Data Commons (GDC)\n\n")

# ------------------------------------------------------------------------------
# SECTION 1: DOWNLOAD TARGET-AML DATA
# ------------------------------------------------------------------------------

cat("SECTION 1: DOWNLOADING TARGET-AML DATA\n")
cat("==============================================================================\n\n")

# Query TARGET-AML RNA-seq data
cat("Querying GDC for TARGET-AML RNA-seq data...\n")

tryCatch({
  query_target <- GDCquery(
    project = "TARGET-AML",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    experimental.strategy = "RNA-Seq"
  )

  cat(sprintf("Found %d samples\n\n", length(query_target$results[[1]]$cases)))

  # Download data (this may take time)
  cat("Downloading TARGET-AML data...\n")
  cat("(This may take 10-30 minutes depending on connection speed)\n\n")

  GDCdownload(query_target, method = "api")

  cat("Preparing expression matrix...\n")
  target_data <- GDCprepare(query_target, summarizedExperiment = TRUE)

  # Extract counts
  target_counts <- assay(target_data, "unstranded")
  target_metadata <- as.data.frame(colData(target_data))

  cat(sprintf("Expression data: %d genes × %d samples\n",
              nrow(target_counts), ncol(target_counts)))

  # Save raw data
  saveRDS(target_counts, "03_Results/18_TARGET_Validation/target_aml_counts_raw.rds")
  saveRDS(target_metadata, "03_Results/18_TARGET_Validation/target_aml_metadata.rds")

  cat("✓ Data downloaded successfully\n\n")

  data_available <- TRUE

}, error = function(e) {
  cat("✗ Error downloading TARGET-AML data:\n")
  cat(sprintf("  %s\n\n", e$message))
  cat("This may be due to:\n")
  cat("  1. Network connectivity issues\n")
  cat("  2. GDC API access restrictions\n")
  cat("  3. Missing TCGAbiolinks dependencies\n\n")
  cat("WORKAROUND: Manual download from GDC Data Portal\n")
  cat("  URL: https://portal.gdc.cancer.gov/projects/TARGET-AML\n\n")

  data_available <<- FALSE
})

if (!data_available) {
  cat("==============================================================================\n")
  cat("TARGET-AML VALIDATION INCOMPLETE\n")
  cat("==============================================================================\n\n")
  cat("Unable to download data automatically.\n")
  cat("Please download manually from GDC and re-run this script.\n\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# SECTION 2: PROCESS AND NORMALIZE DATA
# ------------------------------------------------------------------------------

cat("SECTION 2: PROCESSING TARGET-AML DATA\n")
cat("==============================================================================\n\n")

# Load BeatAML normalization parameters for consistency
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
beataml_genes <- rownames(beataml_expr)

cat("Matching genes with BeatAML...\n")

# Convert Ensembl IDs to gene symbols if needed
if (grepl("^ENSG", rownames(target_counts)[1])) {
  cat("Converting Ensembl IDs to gene symbols...\n")

  # Get gene annotations
  gene_info <- rowData(target_data)
  gene_map <- setNames(gene_info$gene_name, rownames(gene_info))

  # Map to symbols
  rownames(target_counts) <- gene_map[rownames(target_counts)]

  # Remove unmapped genes
  target_counts <- target_counts[!is.na(rownames(target_counts)), ]
  target_counts <- target_counts[rownames(target_counts) != "", ]
}

# Find common genes
common_genes <- intersect(rownames(target_counts), beataml_genes)
cat(sprintf("Common genes: %d / %d BeatAML genes (%.1f%%)\n",
            length(common_genes),
            length(beataml_genes),
            length(common_genes) / length(beataml_genes) * 100))

if (length(common_genes) < 5000) {
  cat("⚠ WARNING: Low gene overlap - may affect classifier performance\n\n")
}

# Subset to common genes
target_counts_matched <- target_counts[common_genes, ]

# Normalize: TPM or log2(CPM+1)
cat("\nNormalizing expression data...\n")

# Calculate library sizes
lib_sizes <- colSums(target_counts_matched)

# CPM normalization
target_cpm <- sweep(target_counts_matched, 2, lib_sizes / 1e6, "/")

# Log transform
target_expr_normalized <- log2(target_cpm + 1)

cat(sprintf("Normalized expression: %d genes × %d samples\n",
            nrow(target_expr_normalized), ncol(target_expr_normalized)))
cat(sprintf("Value range: %.3f to %.3f\n\n",
            min(target_expr_normalized), max(target_expr_normalized)))

# Save normalized data
saveRDS(target_expr_normalized,
        "03_Results/18_TARGET_Validation/target_aml_expression_normalized.rds")

# ------------------------------------------------------------------------------
# SECTION 3: EXTRACT CLINICAL DATA
# ------------------------------------------------------------------------------

cat("SECTION 3: EXTRACTING CLINICAL DATA\n")
cat("==============================================================================\n\n")

# Extract survival and clinical variables
target_clinical <- target_metadata %>%
  dplyr::select(
    sample_id = barcode,
    age = age_at_diagnosis,
    sex = gender,
    vital_status,
    days_to_death,
    days_to_last_follow_up
  ) %>%
  mutate(
    # Calculate OS
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0),
    OS_DAYS = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    OS_MONTHS = OS_DAYS / 30.44,

    # Standardize sex
    SEX = ifelse(sex == "male", "M", "F"),

    # Age in years
    AGE_YEARS = age / 365.25
  ) %>%
  filter(!is.na(OS_DAYS) & OS_DAYS > 0)

cat(sprintf("Clinical data: %d samples with survival information\n", nrow(target_clinical)))
cat(sprintf("Events: %d deaths (%.1f%%)\n",
            sum(target_clinical$OS_STATUS),
            mean(target_clinical$OS_STATUS) * 100))
cat(sprintf("Median follow-up: %.1f months\n", median(target_clinical$OS_MONTHS)))
cat(sprintf("Age range: %.1f - %.1f years (pediatric)\n",
            min(target_clinical$AGE_YEARS), max(target_clinical$AGE_YEARS)))

cat("\nAge distribution:\n")
print(summary(target_clinical$AGE_YEARS))
cat("\n")

# Save clinical data
write.csv(target_clinical,
          "03_Results/18_TARGET_Validation/target_aml_clinical.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 4: APPLY BEATAML CLASSIFIER
# ------------------------------------------------------------------------------

cat("SECTION 4: APPLYING BEATAML-TRAINED CLASSIFIER\n")
cat("==============================================================================\n\n")

# Load BeatAML 50-gene signature and classifier
signature_genes <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")
beataml_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")

cat(sprintf("Loading 50-gene signature (%d genes)\n", nrow(signature_genes)))

# Check gene availability
available_signature_genes <- intersect(signature_genes$gene, rownames(target_expr_normalized))
cat(sprintf("Signature genes available in TARGET: %d / 50 (%.1f%%)\n",
            length(available_signature_genes),
            length(available_signature_genes) / 50 * 100))

if (length(available_signature_genes) < 30) {
  cat("⚠ WARNING: <60% signature genes available - predictions may be unreliable\n\n")
}

# Match samples between expression and clinical data
common_samples <- intersect(colnames(target_expr_normalized), target_clinical$sample_id)
target_expr_matched <- target_expr_normalized[, common_samples]
target_clinical_matched <- target_clinical %>%
  filter(sample_id %in% common_samples)

cat(sprintf("Matched samples: %d\n\n", length(common_samples)))

# Prepare prediction data
target_pred_data <- data.frame(t(target_expr_matched[available_signature_genes, ]))

# For missing genes, impute with BeatAML training mean
missing_genes <- setdiff(signature_genes$gene, available_signature_genes)
if (length(missing_genes) > 0) {
  cat(sprintf("Imputing %d missing genes with BeatAML means...\n", length(missing_genes)))

  beataml_expr_signature <- beataml_expr[intersect(signature_genes$gene, rownames(beataml_expr)), ]

  for (gene in missing_genes) {
    if (gene %in% rownames(beataml_expr_signature)) {
      target_pred_data[[gene]] <- mean(beataml_expr_signature[gene, ])
    }
  }
}

# Reorder columns to match training data
target_pred_data <- target_pred_data[, signature_genes$gene]

# Predict cluster assignments
cat("Predicting cluster assignments...\n")

target_predictions <- predict(beataml_classifier, target_pred_data, type = "prob")
target_cluster <- predict(beataml_classifier, target_pred_data)

# Calculate prediction confidence
target_confidence <- apply(target_predictions, 1, max)

cat("\nPrediction summary:\n")
cat(sprintf("  Cluster 1: %d samples (%.1f%%)\n",
            sum(target_cluster == 1),
            mean(target_cluster == 1) * 100))
cat(sprintf("  Cluster 2: %d samples (%.1f%%)\n",
            sum(target_cluster == 2),
            mean(target_cluster == 2) * 100))
cat(sprintf("  Mean confidence: %.3f\n", mean(target_confidence)))
cat(sprintf("  Low confidence (<0.6): %d samples (%.1f%%)\n",
            sum(target_confidence < 0.6),
            mean(target_confidence < 0.6) * 100))
cat("\n")

# Create results dataframe
target_results <- target_clinical_matched %>%
  mutate(
    cluster = as.numeric(target_cluster),
    confidence = target_confidence
  )

# Save cluster assignments
write.csv(target_results,
          "03_Results/18_TARGET_Validation/target_cluster_assignments.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 5: SURVIVAL VALIDATION
# ------------------------------------------------------------------------------

cat("SECTION 5: SURVIVAL VALIDATION IN TARGET-AML\n")
cat("==============================================================================\n\n")

# Kaplan-Meier analysis
surv_obj <- Surv(time = target_results$OS_MONTHS, event = target_results$OS_STATUS)

# Log-rank test
logrank_test <- survdiff(surv_obj ~ cluster, data = target_results)
logrank_p <- 1 - pchisq(logrank_test$chisq, df = 1)

# Cox regression
cox_model <- coxph(surv_obj ~ cluster, data = target_results)
cox_summary <- summary(cox_model)

hr <- exp(coef(cox_model))
hr_ci <- exp(confint(cox_model))

# Median survival
km_fit <- survfit(surv_obj ~ cluster, data = target_results)
medians <- summary(km_fit)$table[, "median"]

cat("** SURVIVAL ANALYSIS RESULTS **\n\n")
cat(sprintf("Log-rank test: p = %.4f %s\n",
            logrank_p,
            ifelse(logrank_p < 0.05, "***", "")))
cat(sprintf("Cox HR: %.3f (95%% CI: %.3f - %.3f), p = %.4f\n",
            hr, hr_ci[1], hr_ci[2], cox_summary$coefficients[, "Pr(>|z|)"]))
cat("\n")

cat("Median survival by cluster:\n")
cat(sprintf("  Cluster 1: %.1f months (n=%d)\n", medians[1], sum(target_results$cluster == 1)))
cat(sprintf("  Cluster 2: %.1f months (n=%d)\n", medians[2], sum(target_results$cluster == 2)))
cat(sprintf("  Difference: %.1f months\n\n", abs(medians[1] - medians[2])))

# Save survival results
survival_summary <- data.frame(
  cohort = "TARGET-AML",
  n_samples = nrow(target_results),
  n_events = sum(target_results$OS_STATUS),
  cluster1_n = sum(target_results$cluster == 1),
  cluster2_n = sum(target_results$cluster == 2),
  cluster1_median_os = medians[1],
  cluster2_median_os = medians[2],
  hr = hr,
  hr_lower = hr_ci[1],
  hr_upper = hr_ci[2],
  cox_p = cox_summary$coefficients[, "Pr(>|z|)"],
  logrank_p = logrank_p,
  significant = logrank_p < 0.05
)

write.csv(survival_summary,
          "03_Results/18_TARGET_Validation/target_survival_validation.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------
# SECTION 6: COMPARE WITH ADULT AML COHORTS
# ------------------------------------------------------------------------------

cat("SECTION 6: COMPARISON WITH ADULT AML COHORTS\n")
cat("==============================================================================\n\n")

# Load BeatAML and TCGA results
beataml_survival <- read.csv("03_Results/11_Survival_Reanalysis/07_meta_analysis_cohort_results.csv")

# Combine results
all_cohorts <- rbind(
  beataml_survival,
  survival_summary %>%
    dplyr::select(cohort, n_samples, n_events, hr, hr_lower, hr_upper,
                  cox_p, logrank_p, significant)
)

cat("Cross-cohort comparison:\n")
print(all_cohorts %>%
        dplyr::select(cohort, n_samples, n_events, hr, logrank_p, significant))
cat("\n")

# Test for heterogeneity
cat("Heterogeneity testing:\n\n")

# Q statistic
weights <- 1 / ((log(all_cohorts$hr_upper) - log(all_cohorts$hr_lower)) / (2 * 1.96))^2
pooled_log_hr <- sum(log(all_cohorts$hr) * weights) / sum(weights)
Q <- sum(weights * (log(all_cohorts$hr) - pooled_log_hr)^2)
Q_df <- nrow(all_cohorts) - 1
Q_p <- 1 - pchisq(Q, df = Q_df)

cat(sprintf("Cochran's Q: %.3f (df=%d), p = %.4f\n", Q, Q_df, Q_p))

# I² statistic
I2 <- max(0, (Q - Q_df) / Q * 100)
cat(sprintf("I² statistic: %.1f%%\n", I2))

if (Q_p < 0.10) {
  cat("⚠ Significant heterogeneity detected\n")
} else {
  cat("✓ No significant heterogeneity\n")
}
cat("\n")

# Age effect
cat("** AGE COMPARISON **\n\n")
cat("TARGET-AML (pediatric):\n")
cat(sprintf("  Age range: %.1f - %.1f years\n",
            min(target_results$AGE_YEARS), max(target_results$AGE_YEARS)))
cat(sprintf("  Median age: %.1f years\n", median(target_results$AGE_YEARS)))
cat("\n")

cat("Interpretation:\n")
if (survival_summary$significant) {
  cat("✓ Molecular subtypes show prognostic value in PEDIATRIC AML\n")
  cat("  → Effect extends beyond adult populations\n")
  cat("  → Suggests age-independent biological mechanism\n")
} else {
  cat("✗ No significant survival difference in pediatric cohort\n")
  cat("  → May indicate age-specific biology\n")
  cat("  → Or insufficient power (check sample size)\n")
}
cat("\n")

# ------------------------------------------------------------------------------
# SECTION 7: VISUALIZATIONS
# ------------------------------------------------------------------------------

cat("SECTION 7: CREATING VISUALIZATIONS\n")
cat("==============================================================================\n\n")

# Kaplan-Meier plot
pdf("04_Figures/18_TARGET_Validation/target_kaplan_meier.pdf",
    width = 8, height = 7)

plot(km_fit,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Overall Survival (months)",
     ylab = "Survival Probability",
     main = sprintf("TARGET-AML Kaplan-Meier Curves\n(Pediatric AML, n=%d, p=%.4f)",
                    nrow(target_results), logrank_p))

legend("topright",
       legend = c(sprintf("Cluster 1 (n=%d)", sum(target_results$cluster == 1)),
                  sprintf("Cluster 2 (n=%d)", sum(target_results$cluster == 2))),
       col = c("blue", "red"),
       lwd = 2,
       bty = "n")

# Add risk table
grid.text <- sprintf("HR=%.2f (95%% CI: %.2f-%.2f), p=%.4f",
                     hr, hr_ci[1], hr_ci[2], cox_summary$coefficients[, "Pr(>|z|)"])
mtext(grid.text, side = 1, line = 4, cex = 0.9)

dev.off()
cat("✓ Saved: target_kaplan_meier.pdf\n")

# Forest plot comparing all cohorts
pdf("04_Figures/18_TARGET_Validation/forest_plot_all_cohorts.pdf",
    width = 10, height = 6)

par(mar = c(5, 8, 4, 2))

# Plot
y_pos <- seq(nrow(all_cohorts), 1, -1)

plot(all_cohorts$hr, y_pos,
     xlim = c(0.5, 3),
     ylim = c(0.5, nrow(all_cohorts) + 0.5),
     pch = 15,
     cex = 1.5,
     xlab = "Hazard Ratio (95% CI)",
     ylab = "",
     yaxt = "n",
     main = "Prognostic Effect Across Cohorts")

# Add CI lines
segments(all_cohorts$hr_lower, y_pos, all_cohorts$hr_upper, y_pos, lwd = 2)

# Add reference line
abline(v = 1, lty = 2, col = "gray")

# Add cohort labels
axis(2, at = y_pos, labels = all_cohorts$cohort, las = 1)

# Add sample sizes
text(2.8, y_pos, sprintf("n=%d", all_cohorts$n_samples), cex = 0.8)

dev.off()
cat("✓ Saved: forest_plot_all_cohorts.pdf\n\n")

# ------------------------------------------------------------------------------
# SECTION 8: SUMMARY
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("TARGET-AML VALIDATION COMPLETE\n")
cat("==============================================================================\n\n")

cat("COHORT CHARACTERISTICS:\n")
cat(sprintf("  Samples: %d pediatric AML patients\n", nrow(target_results)))
cat(sprintf("  Age: %.1f ± %.1f years (range %.1f-%.1f)\n",
            mean(target_results$AGE_YEARS), sd(target_results$AGE_YEARS),
            min(target_results$AGE_YEARS), max(target_results$AGE_YEARS)))
cat(sprintf("  Events: %d (%.1f%%)\n",
            sum(target_results$OS_STATUS),
            mean(target_results$OS_STATUS) * 100))
cat("\n")

cat("VALIDATION RESULTS:\n")
cat(sprintf("  Log-rank p-value: %.4f %s\n",
            logrank_p,
            ifelse(logrank_p < 0.05, "(SIGNIFICANT)", "(not significant)")))
cat(sprintf("  Hazard ratio: %.2f (95%% CI: %.2f-%.2f)\n",
            hr, hr_ci[1], hr_ci[2]))
cat(sprintf("  Effect direction: %s\n",
            ifelse(hr > 1, "Cluster 2 worse (consistent with adults)",
                   "Cluster 1 worse (opposite to adults)")))
cat("\n")

cat("INTERPRETATION:\n")
if (survival_summary$significant && abs(log(hr) - log(1.35)) < 0.5) {
  cat("✓✓✓ STRONG VALIDATION\n")
  cat("  - Significant survival difference in pediatric AML\n")
  cat("  - Effect size consistent with adult cohorts\n")
  cat("  - Subtypes show age-independent prognostic value\n")
} else if (survival_summary$significant) {
  cat("✓ PARTIAL VALIDATION\n")
  cat("  - Significant survival difference detected\n")
  cat("  - Effect size differs from adult cohorts\n")
  cat("  - May indicate age-specific biology\n")
} else {
  cat("✗ NO VALIDATION\n")
  cat("  - No significant survival difference\n")
  cat("  Possible reasons:\n")
  cat("    1. Insufficient statistical power\n")
  cat("    2. Age-specific biology (pediatric vs adult)\n")
  cat("    3. Treatment differences\n")
  cat("    4. Classifier not transferable to pediatric AML\n")
}
cat("\n")

cat("FILES GENERATED:\n")
cat("  - target_aml_expression_normalized.rds\n")
cat("  - target_aml_clinical.csv\n")
cat("  - target_cluster_assignments.csv\n")
cat("  - target_survival_validation.csv\n")
cat("  - target_kaplan_meier.pdf\n")
cat("  - forest_plot_all_cohorts.pdf\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
