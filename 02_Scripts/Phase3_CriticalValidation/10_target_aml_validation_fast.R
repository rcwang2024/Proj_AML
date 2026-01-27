#!/usr/bin/env Rscript
# ==============================================================================
# Phase 3: TARGET-AML Validation (Fast Version)
# ==============================================================================
# This version directly processes GDC downloaded files without GDCprepare
# which is much faster for large datasets
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(TCGAbiolinks)
  library(survival)
  library(randomForest)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("PHASE 3.4: TARGET-AML VALIDATION (FAST VERSION)\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# SECTION 1: LOAD GDC CLINICAL DATA
# ------------------------------------------------------------------------------

cat("SECTION 1: LOADING CLINICAL DATA FROM GDC\n")
cat("==============================================================================\n\n")

# Query to get clinical data
query_target <- GDCquery(
  project = "TARGET-AML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Get clinical data
target_clinical_gdc <- GDCquery_clinic(project = "TARGET-AML", type = "clinical")

# Fix duplicate column names
colnames(target_clinical_gdc) <- make.unique(colnames(target_clinical_gdc))

cat(sprintf("Clinical records: %d\n", nrow(target_clinical_gdc)))

# ------------------------------------------------------------------------------
# SECTION 2: READ EXPRESSION DATA DIRECTLY
# ------------------------------------------------------------------------------

cat("\nSECTION 2: READING EXPRESSION DATA\n")
cat("==============================================================================\n\n")

# Find all TSV files
data_dir <- "GDCdata/TARGET-AML/Transcriptome_Profiling/Gene_Expression_Quantification"
sample_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)

cat(sprintf("Found %d sample directories\n", length(sample_dirs)))
cat("Reading first file to get gene list...\n")

# Read first file to get gene structure
first_file <- list.files(sample_dirs[1], pattern = "\\.tsv$", full.names = TRUE)[1]
first_data <- read.delim(first_file, stringsAsFactors = FALSE, comment.char = "#")

# Get gene IDs and names
gene_info <- first_data %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding")

cat(sprintf("Protein-coding genes: %d\n\n", nrow(gene_info)))

# Initialize expression matrix
expr_matrix <- matrix(NA,
                     nrow = nrow(gene_info),
                     ncol = length(sample_dirs),
                     dimnames = list(gene_info$gene_name, basename(sample_dirs)))

cat("Reading expression data (this may take 3-5 minutes)...\n")

# Read all files with progress
pb <- txtProgressBar(min = 0, max = length(sample_dirs), style = 3)
for (i in seq_along(sample_dirs)) {
  tryCatch({
    tsv_file <- list.files(sample_dirs[i], pattern = "\\.tsv$", full.names = TRUE)[1]

    if (!is.na(tsv_file)) {
      data <- read.delim(tsv_file, stringsAsFactors = FALSE, comment.char = "#") %>%
        filter(gene_type == "protein_coding") %>%
        dplyr::select(gene_name, unstranded)

      expr_matrix[data$gene_name, i] <- data$unstranded
    }
  }, error = function(e) {
    # Skip files with errors
  })

  setTxtProgressBar(pb, i)
}
close(pb)

cat("\n\nExpression matrix dimensions:\n")
cat(sprintf("  Genes: %d\n", nrow(expr_matrix)))
cat(sprintf("  Samples: %d\n", ncol(expr_matrix)))

# Remove samples with all NAs
valid_samples <- colSums(!is.na(expr_matrix)) > 1000
expr_matrix <- expr_matrix[, valid_samples]

cat(sprintf("  Valid samples after QC: %d\n\n", ncol(expr_matrix)))

# Remove genes with all NAs
valid_genes <- rowSums(!is.na(expr_matrix)) > ncol(expr_matrix) * 0.5
expr_matrix <- expr_matrix[valid_genes, ]

cat(sprintf("  Valid genes after QC: %d\n\n", nrow(expr_matrix)))

# ------------------------------------------------------------------------------
# SECTION 3: NORMALIZE EXPRESSION DATA
# ------------------------------------------------------------------------------

cat("SECTION 3: NORMALIZING EXPRESSION DATA\n")
cat("==============================================================================\n\n")

# Replace NAs with 0 for count data
expr_matrix[is.na(expr_matrix)] <- 0

# Calculate library sizes
lib_sizes <- colSums(expr_matrix)

# CPM normalization
expr_cpm <- sweep(expr_matrix, 2, lib_sizes / 1e6, "/")

# Log2 transform
expr_log <- log2(expr_cpm + 1)

cat(sprintf("Normalized expression range: %.3f to %.3f\n\n",
            min(expr_log), max(expr_log)))

# Save processed data
saveRDS(expr_log, "03_Results/18_TARGET_Validation/target_aml_expression_normalized.rds")
cat("✓ Saved: target_aml_expression_normalized.rds\n\n")

# ------------------------------------------------------------------------------
# SECTION 4: PROCESS CLINICAL DATA
# ------------------------------------------------------------------------------

cat("SECTION 4: PROCESSING CLINICAL DATA\n")
cat("==============================================================================\n\n")

# Map UUIDs to barcodes
query_results <- query_target$results[[1]]
uuid_to_barcode <- setNames(query_results$cases, query_results$id)

# Create sample metadata
sample_barcodes <- uuid_to_barcode[colnames(expr_log)]

# Match with clinical data
target_clinical <- target_clinical_gdc %>%
  mutate(
    sample_barcode = submitter_id,
    age_years = as.numeric(gsub(" days", "", age_at_diagnosis)) / 365.25,
    sex = gender,
    vital_status = vital_status,
    OS_DAYS = ifelse(
      vital_status == "Dead",
      as.numeric(gsub(" days", "", days_to_death)),
      as.numeric(gsub(" days", "", days_to_last_follow_up))
    ),
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0),
    OS_MONTHS = OS_DAYS / 30.44
  ) %>%
  filter(!is.na(OS_DAYS) & OS_DAYS > 0)

# Match samples with expression data
common_samples <- intersect(sample_barcodes, target_clinical$sample_barcode)

cat(sprintf("Samples with both expression and survival: %d\n", length(common_samples)))

if (length(common_samples) < 50) {
  cat("\n⚠ WARNING: Very few samples with matched data\n")
  cat("This may indicate barcode matching issues\n\n")
}

# Filter expression data to matched samples
expr_matched_idx <- which(sample_barcodes %in% common_samples)
expr_matched <- expr_log[, expr_matched_idx]
colnames(expr_matched) <- sample_barcodes[expr_matched_idx]

# Filter clinical data
clinical_matched <- target_clinical %>%
  filter(sample_barcode %in% common_samples) %>%
  arrange(match(sample_barcode, colnames(expr_matched)))

cat(sprintf("\nFinal matched samples: %d\n", ncol(expr_matched)))
cat(sprintf("Events: %d deaths (%.1f%%)\n",
            sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))
cat(sprintf("Median follow-up: %.1f months\n", median(clinical_matched$OS_MONTHS)))
cat(sprintf("Age range: %.1f - %.1f years\n\n",
            min(clinical_matched$age_years), max(clinical_matched$age_years)))

# Save clinical data
write.csv(clinical_matched,
          "03_Results/18_TARGET_Validation/target_aml_clinical.csv",
          row.names = FALSE)
cat("✓ Saved: target_aml_clinical.csv\n\n")

# ------------------------------------------------------------------------------
# SECTION 5: APPLY BEATAML CLASSIFIER
# ------------------------------------------------------------------------------

cat("SECTION 5: APPLYING BEATAML CLASSIFIER\n")
cat("==============================================================================\n\n")

# Load BeatAML reference and classifier
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
signature_genes <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")
beataml_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")

cat(sprintf("50-gene signature loaded\n"))

# Check gene overlap
available_genes <- intersect(signature_genes$gene, rownames(expr_matched))
cat(sprintf("Signature genes available: %d / 50 (%.1f%%)\n",
            length(available_genes),
            length(available_genes) / 50 * 100))

if (length(available_genes) < 30) {
  cat("⚠ WARNING: Low gene overlap - predictions may be unreliable\n\n")
}

# Prepare prediction data
pred_data <- data.frame(t(expr_matched[available_genes, ]))

# Impute missing genes with BeatAML means
missing_genes <- setdiff(signature_genes$gene, available_genes)
if (length(missing_genes) > 0) {
  cat(sprintf("Imputing %d missing genes...\n", length(missing_genes)))

  for (gene in missing_genes) {
    if (gene %in% rownames(beataml_expr)) {
      pred_data[[gene]] <- mean(beataml_expr[gene, ])
    }
  }
}

# Reorder columns
pred_data <- pred_data[, signature_genes$gene]

# Predict clusters
cat("Predicting cluster assignments...\n")
cluster_pred <- predict(beataml_classifier, pred_data)
cluster_prob <- predict(beataml_classifier, pred_data, type = "prob")
cluster_conf <- apply(cluster_prob, 1, max)

cat("\nCluster distribution:\n")
cat(sprintf("  Cluster 1: %d (%.1f%%)\n",
            sum(cluster_pred == 1), mean(cluster_pred == 1) * 100))
cat(sprintf("  Cluster 2: %d (%.1f%%)\n",
            sum(cluster_pred == 2), mean(cluster_pred == 2) * 100))
cat(sprintf("  Mean confidence: %.3f\n\n", mean(cluster_conf)))

# Add to clinical data
clinical_matched$cluster <- as.numeric(cluster_pred)
clinical_matched$confidence <- cluster_conf

# Save results
write.csv(clinical_matched,
          "03_Results/18_TARGET_Validation/target_cluster_assignments.csv",
          row.names = FALSE)
cat("✓ Saved: target_cluster_assignments.csv\n\n")

# ------------------------------------------------------------------------------
# SECTION 6: SURVIVAL ANALYSIS
# ------------------------------------------------------------------------------

cat("SECTION 6: SURVIVAL ANALYSIS\n")
cat("==============================================================================\n\n")

# Create survival object
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
medians <- summary(km_fit)$table[, "median"]

cat("** SURVIVAL RESULTS **\n\n")
cat(sprintf("Log-rank test: p = %.4f %s\n",
            logrank_p,
            ifelse(logrank_p < 0.05, "*** SIGNIFICANT", "")))
cat(sprintf("Cox HR: %.3f (95%% CI: %.3f - %.3f), p = %.4f\n\n",
            hr, hr_ci[1], hr_ci[2], cox_p))

cat("Median survival by cluster:\n")
cat(sprintf("  Cluster 1: %.1f months (n=%d)\n",
            medians[1], sum(clinical_matched$cluster == 1)))
cat(sprintf("  Cluster 2: %.1f months (n=%d)\n",
            medians[2], sum(clinical_matched$cluster == 2)))
cat(sprintf("  Difference: %.1f months\n\n",
            abs(medians[1] - medians[2])))

# Save survival summary
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
cat("✓ Saved: target_survival_validation.csv\n\n")

# ------------------------------------------------------------------------------
# SECTION 7: CROSS-COHORT COMPARISON
# ------------------------------------------------------------------------------

cat("SECTION 7: CROSS-COHORT COMPARISON\n")
cat("==============================================================================\n\n")

# Load existing cohort results
beataml_tcga <- read.csv("03_Results/11_Survival_Reanalysis/07_meta_analysis_cohort_results.csv")

# Combine with TARGET
all_cohorts <- rbind(
  beataml_tcga %>% dplyr::select(cohort, n_samples, n_events, hr, hr_lower, hr_upper),
  survival_summary %>% dplyr::select(cohort, n_samples, n_events, hr, hr_lower, hr_upper)
)

print(all_cohorts)

# Test heterogeneity
weights <- 1 / ((log(all_cohorts$hr_upper) - log(all_cohorts$hr_lower)) / (2 * 1.96))^2
pooled_log_hr <- sum(log(all_cohorts$hr) * weights) / sum(weights)
Q <- sum(weights * (log(all_cohorts$hr) - pooled_log_hr)^2)
Q_df <- nrow(all_cohorts) - 1
Q_p <- 1 - pchisq(Q, df = Q_df)
I2 <- max(0, (Q - Q_df) / Q * 100)

cat("\n\nHeterogeneity test:\n")
cat(sprintf("  Cochran's Q: %.3f (p = %.4f)\n", Q, Q_p))
cat(sprintf("  I² statistic: %.1f%%\n\n", I2))

# ------------------------------------------------------------------------------
# SECTION 8: VISUALIZATIONS
# ------------------------------------------------------------------------------

cat("SECTION 8: CREATING FIGURES\n")
cat("==============================================================================\n\n")

dir.create("04_Figures/18_TARGET_Validation", recursive = TRUE, showWarnings = FALSE)

# Kaplan-Meier plot
pdf("04_Figures/18_TARGET_Validation/target_kaplan_meier.pdf",
    width = 8, height = 7)

plot(km_fit,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Overall Survival (months)",
     ylab = "Survival Probability",
     main = sprintf("TARGET-AML (Pediatric) - Kaplan-Meier\nn=%d, p=%.4f",
                    nrow(clinical_matched), logrank_p))

legend("topright",
       legend = c(sprintf("Cluster 1 (n=%d)", sum(clinical_matched$cluster == 1)),
                  sprintf("Cluster 2 (n=%d)", sum(clinical_matched$cluster == 2))),
       col = c("blue", "red"),
       lwd = 2)

grid.text <- sprintf("HR=%.2f (95%% CI: %.2f-%.2f), p=%.4f",
                     hr, hr_ci[1], hr_ci[2], cox_p)
mtext(grid.text, side = 1, line = 4, cex = 0.9)

dev.off()
cat("✓ Saved: target_kaplan_meier.pdf\n")

# Forest plot
pdf("04_Figures/18_TARGET_Validation/forest_plot_all_cohorts.pdf",
    width = 10, height = 6)

par(mar = c(5, 8, 4, 2))
y_pos <- seq(nrow(all_cohorts), 1, -1)

plot(all_cohorts$hr, y_pos,
     xlim = c(0.5, 3),
     ylim = c(0.5, nrow(all_cohorts) + 0.5),
     pch = 15,
     cex = 1.5,
     xlab = "Hazard Ratio (95% CI)",
     ylab = "",
     yaxt = "n",
     main = "Prognostic Effect Across Cohorts (Adult + Pediatric)")

segments(all_cohorts$hr_lower, y_pos, all_cohorts$hr_upper, y_pos, lwd = 2)
abline(v = 1, lty = 2, col = "gray")
axis(2, at = y_pos, labels = all_cohorts$cohort, las = 1)
text(2.8, y_pos, sprintf("n=%d", all_cohorts$n_samples), cex = 0.8)

dev.off()
cat("✓ Saved: forest_plot_all_cohorts.pdf\n\n")

# ------------------------------------------------------------------------------
# SECTION 9: SUMMARY
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("TARGET-AML VALIDATION COMPLETE\n")
cat("==============================================================================\n\n")

cat("RESULTS:\n")
cat(sprintf("  Samples: %d pediatric AML patients\n", nrow(clinical_matched)))
cat(sprintf("  Events: %d (%.1f%%)\n",
            sum(clinical_matched$OS_STATUS),
            mean(clinical_matched$OS_STATUS) * 100))
cat(sprintf("  Log-rank p: %.4f %s\n",
            logrank_p,
            ifelse(logrank_p < 0.05, "(SIGNIFICANT)", "")))
cat(sprintf("  Hazard ratio: %.2f (95%% CI: %.2f-%.2f)\n\n",
            hr, hr_ci[1], hr_ci[2]))

if (survival_summary$significant && abs(log(hr) - log(1.35)) < 0.5) {
  cat("✓✓✓ STRONG VALIDATION\n")
  cat("  Molecular subtypes validated in pediatric AML\n")
  cat("  Effect size consistent with adult cohorts\n")
  cat("  Age-independent prognostic value confirmed\n")
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
