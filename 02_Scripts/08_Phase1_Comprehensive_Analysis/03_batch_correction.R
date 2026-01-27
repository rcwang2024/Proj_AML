#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Task 1: Apply Batch Correction
# ==============================================================================
# Objective: Apply ComBat batch correction to remove centerID batch effects
# Based on: PROJECT_SUMMARY_FOR_CONTINUATION.md recommendations
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(sva)  # For ComBat
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

# Create output directories
dir.create("03_Results/04_Batch_Corrected", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/02_Batch_Correction", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("TASK 1: BATCH CORRECTION OF EXPRESSION DATA\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# Step 1: Load Expression Data
# ------------------------------------------------------------------------------

cat("Step 1: Loading expression data...\n")
expression_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
                        header = TRUE, data.table = FALSE)

# Identify which columns are metadata vs expression data
# First column is gene ID, check which columns have numeric data
first_row <- expression_raw[1, ]
numeric_cols <- sapply(first_row, function(x) !is.na(suppressWarnings(as.numeric(x))))

# Find where numeric columns start
numeric_start <- which(numeric_cols)[1]

cat(sprintf("  First %d columns appear to be metadata\n", numeric_start - 1))

# Extract gene info and expression matrix
gene_info <- expression_raw[, 1:(numeric_start-1)]
colnames(gene_info)[1] <- "stable_id"

expr_data_cols <- numeric_start:ncol(expression_raw)
expr_matrix <- as.matrix(expression_raw[, expr_data_cols])
mode(expr_matrix) <- "numeric"  # Ensure numeric
rownames(expr_matrix) <- expression_raw[[1]]

cat(sprintf("  Expression matrix: %d genes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("  Value range: %.3f to %.3f\n", min(expr_matrix, na.rm=T), max(expr_matrix, na.rm=T)))

# Data appears to be log2-transformed based on range (from investigation)
cat("  Data type: log2-transformed (likely log2(TPM+1))\n\n")

# ------------------------------------------------------------------------------
# Step 2: Load Clinical Data and Extract Batch Variable
# ------------------------------------------------------------------------------

cat("Step 2: Loading clinical data for batch variable...\n")
clinical_raw <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
clinical_raw <- as.data.frame(clinical_raw)

cat(sprintf("  Clinical data: %d samples\n", nrow(clinical_raw)))

# Check for centerID (sequencing center)
if (!"sequencingCenter" %in% colnames(clinical_raw)) {
  cat("  ⚠ WARNING: 'sequencingCenter' column not found\n")
  cat("  Available columns with 'center' in name:\n")
  center_cols <- grep("center", colnames(clinical_raw), ignore.case = TRUE, value = TRUE)
  cat(paste("   ", center_cols, collapse = "\n"), "\n")

  # Use centerID if available
  if ("centerID" %in% colnames(clinical_raw)) {
    batch_col <- "centerID"
    cat(sprintf("  Using '%s' as batch variable\n\n", batch_col))
  } else {
    stop("Cannot find batch variable (centerID or sequencingCenter)")
  }
} else {
  batch_col <- "sequencingCenter"
  cat(sprintf("  Using '%s' as batch variable\n\n", batch_col))
}

# Map batch info to expression samples
expr_sample_ids <- colnames(expr_matrix)
clinical_rna <- clinical_raw[!is.na(clinical_raw$dbgap_rnaseq_sample), ]

# Create batch mapping
batch_mapping <- data.frame(
  sample_id = expr_sample_ids,
  stringsAsFactors = FALSE
)

batch_mapping <- merge(batch_mapping,
                       clinical_rna[, c("dbgap_rnaseq_sample", batch_col)],
                       by.x = "sample_id",
                       by.y = "dbgap_rnaseq_sample",
                       all.x = TRUE)

colnames(batch_mapping)[2] <- "batch"

# Summary of batch variable
cat("Batch variable distribution:\n")
print(table(batch_mapping$batch, useNA = "always"))
cat("\n")

samples_with_batch <- sum(!is.na(batch_mapping$batch))
cat(sprintf("Samples with batch info: %d / %d (%.1f%%)\n\n",
            samples_with_batch, nrow(batch_mapping),
            samples_with_batch/nrow(batch_mapping)*100))

# ------------------------------------------------------------------------------
# Step 3: Pre-Batch Correction PCA
# ------------------------------------------------------------------------------

cat("Step 3: Running PCA before batch correction...\n")

# Keep only samples with batch info for fair comparison
samples_to_keep <- batch_mapping$sample_id[!is.na(batch_mapping$batch)]
expr_subset <- expr_matrix[, samples_to_keep]
batch_subset <- batch_mapping$batch[!is.na(batch_mapping$batch)]

cat(sprintf("  Samples for analysis: %d\n", ncol(expr_subset)))

# PCA
pca_before <- prcomp(t(expr_subset), scale. = TRUE, center = TRUE)
pca_df_before <- data.frame(
  sample_id = colnames(expr_subset),
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  PC3 = pca_before$x[, 3],
  PC4 = pca_before$x[, 4],
  PC5 = pca_before$x[, 5],
  batch = batch_subset,
  stringsAsFactors = FALSE
)

# Variance explained
var_explained_before <- summary(pca_before)$importance[2, 1:10] * 100

cat("\nVariance explained by PCs (before correction):\n")
for (i in 1:5) {
  cat(sprintf("  PC%d: %.2f%%\n", i, var_explained_before[i]))
}
cat("\n")

# Test for batch effect using ANOVA
aov_pc1_before <- aov(PC1 ~ batch, data = pca_df_before)
aov_pc2_before <- aov(PC2 ~ batch, data = pca_df_before)

pval_pc1_before <- summary(aov_pc1_before)[[1]]$`Pr(>F)`[1]
pval_pc2_before <- summary(aov_pc2_before)[[1]]$`Pr(>F)`[1]

cat("Statistical test for batch effects (ANOVA):\n")
cat(sprintf("  PC1 ~ batch: F = %.2f, p = %.2e\n",
            summary(aov_pc1_before)[[1]]$`F value`[1], pval_pc1_before))
cat(sprintf("  PC2 ~ batch: F = %.2f, p = %.2e\n",
            summary(aov_pc2_before)[[1]]$`F value`[1], pval_pc2_before))

if (pval_pc1_before < 0.05 || pval_pc2_before < 0.05) {
  cat("\n  ⚠ SIGNIFICANT BATCH EFFECT DETECTED - Correction needed\n\n")
  apply_correction <- TRUE
} else {
  cat("\n  ✓ No significant batch effect detected\n\n")
  apply_correction <- FALSE
}

# Plot PCA before correction
pdf("04_Figures/02_Batch_Correction/PCA_before_batch_correction.pdf",
    width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# PC1 vs PC2 colored by batch
batch_colors <- rainbow(length(unique(batch_subset)))
names(batch_colors) <- unique(batch_subset)

plot(pca_df_before$PC1, pca_df_before$PC2,
     xlab = paste0("PC1 (", round(var_explained_before[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained_before[2], 1), "%)"),
     main = "PCA Before Batch Correction\n(Colored by Batch)",
     pch = 19, col = batch_colors[pca_df_before$batch], cex = 0.8)
legend("topright", legend = names(batch_colors), col = batch_colors,
       pch = 19, cex = 0.7, title = "Batch")

# PC3 vs PC4
plot(pca_df_before$PC3, pca_df_before$PC4,
     xlab = paste0("PC3 (", round(var_explained_before[3], 1), "%)"),
     ylab = paste0("PC4 (", round(var_explained_before[4], 1), "%)"),
     main = "PC3 vs PC4 Before Correction",
     pch = 19, col = batch_colors[pca_df_before$batch], cex = 0.8)

# Boxplot of PC1 by batch
boxplot(PC1 ~ batch, data = pca_df_before,
        main = "PC1 Distribution by Batch",
        xlab = "Batch", ylab = "PC1",
        col = batch_colors, las = 2)

# Boxplot of PC2 by batch
boxplot(PC2 ~ batch, data = pca_df_before,
        main = "PC2 Distribution by Batch",
        xlab = "Batch", ylab = "PC2",
        col = batch_colors, las = 2)

dev.off()

cat("✓ Saved: 04_Figures/02_Batch_Correction/PCA_before_batch_correction.pdf\n\n")

# ------------------------------------------------------------------------------
# Step 4: Apply ComBat Batch Correction
# ------------------------------------------------------------------------------

if (apply_correction) {
  cat("Step 4: Applying ComBat batch correction...\n")

  cat(sprintf("  Samples: %d\n", ncol(expr_subset)))
  cat(sprintf("  Genes: %d\n", nrow(expr_subset)))
  cat(sprintf("  Batches: %d unique levels\n", length(unique(batch_subset))))

  # Apply ComBat
  # Note: mod=NULL means we're not preserving any biological covariates
  # In a real analysis, you might want to preserve age, sex, etc.
  cat("\n  Running ComBat (this may take a few minutes)...\n")

  expr_corrected <- ComBat(
    dat = expr_subset,
    batch = batch_subset,
    mod = NULL,  # No covariates to preserve
    par.prior = TRUE,  # Use parametric adjustments
    prior.plots = FALSE
  )

  cat("  ✓ ComBat completed successfully\n\n")

  # ------------------------------------------------------------------------------
  # Step 5: Post-Correction PCA
  # ------------------------------------------------------------------------------

  cat("Step 5: Running PCA after batch correction...\n")

  pca_after <- prcomp(t(expr_corrected), scale. = TRUE, center = TRUE)
  pca_df_after <- data.frame(
    sample_id = colnames(expr_corrected),
    PC1 = pca_after$x[, 1],
    PC2 = pca_after$x[, 2],
    PC3 = pca_after$x[, 3],
    PC4 = pca_after$x[, 4],
    PC5 = pca_after$x[, 5],
    batch = batch_subset,
    stringsAsFactors = FALSE
  )

  var_explained_after <- summary(pca_after)$importance[2, 1:10] * 100

  cat("\nVariance explained by PCs (after correction):\n")
  for (i in 1:5) {
    cat(sprintf("  PC%d: %.2f%%\n", i, var_explained_after[i]))
  }
  cat("\n")

  # Test batch effect after correction
  aov_pc1_after <- aov(PC1 ~ batch, data = pca_df_after)
  aov_pc2_after <- aov(PC2 ~ batch, data = pca_df_after)

  pval_pc1_after <- summary(aov_pc1_after)[[1]]$`Pr(>F)`[1]
  pval_pc2_after <- summary(aov_pc2_after)[[1]]$`Pr(>F)`[1]

  cat("Statistical test after correction:\n")
  cat(sprintf("  PC1 ~ batch: F = %.2f, p = %.2e\n",
              summary(aov_pc1_after)[[1]]$`F value`[1], pval_pc1_after))
  cat(sprintf("  PC2 ~ batch: F = %.2f, p = %.2e\n",
              summary(aov_pc2_after)[[1]]$`F value`[1], pval_pc2_after))

  if (pval_pc1_after < 0.05 || pval_pc2_after < 0.05) {
    cat("\n  ⚠ WARNING: Some batch effect remains\n\n")
  } else {
    cat("\n  ✓ Batch effect successfully removed\n\n")
  }

  # Plot PCA after correction
  pdf("04_Figures/02_Batch_Correction/PCA_after_batch_correction.pdf",
      width = 12, height = 10)

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

  plot(pca_df_after$PC1, pca_df_after$PC2,
       xlab = paste0("PC1 (", round(var_explained_after[1], 1), "%)"),
       ylab = paste0("PC2 (", round(var_explained_after[2], 1), "%)"),
       main = "PCA After Batch Correction\n(Colored by Batch)",
       pch = 19, col = batch_colors[pca_df_after$batch], cex = 0.8)
  legend("topright", legend = names(batch_colors), col = batch_colors,
         pch = 19, cex = 0.7, title = "Batch")

  plot(pca_df_after$PC3, pca_df_after$PC4,
       xlab = paste0("PC3 (", round(var_explained_after[3], 1), "%)"),
       ylab = paste0("PC4 (", round(var_explained_after[4], 1), "%)"),
       main = "PC3 vs PC4 After Correction",
       pch = 19, col = batch_colors[pca_df_after$batch], cex = 0.8)

  boxplot(PC1 ~ batch, data = pca_df_after,
          main = "PC1 Distribution After Correction",
          xlab = "Batch", ylab = "PC1",
          col = batch_colors, las = 2)

  boxplot(PC2 ~ batch, data = pca_df_after,
          main = "PC2 Distribution After Correction",
          xlab = "Batch", ylab = "PC2",
          col = batch_colors, las = 2)

  dev.off()

  cat("✓ Saved: 04_Figures/02_Batch_Correction/PCA_after_batch_correction.pdf\n\n")

  # ------------------------------------------------------------------------------
  # Step 6: Save Batch-Corrected Data
  # ------------------------------------------------------------------------------

  cat("Step 6: Saving batch-corrected expression data...\n")

  # Reconstruct full expression matrix with batch-corrected values
  expr_final <- expr_matrix
  expr_final[, samples_to_keep] <- expr_corrected

  # Save as text file (same format as input)
  expr_output <- cbind(gene_info, expr_final)
  expr_output <- as.data.frame(expr_output, stringsAsFactors = FALSE)

  fwrite(expr_output,
         "03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.txt",
         sep = "\t", quote = FALSE, row.names = FALSE)

  cat("✓ Saved: 03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.txt\n")

  # Also save as RDS for faster loading
  saveRDS(expr_final,
          "03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

  cat("✓ Saved: 03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds\n\n")

  # Save summary report
  summary_report <- data.frame(
    Metric = c("Total Samples", "Samples with Batch Info", "Samples Corrected",
               "Number of Batches", "Genes",
               "PC1 ~ Batch p-value (before)", "PC1 ~ Batch p-value (after)",
               "PC2 ~ Batch p-value (before)", "PC2 ~ Batch p-value (after)",
               "Batch Effect Removed"),
    Value = c(
      ncol(expr_matrix),
      samples_with_batch,
      ncol(expr_corrected),
      length(unique(batch_subset)),
      nrow(expr_matrix),
      sprintf("%.2e", pval_pc1_before),
      sprintf("%.2e", pval_pc1_after),
      sprintf("%.2e", pval_pc2_before),
      sprintf("%.2e", pval_pc2_after),
      ifelse(pval_pc1_after > 0.05 && pval_pc2_after > 0.05, "YES", "PARTIAL")
    ),
    stringsAsFactors = FALSE
  )

  write.csv(summary_report,
            "03_Results/04_Batch_Corrected/batch_correction_summary.csv",
            row.names = FALSE)

  cat("✓ Saved: 03_Results/04_Batch_Corrected/batch_correction_summary.csv\n\n")

} else {
  cat("Step 4: Batch correction not needed (no significant batch effect)\n")
  cat("  Saving original expression data to output directory...\n\n")

  expr_output <- expression_raw
  fwrite(expr_output,
         "03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.txt",
         sep = "\t", quote = FALSE, row.names = FALSE)

  saveRDS(expr_matrix,
          "03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
}

# ------------------------------------------------------------------------------
# Final Summary
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("BATCH CORRECTION COMPLETE\n")
cat("==============================================================================\n\n")

cat("Summary:\n")
if (apply_correction) {
  cat(sprintf("  ✓ Batch effects detected (PC1 p=%.2e, PC2 p=%.2e)\n",
              pval_pc1_before, pval_pc2_before))
  cat("  ✓ ComBat batch correction applied\n")
  cat(sprintf("  ✓ Post-correction: PC1 p=%.2e, PC2 p=%.2e\n",
              pval_pc1_after, pval_pc2_after))
  cat(sprintf("  ✓ Batch-corrected data saved (%d genes × %d samples)\n",
              nrow(expr_matrix), ncol(expr_matrix)))
} else {
  cat("  ✓ No significant batch effects detected\n")
  cat("  ✓ Original data copied to output directory\n")
}

cat("\nOutput files:\n")
cat("  - 03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.txt\n")
cat("  - 03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds\n")
cat("  - 03_Results/04_Batch_Corrected/batch_correction_summary.csv\n")
cat("  - 04_Figures/02_Batch_Correction/PCA_before_batch_correction.pdf\n")
if (apply_correction) {
  cat("  - 04_Figures/02_Batch_Correction/PCA_after_batch_correction.pdf\n")
}

cat("\nNext steps:\n")
cat("  → Review batch correction figures\n")
cat("  → Review outlier samples (7 identified)\n")
cat("  → Proceed with molecular subtyping using batch-corrected data\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
