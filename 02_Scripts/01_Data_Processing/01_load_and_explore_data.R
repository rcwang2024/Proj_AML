#!/usr/bin/env Rscript
# ============================================================================
# Beat AML Multi-Omics Data Loading and Exploration
# ============================================================================
# Purpose: Load and explore all Beat AML data types
# Author: AML Research Team
# Date: 2025-10-02
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(data.table)
  library(ComplexHeatmap)
  library(pheatmap)
})

# Set working directory
setwd("D:/Projects/Project_AML")

# Create output directories
dir.create("03_Results/01_Processed_Data", recursive = TRUE, showWarnings = FALSE)
dir.create("03_Results/02_QC_Reports", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/01_QC_Figures", recursive = TRUE, showWarnings = FALSE)

cat("====================================================================\n")
cat("BEAT AML MULTI-OMICS DATA LOADING\n")
cat("====================================================================\n\n")

# ============================================================================
# 1. LOAD EXPRESSION DATA
# ============================================================================
cat("1. Loading Expression Data...\n")

expr_data <- fread("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
                   data.table = FALSE)

cat("   - Dimensions:", nrow(expr_data), "genes x", ncol(expr_data)-4, "samples\n")
cat("   - Columns:", paste(colnames(expr_data)[1:10], collapse = ", "), "...\n")

# Extract gene info and expression matrix
gene_info <- expr_data[, 1:4]
colnames(gene_info) <- c("ensembl_id", "gene_symbol", "description", "biotype")

expr_matrix <- expr_data[, -c(1:4)]
rownames(expr_matrix) <- gene_info$gene_symbol

cat("   - Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
cat("   - Sample IDs format:", paste(head(colnames(expr_matrix), 3), collapse = ", "), "\n\n")

# ============================================================================
# 2. LOAD CLINICAL DATA
# ============================================================================
cat("2. Loading Clinical Data...\n")

clinical_data <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")

cat("   - Dimensions:", nrow(clinical_data), "patients x", ncol(clinical_data), "variables\n")
cat("   - Key variables:", paste(colnames(clinical_data)[1:10], collapse = ", "), "...\n")

# Check for survival data
survival_cols <- grep("survival|death|vital|days|age|gender|sex",
                      colnames(clinical_data),
                      ignore.case = TRUE, value = TRUE)
cat("   - Survival-related columns:", paste(survival_cols[1:min(5, length(survival_cols))], collapse = ", "), "...\n\n")

# ============================================================================
# 3. LOAD MUTATION DATA
# ============================================================================
cat("3. Loading Mutation Data...\n")

mutation_data <- fread("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
                       data.table = FALSE)

cat("   - Dimensions:", nrow(mutation_data), "variants x", ncol(mutation_data), "columns\n")
cat("   - Columns:", paste(colnames(mutation_data)[1:10], collapse = ", "), "...\n")

# Count samples and genes with mutations
n_mut_samples <- length(unique(mutation_data$dbgap_sample_id))
n_mut_genes <- length(unique(mutation_data$gene))

cat("   - Unique samples:", n_mut_samples, "\n")
cat("   - Unique genes:", n_mut_genes, "\n")

# Top mutated genes
top_genes <- mutation_data %>%
  count(gene, sort = TRUE) %>%
  head(10)

cat("   - Top mutated genes:\n")
print(top_genes)
cat("\n")

# ============================================================================
# 4. LOAD DRUG RESPONSE DATA
# ============================================================================
cat("4. Loading Drug Response Data (AUC)...\n")

drug_auc <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt",
                  data.table = FALSE)

cat("   - Dimensions:", nrow(drug_auc), "measurements x", ncol(drug_auc), "columns\n")
cat("   - Columns:", paste(colnames(drug_auc)[1:15], collapse = ", "), "...\n")

# Count drugs and samples
n_drugs <- length(unique(drug_auc$inhibitor))
n_drug_samples <- length(unique(drug_auc$dbgap_subject_id))

cat("   - Unique drugs:", n_drugs, "\n")
cat("   - Unique samples with drug data:", n_drug_samples, "\n")

# Top drugs tested
top_drugs <- drug_auc %>%
  count(inhibitor, sort = TRUE) %>%
  head(10)

cat("   - Top tested drugs:\n")
print(top_drugs)
cat("\n")

# ============================================================================
# 5. SAMPLE OVERLAP ANALYSIS
# ============================================================================
cat("5. Analyzing Sample Overlap Across Data Types...\n\n")

# Extract sample IDs (need to standardize formats)
expr_samples <- colnames(expr_matrix)
clinical_samples <- clinical_data$dbgap_rnaseq_sample
mut_samples <- unique(mutation_data$dbgap_sample_id)
drug_samples <- unique(drug_auc$dbgap_rnaseq_sample)

# Create overlap summary
overlap_summary <- data.frame(
  Data_Type = c("Expression", "Clinical", "Mutation", "Drug_Response"),
  N_Samples = c(
    length(expr_samples),
    sum(!is.na(clinical_samples)),
    length(mut_samples),
    sum(!is.na(drug_samples))
  )
)

cat("Sample counts by data type:\n")
print(overlap_summary)
cat("\n")

# Find complete cases (all 4 data types)
# Note: Need to match sample IDs carefully
all_samples <- Reduce(intersect, list(
  expr_samples,
  clinical_samples[!is.na(clinical_samples)],
  mut_samples,
  drug_samples[!is.na(drug_samples)]
))

cat("Samples with ALL 4 data types:", length(all_samples), "\n")

# Expression + Clinical + Mutation
expr_clin_mut <- Reduce(intersect, list(
  expr_samples,
  clinical_samples[!is.na(clinical_samples)],
  mut_samples
))
cat("Samples with Expression + Clinical + Mutation:", length(expr_clin_mut), "\n")

# Expression + Clinical + Drug
expr_clin_drug <- Reduce(intersect, list(
  expr_samples,
  clinical_samples[!is.na(clinical_samples)],
  drug_samples[!is.na(drug_samples)]
))
cat("Samples with Expression + Clinical + Drug:", length(expr_clin_drug), "\n")

# Expression + Clinical
expr_clin <- intersect(expr_samples, clinical_samples[!is.na(clinical_samples)])
cat("Samples with Expression + Clinical:", length(expr_clin), "\n\n")

# ============================================================================
# 6. BASIC QUALITY METRICS
# ============================================================================
cat("6. Computing Basic Quality Metrics...\n\n")

# Expression QC
cat("Expression QC:\n")
cat("   - Median expression per sample:", median(colMeans(expr_matrix, na.rm = TRUE)), "\n")
cat("   - % genes with mean > 0:",
    round(100 * mean(rowMeans(expr_matrix, na.rm = TRUE) > 0), 2), "%\n")
cat("   - % missing values:", round(100 * mean(is.na(expr_matrix)), 4), "%\n\n")

# Mutation QC
cat("Mutation QC:\n")
vaf_summary <- summary(mutation_data$t_vaf)
cat("   - Tumor VAF distribution:\n")
print(vaf_summary)
cat("\n")

# Variant types
var_types <- mutation_data %>%
  count(variant_classification, sort = TRUE) %>%
  head(10)
cat("   - Top variant types:\n")
print(var_types)
cat("\n")

# Drug response QC
cat("Drug Response QC:\n")
auc_summary <- summary(drug_auc$auc, na.rm = TRUE)
cat("   - AUC distribution:\n")
print(auc_summary)
cat("\n")

# ============================================================================
# 7. SAVE PROCESSED DATA
# ============================================================================
cat("7. Saving Processed Data...\n\n")

# Save gene info
write.csv(gene_info,
          "03_Results/01_Processed_Data/gene_annotations.csv",
          row.names = FALSE)

# Save expression matrix
saveRDS(expr_matrix,
        "03_Results/01_Processed_Data/expression_matrix.rds")

# Save clinical data
saveRDS(clinical_data,
        "03_Results/01_Processed_Data/clinical_data.rds")

# Save mutation data
saveRDS(mutation_data,
        "03_Results/01_Processed_Data/mutation_data.rds")

# Save drug response data
saveRDS(drug_auc,
        "03_Results/01_Processed_Data/drug_response_auc.rds")

# Save overlap summary
write.csv(overlap_summary,
          "03_Results/02_QC_Reports/sample_overlap_summary.csv",
          row.names = FALSE)

cat("Processed data saved to: 03_Results/01_Processed_Data/\n")
cat("QC reports saved to: 03_Results/02_QC_Reports/\n\n")

# ============================================================================
# 8. CREATE SUMMARY REPORT
# ============================================================================
cat("8. Creating Summary Report...\n\n")

summary_report <- paste0(
  "====================================================================\n",
  "BEAT AML MULTI-OMICS DATA SUMMARY\n",
  "====================================================================\n\n",
  "DATA LOADED:\n",
  "1. Expression: ", nrow(expr_matrix), " genes × ", ncol(expr_matrix), " samples\n",
  "2. Clinical: ", nrow(clinical_data), " patients × ", ncol(clinical_data), " variables\n",
  "3. Mutation: ", nrow(mutation_data), " variants (", n_mut_samples, " samples, ", n_mut_genes, " genes)\n",
  "4. Drug Response: ", nrow(drug_auc), " measurements (", n_drug_samples, " samples, ", n_drugs, " drugs)\n\n",
  "SAMPLE OVERLAP:\n",
  "- All 4 data types: ", length(all_samples), " samples\n",
  "- Expression + Clinical + Mutation: ", length(expr_clin_mut), " samples\n",
  "- Expression + Clinical + Drug: ", length(expr_clin_drug), " samples\n",
  "- Expression + Clinical: ", length(expr_clin), " samples\n\n",
  "TOP MUTATED GENES:\n",
  paste(capture.output(print(top_genes)), collapse = "\n"), "\n\n",
  "TOP TESTED DRUGS:\n",
  paste(capture.output(print(top_drugs)), collapse = "\n"), "\n\n",
  "====================================================================\n"
)

writeLines(summary_report, "03_Results/02_QC_Reports/data_loading_summary.txt")

cat(summary_report)

cat("\n✓ Data loading complete!\n")
cat("✓ Summary saved to: 03_Results/02_QC_Reports/data_loading_summary.txt\n\n")

# ============================================================================
# 9. SESSION INFO
# ============================================================================
sink("03_Results/02_QC_Reports/session_info.txt")
cat("R Session Information\n")
cat("======================\n\n")
print(sessionInfo())
sink()

cat("Session info saved to: 03_Results/02_QC_Reports/session_info.txt\n")
cat("\n====================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("====================================================================\n")
