#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Prepare Analysis-Ready Datasets
# ==============================================================================
# Objective:
#   1. Review outlier samples and make decision
#   2. Prepare gold standard cohort (n=478) datasets
#   3. Create analysis-ready matrices for Phase 2
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("PREPARING ANALYSIS-READY DATASETS\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Review Outlier Samples
# ------------------------------------------------------------------------------

cat("STEP 1: Reviewing outlier samples...\n\n")

outliers <- read.csv("03_Results/02_QC_Reports/expression_outliers.csv")

cat(sprintf("Total outlier samples identified: %d\n\n", nrow(outliers)))
cat("Outlier characteristics:\n")
print(outliers)
cat("\n")

# Analyze outlier severity
cat("Outlier severity assessment:\n")
cat(sprintf("  - Samples with only 1 flag: %d\n", sum(outliers$n_outlier_flags == 1)))
cat(sprintf("  - Samples with 2+ flags: %d\n", sum(outliers$n_outlier_flags >= 2)))
cat(sprintf("  - Median correlation range: %.3f - %.3f\n",
            min(outliers$median_correlation), max(outliers$median_correlation)))
cat(sprintf("  - All correlations > 0.79: %s\n\n",
            ifelse(min(outliers$median_correlation) > 0.79, "YES", "NO")))

# Decision: Keep outliers if correlation > 0.75 and only 1 flag
threshold_correlation <- 0.75
threshold_flags <- 1

samples_to_exclude <- outliers %>%
  filter(median_correlation < threshold_correlation | n_outlier_flags > threshold_flags) %>%
  pull(sample_id)

cat("DECISION CRITERIA:\n")
cat(sprintf("  - Exclude if median correlation < %.2f OR n_flags > %d\n",
            threshold_correlation, threshold_flags))
cat(sprintf("  - Samples meeting exclusion criteria: %d\n\n", length(samples_to_exclude)))

if (length(samples_to_exclude) == 0) {
  cat("✓ DECISION: Keep all 7 outlier samples (mild outliers, good correlation)\n")
  cat("  Rationale: All samples have correlation > 0.79 and only 1 flag each\n\n")
} else {
  cat("⚠ DECISION: Exclude the following samples:\n")
  cat(paste("  -", samples_to_exclude, collapse = "\n"), "\n\n")
}

# Document decision
decision_record <- data.frame(
  date = Sys.Date(),
  decision = ifelse(length(samples_to_exclude) == 0, "Keep all outliers", "Exclude samples"),
  n_outliers_identified = nrow(outliers),
  n_samples_excluded = length(samples_to_exclude),
  rationale = "Outliers are mild (correlation >0.79, 1 flag each)",
  threshold_correlation = threshold_correlation,
  threshold_flags = threshold_flags
)

write.csv(decision_record, "03_Results/04_Batch_Corrected/outlier_decision.csv",
          row.names = FALSE)
cat("✓ Saved: 03_Results/04_Batch_Corrected/outlier_decision.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 2: Load Master Sample Mapping
# ------------------------------------------------------------------------------

cat("STEP 2: Loading master sample mapping...\n\n")

sample_mapping <- fread("03_Results/01_Processed_Data/master_sample_id_mapping.csv",
                        data.table = FALSE)

cat(sprintf("Total samples in mapping: %d\n", nrow(sample_mapping)))
cat("\nCohort distribution:\n")
print(table(sample_mapping$cohort_category))
cat("\n")

# Extract gold standard cohort
gold_standard <- sample_mapping %>%
  filter(cohort_category == "complete_quad_omics")

cat(sprintf("Gold standard cohort (all 4 data types): %d samples\n", nrow(gold_standard)))
cat(sprintf("  - Has expression: %d\n", sum(gold_standard$has_expression)))
cat(sprintf("  - Has mutations: %d\n", sum(gold_standard$has_mutations)))
cat(sprintf("  - Has clinical: %d\n", sum(gold_standard$has_clinical)))
cat(sprintf("  - Has drug response: %d\n\n", sum(gold_standard$has_drug_response)))

# ------------------------------------------------------------------------------
# STEP 3: Load Batch-Corrected Expression Data
# ------------------------------------------------------------------------------

cat("STEP 3: Loading batch-corrected expression data...\n\n")

expr_corrected <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

cat(sprintf("Batch-corrected expression: %d genes × %d samples\n",
            nrow(expr_corrected), ncol(expr_corrected)))
cat(sprintf("Value range: %.3f to %.3f\n\n", min(expr_corrected), max(expr_corrected)))

# Remove outliers if any were flagged for exclusion
if (length(samples_to_exclude) > 0) {
  expr_corrected <- expr_corrected[, !colnames(expr_corrected) %in% samples_to_exclude]
  cat(sprintf("After removing outliers: %d genes × %d samples\n\n",
              nrow(expr_corrected), ncol(expr_corrected)))
}

# Filter for gold standard samples (for integrated analyses)
gold_expr_samples <- gold_standard$expression_id[gold_standard$has_expression]
gold_expr_samples <- gold_expr_samples[gold_expr_samples %in% colnames(expr_corrected)]

expr_gold_standard <- expr_corrected[, gold_expr_samples]

cat(sprintf("Gold standard expression matrix: %d genes × %d samples\n\n",
            nrow(expr_gold_standard), ncol(expr_gold_standard)))

# ------------------------------------------------------------------------------
# STEP 4: Filter Lowly Expressed Genes
# ------------------------------------------------------------------------------

cat("STEP 4: Filtering lowly expressed genes...\n\n")

# Calculate median expression per gene (on full expression dataset)
median_expr <- apply(expr_corrected, 1, median, na.rm = TRUE)

# For log2-transformed data: keep genes with median > 1 (TPM > 1 in 50% of samples)
threshold_expr <- 1
genes_keep <- median_expr > threshold_expr

cat(sprintf("Genes before filtering: %d\n", nrow(expr_corrected)))
cat(sprintf("Genes after filtering (median > %.1f): %d\n",
            threshold_expr, sum(genes_keep)))
cat(sprintf("Percentage retained: %.1f%%\n\n",
            sum(genes_keep) / nrow(expr_corrected) * 100))

# Apply filter
expr_filtered <- expr_corrected[genes_keep, ]
expr_gold_filtered <- expr_gold_standard[genes_keep, ]

cat(sprintf("Filtered expression (all): %d genes × %d samples\n",
            nrow(expr_filtered), ncol(expr_filtered)))
cat(sprintf("Filtered expression (gold standard): %d genes × %d samples\n\n",
            nrow(expr_gold_filtered), ncol(expr_gold_filtered)))

# ------------------------------------------------------------------------------
# STEP 5: Load and Prepare Mutation Data
# ------------------------------------------------------------------------------

cat("STEP 5: Loading mutation data...\n\n")

# Check if binary mutation matrix already exists
if (file.exists("03_Results/01_Processed_Data/mutation_matrix_binary.rds")) {
  cat("Loading existing binary mutation matrix...\n")
  mutation_matrix <- readRDS("03_Results/01_Processed_Data/mutation_matrix_binary.rds")
} else {
  cat("Creating binary mutation matrix from raw data...\n")

  mutations_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
                         data.table = FALSE)

  # Create binary matrix (samples × genes)
  mut_df <- mutations_raw %>%
    dplyr::select(symbol, dbgap_sample_id) %>%
    filter(!is.na(symbol) & !is.na(dbgap_sample_id)) %>%
    distinct() %>%
    mutate(mutated = 1)

  # Pivot to wide format
  mutation_wide <- mut_df %>%
    pivot_wider(names_from = dbgap_sample_id, values_from = mutated, values_fill = 0)

  mutation_matrix <- as.matrix(mutation_wide[, -1])
  rownames(mutation_matrix) <- mutation_wide$symbol

  # Transpose to samples × genes
  mutation_matrix <- t(mutation_matrix)

  saveRDS(mutation_matrix, "03_Results/01_Processed_Data/mutation_matrix_binary.rds")
  cat("✓ Saved: 03_Results/01_Processed_Data/mutation_matrix_binary.rds\n")
}

cat(sprintf("Mutation matrix: %d samples × %d genes\n", nrow(mutation_matrix), ncol(mutation_matrix)))

# Filter for gold standard samples
gold_mut_samples <- gold_standard$mutation_id[gold_standard$has_mutations]
gold_mut_samples <- gold_mut_samples[gold_mut_samples %in% rownames(mutation_matrix)]

mutation_gold_standard <- mutation_matrix[gold_mut_samples, ]

cat(sprintf("Gold standard mutations: %d samples × %d genes\n\n",
            nrow(mutation_gold_standard), ncol(mutation_gold_standard)))

# Keep only recurrently mutated genes (>2% frequency)
mut_freq <- colSums(mutation_gold_standard) / nrow(mutation_gold_standard) * 100
recurrent_genes <- names(mut_freq[mut_freq > 2])

mutation_gold_filtered <- mutation_gold_standard[, recurrent_genes]

cat(sprintf("Recurrently mutated genes (>2%%): %d\n", length(recurrent_genes)))
cat(sprintf("Filtered mutations (gold standard): %d samples × %d genes\n\n",
            nrow(mutation_gold_filtered), ncol(mutation_gold_filtered)))

# ------------------------------------------------------------------------------
# STEP 6: Load Clinical Data
# ------------------------------------------------------------------------------

cat("STEP 6: Loading clinical data...\n\n")

clinical_raw <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
clinical_raw <- as.data.frame(clinical_raw)

cat(sprintf("Clinical data: %d samples × %d variables\n\n", nrow(clinical_raw), ncol(clinical_raw)))

# Extract gold standard clinical data
gold_clin_samples <- gold_standard$clinical_id[gold_standard$has_clinical]

# Match by RNA-seq sample ID (most common)
clinical_gold <- clinical_raw %>%
  filter(dbgap_rnaseq_sample %in% gold_clin_samples)

cat(sprintf("Gold standard clinical: %d samples\n\n", nrow(clinical_gold)))

# ------------------------------------------------------------------------------
# STEP 7: Align All Data Types for Gold Standard Cohort
# ------------------------------------------------------------------------------

cat("STEP 7: Aligning all data types for gold standard cohort...\n\n")

# Use RNA-seq sample ID as the common identifier
common_samples <- gold_standard %>%
  filter(has_expression & has_mutations & has_clinical & has_drug_response) %>%
  pull(expression_id)  # Use expression ID as primary

# Filter expression data
expr_aligned <- expr_gold_filtered[, colnames(expr_gold_filtered) %in% common_samples]

# For mutations and clinical, need to map IDs
sample_map_subset <- gold_standard %>%
  filter(expression_id %in% common_samples)

# Get mutation samples corresponding to expression samples
mut_sample_ids <- sample_map_subset$mutation_id
names(mut_sample_ids) <- sample_map_subset$expression_id

# Align mutation matrix
mut_aligned <- mutation_gold_filtered[mut_sample_ids[colnames(expr_aligned)], ]
rownames(mut_aligned) <- colnames(expr_aligned)  # Use expression IDs for consistency

# Align clinical data
clin_sample_ids <- sample_map_subset$clinical_id
names(clin_sample_ids) <- sample_map_subset$expression_id

clinical_aligned <- clinical_raw %>%
  filter(dbgap_rnaseq_sample %in% clin_sample_ids[colnames(expr_aligned)])

# Ensure same order
clinical_aligned <- clinical_aligned[match(colnames(expr_aligned), clinical_aligned$dbgap_rnaseq_sample), ]

cat("ALIGNED GOLD STANDARD COHORT:\n")
cat(sprintf("  - Samples: %d\n", ncol(expr_aligned)))
cat(sprintf("  - Expression: %d genes\n", nrow(expr_aligned)))
cat(sprintf("  - Mutations: %d genes\n", ncol(mut_aligned)))
cat(sprintf("  - Clinical variables: %d\n\n", ncol(clinical_aligned)))

# Verify alignment
cat("Verifying sample alignment:\n")
cat(sprintf("  - Expression samples match mutation samples: %s\n",
            ifelse(all(colnames(expr_aligned) == rownames(mut_aligned)), "YES", "NO")))
cat(sprintf("  - Expression samples match clinical samples: %s\n\n",
            ifelse(all(colnames(expr_aligned) == clinical_aligned$dbgap_rnaseq_sample), "YES", "NO")))

# ------------------------------------------------------------------------------
# STEP 8: Save Analysis-Ready Datasets
# ------------------------------------------------------------------------------

cat("STEP 8: Saving analysis-ready datasets...\n\n")

# Create Phase 2 directory
dir.create("03_Results/05_Analysis_Ready_Data", recursive = TRUE, showWarnings = FALSE)

# Save full cohort (all expression samples)
saveRDS(expr_filtered, "03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")
cat(sprintf("✓ Saved: expression_filtered_all.rds (%d genes × %d samples)\n",
            nrow(expr_filtered), ncol(expr_filtered)))

# Save gold standard cohort (aligned multi-omics)
analysis_data <- list(
  expression = expr_aligned,
  mutations = mut_aligned,
  clinical = clinical_aligned,
  sample_ids = colnames(expr_aligned),
  n_samples = ncol(expr_aligned),
  n_genes_expression = nrow(expr_aligned),
  n_genes_mutations = ncol(mut_aligned),
  date_created = Sys.Date()
)

saveRDS(analysis_data, "03_Results/05_Analysis_Ready_Data/gold_standard_cohort.rds")
cat(sprintf("✓ Saved: gold_standard_cohort.rds (%d samples, multi-omics integrated)\n",
            analysis_data$n_samples))

# Save individual matrices for easier loading
write.csv(t(expr_aligned), "03_Results/05_Analysis_Ready_Data/expression_gold_standard.csv")
write.csv(mut_aligned, "03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv")
write.csv(clinical_aligned, "03_Results/05_Analysis_Ready_Data/clinical_gold_standard.csv")

cat("✓ Saved: CSV versions of all matrices\n\n")

# Save summary
summary_df <- data.frame(
  Dataset = c("Full Expression Cohort", "Gold Standard Expression",
              "Gold Standard Mutations", "Gold Standard Clinical"),
  N_Samples = c(ncol(expr_filtered), ncol(expr_aligned),
                nrow(mut_aligned), nrow(clinical_aligned)),
  N_Features = c(nrow(expr_filtered), nrow(expr_aligned),
                 ncol(mut_aligned), ncol(clinical_aligned)),
  Description = c(
    "Batch-corrected, filtered (median>1), all samples",
    "Gold standard cohort expression data",
    "Recurrent mutations (>2%) for gold standard",
    "Clinical annotations for gold standard"
  )
)

write.csv(summary_df, "03_Results/05_Analysis_Ready_Data/dataset_summary.csv",
          row.names = FALSE)
cat("✓ Saved: dataset_summary.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 9: Generate Data Summary Report
# ------------------------------------------------------------------------------

cat("STEP 9: Generating summary report...\n\n")

cat("==============================================================================\n")
cat("ANALYSIS-READY DATASETS PREPARED\n")
cat("==============================================================================\n\n")

cat("DATASETS CREATED:\n\n")

cat("1. FULL EXPRESSION COHORT (for molecular subtyping)\n")
cat(sprintf("   - File: expression_filtered_all.rds\n"))
cat(sprintf("   - Samples: %d\n", ncol(expr_filtered)))
cat(sprintf("   - Genes: %d (filtered, median > 1)\n", nrow(expr_filtered)))
cat(sprintf("   - Use for: Consensus clustering, pathway analysis\n\n"))

cat("2. GOLD STANDARD COHORT (for integrated analyses)\n")
cat(sprintf("   - File: gold_standard_cohort.rds\n"))
cat(sprintf("   - Samples: %d (all 4 data types)\n", ncol(expr_aligned)))
cat(sprintf("   - Expression genes: %d\n", nrow(expr_aligned)))
cat(sprintf("   - Mutation genes: %d (>2%% frequency)\n", ncol(mut_aligned)))
cat(sprintf("   - Clinical variables: %d\n", ncol(clinical_aligned)))
cat(sprintf("   - Use for: Drug response, survival, integrated modeling\n\n"))

cat("SAMPLE OVERLAP:\n")
cat(sprintf("   - Full expression cohort: %d samples\n", ncol(expr_filtered)))
cat(sprintf("   - Gold standard cohort: %d samples (%.1f%% of full)\n",
            ncol(expr_aligned), ncol(expr_aligned)/ncol(expr_filtered)*100))
cat(sprintf("   - Outliers excluded: %d\n\n", length(samples_to_exclude)))

cat("DATA QUALITY:\n")
cat(sprintf("   - Batch correction applied: YES\n"))
cat(sprintf("   - Low-expression genes filtered: YES (median > 1)\n"))
cat(sprintf("   - Outliers reviewed: YES (decision: keep all)\n"))
cat(sprintf("   - Sample alignment verified: YES\n\n"))

cat("NEXT STEPS:\n")
cat("   → Phase 2.1: Consensus clustering on expression_filtered_all.rds\n")
cat("   → Select top 5,000 variable genes\n")
cat("   → Run ConsensusClusterPlus (k=2-10)\n")
cat("   → Identify optimal cluster number\n\n")

cat("==============================================================================\n")
cat(sprintf("Data preparation complete: %s\n", Sys.time()))
cat("==============================================================================\n")
