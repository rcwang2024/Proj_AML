#!/usr/bin/env Rscript
# Phase 4 Part 1: Sample Attrition Analysis
# Purpose: Document every step of sample filtering to explain n=459 in multivariate

suppressPackageStartupMessages({
  library(dplyr)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("SAMPLE ATTRITION ANALYSIS\n")
cat("==============================================================================\n\n")

# Create output directories
dir.create("03_Results/21_Manuscript_Prep", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/20_Manuscript_Prep", recursive = TRUE, showWarnings = FALSE)

# Load all datasets
cat("Loading datasets...\n")

# Expression data
expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
cat(sprintf("Expression data: %d genes x %d samples\n", nrow(expr_data), ncol(expr_data)))

# Cluster assignments
cluster_assignments <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
cat(sprintf("Cluster assignments: %d samples\n", nrow(cluster_assignments)))

# Clinical data
clinical_data <- read.csv("03_Results/01_Processed_Data/beataml_clinical_processed.csv")
cat(sprintf("Clinical data: %d samples\n", nrow(clinical_data)))

# Survival data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
cat(sprintf("Survival data: %d samples\n", nrow(survival_data)))

# Mutation data
mutation_matrix <- readRDS("03_Results/10_Mutations/mutation_matrix.rds")
cat(sprintf("Mutation data: %d samples\n", nrow(mutation_matrix)))

cat("\n")

# Step-by-step tracking
cat("=== STEP-BY-STEP SAMPLE TRACKING ===\n\n")

# Step 1: Expression data
step1_n <- ncol(expr_data)
cat(sprintf("Step 1 - Expression data: %d samples\n", step1_n))

# Step 2: With cluster assignments
step2_n <- nrow(cluster_assignments)
cat(sprintf("Step 2 - With cluster assignments: %d samples\n", step2_n))

# Step 3: With clinical data
step3_n <- nrow(clinical_data)
cat(sprintf("Step 3 - With clinical data: %d samples\n", step3_n))

# Step 4: With complete survival data
survival_complete <- survival_data %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS))
step4_n <- nrow(survival_complete)
cat(sprintf("Step 4 - With survival data (no NA): %d samples\n", step4_n))

# Step 5: Survival + Clusters
survival_cluster <- survival_data %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS) & !is.na(cluster_assignment))
step5_n <- nrow(survival_cluster)
cat(sprintf("Step 5 - Survival + Clusters: %d samples (used for univariate Cox)\n", step5_n))

# Step 6: With mutation data
multivar_data <- survival_cluster %>%
  inner_join(mutation_matrix, by="lab_id")
step6_n <- nrow(multivar_data)
cat(sprintf("Step 6 - With mutation data merged: %d samples\n", step6_n))

# Step 7: Complete for multivariate (no missing in key variables)
multivar_complete <- multivar_data %>%
  filter(
    !is.na(NPM1), !is.na(TP53), !is.na(RUNX1), !is.na(DNMT3A),
    !is.na(FLT3), !is.na(TET2), !is.na(ASXL1),
    !is.na(age_at_diagnosis), !is.na(sex)
  )
step7_n <- nrow(multivar_complete)
cat(sprintf("Step 7 - Complete for multivariate (no missing values): %d samples\n", step7_n))

# Create attrition table
attrition_table <- data.frame(
  Step = c(
    "1. Expression data",
    "2. With cluster assignment",
    "3. With clinical data",
    "4. With survival data",
    "5. Survival + Clusters",
    "6. With mutation data",
    "7. Complete for multivariate"
  ),
  N = c(step1_n, step2_n, step3_n, step4_n, step5_n, step6_n, step7_n),
  Description = c(
    "Raw expression matrix",
    "Samples successfully clustered",
    "Samples with clinical information",
    "Complete survival outcomes (no NA)",
    "Complete for survival analysis",
    "Expression + survival + mutations merged",
    "No missing values in key variables"
  )
)

attrition_table$N_excluded <- c(0, -diff(attrition_table$N))
attrition_table$Percent_retained <- round(attrition_table$N / attrition_table$N[1] * 100, 1)

cat("\n=== ATTRITION TABLE ===\n")
print(attrition_table)

write.csv(attrition_table,
          "03_Results/21_Manuscript_Prep/sample_attrition_table.csv",
          row.names = FALSE)

# Analyze reasons for exclusion
cat("\n=== REASONS FOR MULTIVARIATE EXCLUSION ===\n")

# Check missing mutations
missing_mutations <- multivar_data %>%
  summarise(
    Missing_NPM1 = sum(is.na(NPM1)),
    Missing_TP53 = sum(is.na(TP53)),
    Missing_RUNX1 = sum(is.na(RUNX1)),
    Missing_DNMT3A = sum(is.na(DNMT3A)),
    Missing_FLT3 = sum(is.na(FLT3)),
    Missing_TET2 = sum(is.na(TET2)),
    Missing_ASXL1 = sum(is.na(ASXL1)),
    Missing_any_mutation = sum(
      is.na(NPM1) | is.na(TP53) | is.na(RUNX1) | is.na(TET2) |
      is.na(FLT3) | is.na(DNMT3A) | is.na(ASXL1)
    )
  )

cat("\nMissing mutation data:\n")
print(t(missing_mutations))

# Check missing clinical
missing_clinical <- multivar_data %>%
  summarise(
    Missing_AGE = sum(is.na(age_at_diagnosis)),
    Missing_SEX = sum(is.na(sex))
  )

cat("\nMissing clinical data:\n")
print(t(missing_clinical))

# Total excluded
total_excluded <- step5_n - step7_n
cat(sprintf("\nTotal excluded from multivariate: %d samples (%.1f%%)\n",
            total_excluded,
            total_excluded / step5_n * 100))

# Save detailed exclusion
exclusion_reasons <- data.frame(
  Reason = c("Missing mutation data", "Missing age", "Missing sex", "Total excluded", "Final n"),
  N = c(
    missing_mutations$Missing_any_mutation,
    missing_clinical$Missing_AGE,
    missing_clinical$Missing_SEX,
    total_excluded,
    step7_n
  )
)

write.csv(exclusion_reasons,
          "03_Results/21_Manuscript_Prep/multivariate_exclusion_reasons.csv",
          row.names = FALSE)

cat("\n✅ Sample attrition analysis complete\n")
cat("\nFiles generated:\n")
cat("  - 03_Results/21_Manuscript_Prep/sample_attrition_table.csv\n")
cat("  - 03_Results/21_Manuscript_Prep/multivariate_exclusion_reasons.csv\n")
