#!/usr/bin/env Rscript
# Phase 4 Part 1: Sample Attrition Analysis (Simplified)
# Purpose: Track sample sizes from expression to multivariate analysis

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("SAMPLE ATTRITION ANALYSIS\n")
cat("==============================================================================\n\n")

# Create output directory
dir.create("03_Results/21_Manuscript_Prep", recursive = TRUE, showWarnings = FALSE)

# Step 1: Expression data
expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
n_expr <- ncol(expr_data)
cat(sprintf("Step 1 - Expression data: %d samples\n", n_expr))

# Step 2: Cluster assignments
cluster_data <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
n_cluster <- nrow(cluster_data)
cat(sprintf("Step 2 - With cluster assignments: %d samples\n", n_cluster))

# Step 3: Survival data with clusters
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
n_survival <- nrow(survival_data)
cat(sprintf("Step 3 - Merged with survival data: %d samples\n", n_survival))

# Step 4: Complete survival data (no NA, OS_months > 0)
survival_complete <- survival_data %>%
  filter(!is.na(OS_months) & OS_months > 0 & !is.na(OS_event))
n_survival_complete <- nrow(survival_complete)
cat(sprintf("Step 4 - Complete survival data: %d samples\n", n_survival_complete))

# Step 5: Load mutation data
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)
n_mutations <- nrow(mutations)
cat(sprintf("Step 5 - Mutation data available: %d samples\n", n_mutations))

# Step 6: Merge survival + mutations
merged_data <- survival_complete %>%
  left_join(
    mutations %>% rownames_to_column("sample_id"),
    by = "sample_id"
  )
n_merged <- nrow(merged_data)
cat(sprintf("Step 6 - Survival + mutations merged: %d samples\n", n_merged))

# Step 7: Clean and filter for multivariate
# Select key prognostic mutations
key_mutations <- c("FLT3", "NPM1", "DNMT3A", "TP53", "IDH1", "IDH2",
                   "RUNX1", "ASXL1", "TET2", "NRAS", "KRAS")
available_mutations <- key_mutations[key_mutations %in% colnames(merged_data)]

analysis_data <- merged_data %>%
  rename(
    OS_MONTHS = OS_months,
    OS_STATUS = OS_event,
    cluster_assignment = cluster,
    AGE = age,
    SEX = sex
  ) %>%
  filter(!is.na(AGE) & !is.na(SEX) & !is.na(cluster_assignment))

n_clinical_complete <- nrow(analysis_data)
cat(sprintf("Step 7 - Complete clinical data: %d samples\n", n_clinical_complete))

# Step 8: Complete for multivariate (no missing mutations)
# Check how many complete cases for each mutation
cat("\n=== MISSING DATA BY MUTATION ===\n")
for (mut in available_mutations) {
  if (mut %in% colnames(analysis_data)) {
    n_missing <- sum(is.na(analysis_data[[mut]]))
    pct_missing <- n_missing / nrow(analysis_data) * 100
    cat(sprintf("  %s: %d missing (%.1f%%)\n", mut, n_missing, pct_missing))
  }
}

# Select mutations used in actual multivariate model
# From Phase 3, these were: TP53, TET2, RUNX1, ASXL1
model_mutations <- c("TP53", "TET2", "RUNX1", "ASXL1")

multivar_complete <- analysis_data
for (mut in model_mutations) {
  if (mut %in% colnames(multivar_complete)) {
    multivar_complete <- multivar_complete %>% filter(!is.na(!!sym(mut)))
  }
}

n_multivar <- nrow(multivar_complete)
cat(sprintf("\nStep 8 - Complete for multivariate (no missing key mutations): %d samples\n", n_multivar))

# Create attrition table
attrition_table <- data.frame(
  Step = 1:8,
  Stage = c(
    "Expression data",
    "With cluster assignment",
    "With survival data",
    "Complete survival (no NA)",
    "Mutation data available",
    "Survival + mutations merged",
    "Complete clinical (age, sex)",
    "Complete for multivariate"
  ),
  N = c(n_expr, n_cluster, n_survival, n_survival_complete,
        n_mutations, n_merged, n_clinical_complete, n_multivar),
  Description = c(
    "Batch-corrected expression matrix",
    "Samples successfully clustered (k=2)",
    "Merged with clinical/survival data",
    "Non-missing OS_months > 0",
    "Samples with mutation profiling",
    "Combined survival + mutation data",
    "Non-missing age and sex",
    "No missing values in key mutations (TP53, TET2, RUNX1, ASXL1)"
  )
)

# Calculate exclusions
attrition_table$N_excluded <- c(0, -diff(attrition_table$N))
attrition_table$Percent_retained <- round(attrition_table$N / attrition_table$N[1] * 100, 1)

cat("\n=== ATTRITION TABLE ===\n\n")
print(attrition_table, row.names = FALSE)

write.csv(attrition_table,
          "03_Results/21_Manuscript_Prep/sample_attrition_table.csv",
          row.names = FALSE)

# Summary statistics
cat("\n=== SUMMARY ===\n\n")
cat(sprintf("Starting samples: %d\n", n_expr))
cat(sprintf("Final multivariate n: %d\n", n_multivar))
cat(sprintf("Total excluded: %d (%.1f%%)\n",
            n_expr - n_multivar,
            (n_expr - n_multivar) / n_expr * 100))

# Key exclusion reasons
cat("\n=== KEY EXCLUSION REASONS ===\n\n")
cat(sprintf("Survival data unavailable/invalid: %d\n", n_cluster - n_survival_complete))
cat(sprintf("Mutation data unavailable: %d\n", n_survival_complete - n_merged))
cat(sprintf("Missing clinical variables: %d\n", n_merged - n_clinical_complete))
cat(sprintf("Missing key mutations: %d\n", n_clinical_complete - n_multivar))

cat("\n✅ Sample attrition analysis complete\n")
cat("   File: 03_Results/21_Manuscript_Prep/sample_attrition_table.csv\n")
