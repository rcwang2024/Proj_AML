#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Multi-Omics Integration Project
# Phase 1: Comprehensive Data Foundation & Quality Control
# ==============================================================================
# Author: AML Research Analysis Pipeline
# Date: 2025-10-04
# Description: Publication-quality data loading, QC, preprocessing, integration
# Based on: 2024-2025 AML research best practices
# ==============================================================================

# Suppress startup messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(sva)          # For ComBat batch correction
  library(pheatmap)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# Set options for reproducibility
set.seed(42)
options(stringsAsFactors = FALSE)

# Set working directory
setwd("D:/Projects/Project_AML")

# Create organized output directory structure
dirs_to_create <- c(
  "data/processed",
  "data/metadata",
  "results/qc",
  "results/tables",
  "results/figures",
  "results/phase1_outputs"
)

for (dir in dirs_to_create) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Function to print section headers
print_header <- function(text, level = 1) {
  width <- 80
  if (level == 1) {
    cat(paste0("\n", strrep("=", width), "\n"))
    cat(paste0(text, "\n"))
    cat(paste0(strrep("=", width), "\n\n"))
  } else {
    cat(paste0("\n", strrep("-", width), "\n"))
    cat(paste0(text, "\n"))
    cat(paste0(strrep("-", width), "\n\n"))
  }
}

# Function to save results with timestamp
save_result <- function(obj, filename) {
  if (grepl("\\.rds$", filename)) {
    saveRDS(obj, filename)
  } else if (grepl("\\.csv$", filename)) {
    write.csv(obj, filename, row.names = FALSE)
  } else {
    save(obj, file = filename)
  }
  cat(paste0("✓ Saved: ", filename, "\n"))
}

# ==============================================================================
# PHASE 1.1: DATA LOADING AND INITIAL EXAMINATION
# ==============================================================================

print_header("PHASE 1.1: DATA LOADING AND INITIAL EXAMINATION")

cat("Loading BeatAML data files...\n\n")

# 1. Expression Data
cat("1. Loading expression data...\n")
expr_start <- Sys.time()
expression_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
                        header = TRUE, data.table = FALSE)
expr_time <- difftime(Sys.time(), expr_start, units = "secs")

cat(sprintf("   Dimensions: %d genes × %d samples\n",
            nrow(expression_raw), ncol(expression_raw) - 1))
cat(sprintf("   Load time: %.1f seconds\n", expr_time))
cat(sprintf("   File size: %.1f MB\n",
            file.info("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt")$size / 1e6))
cat(sprintf("   Gene ID column: '%s'\n\n", colnames(expression_raw)[1]))

# 2. Mutation Data
cat("2. Loading mutation data...\n")
mutations_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
                       header = TRUE, data.table = FALSE)

cat(sprintf("   Total mutation calls: %d\n", nrow(mutations_raw)))
cat(sprintf("   Columns: %d\n", ncol(mutations_raw)))
cat("   Column names:\n")
cat(paste0("      ", paste(colnames(mutations_raw), collapse = ", "), "\n\n"))

# 3. Clinical Data
cat("3. Loading clinical data...\n")
clinical_raw <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
clinical_raw <- as.data.frame(clinical_raw)

cat(sprintf("   Samples: %d\n", nrow(clinical_raw)))
cat(sprintf("   Variables: %d\n", ncol(clinical_raw)))
cat("   First 15 variables:\n")
cat(paste0("      ", paste(colnames(clinical_raw)[1:min(15, ncol(clinical_raw))],
                           collapse = ", "), "\n\n"))

# 4. Drug Response Data
cat("4. Loading drug response data...\n")
drug_response_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt",
                           header = TRUE, data.table = FALSE)

cat(sprintf("   Samples: %d\n", nrow(drug_response_raw)))
cat(sprintf("   Drugs tested: %d\n", ncol(drug_response_raw) - 1))
cat(sprintf("   First column (sample ID): '%s'\n\n", colnames(drug_response_raw)[1]))

# 5. Drug Family Annotations
cat("5. Loading drug family annotations...\n")
drug_families <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_drug_families.xlsx")
drug_families <- as.data.frame(drug_families)

cat(sprintf("   Annotated drugs: %d\n", nrow(drug_families)))
cat(sprintf("   Annotation columns: %s\n\n",
            paste(colnames(drug_families), collapse = ", ")))

# Save initial data summary
initial_summary <- data.frame(
  Data_Type = c("Expression", "Mutations", "Clinical", "Drug_Response", "Drug_Annotations"),
  N_Rows = c(nrow(expression_raw), nrow(mutations_raw), nrow(clinical_raw),
             nrow(drug_response_raw), nrow(drug_families)),
  N_Cols = c(ncol(expression_raw), ncol(mutations_raw), ncol(clinical_raw),
             ncol(drug_response_raw), ncol(drug_families)),
  File_Size_MB = c(
    file.info("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt")$size / 1e6,
    file.info("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt")$size / 1e6,
    file.info("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")$size / 1e6,
    file.info("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt")$size / 1e6,
    file.info("01_Data/BeatAML_Downloaded_Data/beataml_drug_families.xlsx")$size / 1e6
  )
)

save_result(initial_summary, "results/tables/initial_data_summary.csv")

# ==============================================================================
# PHASE 1.2: MASTER SAMPLE ID MAPPING
# ==============================================================================

print_header("PHASE 1.2: MASTER SAMPLE ID MAPPING")

cat("Extracting sample IDs from each dataset...\n\n")

# Extract sample IDs - need to identify the correct column names
sample_ids_expression <- colnames(expression_raw)[-1]  # All columns except gene names

# For mutations - find the sample ID column
mut_sample_col <- grep("sample|dbgap", colnames(mutations_raw),
                       ignore.case = TRUE, value = TRUE)[1]
if (is.na(mut_sample_col)) {
  cat("⚠ Could not auto-detect mutation sample column, using column 2\n")
  mut_sample_col <- colnames(mutations_raw)[2]
}
sample_ids_mutations <- unique(mutations_raw[[mut_sample_col]])

# For clinical - find the sample ID column
clin_sample_col <- grep("sample|dbgap|id", colnames(clinical_raw),
                        ignore.case = TRUE, value = TRUE)[1]
if (is.na(clin_sample_col)) {
  cat("⚠ Could not auto-detect clinical sample column, using column 1\n")
  clin_sample_col <- colnames(clinical_raw)[1]
}
sample_ids_clinical <- clinical_raw[[clin_sample_col]]

# For drug response - first column should be sample ID
drug_sample_col <- colnames(drug_response_raw)[1]
sample_ids_drugs <- drug_response_raw[[drug_sample_col]]

cat(sprintf("Expression samples: %d (column source: expression matrix columns)\n",
            length(sample_ids_expression)))
cat(sprintf("Mutation samples: %d (column: %s)\n",
            length(sample_ids_mutations), mut_sample_col))
cat(sprintf("Clinical samples: %d (column: %s)\n",
            length(sample_ids_clinical), clin_sample_col))
cat(sprintf("Drug response samples: %d (column: %s)\n\n",
            length(sample_ids_drugs), drug_sample_col))

# Create comprehensive sample manifest
all_samples <- unique(c(sample_ids_expression, sample_ids_mutations,
                        sample_ids_clinical, sample_ids_drugs))

cat(sprintf("Total unique samples across all datasets: %d\n\n", length(all_samples)))

# Build overlap matrix
overlap_data <- data.frame(
  sample_id = all_samples,
  has_expression = all_samples %in% sample_ids_expression,
  has_mutations = all_samples %in% sample_ids_mutations,
  has_clinical = all_samples %in% sample_ids_clinical,
  has_drugs = all_samples %in% sample_ids_drugs
)

# Add cohort classification
overlap_data$n_datatypes <- rowSums(overlap_data[, c("has_expression", "has_mutations",
                                                      "has_clinical", "has_drugs")])

overlap_data$cohort_category <- case_when(
  overlap_data$n_datatypes == 4 ~ "Complete_All4",
  overlap_data$n_datatypes == 3 ~ "Partial_3types",
  overlap_data$n_datatypes == 2 ~ "Partial_2types",
  overlap_data$n_datatypes == 1 ~ "Single_type",
  TRUE ~ "Unknown"
)

# Summary statistics
print_header("DATA AVAILABILITY SUMMARY", level = 2)

data_avail_summary <- data.frame(
  Data_Type = c("Expression", "Mutations", "Clinical", "Drug_Response"),
  N_Samples = c(sum(overlap_data$has_expression),
                sum(overlap_data$has_mutations),
                sum(overlap_data$has_clinical),
                sum(overlap_data$has_drugs)),
  Percent_of_Total = sprintf("%.1f%%", c(
    sum(overlap_data$has_expression) / nrow(overlap_data) * 100,
    sum(overlap_data$has_mutations) / nrow(overlap_data) * 100,
    sum(overlap_data$has_clinical) / nrow(overlap_data) * 100,
    sum(overlap_data$has_drugs) / nrow(overlap_data) * 100
  ))
)

print(data_avail_summary)
cat("\n")

# Key overlaps
overlap_expr_mut <- sum(overlap_data$has_expression & overlap_data$has_mutations)
overlap_expr_mut_clin <- sum(overlap_data$has_expression & overlap_data$has_mutations &
                              overlap_data$has_clinical)
overlap_all4 <- sum(overlap_data$n_datatypes == 4)

print_header("KEY SAMPLE OVERLAPS", level = 2)
cat(sprintf("Expression + Mutations: %d samples (%.1f%%)\n",
            overlap_expr_mut, overlap_expr_mut / nrow(overlap_data) * 100))
cat(sprintf("Expression + Mutations + Clinical: %d samples (%.1f%%)\n",
            overlap_expr_mut_clin, overlap_expr_mut_clin / nrow(overlap_data) * 100))
cat(sprintf("ALL 4 data types (GOLD STANDARD): %d samples (%.1f%%)\n\n",
            overlap_all4, overlap_all4 / nrow(overlap_data) * 100))

# Cohort breakdown
cat("Cohort distribution by number of data types:\n")
cohort_table <- table(overlap_data$cohort_category)
print(cohort_table)
cat("\n")

# Save master sample manifest
save_result(overlap_data, "data/processed/master_sample_manifest.csv")

# Create a Venn diagram-style summary
overlap_combinations <- overlap_data %>%
  group_by(has_expression, has_mutations, has_clinical, has_drugs) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(desc(count))

save_result(overlap_combinations, "results/tables/sample_overlap_combinations.csv")

cat("\nTop 10 sample overlap combinations:\n")
print(head(overlap_combinations, 10))
cat("\n")

# ==============================================================================
# SAVE CHECKPOINT
# ==============================================================================

cat("\n")
print_header("CHECKPOINT: PHASE 1.1-1.2 COMPLETE", level = 2)
cat("Saved files:\n")
cat("  - results/tables/initial_data_summary.csv\n")
cat("  - data/processed/master_sample_manifest.csv\n")
cat("  - results/tables/sample_overlap_combinations.csv\n\n")

cat("Next steps:\n")
cat("  → Phase 1.3: Generate comprehensive QC report\n")
cat("  → Phase 1.4: Expression data preprocessing\n\n")

# Save workspace for continuation
save.image("data/processed/phase1_checkpoint_1.RData")
cat("✓ Workspace saved: data/processed/phase1_checkpoint_1.RData\n\n")

cat("=" * 80, "\n")
cat("Phase 1.1-1.2 execution complete!\n")
cat("Run time: ", format(Sys.time()), "\n")
cat("=" * 80, "\n\n")
