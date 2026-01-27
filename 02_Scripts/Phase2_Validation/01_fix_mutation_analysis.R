# TASK 1.1: Diagnose Sample ID Mismatch
# Fix mutation analysis - highest priority

library(tidyverse)
library(data.table)

cat("=== TASK 1: FIX MUTATION ANALYSIS ===\n\n")

# 1.1 Diagnose Sample ID Mismatch
cat("### 1.1 Diagnosing Sample ID Mismatch ###\n\n")

# Load expression data (batch-corrected)
if (file.exists("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")) {
  expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
} else {
  cat("ERROR: Batch-corrected expression data not found.\n")
  cat("Trying alternative path...\n")
  expr_data <- readRDS("03_Results/04_Batch_Corrected_Data/beataml_expression_batchcorrected.rds")
}

expr_sample_ids <- colnames(expr_data)
cat("Expression data loaded:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n\n")

# Load mutation data
mutation_files <- c(
  "01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
  "01_Data/beataml_mutations.txt",
  "01_Data/mutations.txt"
)

mutation_data <- NULL
mutation_file_used <- NULL

for (file in mutation_files) {
  if (file.exists(file)) {
    cat("Found mutation file:", file, "\n")
    mutation_data <- fread(file)
    mutation_file_used <- file
    break
  }
}

if (is.null(mutation_data)) {
  cat("ERROR: No mutation file found.\n")
  cat("Expression sample IDs look like:\n")
  print(head(expr_sample_ids, 20))
  cat("\nPlease provide mutation data file path.\n")
  quit(save = "no", status = 1)
}

cat("Mutation data loaded:", nrow(mutation_data), "rows\n")
cat("Mutation data columns:\n")
print(colnames(mutation_data))
cat("\n")

# Identify sample ID column
possible_id_cols <- c("Tumor_Sample_Barcode", "sample_id", "Sample",
                      "patient_id", "case_id", "submitter_id", "dbgap_sample_id")

id_col <- NULL
for (col in possible_id_cols) {
  if (col %in% colnames(mutation_data)) {
    id_col <- col
    cat("Found sample ID column:", col, "\n")
    break
  }
}

if (is.null(id_col)) {
  cat("ERROR: Cannot identify sample ID column in mutation data\n")
  cat("Available columns:\n")
  print(colnames(mutation_data))
  quit(save = "no", status = 1)
}

mut_sample_ids <- unique(mutation_data[[id_col]])

cat("\n=== SAMPLE ID DIAGNOSIS ===\n")
cat("Expression samples (n):", length(expr_sample_ids), "\n")
cat("Mutation samples (n):", length(mut_sample_ids), "\n\n")

cat("Expression ID examples (first 15):\n")
print(head(expr_sample_ids, 15))

cat("\nMutation ID examples (first 15):\n")
print(head(mut_sample_ids, 15))

# Try direct matching
direct_match <- intersect(expr_sample_ids, mut_sample_ids)
cat("\nDirect matches:", length(direct_match), "\n\n")

if (length(direct_match) == 0) {
  cat("⚠️ NO DIRECT MATCHES - Need to harmonize IDs\n\n")

  # Try common transformations
  cat("Trying common transformations...\n\n")

  # Expression IDs look like they have 'R' suffix (RNA)
  # Mutation IDs might have 'D' suffix (DNA)
  # Common pattern: remove suffix and match on base ID

  # For expression: extract base ID (before R/D suffix)
  expr_clean <- str_replace(expr_sample_ids, "R$", "")  # Remove trailing R
  expr_clean <- str_replace(expr_clean, "D$", "")      # Remove trailing D if any

  # For mutations: extract base ID
  mut_clean <- str_replace(mut_sample_ids, "R$", "")
  mut_clean <- str_replace(mut_clean, "D$", "")

  cat("After removing R/D suffixes:\n")
  cat("Expression IDs (cleaned, first 10):\n")
  print(head(expr_clean, 10))
  cat("\nMutation IDs (cleaned, first 10):\n")
  print(head(mut_clean, 10))

  cleaned_match <- intersect(expr_clean, mut_clean)
  cat("\nMatches after cleaning:", length(cleaned_match), "\n\n")

  if (length(cleaned_match) > 0) {
    cat("✓ SUCCESS: Found", length(cleaned_match), "matches after cleaning\n\n")

    # Create mapping table
    expr_mapping <- data.frame(
      original_id = expr_sample_ids,
      cleaned_id = expr_clean,
      stringsAsFactors = FALSE
    )

    mut_mapping <- data.frame(
      mutation_id = mut_sample_ids,
      cleaned_id = mut_clean,
      stringsAsFactors = FALSE
    )

    # Join to get matched pairs
    id_mapping <- inner_join(expr_mapping, mut_mapping, by = "cleaned_id")

    cat("Final matched samples:", nrow(id_mapping), "\n")
    cat("\nSample of mapping (first 10 rows):\n")
    print(head(id_mapping, 10))

    # Save mapping
    write.csv(id_mapping, "03_Results/10_Mutations/sample_id_mapping.csv", row.names = FALSE)
    cat("\n✓ Saved ID mapping to: 03_Results/10_Mutations/sample_id_mapping.csv\n")

    # Also save to root 03_Results for easier access
    write.csv(id_mapping, "03_Results/sample_id_mapping.csv", row.names = FALSE)
    cat("✓ Also saved to: 03_Results/sample_id_mapping.csv\n\n")

  } else {
    cat("\n⚠️ STILL NO MATCHES - Manual intervention required\n")
    cat("Saving sample IDs for manual inspection...\n")

    write.csv(data.frame(expr_ids = expr_sample_ids),
              "03_Results/10_Mutations/expression_sample_ids.csv", row.names = FALSE)
    write.csv(data.frame(mut_ids = mut_sample_ids),
              "03_Results/10_Mutations/mutation_sample_ids.csv", row.names = FALSE)

    cat("Saved expression and mutation IDs for inspection\n")
  }
} else {
  cat("✓ SUCCESS: Found", length(direct_match), "direct matches\n")

  # Create simple mapping
  id_mapping <- data.frame(
    original_id = direct_match,
    cleaned_id = direct_match,
    mutation_id = direct_match,
    stringsAsFactors = FALSE
  )

  write.csv(id_mapping, "03_Results/10_Mutations/sample_id_mapping.csv", row.names = FALSE)
  write.csv(id_mapping, "03_Results/sample_id_mapping.csv", row.names = FALSE)
  cat("✓ Saved ID mapping\n")
}

cat("\n### Task 1.1 COMPLETE ###\n")
