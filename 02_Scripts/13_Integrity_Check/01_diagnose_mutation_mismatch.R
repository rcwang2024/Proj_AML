#!/usr/bin/env Rscript
# TASK 1.1: Diagnose Mutation Analysis Sample ID Mismatch
# Date: 2025-10-09

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("═══════════════════════════════════════════════════════════════\n")
cat("  TASK 1.1: DIAGNOSING MUTATION SAMPLE ID MISMATCH\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load Expression Data Sample IDs
# ------------------------------------------------------------------------------
cat("STEP 1: Loading expression data sample IDs...\n")

expr_file <- "03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds"
if (!file.exists(expr_file)) {
  stop("❌ Expression data not found: ", expr_file)
}

expr_data <- readRDS(expr_file)
expr_samples <- colnames(expr_data)
cat(sprintf("✓ Expression data: %d genes × %d samples\n", nrow(expr_data), ncol(expr_data)))
cat(sprintf("✓ Sample ID format (first 5): %s\n", paste(head(expr_samples, 5), collapse=", ")))
cat("\n")

# ------------------------------------------------------------------------------
# STEP 2: Search for Mutation Data Files
# ------------------------------------------------------------------------------
cat("STEP 2: Searching for mutation data files...\n")

mutation_file_candidates <- c(
  "01_Data/BeatAML_Downloaded_Data/beataml_wes_waveA_mutations.txt",
  "01_Data/BeatAML_Downloaded_Data/beataml_wes_waveB_mutations.txt",
  "01_Data/BeatAML_Downloaded_Data/beataml_wes_combined_mutations.txt",
  "01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
  "01_Data/BeatAML_Downloaded_Data/beataml_variants.txt"
)

mutation_files_found <- list()
for (f in mutation_file_candidates) {
  if (file.exists(f)) {
    mutation_files_found[[f]] <- f
    cat(sprintf("  ✓ Found: %s\n", f))
  }
}

if (length(mutation_files_found) == 0) {
  cat("\n❌ No mutation files found in expected locations.\n")
  cat("Searching entire 01_Data directory...\n\n")

  # Search for any files with "mutation" or "variant" in name
  all_data_files <- list.files("01_Data", pattern = ".*", recursive = TRUE, full.names = TRUE)
  mutation_pattern <- grepl("mutation|variant|maf|vcf|wes", all_data_files, ignore.case = TRUE)
  potential_files <- all_data_files[mutation_pattern]

  if (length(potential_files) > 0) {
    cat("Found potential mutation files:\n")
    for (f in potential_files) {
      cat(sprintf("  - %s\n", f))
    }
  } else {
    cat("❌ No mutation-related files found in 01_Data directory.\n")
  }
  cat("\n")
} else {
  cat(sprintf("\n✓ Found %d mutation file(s)\n\n", length(mutation_files_found)))
}

# ------------------------------------------------------------------------------
# STEP 3: Load and Inspect Each Mutation File
# ------------------------------------------------------------------------------
cat("STEP 3: Loading and inspecting mutation files...\n")

mutation_data_list <- list()

for (mfile in names(mutation_files_found)) {
  cat(sprintf("\n--- Analyzing: %s ---\n", basename(mfile)))

  tryCatch({
    # Try reading as tab-delimited
    mut_data <- read.delim(mfile, header = TRUE, stringsAsFactors = FALSE, nrows = 10)

    cat(sprintf("✓ Dimensions (first 10 rows): %d rows × %d columns\n", nrow(mut_data), ncol(mut_data)))
    cat("✓ Column names:\n")
    print(colnames(mut_data))
    cat("\n")

    # Look for sample ID columns
    sample_id_cols <- colnames(mut_data)[grepl("sample|patient|barcode|lab|dbgap|specimen",
                                                 colnames(mut_data), ignore.case = TRUE)]

    if (length(sample_id_cols) > 0) {
      cat("✓ Potential sample ID columns:\n")
      for (col in sample_id_cols) {
        cat(sprintf("  - %s (first 5 values): %s\n",
                    col,
                    paste(head(mut_data[[col]], 5), collapse=", ")))
      }
    } else {
      cat("⚠ No obvious sample ID column found\n")
    }

    # Load full file
    mut_data_full <- read.delim(mfile, header = TRUE, stringsAsFactors = FALSE)
    mutation_data_list[[mfile]] <- mut_data_full

  }, error = function(e) {
    cat(sprintf("❌ Error reading file: %s\n", e$message))
  })
}

cat("\n")

# ------------------------------------------------------------------------------
# STEP 4: Attempt to Match Sample IDs
# ------------------------------------------------------------------------------
cat("STEP 4: Attempting to match sample IDs across datasets...\n\n")

for (mfile in names(mutation_data_list)) {
  mut_data <- mutation_data_list[[mfile]]
  cat(sprintf("--- %s ---\n", basename(mfile)))

  # Identify sample ID columns
  sample_id_cols <- colnames(mut_data)[grepl("sample|patient|barcode|lab|dbgap|specimen",
                                               colnames(mut_data), ignore.case = TRUE)]

  if (length(sample_id_cols) == 0) {
    cat("⚠ No sample ID column identified, skipping\n\n")
    next
  }

  for (col in sample_id_cols) {
    mut_samples <- unique(mut_data[[col]])
    cat(sprintf("\nTesting column: %s (%d unique samples)\n", col, length(mut_samples)))

    # Direct match
    direct_match <- intersect(expr_samples, mut_samples)
    cat(sprintf("  Direct match: %d samples\n", length(direct_match)))

    # Try common transformations
    # 1. Remove hyphens
    mut_samples_no_hyphen <- gsub("-", "_", mut_samples)
    match_no_hyphen <- intersect(expr_samples, mut_samples_no_hyphen)
    cat(sprintf("  Match (no hyphens): %d samples\n", length(match_no_hyphen)))

    # 2. Extract prefix (e.g., "11-00001" -> "11-00001")
    mut_samples_prefix <- sapply(strsplit(as.character(mut_samples), "_"), `[`, 1)
    match_prefix <- intersect(expr_samples, mut_samples_prefix)
    cat(sprintf("  Match (prefix only): %d samples\n", length(match_prefix)))

    # 3. Add "X" prefix (common in R when column names start with number)
    mut_samples_X <- paste0("X", mut_samples)
    match_X <- intersect(expr_samples, mut_samples_X)
    cat(sprintf("  Match (with X prefix): %d samples\n", length(match_X)))

    # 4. Check if expression IDs are dbgap_rnaseq_sample format
    # Try to load clinical data for mapping
    clinical_file <- "01_Data/BeatAML_Downloaded_Data/beataml_waves1to4_clinical.txt"
    if (file.exists(clinical_file)) {
      clinical <- read.delim(clinical_file, stringsAsFactors = FALSE)

      if ("dbgap_rnaseq_sample" %in% colnames(clinical) &&
          "labId" %in% colnames(clinical)) {

        # Create mapping
        id_map <- clinical %>%
          select(dbgap_rnaseq_sample, labId) %>%
          filter(!is.na(dbgap_rnaseq_sample) & !is.na(labId))

        # Try matching mutation samples to labId
        mut_in_map <- mut_samples[mut_samples %in% id_map$labId]
        if (length(mut_in_map) > 0) {
          mapped_expr_ids <- id_map %>%
            filter(labId %in% mut_in_map) %>%
            filter(dbgap_rnaseq_sample %in% expr_samples)

          cat(sprintf("  Match (via labId mapping): %d samples\n", nrow(mapped_expr_ids)))

          if (nrow(mapped_expr_ids) > 0) {
            cat("\n✓✓✓ SUCCESS: Found matching samples via labId mapping! ✓✓✓\n")
            cat(sprintf("Matched samples (first 10): %s\n",
                        paste(head(mapped_expr_ids$dbgap_rnaseq_sample, 10), collapse=", ")))

            # Save mapping
            output_file <- "03_Results/mutation_sample_id_mapping.csv"
            write.csv(mapped_expr_ids, output_file, row.names = FALSE)
            cat(sprintf("\n✓ Saved ID mapping to: %s\n", output_file))

            # Save matched sample list
            matched_samples <- data.frame(
              dbgap_rnaseq_sample = mapped_expr_ids$dbgap_rnaseq_sample,
              labId = mapped_expr_ids$labId
            )
            write.csv(matched_samples,
                      "03_Results/matched_mutation_samples.csv",
                      row.names = FALSE)
          }
        }
      }
    }
  }
  cat("\n")
}

# ------------------------------------------------------------------------------
# STEP 5: Summary Report
# ------------------------------------------------------------------------------
cat("═══════════════════════════════════════════════════════════════\n")
cat("  DIAGNOSTIC SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat(sprintf("Expression samples: %d\n", length(expr_samples)))
cat(sprintf("Expression ID format: %s\n",
            ifelse(grepl("^dbgap", expr_samples[1]), "dbgap_rnaseq_sample", "unknown")))

cat(sprintf("\nMutation files found: %d\n", length(mutation_data_list)))

if (file.exists("03_Results/mutation_sample_id_mapping.csv")) {
  mapping <- read.csv("03_Results/mutation_sample_id_mapping.csv")
  cat(sprintf("\n✓✓✓ SOLUTION FOUND ✓✓✓\n"))
  cat(sprintf("Successfully matched %d samples between mutation and expression data\n",
              nrow(mapping)))
  cat("\nNext steps:\n")
  cat("  1. Run TASK 1.2 to create mutation matrix\n")
  cat("  2. Run TASK 1.3 to test mutation enrichment by subtype\n")
} else {
  cat("\n❌ No matching samples found\n")
  cat("Possible reasons:\n")
  cat("  1. Mutation data uses different sample identifiers\n")
  cat("  2. No overlap between expression and mutation cohorts\n")
  cat("  3. Mutation files not in expected format\n")
  cat("\nRecommendations:\n")
  cat("  1. Check BeatAML documentation for sample ID format\n")
  cat("  2. Verify mutation files are present in 01_Data directory\n")
  cat("  3. Consider requesting mutation data if missing\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("  DIAGNOSTIC COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n")
