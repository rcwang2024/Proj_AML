################################################################################
# Task 1.3: Create Sample Inventory and Overlap Analysis
################################################################################
# This script creates a comprehensive inventory of samples across all datasets
# and analyzes overlap for multi-omics integration.
#
# Outputs:
# - Sample inventory table
# - Overlap analysis
# - Integration-ready sample list
#
# Author: Data Analysis Pipeline
# Date: 2025-10-02
# Project: AML Multi-Omics Integration
################################################################################

library(readr)
library(readxl)
library(dplyr)
library(tibble)

#' Extract sample IDs from a dataframe
#' @param df Data frame
#' @param dataset_name Name of dataset
#' @return List with sample_ids, id_column, id_format
extract_sample_ids <- function(df, dataset_name) {

  if (dataset_name == "expression") {
    # Sample IDs are column names (skip first 4 annotation columns)
    sample_ids <- colnames(df)[5:ncol(df)]
    sample_ids <- sample_ids[grepl("^BA", sample_ids)]
    id_column <- "column_names"
    id_format <- "BA####R"

  } else if (dataset_name == "drug_auc") {
    sample_ids <- unique(df$dbgap_rnaseq_sample[!is.na(df$dbgap_rnaseq_sample)])
    id_column <- "dbgap_rnaseq_sample"
    id_format <- "BA####R"

  } else if (dataset_name == "clinical") {
    if ("dbgap_rnaseq_sample" %in% colnames(df)) {
      sample_ids <- unique(df$dbgap_rnaseq_sample[!is.na(df$dbgap_rnaseq_sample)])
      id_column <- "dbgap_rnaseq_sample"
      id_format <- "BA####R"
    } else if ("dbgap_subject_id" %in% colnames(df)) {
      sample_ids <- unique(df$dbgap_subject_id[!is.na(df$dbgap_subject_id)])
      id_column <- "dbgap_subject_id"
      id_format <- "numeric"
    } else {
      sample_ids <- character(0)
      id_column <- "unknown"
      id_format <- "unknown"
    }

  } else if (dataset_name == "mutations") {
    if ("dbgap_sample_id" %in% colnames(df)) {
      sample_ids <- unique(df$dbgap_sample_id[!is.na(df$dbgap_sample_id)])
      id_column <- "dbgap_sample_id"
      id_format <- "BA####D"
    } else {
      sample_ids <- character(0)
      id_column <- "unknown"
      id_format <- "unknown"
    }

  } else if (dataset_name == "raw_inhibitor") {
    sample_ids <- unique(df$dbgap_rnaseq_sample[!is.na(df$dbgap_rnaseq_sample)])
    id_column <- "dbgap_rnaseq_sample"
    id_format <- "BA####R"

  } else {
    sample_ids <- character(0)
    id_column <- "unknown"
    id_format <- "unknown"
  }

  return(list(
    sample_ids = sample_ids,
    id_column = id_column,
    id_format = id_format
  ))
}

#' Load all datasets
#' @param data_dir Path to data directory
#' @return List of datasets
load_datasets <- function(data_dir) {

  cat("Loading datasets...\n")
  cat(strrep("=", 80), "\n")

  datasets <- list()

  # 1. Expression data
  cat("\n1. Loading expression data...\n")
  tryCatch({
    expr <- read_delim(file.path(data_dir, "beataml_expression.txt"),
                      delim = "\t", show_col_types = FALSE)
    datasets$expression <- expr
    cat(sprintf("   ✓ Loaded: %d x %d\n", nrow(expr), ncol(expr)))
  }, error = function(e) {
    cat("   ✗ Error:", e$message, "\n")
    datasets$expression <- NULL
  })

  # 2. Drug AUC data
  cat("\n2. Loading drug AUC data...\n")
  tryCatch({
    drug_auc <- read_delim(file.path(data_dir, "beataml_drug_auc.txt"),
                          delim = "\t", show_col_types = FALSE)
    datasets$drug_auc <- drug_auc
    cat(sprintf("   ✓ Loaded: %d x %d\n", nrow(drug_auc), ncol(drug_auc)))
  }, error = function(e) {
    cat("   ✗ Error:", e$message, "\n")
    datasets$drug_auc <- NULL
  })

  # 3. Clinical data
  cat("\n3. Loading clinical data...\n")
  tryCatch({
    clinical <- read_excel(file.path(data_dir, "beataml_clinical.xlsx"))
    datasets$clinical <- clinical
    cat(sprintf("   ✓ Loaded: %d x %d\n", nrow(clinical), ncol(clinical)))
  }, error = function(e) {
    cat("   ✗ Error:", e$message, "\n")
    datasets$clinical <- NULL
  })

  # 4. Mutations data
  cat("\n4. Loading mutations data...\n")
  tryCatch({
    mutations <- read_delim(file.path(data_dir, "beataml_mutations.txt"),
                           delim = "\t", show_col_types = FALSE)
    datasets$mutations <- mutations
    cat(sprintf("   ✓ Loaded: %d x %d\n", nrow(mutations), ncol(mutations)))
  }, error = function(e) {
    cat("   ✗ Error:", e$message, "\n")
    datasets$mutations <- NULL
  })

  # 5. Raw inhibitor data
  cat("\n5. Loading raw inhibitor data...\n")
  tryCatch({
    raw_inhib <- read_delim(file.path(data_dir, "beataml_raw_inhibitor.txt"),
                           delim = "\t", show_col_types = FALSE)
    datasets$raw_inhibitor <- raw_inhib
    cat(sprintf("   ✓ Loaded: %d x %d\n", nrow(raw_inhib), ncol(raw_inhib)))
  }, error = function(e) {
    cat("   ✗ Error:", e$message, "\n")
    datasets$raw_inhibitor <- NULL
  })

  cat("\n", strrep("=", 80), "\n")
  return(datasets)
}

#' Create sample inventory
#' @param datasets List of datasets
#' @return List with inventory dataframe and sample_sets
create_sample_inventory <- function(datasets) {

  cat("\nCreating sample inventory...\n")
  cat(strrep("=", 80), "\n")

  inventory <- data.frame(
    data_type = character(),
    n_samples = integer(),
    sample_id_column = character(),
    id_format = character(),
    date_analyzed = character(),
    stringsAsFactors = FALSE
  )

  sample_sets <- list()

  for (dataset_name in names(datasets)) {
    df <- datasets[[dataset_name]]

    if (is.null(df)) {
      cat(sprintf("\n✗ Skipping %s (not loaded)\n", dataset_name))
      next
    }

    cat(sprintf("\nProcessing %s...\n", dataset_name))

    # Extract sample IDs
    result <- extract_sample_ids(df, dataset_name)
    sample_ids <- result$sample_ids
    id_column <- result$id_column
    id_format <- result$id_format

    n_samples <- length(sample_ids)

    cat(sprintf("  ID Column: %s\n", id_column))
    cat(sprintf("  ID Format: %s\n", id_format))
    cat(sprintf("  N Samples: %s\n", format(n_samples, big.mark = ",")))

    # Store sample set
    sample_sets[[dataset_name]] <- sample_ids

    # Add to inventory
    inventory <- rbind(inventory, data.frame(
      data_type = dataset_name,
      n_samples = n_samples,
      sample_id_column = id_column,
      id_format = id_format,
      date_analyzed = format(Sys.Date(), "%Y-%m-%d"),
      stringsAsFactors = FALSE
    ))
  }

  cat("\n", strrep("=", 80), "\n")

  return(list(
    inventory = inventory,
    sample_sets = sample_sets
  ))
}

#' Analyze overlap between datasets
#' @param sample_sets List of sample sets
#' @return Core samples
analyze_overlap <- function(sample_sets) {

  cat("\nAnalyzing sample overlap...\n")
  cat(strrep("=", 80), "\n")

  # Separate RNA and DNA samples
  rna_samples <- list()
  dna_samples <- list()

  for (name in names(sample_sets)) {
    samples <- sample_sets[[name]]
    if (length(samples) > 0) {
      if (is.character(samples) && grepl("R", samples[1])) {
        rna_samples[[name]] <- samples
      } else if (is.character(samples) && grepl("D", samples[1])) {
        dna_samples[[name]] <- samples
      }
    }
  }

  cat("\nRNA-seq based datasets:\n")
  for (name in names(rna_samples)) {
    cat(sprintf("  %s: %d samples\n", name, length(rna_samples[[name]])))
  }

  cat("\nDNA-seq based datasets:\n")
  for (name in names(dna_samples)) {
    cat(sprintf("  %s: %d samples\n", name, length(dna_samples[[name]])))
  }

  # Calculate overlaps for RNA samples
  if (length(rna_samples) >= 2) {
    cat("\n", strrep("-", 80), "\n")
    cat("RNA Sample Overlaps:\n")
    cat(strrep("-", 80), "\n")

    names_list <- names(rna_samples)
    for (i in 1:(length(names_list)-1)) {
      for (j in (i+1):length(names_list)) {
        name1 <- names_list[i]
        name2 <- names_list[j]

        overlap <- intersect(rna_samples[[name1]], rna_samples[[name2]])
        only1 <- setdiff(rna_samples[[name1]], rna_samples[[name2]])
        only2 <- setdiff(rna_samples[[name2]], rna_samples[[name1]])

        cat(sprintf("\n%s ∩ %s:\n", name1, name2))
        cat(sprintf("  Overlap: %d samples\n", length(overlap)))
        cat(sprintf("  Only in %s: %d samples\n", name1, length(only1)))
        cat(sprintf("  Only in %s: %d samples\n", name2, length(only2)))
      }
    }
  }

  # Find core multi-omics samples
  core_samples <- character(0)

  if ("expression" %in% names(rna_samples) && "drug_auc" %in% names(rna_samples)) {
    core_samples <- intersect(rna_samples$expression, rna_samples$drug_auc)

    if ("clinical" %in% names(rna_samples)) {
      core_samples <- intersect(core_samples, rna_samples$clinical)
    }

    cat("\n", strrep("=", 80), "\n")
    cat("CORE MULTI-OMICS SAMPLES\n")
    cat(strrep("=", 80), "\n")
    cat(sprintf("Samples with Expression + Drug Response: %d\n", length(core_samples)))
  }

  return(core_samples)
}

#' Create integration table
#' @param datasets List of datasets
#' @param core_samples Vector of core sample IDs
#' @return Integration dataframe
create_integration_table <- function(datasets, core_samples) {

  cat("\nCreating integration table...\n")
  cat(strrep("=", 80), "\n")

  integration_data <- data.frame(
    sample_id = core_samples,
    stringsAsFactors = FALSE
  )

  # Get sample lists
  expr_samples <- extract_sample_ids(datasets$expression, "expression")$sample_ids
  drug_samples <- extract_sample_ids(datasets$drug_auc, "drug_auc")$sample_ids

  integration_data$has_expression <- integration_data$sample_id %in% expr_samples
  integration_data$has_drug_response <- integration_data$sample_id %in% drug_samples

  if (!is.null(datasets$clinical)) {
    clin_samples <- extract_sample_ids(datasets$clinical, "clinical")$sample_ids
    integration_data$has_clinical <- integration_data$sample_id %in% clin_samples
  } else {
    integration_data$has_clinical <- FALSE
  }

  # Check for mutations (map RNA -> DNA)
  if (!is.null(datasets$mutations)) {
    mut_samples <- extract_sample_ids(datasets$mutations, "mutations")$sample_ids
    integration_data$dna_sample_id <- gsub("R$", "D", integration_data$sample_id)
    integration_data$has_mutations <- integration_data$dna_sample_id %in% mut_samples
  } else {
    integration_data$has_mutations <- FALSE
  }

  cat(sprintf("\nIntegration table created: %d samples\n", nrow(integration_data)))

  cat("\nData availability:\n")
  for (col in c("has_expression", "has_drug_response", "has_clinical", "has_mutations")) {
    if (col %in% colnames(integration_data)) {
      n <- sum(integration_data[[col]], na.rm = TRUE)
      pct <- (n / nrow(integration_data)) * 100
      cat(sprintf("  %s: %d (%.1f%%)\n",
                  gsub("has_", "", col), n, pct))
    }
  }

  # Count complete cases
  complete <- integration_data[integration_data$has_expression &
                               integration_data$has_drug_response &
                               integration_data$has_mutations, ]
  cat(sprintf("\nComplete multi-omics (Expression + Drug + Mutations): %d samples\n",
              nrow(complete)))

  return(integration_data)
}

#' Main function
main <- function() {

  # Set paths
  script_dir <- dirname(sys.frame(1)$ofile)
  project_root <- dirname(dirname(script_dir))
  data_dir <- file.path(project_root, "01_Data", "BeatAML_Downloaded_Data")
  output_dir <- file.path(project_root, "03_Results", "01_Processed_Data")

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat(strrep("=", 80), "\n")
  cat("BeatAML Sample Inventory and Overlap Analysis\n")
  cat(strrep("=", 80), "\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Data Directory:", data_dir, "\n")
  cat("Output Directory:", output_dir, "\n")
  cat(strrep("=", 80), "\n")

  # Load datasets
  datasets <- load_datasets(data_dir)

  # Create inventory
  result <- create_sample_inventory(datasets)
  inventory_df <- result$inventory
  sample_sets <- result$sample_sets

  # Display inventory
  cat("\n", strrep("=", 80), "\n")
  cat("SAMPLE INVENTORY\n")
  cat(strrep("=", 80), "\n")
  print(inventory_df, row.names = FALSE)

  # Save inventory
  inventory_file <- file.path(output_dir, "sample_inventory.csv")
  write.csv(inventory_df, inventory_file, row.names = FALSE)
  cat(sprintf("\n✓ Saved inventory to: %s\n", inventory_file))

  # Analyze overlap
  core_samples <- analyze_overlap(sample_sets)

  # Create integration table
  if (length(core_samples) > 0) {
    integration_df <- create_integration_table(datasets, core_samples)

    # Save integration table
    integration_file <- file.path(output_dir, "sample_integration_table.csv")
    write.csv(integration_df, integration_file, row.names = FALSE)
    cat(sprintf("\n✓ Saved integration table to: %s\n", integration_file))

    # Save core sample list
    core_file <- file.path(output_dir, "core_multi_omics_samples.txt")
    writeLines(c(
      "# Core Multi-Omics Samples (Expression + Drug Response)",
      sprintf("# Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf("# Total samples: %d", length(core_samples)),
      "",
      sort(core_samples)
    ), core_file)
    cat(sprintf("✓ Saved core sample list to: %s\n", core_file))
  }

  cat("\n", strrep("=", 80), "\n")
  cat("✓ Sample inventory analysis complete!\n")
  cat(strrep("=", 80), "\n")
}

# Run if not in interactive mode
if (!interactive()) {
  main()
} else {
  cat("Script loaded. Run main() to execute.\n")
}
