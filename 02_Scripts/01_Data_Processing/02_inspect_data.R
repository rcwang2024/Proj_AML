################################################################################
# Task 1.2: Inspect BeatAML Data Files in Detail
################################################################################
# This script provides comprehensive inspection of each downloaded BeatAML file:
# - File format and structure
# - Dimensions (rows x columns)
# - Column names and data types
# - First 10 rows preview
# - Basic statistics
# - Missing data analysis
#
# Author: Data Analysis Pipeline
# Date: 2025-10-02
# Project: AML Multi-Omics Integration
################################################################################

library(readr)
library(readxl)
library(dplyr)
library(tibble)

#' Format bytes to human readable
#' @param size Numeric bytes
#' @return Character string
format_bytes <- function(size) {
  units <- c("B", "KB", "MB", "GB", "TB")
  unit_idx <- 1

  while (size >= 1024 && unit_idx < length(units)) {
    size <- size / 1024
    unit_idx <- unit_idx + 1
  }

  return(sprintf("%.2f %s", size, units[unit_idx]))
}

#' Get basic statistics for a column
#' @param df Data frame
#' @param col Column name
#' @return List of statistics
get_basic_stats <- function(df, col) {
  stats <- list()

  tryCatch({
    if (is.numeric(df[[col]])) {
      stats$min <- min(df[[col]], na.rm = TRUE)
      stats$max <- max(df[[col]], na.rm = TRUE)
      stats$mean <- mean(df[[col]], na.rm = TRUE)
      stats$median <- median(df[[col]], na.rm = TRUE)
    } else {
      stats$unique <- length(unique(df[[col]]))
      top_val <- names(sort(table(df[[col]]), decreasing = TRUE)[1])
      if (length(top_val) > 0) {
        stats$top_value <- top_val
        stats$top_count <- sum(df[[col]] == top_val, na.rm = TRUE)
      }
    }
  }, error = function(e) {})

  return(stats)
}

#' Inspect a single data file
#' @param filepath Character path to file
#' @param filename Character filename
#' @param output_con File connection for output
inspect_file <- function(filepath, filename, output_con) {

  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("INSPECTING:", filename, "\n")
  cat(strrep("=", 80), "\n")

  writeLines(paste0("\n", strrep("=", 80)), output_con)
  writeLines(paste("File:", filename), output_con)
  writeLines(paste0(strrep("=", 80), "\n"), output_con)

  # File size
  file_size <- file.info(filepath)$size
  cat("File Size:", format_bytes(file_size), "\n")
  writeLines(paste("File Size:", format_bytes(file_size)), output_con)

  # Determine file type and read
  df <- NULL

  tryCatch({
    if (grepl("\\.xlsx$", filename)) {
      cat("File Format: Excel (.xlsx)\n")
      writeLines("File Format: Excel (.xlsx)", output_con)
      df <- read_excel(filepath)

    } else if (grepl("\\.txt$", filename)) {
      # Try to detect delimiter
      first_line <- readLines(filepath, n = 1)

      if (grepl("\t", first_line)) {
        delimiter <- "\t"
        format_name <- "Tab-delimited text"
      } else if (grepl(",", first_line)) {
        delimiter <- ","
        format_name <- "Comma-separated text"
      } else {
        delimiter <- "\t"  # default
        format_name <- "Text file"
      }

      cat("File Format:", format_name, "\n")
      writeLines(paste("File Format:", format_name), output_con)

      df <- read_delim(filepath, delim = delimiter, show_col_types = FALSE)

    } else {
      cat("File Format: Unknown\n")
      writeLines("File Format: Unknown", output_con)
      return(NULL)
    }

  }, error = function(e) {
    cat("✗ ERROR reading file:", e$message, "\n")
    writeLines(paste("ERROR reading file:", e$message), output_con)
    return(NULL)
  })

  if (is.null(df)) return(NULL)

  # Dimensions
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  cat(sprintf("Dimensions: %s rows × %s columns\n",
              format(n_rows, big.mark = ","),
              format(n_cols, big.mark = ",")))
  writeLines(sprintf("Dimensions: %s rows × %s columns\n",
                     format(n_rows, big.mark = ","),
                     format(n_cols, big.mark = ",")), output_con)

  # Column information
  cat("\nColumn Information:\n")
  cat(strrep("-", 80), "\n")
  writeLines("Column Information:", output_con)
  writeLines(strrep("-", 80), output_con)

  col_info <- data.frame(
    Column = character(),
    Type = character(),
    Missing = integer(),
    Missing_Pct = character(),
    stringsAsFactors = FALSE
  )

  for (col in colnames(df)) {
    dtype <- class(df[[col]])[1]
    n_missing <- sum(is.na(df[[col]]))
    pct_missing <- (n_missing / n_rows) * 100

    col_info <- rbind(col_info, data.frame(
      Column = col,
      Type = dtype,
      Missing = n_missing,
      Missing_Pct = sprintf("%.1f%%", pct_missing),
      stringsAsFactors = FALSE
    ))

    cat("  ", col, "\n", sep = "")
    cat("    Type:", dtype, "\n")
    cat(sprintf("    Missing: %s (%.1f%%)\n",
                format(n_missing, big.mark = ","), pct_missing))

    # Basic stats
    stats <- get_basic_stats(df, col)
    if (!is.null(stats$min)) {
      cat(sprintf("    Range: [%.2f, %.2f]\n", stats$min, stats$max))
      cat(sprintf("    Mean: %.2f, Median: %.2f\n", stats$mean, stats$median))
    } else if (!is.null(stats$unique)) {
      cat(sprintf("    Unique values: %s\n", format(stats$unique, big.mark = ",")))
      if (!is.null(stats$top_value)) {
        cat(sprintf("    Most common: '%s' (%s times)\n",
                    stats$top_value, format(stats$top_count, big.mark = ",")))
      }
    }
  }

  # Save column info
  writeLines(capture.output(print(col_info, row.names = FALSE)), output_con)
  writeLines("", output_con)

  # First 10 rows preview
  cat("\nFirst 10 Rows Preview:\n")
  cat(strrep("-", 80), "\n")
  preview <- head(df, 10)
  print(preview)

  writeLines("First 10 Rows Preview:", output_con)
  writeLines(strrep("-", 80), output_con)
  writeLines(capture.output(print(preview)), output_con)
  writeLines("", output_con)

  # Overall missing data
  total_cells <- n_rows * n_cols
  total_missing <- sum(is.na(df))
  pct_missing <- (total_missing / total_cells) * 100

  cat(sprintf("\nOverall Missing Data: %s / %s (%.2f%%)\n",
              format(total_missing, big.mark = ","),
              format(total_cells, big.mark = ","),
              pct_missing))
  writeLines(sprintf("Overall Missing Data: %s / %s (%.2f%%)",
                     format(total_missing, big.mark = ","),
                     format(total_cells, big.mark = ","),
                     pct_missing), output_con)

  # Data quality issues
  cat("\nData Quality Check:\n")
  writeLines("\nData Quality Check:", output_con)

  # Check for duplicate rows
  n_duplicates <- sum(duplicated(df))
  if (n_duplicates > 0) {
    msg <- sprintf("  ⚠ %s duplicate rows found", format(n_duplicates, big.mark = ","))
    cat(msg, "\n")
    writeLines(msg, output_con)
  } else {
    cat("  ✓ No duplicate rows\n")
    writeLines("  No duplicate rows", output_con)
  }

  # Check for columns with all missing data
  all_missing_cols <- colnames(df)[colSums(is.na(df)) == n_rows]
  if (length(all_missing_cols) > 0) {
    msg <- sprintf("  ⚠ %d columns with all missing data: %s",
                   length(all_missing_cols),
                   paste(all_missing_cols, collapse = ", "))
    cat(msg, "\n")
    writeLines(msg, output_con)
  } else {
    cat("  ✓ No columns with all missing data\n")
    writeLines("  No columns with all missing data", output_con)
  }

  # Check for columns with >50% missing
  high_missing_cols <- colnames(df)[colSums(is.na(df)) / n_rows > 0.5]
  if (length(high_missing_cols) > 0) {
    cat(sprintf("  ⚠ %d columns with >50%% missing data:\n", length(high_missing_cols)))
    writeLines(sprintf("  Columns with >50%% missing data:", length(high_missing_cols)),
               output_con)
    for (col in high_missing_cols) {
      pct <- (sum(is.na(df[[col]])) / n_rows) * 100
      msg <- sprintf("      - %s (%.1f%% missing)", col, pct)
      cat(msg, "\n")
      writeLines(paste0("    - ", col, sprintf(" (%.1f%% missing)", pct)), output_con)
    }
  } else {
    cat("  ✓ No columns with >50% missing data\n")
    writeLines("  No columns with >50% missing data", output_con)
  }

  cat("\n", strrep("=", 80), "\n", sep = "")
  writeLines("", output_con)
}

#' Main inspection function
main <- function() {

  # Set paths
  script_dir <- dirname(sys.frame(1)$ofile)
  project_root <- dirname(dirname(script_dir))
  data_dir <- file.path(project_root, "01_Data", "BeatAML_Downloaded_Data")
  output_dir <- file.path(project_root, "03_Results", "02_QC_Reports")

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_file_path <- file.path(output_dir, "data_inspection_summary.txt")

  # Files to inspect
  files_to_inspect <- c(
    "beataml_expression.txt",
    "beataml_drug_auc.txt",
    "beataml_clinical.xlsx",
    "beataml_mutations.txt",
    "beataml_raw_inhibitor.txt",
    "beataml_drug_families.xlsx"
  )

  cat(strrep("=", 80), "\n")
  cat("BeatAML Data Files - Detailed Inspection\n")
  cat(strrep("=", 80), "\n")
  cat("Inspection Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Output File:", output_file_path, "\n")
  cat(strrep("=", 80), "\n")

  # Open output file
  output_con <- file(output_file_path, "w")

  writeLines("BeatAML Data Files - Detailed Inspection", output_con)
  writeLines(strrep("=", 80), output_con)
  writeLines(paste("Inspection Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), output_con)
  writeLines(strrep("=", 80), output_con)

  for (filename in files_to_inspect) {
    filepath <- file.path(data_dir, filename)
    if (file.exists(filepath)) {
      inspect_file(filepath, filename, output_con)
    } else {
      cat("\n✗ File not found:", filename, "\n")
      writeLines(paste("\nFile not found:", filename), output_con)
    }
  }

  close(output_con)

  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("✓ Inspection complete! Report saved to:\n")
  cat(" ", output_file_path, "\n")
  cat(strrep("=", 80), "\n")
}

# Run if not in interactive mode
if (!interactive()) {
  main()
} else {
  cat("Script loaded. Run main() to execute.\n")
}
