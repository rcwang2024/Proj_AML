################################################################################
# Task 1.1: Verify Downloaded BeatAML Data Files
################################################################################
# This script verifies that all expected BeatAML data files are present,
# checks their integrity, and reports file sizes and basic properties.
#
# Author: Data Analysis Pipeline
# Date: 2025-10-02
# Project: AML Multi-Omics Integration
################################################################################

library(tools)  # for md5sum

# Expected files with their approximate sizes (in MB)
EXPECTED_FILES <- list(
  beataml_expression.txt = 269,
  beataml_drug_auc.txt = 19,
  beataml_clinical.xlsx = 0.5,
  beataml_mutations.txt = 3.5,
  beataml_raw_inhibitor.txt = 48,
  beataml_drug_families.xlsx = 0.1
)

#' Get file size in megabytes
#' @param filepath Character string of file path
#' @return Numeric size in MB or NA
get_file_size_mb <- function(filepath) {
  tryCatch({
    size_bytes <- file.info(filepath)$size
    size_mb <- size_bytes / (1024 * 1024)
    return(size_mb)
  }, error = function(e) {
    return(NA)
  })
}

#' Check if file can be opened and read
#' @param filepath Character string of file path
#' @return Logical TRUE if readable
check_file_readable <- function(filepath) {
  tryCatch({
    # Try to open and read first line
    con <- file(filepath, "rb")
    readBin(con, "raw", n = 1024)
    close(con)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

#' Calculate MD5 checksum for file integrity
#' @param filepath Character string of file path
#' @return Character MD5 hash or NA
calculate_md5 <- function(filepath) {
  tryCatch({
    md5 <- md5sum(filepath)
    return(as.character(md5))
  }, error = function(e) {
    return(NA)
  })
}

#' Main verification function
#' @param data_dir Character string of data directory path
#' @param output_log Character string of output log file path
#' @return Logical TRUE if all checks passed
verify_files <- function(data_dir, output_log) {

  cat(strrep("=", 80), "\n")
  cat("BeatAML Data Files Verification\n")
  cat(strrep("=", 80), "\n")
  cat("Data Directory:", data_dir, "\n")
  cat("Verification Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat(strrep("=", 80), "\n\n")

  results <- list()
  all_passed <- TRUE

  for (filename in names(EXPECTED_FILES)) {
    expected_size_mb <- EXPECTED_FILES[[filename]]
    filepath <- file.path(data_dir, filename)

    cat("Checking:", filename, "\n")
    cat(strrep("-", 60), "\n")

    # Check existence
    exists <- file.exists(filepath)
    cat("  Exists:", ifelse(exists, "✓ YES", "✗ NO"), "\n")

    if (!exists) {
      results[[filename]] <- list(
        filename = filename,
        exists = FALSE,
        size_mb = NA,
        readable = FALSE,
        status = "MISSING"
      )
      all_passed <- FALSE
      cat("  Status: MISSING\n\n")
      next
    }

    # Check size
    actual_size_mb <- get_file_size_mb(filepath)
    cat("  Expected Size: ~", expected_size_mb, "MB\n")
    cat("  Actual Size:", sprintf("%.2f MB", actual_size_mb), "\n")

    # Size tolerance check (±20%)
    size_ok <- (actual_size_mb >= expected_size_mb * 0.8 &
                actual_size_mb <= expected_size_mb * 1.2)
    cat("  Size Check:", ifelse(size_ok, "✓ PASS", "⚠ WARNING"), "\n")

    # Check readability
    readable <- check_file_readable(filepath)
    cat("  Readable:", ifelse(readable, "✓ YES", "✗ NO"), "\n")

    # Calculate checksum (for smaller files only)
    checksum <- NA
    if (actual_size_mb < 50) {  # Only for files < 50MB
      checksum <- calculate_md5(filepath)
      if (!is.na(checksum)) {
        cat("  MD5 Checksum:", substr(checksum, 1, 16), "...\n")
      } else {
        cat("  MD5 Checksum: Failed\n")
      }
    }

    # Overall status
    if (readable && size_ok) {
      status <- "OK"
      cat("  Status: ✓ OK\n")
    } else if (readable) {
      status <- "WARNING"
      cat("  Status: ⚠ WARNING (size mismatch)\n")
      all_passed <- FALSE
    } else {
      status <- "ERROR"
      cat("  Status: ✗ ERROR (not readable)\n")
      all_passed <- FALSE
    }

    results[[filename]] <- list(
      filename = filename,
      exists = exists,
      size_mb = actual_size_mb,
      readable = readable,
      checksum = checksum,
      status = status
    )

    cat("\n")
  }

  # Summary
  cat(strrep("=", 80), "\n")
  cat("VERIFICATION SUMMARY\n")
  cat(strrep("=", 80), "\n")

  ok_count <- sum(sapply(results, function(r) r$status == "OK"))
  warning_count <- sum(sapply(results, function(r) r$status == "WARNING"))
  error_count <- sum(sapply(results, function(r) r$status %in% c("ERROR", "MISSING")))

  cat("Total Files Checked:", length(EXPECTED_FILES), "\n")
  cat("  ✓ OK:", ok_count, "\n")
  cat("  ⚠ WARNING:", warning_count, "\n")
  cat("  ✗ ERROR/MISSING:", error_count, "\n\n")

  if (all_passed) {
    cat("Overall Status: ✓ ALL CHECKS PASSED\n")
  } else {
    cat("Overall Status: ⚠ SOME ISSUES DETECTED\n")
  }

  cat(strrep("=", 80), "\n")

  # Write log file
  log_con <- file(output_log, "w")

  writeLines("BeatAML Data Files Verification Log", log_con)
  writeLines(strrep("=", 80), log_con)
  writeLines(paste("Verification Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_con)
  writeLines(paste("Data Directory:", data_dir), log_con)
  writeLines(strrep("=", 80), log_con)
  writeLines("", log_con)

  for (result in results) {
    writeLines(paste("File:", result$filename), log_con)
    writeLines(paste("  Exists:", result$exists), log_con)
    writeLines(paste("  Size (MB):",
                    ifelse(is.na(result$size_mb), "N/A", sprintf("%.2f", result$size_mb))),
              log_con)
    writeLines(paste("  Readable:", result$readable), log_con)
    if (!is.na(result$checksum)) {
      writeLines(paste("  MD5:", result$checksum), log_con)
    }
    writeLines(paste("  Status:", result$status), log_con)
    writeLines("", log_con)
  }

  writeLines("\nSummary:", log_con)
  writeLines(paste("  Total Files:", length(EXPECTED_FILES)), log_con)
  writeLines(paste("  OK:", ok_count), log_con)
  writeLines(paste("  WARNING:", warning_count), log_con)
  writeLines(paste("  ERROR/MISSING:", error_count), log_con)
  writeLines(paste("  Overall:", ifelse(all_passed, "PASSED", "ISSUES DETECTED")), log_con)

  close(log_con)

  cat("\nLog file saved to:", output_log, "\n")

  return(all_passed)
}

# Main execution
if (!interactive()) {
  # Set paths
  script_dir <- dirname(sys.frame(1)$ofile)
  project_root <- dirname(dirname(script_dir))
  data_dir <- file.path(project_root, "01_Data", "BeatAML_Downloaded_Data")
  output_log <- file.path(project_root, "06_Documentation", "Data_Analysis_Log.txt")

  # Create output directory if needed
  dir.create(dirname(output_log), recursive = TRUE, showWarnings = FALSE)

  # Run verification
  success <- verify_files(data_dir, output_log)

  # Exit with appropriate code
  if (!success) {
    quit(status = 1)
  }
} else {
  # For interactive use
  cat("Script loaded. Run verify_files(data_dir, output_log) to execute.\n")
  cat("Example:\n")
  cat("  data_dir <- 'D:/Projects/Project_AML/01_Data/BeatAML_Downloaded_Data'\n")
  cat("  output_log <- 'D:/Projects/Project_AML/06_Documentation/Data_Analysis_Log.txt'\n")
  cat("  verify_files(data_dir, output_log)\n")
}
