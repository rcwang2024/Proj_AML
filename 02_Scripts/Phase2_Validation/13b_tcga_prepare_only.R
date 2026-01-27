# TASK: Prepare Already-Downloaded TCGA-LAML Data
# The download completed successfully, now we just need to prepare and normalize

library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)

cat("=== TCGA-LAML DATA PREPARATION (FROM DOWNLOADED FILES) ===\n\n")

# ============================================================================
# 1. LOAD QUERY (reconstruct without re-downloading)
# ============================================================================

cat("=== STEP 1: RECONSTRUCTING QUERY ===\n\n")

# Create output directory
dir.create("01_Data/TCGA_LAML", showWarnings = FALSE, recursive = TRUE)

# Query TCGA-LAML project (but don't download - files already exist)
cat("Querying TCGA-LAML project...\n")

query_tcga <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

cat("Found", length(query_tcga$results[[1]]$file_id), "files\n")
cat("Files already downloaded in GDCdata/ directory\n\n")

# ============================================================================
# 2. PREPARE EXPRESSION MATRIX
# ============================================================================

cat("=== STEP 2: PREPARING EXPRESSION MATRIX ===\n\n")

cat("Loading data from downloaded files...\n")
cat("This may take 5-10 minutes...\n")

tryCatch({
  # Prepare expression matrix from downloaded data
  tcga_data <- GDCprepare(
    query = query_tcga,
    save = TRUE,
    save.filename = "01_Data/TCGA_LAML/TCGA_LAML_RNAseq.rda"
  )

  cat("✓ Data preparation complete\n\n")

}, error = function(e) {
  cat("✗ Data preparation failed:", e$message, "\n")
  cat("Error details:\n")
  print(e)
  quit(save = "no", status = 1)
})

# ============================================================================
# 3. EXTRACT AND NORMALIZE EXPRESSION DATA
# ============================================================================

cat("=== STEP 3: EXTRACTING AND NORMALIZING EXPRESSION ===\n\n")

# Extract counts
counts <- assay(tcga_data, "unstranded")
cat("Raw counts matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Get gene symbols
genes <- rowData(tcga_data)
rownames(counts) <- genes$gene_name

# Remove genes without symbols
counts <- counts[!is.na(rownames(counts)) & rownames(counts) != "", ]
cat("After removing genes without symbols:", nrow(counts), "genes\n")

# Aggregate duplicate gene symbols (take mean)
cat("Aggregating duplicate gene symbols...\n")
counts_aggregated <- aggregate(counts,
                               by = list(gene = rownames(counts)),
                               FUN = mean)
rownames(counts_aggregated) <- counts_aggregated$gene
counts_aggregated$gene <- NULL

cat("After aggregation:", nrow(counts_aggregated), "genes\n\n")

# Normalize with DESeq2
cat("Normalizing with DESeq2...\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_aggregated),
  colData = data.frame(sample = colnames(counts_aggregated)),
  design = ~ 1
)

# Size factor normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

cat("Normalized counts:", nrow(normalized_counts), "genes x", ncol(normalized_counts), "samples\n")

# Log2 transform (log2(x + 1))
expr_log2 <- log2(normalized_counts + 1)

cat("✓ Normalization complete\n\n")

# ============================================================================
# 4. DOWNLOAD AND PREPARE CLINICAL DATA
# ============================================================================

cat("=== STEP 4: DOWNLOADING CLINICAL DATA ===\n\n")

# Query clinical data
cat("Querying clinical data...\n")
clinical_query <- GDCquery_clinic(project = "TCGA-LAML", type = "clinical")

cat("Clinical data:", nrow(clinical_query), "samples\n")
cat("Clinical variables:", ncol(clinical_query), "\n")

# Extract key variables
tcga_clinical <- clinical_query %>%
  select(
    sample_id = submitter_id,
    age = age_at_diagnosis,
    sex = gender,
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    race,
    ethnicity
  ) %>%
  mutate(
    # Convert age from days to years
    age = age / 365.25,

    # Create survival variables
    OS_event = ifelse(vital_status == "Dead", 1, 0),
    OS_days = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS_months = OS_days / 30.44  # Average days per month
  )

cat("\nSurvival data:\n")
cat("  Samples with survival:", sum(!is.na(tcga_clinical$OS_days)), "\n")
cat("  Deaths:", sum(tcga_clinical$OS_event == 1, na.rm = TRUE), "\n")
cat("  Censored:", sum(tcga_clinical$OS_event == 0, na.rm = TRUE), "\n\n")

# ============================================================================
# 5. MATCH SAMPLES
# ============================================================================

cat("=== STEP 5: MATCHING EXPRESSION AND CLINICAL DATA ===\n\n")

# TCGA barcodes: TCGA-XX-XXXX-XXA-XXX
# Clinical uses: TCGA-XX-XXXX
# Expression uses full barcode

# Extract patient IDs from expression barcodes
expr_barcodes <- colnames(expr_log2)
patient_ids <- substr(expr_barcodes, 1, 12)  # First 12 characters = patient ID

cat("Expression barcodes:", length(expr_barcodes), "\n")
cat("Unique patients:", length(unique(patient_ids)), "\n")

# Match with clinical
common_patients <- intersect(patient_ids, tcga_clinical$sample_id)
cat("Matched patients:", length(common_patients), "\n\n")

# For patients with multiple samples, take the first (primary tumor)
matched_samples <- sapply(common_patients, function(pid) {
  samples <- expr_barcodes[patient_ids == pid]
  samples[1]  # Take first sample
})

# Subset expression data
expr_matched <- expr_log2[, matched_samples]
colnames(expr_matched) <- names(matched_samples)  # Use patient IDs

# Subset clinical data
clinical_matched <- tcga_clinical %>%
  filter(sample_id %in% common_patients)

cat("Final matched data:\n")
cat("  Expression:", nrow(expr_matched), "genes x", ncol(expr_matched), "samples\n")
cat("  Clinical:", nrow(clinical_matched), "samples\n\n")

# ============================================================================
# 6. SAVE PROCESSED DATA
# ============================================================================

cat("=== STEP 6: SAVING PROCESSED DATA ===\n\n")

# Save expression data
saveRDS(expr_matched, "01_Data/TCGA_LAML/tcga_laml_expression_normalized.rds")
cat("✓ Saved: tcga_laml_expression_normalized.rds\n")

# Save clinical data
write.csv(clinical_matched,
          "01_Data/TCGA_LAML/tcga_laml_clinical.csv",
          row.names = FALSE)
cat("✓ Saved: tcga_laml_clinical.csv\n\n")

# Create summary
summary_data <- data.frame(
  metric = c("Total samples", "Genes", "Deaths", "Median age",
             "Male %", "Median follow-up (months)"),
  value = c(
    ncol(expr_matched),
    nrow(expr_matched),
    sum(clinical_matched$OS_event, na.rm = TRUE),
    round(median(clinical_matched$age, na.rm = TRUE), 1),
    round(sum(clinical_matched$sex == "male", na.rm = TRUE) / nrow(clinical_matched) * 100, 1),
    round(median(clinical_matched$OS_months, na.rm = TRUE), 1)
  )
)

write.csv(summary_data, "01_Data/TCGA_LAML/tcga_laml_summary.csv", row.names = FALSE)
cat("✓ Saved: tcga_laml_summary.csv\n\n")

# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("TCGA-LAML DATA PREPARATION COMPLETE\n\n")

cat("Downloaded and processed:\n")
cat("  Expression data:", ncol(expr_matched), "samples,", nrow(expr_matched), "genes\n")
cat("  Clinical data:", nrow(clinical_matched), "samples\n")
cat("  Deaths:", sum(clinical_matched$OS_event, na.rm = TRUE), "\n")
cat("  Median age:", round(median(clinical_matched$age, na.rm = TRUE), 1), "years\n")
cat("  Male:", sum(clinical_matched$sex == "male", na.rm = TRUE), "\n")
cat("  Female:", sum(clinical_matched$sex == "female", na.rm = TRUE), "\n\n")

cat("Files saved in: 01_Data/TCGA_LAML/\n")
cat("  - tcga_laml_expression_normalized.rds\n")
cat("  - tcga_laml_clinical.csv\n")
cat("  - tcga_laml_summary.csv\n\n")

cat("Next step: Apply BeatAML classifier to TCGA data\n")
cat("Run script: 14_tcga_apply_classifier.R\n\n")

cat("### TCGA Preparation Task COMPLETE ###\n")
