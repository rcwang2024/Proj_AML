# TASK: Manual TCGA-LAML Data Preparation
# Bypass TC GABiolinks preparation issues by manually loading files

library(tidyverse)
library(TCGAbiolinks)
library(DESeq2)

cat("=== TCGA-LAML MANUAL DATA PREPARATION ===\n\n")

# ============================================================================
# 1. MANUALLY LOAD EXPRESSION FILES
# ============================================================================

cat("=== STEP 1: LOADING EXPRESSION FILES MANUALLY ===\n\n")

# Create output directory
dir.create("01_Data/TCGA_LAML", showWarnings = FALSE, recursive = TRUE)

# Find all count files
count_files <- list.files("GDCdata/TCGA-LAML/",
                          pattern = "*rna_seq.aug_2014.txt.gz$",
                          recursive = TRUE,
                          full.names = TRUE)

if (length(count_files) == 0) {
  # Try alternative pattern for STAR counts
  count_files <- list.files("GDCdata/TCGA-LAML/",
                            pattern = "*.tsv$",
                            recursive = TRUE,
                            full.names = TRUE)
}

cat("Found", length(count_files), "expression files\n")

if (length(count_files) == 0) {
  cat("ERROR: No expression files found\n")
  cat("Looking for files in GDCdata/TCGA-LAML/\n")
  system("find GDCdata/TCGA-LAML/ -type f | head -5")
  quit(save = "no", status = 1)
}

# Load first file to get gene names
cat("Reading first file to get gene names...\n")
first_file <- read.table(count_files[1], header = TRUE, sep = "\t", comment.char = "#", fill = TRUE)

# Remove N_* rows (unmapped, multimapping, etc.)
gene_rows <- !grepl("^N_", first_file[, 1])
first_file <- first_file[gene_rows, ]

gene_names <- first_file$gene_id
cat("Found", length(gene_names), "genes\n\n")

# Initialize matrix and sample IDs
expr_matrix <- matrix(0, nrow = length(gene_names), ncol = length(count_files))
rownames(expr_matrix) <- gene_names
sample_ids <- character(length(count_files))

# Load all files
cat("Loading all expression files...\n")
for (i in 1:length(count_files)) {
  if (i %% 20 == 0) cat("  Loaded", i, "of", length(count_files), "files\n")

  # Read file
  data <- read.table(count_files[i], header = TRUE, sep = "\t", comment.char = "#", fill = TRUE)

  # Remove N_* rows
  data <- data[!grepl("^N_", data[, 1]), ]

  # Extract sample ID from file path
  sample_ids[i] <- basename(dirname(count_files[i]))

  # Get unstranded counts (column 4: gene_id, gene_name, gene_type, unstranded)
  expr_matrix[, i] <- data$unstranded
}

# Set column names after loading all files
colnames(expr_matrix) <- sample_ids

cat("✓ Loaded all expression files\n")
cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n\n")

# ============================================================================
# 2. PROCESS GENE NAMES AND COUNTS
# ============================================================================

cat("=== STEP 2: PROCESSING GENE NAMES ===\n\n")

# Gene names are Ensembl IDs with version, extract gene symbols
# The first column has format: ENSG00000000003.14
gene_ids <- sapply(strsplit(as.character(gene_names), "\\."), `[`, 1)
rownames(expr_matrix) <- gene_ids

# Remove genes with no/low expression
keep_genes <- rowSums(expr_matrix) > 0
expr_matrix <- expr_matrix[keep_genes, ]
cat("After removing zero-expression genes:", nrow(expr_matrix), "genes\n\n")

# ============================================================================
# 3. NORMALIZE WITH DESEQ2
# ============================================================================

cat("=== STEP 3: NORMALIZING WITH DESEQ2 ===\n\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(expr_matrix),
  colData = data.frame(sample = colnames(expr_matrix)),
  design = ~ 1
)

# Size factor normalization
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

cat("Normalized counts:", nrow(normalized_counts), "genes x", ncol(normalized_counts), "samples\n")

# Log2 transform
expr_log2 <- log2(normalized_counts + 1)

cat("✓ Normalization complete\n\n")

# ============================================================================
# 4. DOWNLOAD CLINICAL DATA (Direct GDC API)
# ============================================================================

cat("=== STEP 4: DOWNLOADING CLINICAL DATA (bypassing TCGAbiolinks) ===\n\n")

# Download clinical data directly from GDC API
cat("Downloading clinical data from GDC API...\n")

library(httr)
library(jsonlite)

# Query GDC for TCGA-LAML clinical data
gdc_url <- "https://api.gdc.cancer.gov/cases"
fields <- c("submitter_id", "diagnoses.age_at_diagnosis", "demographic.gender",
            "demographic.vital_status", "diagnoses.days_to_death",
            "demographic.days_to_death",
            "diagnoses.days_to_last_follow_up")

response <- GET(
  url = gdc_url,
  query = list(
    filters = toJSON(list(
      op = "=",
      content = list(
        field = "project.project_id",
        value = "TCGA-LAML"
      )
    ), auto_unbox = TRUE),
    fields = paste(fields, collapse = ","),
    format = "JSON",
    size = 500
  ),
  add_headers(Accept = "application/json")
)

# Check response status
if (status_code(response) != 200) {
  cat("ERROR: GDC API returned status", status_code(response), "\n")
  cat("Response:", content(response, "text"), "\n")
  quit(save = "no", status = 1)
}

gdc_data <- content(response, "text", encoding = "UTF-8") %>% fromJSON()
cases <- gdc_data$data$hits

cat("Downloaded", nrow(cases), "cases from GDC\n")

# Process clinical data
tcga_clinical <- cases %>%
  mutate(
    sample_id = submitter_id,
    age = if(!is.null(diagnoses)) {
      sapply(diagnoses, function(x) {
        ifelse(length(x$age_at_diagnosis) > 0, x$age_at_diagnosis[1] / 365.25, NA)
      })
    } else NA,
    sex = demographic$gender,
    vital_status = demographic$vital_status,
    days_to_death_demo = demographic$days_to_death,
    days_to_death_diag = if(!is.null(diagnoses)) {
      sapply(diagnoses, function(x) {
        ifelse(length(x$days_to_death) > 0, x$days_to_death[1], NA)
      })
    } else NA,
    OS_event = ifelse(!is.na(vital_status) & vital_status %in% c("Dead", "dead"), 1, 0),
    OS_days = ifelse(!is.na(days_to_death_demo), days_to_death_demo,
                     ifelse(!is.na(days_to_death_diag), days_to_death_diag, NA)),
    OS_months = OS_days / 30.44
  ) %>%
  select(sample_id, age, sex, vital_status, OS_event, OS_days, OS_months)

cat("\nSurvival data:\n")
cat("  Deaths:", sum(tcga_clinical$OS_event == 1, na.rm = TRUE), "\n")
cat("  Censored:", sum(tcga_clinical$OS_event == 0, na.rm = TRUE), "\n\n")

# ============================================================================
# 5. MAP FILE UUIDs TO TCGA BARCODES
# ============================================================================

cat("=== STEP 5: MAPPING UUIDs TO TCGA BARCODES ===\n\n")

# Get file UUIDs from expression data
file_uuids <- colnames(expr_log2)
cat("Expression files:", length(file_uuids), "\n")

# Query GDC to map file UUIDs to sample barcodes
cat("Querying GDC for file-to-sample mapping...\n")
files_url <- "https://api.gdc.cancer.gov/files"

# Query in batches to avoid URL length limits
batch_size <- 100
uuid_to_barcode <- data.frame()

for (i in seq(1, length(file_uuids), by = batch_size)) {
  batch_uuids <- file_uuids[i:min(i+batch_size-1, length(file_uuids))]

  file_response <- GET(
    url = files_url,
    query = list(
      filters = toJSON(list(
        op = "in",
        content = list(
          field = "file_id",
          value = batch_uuids
        )
      ), auto_unbox = TRUE),
      fields = "file_id,cases.submitter_id",
      format = "JSON",
      size = batch_size
    ),
    add_headers(Accept = "application/json")
  )

  file_data <- content(file_response, "text", encoding = "UTF-8") %>% fromJSON()
  file_hits <- file_data$data$hits

  if (nrow(file_hits) > 0) {
    batch_mapping <- file_hits %>%
      mutate(
        file_uuid = file_id,
        tcga_barcode = sapply(cases, function(x) {
          if (length(x$submitter_id) > 0) x$submitter_id[1] else NA
        })
      ) %>%
      select(file_uuid, tcga_barcode)

    uuid_to_barcode <- rbind(uuid_to_barcode, batch_mapping)
  }
}

cat("Mapped", nrow(uuid_to_barcode), "files to TCGA barcodes\n\n")

# ============================================================================
# 6. MATCH SAMPLES
# ============================================================================

cat("=== STEP 6: MATCHING SAMPLES ===\n\n")

# Map expression UUIDs to TCGA barcodes
expr_mapping <- uuid_to_barcode %>%
  filter(file_uuid %in% file_uuids)

cat("Expression files with barcodes:", nrow(expr_mapping), "\n")

# Extract patient IDs (first 12 characters of TCGA barcode)
expr_mapping <- expr_mapping %>%
  mutate(patient_id = substr(tcga_barcode, 1, 12))

cat("Unique patients in expression:", length(unique(expr_mapping$patient_id)), "\n")

# Match with clinical
common_patients <- intersect(expr_mapping$patient_id, tcga_clinical$sample_id)
cat("Matched patients:", length(common_patients), "\n\n")

# For patients with multiple samples, take the first
matched_uuids <- sapply(common_patients, function(pid) {
  uuids <- expr_mapping$file_uuid[expr_mapping$patient_id == pid]
  uuids[1]
})

# Subset expression data
expr_matched <- expr_log2[, matched_uuids]
colnames(expr_matched) <- common_patients

# Subset clinical data
clinical_matched <- tcga_clinical %>%
  filter(sample_id %in% common_patients)

cat("Final matched data:\n")
cat("  Expression:", nrow(expr_matched), "genes x", ncol(expr_matched), "samples\n")
cat("  Clinical:", nrow(clinical_matched), "samples\n\n")

# ============================================================================
# 7. SAVE DATA
# ============================================================================

cat("=== STEP 7: SAVING DATA ===\n\n")

saveRDS(expr_matched, "01_Data/TCGA_LAML/tcga_laml_expression_normalized.rds")
cat("✓ Saved: tcga_laml_expression_normalized.rds\n")

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
# 8. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("TCGA-LAML DATA PREPARATION COMPLETE\n\n")

cat("Processed:\n")
cat("  Expression data:", ncol(expr_matched), "samples,", nrow(expr_matched), "genes\n")
cat("  Clinical data:", nrow(clinical_matched), "samples\n")
cat("  Deaths:", sum(clinical_matched$OS_event, na.rm = TRUE), "\n")
cat("  Median age:", round(median(clinical_matched$age, na.rm = TRUE), 1), "years\n\n")

cat("Next step: Apply BeatAML classifier\n")
cat("Run: Rscript 14_tcga_apply_classifier.R\n\n")

cat("### TCGA Preparation COMPLETE ###\n")
