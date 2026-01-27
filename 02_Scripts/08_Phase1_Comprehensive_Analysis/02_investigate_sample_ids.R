#!/usr/bin/env Rscript
# Investigate sample ID formats to fix matching issue

library(data.table)
library(readxl)

setwd("D:/Projects/Project_AML")

# Load data
cat("Loading data files...\n")
expression_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
                        header = TRUE, data.table = FALSE, nrows = 5)
mutations_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
                       header = TRUE, data.table = FALSE)
clinical_raw <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
clinical_raw <- as.data.frame(clinical_raw)

# Expression IDs
expr_ids <- colnames(expression_raw)[-1]
cat("\n=== EXPRESSION SAMPLE IDs ===\n")
cat("First 10 IDs:\n")
cat(paste(head(expr_ids, 10), collapse = "\n"), "\n")
cat(sprintf("\nTotal: %d samples\n", length(expr_ids)))

# Clinical RNA-seq IDs
cat("\n=== CLINICAL RNA-SEQ IDs (dbgap_rnaseq_sample) ===\n")
clin_rna <- clinical_raw$dbgap_rnaseq_sample[!is.na(clinical_raw$dbgap_rnaseq_sample)]
cat("First 10 IDs:\n")
cat(paste(head(clin_rna, 10), collapse = "\n"), "\n")
cat(sprintf("\nTotal: %d samples\n", length(clin_rna)))

# Check matching
matches_rna <- sum(expr_ids %in% clin_rna)
cat(sprintf("\n✓ Expression IDs matching clinical RNA IDs: %d / %d (%.1f%%)\n",
            matches_rna, length(expr_ids), matches_rna/length(expr_ids)*100))

# Mutation IDs
cat("\n=== MUTATION SAMPLE IDs (dbgap_sample_id) ===\n")
mut_ids <- unique(mutations_raw$dbgap_sample_id)
cat("First 10 IDs:\n")
cat(paste(head(mut_ids, 10), collapse = "\n"), "\n")
cat(sprintf("\nTotal: %d samples\n", length(mut_ids)))

# Clinical DNA-seq IDs
cat("\n=== CLINICAL DNA-SEQ IDs (dbgap_dnaseq_sample) ===\n")
clin_dna <- clinical_raw$dbgap_dnaseq_sample[!is.na(clinical_raw$dbgap_dnaseq_sample)]
cat("First 10 IDs:\n")
cat(paste(head(clin_dna, 10), collapse = "\n"), "\n")
cat(sprintf("\nTotal: %d samples\n", length(clin_dna)))

# Check matching
matches_dna <- sum(mut_ids %in% clin_dna)
cat(sprintf("\n✓ Mutation IDs matching clinical DNA IDs: %d / %d (%.1f%%)\n",
            matches_dna, length(mut_ids), matches_dna/length(mut_ids)*100))

# Drug response - check column
cat("\n=== DRUG RESPONSE SAMPLE IDs ===\n")
drug_response_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt",
                           header = TRUE, data.table = FALSE, nrows = 10)
cat("Columns in drug response file:\n")
cat(paste(colnames(drug_response_raw), collapse = ", "), "\n")
cat("\nFirst 10 rows of first column:\n")
print(head(drug_response_raw[,1], 10))

# Clinical subject IDs
cat("\n=== CLINICAL SUBJECT IDs (dbgap_subject_id) ===\n")
clin_subj <- clinical_raw$dbgap_subject_id[!is.na(clinical_raw$dbgap_subject_id)]
cat("First 10 IDs:\n")
cat(paste(head(clin_subj, 10), collapse = "\n"), "\n")
cat(sprintf("\nTotal: %d samples\n", length(clin_subj)))

cat("\n=== SUMMARY ===\n")
cat(sprintf("Expression samples: %d\n", length(expr_ids)))
cat(sprintf("Mutation samples: %d\n", length(mut_ids)))
cat(sprintf("Clinical RNA samples: %d\n", length(clin_rna)))
cat(sprintf("Clinical DNA samples: %d\n", length(clin_dna)))
cat(sprintf("Clinical subject IDs: %d\n", length(clin_subj)))
cat(sprintf("\nExpression-Clinical RNA match rate: %.1f%%\n",
            matches_rna/length(expr_ids)*100))
cat(sprintf("Mutation-Clinical DNA match rate: %.1f%%\n",
            matches_dna/length(mut_ids)*100))
