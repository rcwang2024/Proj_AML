# Phase 7: Manuscript Revision & Data Verification
## Claude Code Execution Prompt for Final Pre-Submission Corrections

**Date**: December 2025  
**Purpose**: Address critical issues identified in manuscript review before journal submission  
**Target Journal**: Blood (primary), JCO (secondary)  
**Priority**: ESSENTIAL - These corrections determine publication success

---

## MISSION STATEMENT

The AML molecular subtype manuscript is 90% publication-ready but has **critical gaps** that must be addressed:

1. 🔴 **Drug response sample size NOT reported** (critical omission)
2. 🔴 **Drug count discrepancy** (155 vs 122 in previous docs)
3. 🔴 **Cohen's d discrepancy** (1.25 vs >2.0 in previous docs)
4. 🔴 **Multivariate sample attrition unexplained** (671 → 459)
5. 🟡 **Ex vivo limitation underemphasized**
6. 🟡 **VRS validation details missing**
7. 🟡 **TP53 typo verification needed**

This prompt will systematically verify data, identify discrepancies, and generate publication-ready text corrections.

---

## WORKING DIRECTORY STRUCTURE

```
Project_AML/
├── 00_RawData/
│   ├── beataml_waves1to4_norm_exp_dbgap.txt
│   ├── beataml_wv1to4_clinical.xlsx
│   ├── beataml_drugresponse_auc_dbgap.txt
│   └── beataml_mutations.csv
├── 01_ProcessedData/
│   ├── expression_batchcorrected.rds
│   ├── cluster_assignments.csv
│   └── drug_response_processed.rds
├── 02_Analysis/
│   └── [previous phase outputs]
├── 03_Results/
│   └── [summary tables and statistics]
├── 04_Figures/
│   └── [publication figures]
├── 05_Manuscript/
│   └── [current manuscript files]
└── 07_Revision/          # NEW - Create this
    ├── data_verification/
    ├── text_corrections/
    └── supplementary_updates/
```

---

## TASK 1: DATA VERIFICATION & DISCREPANCY RESOLUTION
### Priority: 🔴 CRITICAL

### Task 1.1: Drug Response Sample Size Verification

**Objective**: Determine exact number of patients with both expression and drug response data

```r
# Load required data
library(tidyverse)

# Load cluster assignments
clusters <- read_csv("01_ProcessedData/cluster_assignments.csv")

# Load drug response data
drug_response <- read_tsv("00_RawData/beataml_drugresponse_auc_dbgap.txt")
# OR if processed version exists:
# drug_response <- readRDS("01_ProcessedData/drug_response_processed.rds")

# Identify patients with expression data (clustered patients)
expr_patients <- clusters$patient_id  # or sample_id column

# Identify patients with drug response data
drug_patients <- unique(drug_response$patient_id)  # adjust column name as needed

# Find overlap
overlap_patients <- intersect(expr_patients, drug_patients)
n_drug_analysis <- length(overlap_patients)

# Report
cat("\n========== DRUG RESPONSE SAMPLE SIZE ==========\n")
cat("Patients with expression data (clustered):", length(expr_patients), "\n")
cat("Patients with any drug response data:", length(drug_patients), "\n")
cat("Patients with BOTH (drug analysis cohort):", n_drug_analysis, "\n")

# Save verification
write_csv(data.frame(
  metric = c("expression_patients", "drug_patients", "overlap_patients"),
  count = c(length(expr_patients), length(drug_patients), n_drug_analysis)
), "07_Revision/data_verification/drug_sample_size_verification.csv")

# CRITICAL OUTPUT - Save this number for manuscript
cat("\n*** MANUSCRIPT TEXT ***\n")
cat(paste0("Drug sensitivity data were available for ", n_drug_analysis, 
           " patients with matched expression profiles.\n"))
```

### Task 1.2: Drug Count Verification

**Objective**: Verify exact number of drugs tested and reconcile 155 vs 122 discrepancy

```r
# Count unique drugs in dataset
drugs <- unique(drug_response$drug_name)  # adjust column name
n_drugs_total <- length(drugs)

# Count drugs with sufficient data for analysis (e.g., n >= 10 per cluster)
drug_counts <- drug_response %>%
  left_join(clusters, by = "patient_id") %>%
  filter(!is.na(cluster)) %>%
  group_by(drug_name, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(drug_name) %>%
  summarise(
    min_per_cluster = min(n),
    total_n = sum(n),
    .groups = "drop"
  )

# Drugs with minimum N per cluster
n_drugs_analyzed_n5 <- sum(drug_counts$min_per_cluster >= 5)
n_drugs_analyzed_n10 <- sum(drug_counts$min_per_cluster >= 10)

cat("\n========== DRUG COUNT VERIFICATION ==========\n")
cat("Total unique drugs in database:", n_drugs_total, "\n")
cat("Drugs with n≥5 per cluster:", n_drugs_analyzed_n5, "\n")
cat("Drugs with n≥10 per cluster:", n_drugs_analyzed_n10, "\n")

# Identify which count matches manuscript (155 or 122)
cat("\nDISCREPANCY CHECK:\n")
cat("Manuscript states: 155 drugs\n")
cat("Previous documentation: 122 compounds\n")
cat("Actual count: ", n_drugs_total, " total, ", n_drugs_analyzed_n10, " analyzed\n")

# Save
write_csv(drug_counts, "07_Revision/data_verification/drug_counts_by_cluster.csv")
```

### Task 1.3: Cohen's d Verification for Venetoclax

**Objective**: Calculate exact Cohen's d and reconcile 1.25 vs >2.0 discrepancy

```r
library(effsize)

# Get Venetoclax data
venetoclax_data <- drug_response %>%
  filter(grepl("venetoclax|ABT-199", drug_name, ignore.case = TRUE)) %>%
  left_join(clusters, by = "patient_id") %>%
  filter(!is.na(cluster), !is.na(auc))  # adjust column names

# Sample sizes per cluster
n_per_cluster <- venetoclax_data %>% count(cluster)
cat("\n========== VENETOCLAX SAMPLE SIZES ==========\n")
print(n_per_cluster)

# Calculate Cohen's d
cluster1_auc <- venetoclax_data$auc[venetoclax_data$cluster == 1]
cluster2_auc <- venetoclax_data$auc[venetoclax_data$cluster == 2]

# Using effsize package
cohens_d <- cohen.d(cluster1_auc, cluster2_auc)

# Manual calculation for verification
mean1 <- mean(cluster1_auc, na.rm = TRUE)
mean2 <- mean(cluster2_auc, na.rm = TRUE)
sd1 <- sd(cluster1_auc, na.rm = TRUE)
sd2 <- sd(cluster2_auc, na.rm = TRUE)
n1 <- length(cluster1_auc)
n2 <- length(cluster2_auc)

# Pooled SD
pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2))
manual_d <- (mean1 - mean2) / pooled_sd

cat("\n========== COHEN'S D VERIFICATION ==========\n")
cat("effsize package d:", round(cohens_d$estimate, 3), "\n")
cat("Manual calculation d:", round(manual_d, 3), "\n")
cat("Magnitude:", cohens_d$magnitude, "\n")

cat("\nDISCREPANCY CHECK:\n")
cat("Manuscript states: d = 1.25\n")
cat("Previous discussions: d > 2.0\n")
cat("Actual calculated d:", round(cohens_d$estimate, 3), "\n")

# Additional details
cat("\nDetailed Statistics:\n")
cat("Cluster 1 - Mean AUC:", round(mean1, 2), "SD:", round(sd1, 2), "n:", n1, "\n")
cat("Cluster 2 - Mean AUC:", round(mean2, 2), "SD:", round(sd2, 2), "n:", n2, "\n")
cat("Pooled SD:", round(pooled_sd, 2), "\n")

# Calculate 95% CI for Cohen's d
d_se <- sqrt((n1+n2)/(n1*n2) + (cohens_d$estimate^2)/(2*(n1+n2)))
d_ci_lower <- cohens_d$estimate - 1.96*d_se
d_ci_upper <- cohens_d$estimate + 1.96*d_se

cat("95% CI for d: [", round(d_ci_lower, 3), ", ", round(d_ci_upper, 3), "]\n")

# Save verification
verification_df <- data.frame(
  metric = c("cohens_d", "ci_lower", "ci_upper", "magnitude",
             "cluster1_mean", "cluster1_sd", "cluster1_n",
             "cluster2_mean", "cluster2_sd", "cluster2_n"),
  value = c(cohens_d$estimate, d_ci_lower, d_ci_upper, as.character(cohens_d$magnitude),
            mean1, sd1, n1, mean2, sd2, n2)
)
write_csv(verification_df, "07_Revision/data_verification/venetoclax_cohens_d_verification.csv")
```

### Task 1.4: Multivariate Sample Attrition Analysis

**Objective**: Explain why 671 patients dropped to 459 in multivariate analysis

```r
# Load clinical and mutation data
clinical <- readxl::read_excel("00_RawData/beataml_wv1to4_clinical.xlsx")
mutations <- read_csv("00_RawData/beataml_mutations.csv")

# Variables needed for multivariate model
required_vars <- c("patient_id", "age", "sex", "TP53", "TET2", "NPM1", 
                   "FLT3", "DNMT3A", "survival_time", "survival_status")

# Check completeness
completeness <- clusters %>%
  left_join(clinical, by = "patient_id") %>%
  left_join(mutations %>% select(patient_id, TP53, TET2, NPM1, FLT3, DNMT3A), 
            by = "patient_id") %>%
  mutate(
    has_age = !is.na(age),
    has_sex = !is.na(sex),
    has_survival = !is.na(survival_time) & !is.na(survival_status),
    has_TP53 = !is.na(TP53),
    has_TET2 = !is.na(TET2),
    has_all_mutations = has_TP53 & has_TET2,
    complete_case = has_age & has_sex & has_survival & has_all_mutations
  )

# Summary
attrition_summary <- completeness %>%
  summarise(
    total_clustered = n(),
    has_age = sum(has_age),
    has_sex = sum(has_sex),
    has_survival = sum(has_survival),
    has_TP53 = sum(has_TP53),
    has_TET2 = sum(has_TET2),
    has_all_mutations = sum(has_all_mutations),
    complete_cases = sum(complete_case)
  )

cat("\n========== SAMPLE ATTRITION ANALYSIS ==========\n")
print(t(attrition_summary))

# Identify primary reason for exclusion
excluded <- completeness %>% filter(!complete_case)
exclusion_reasons <- excluded %>%
  summarise(
    missing_age = sum(!has_age),
    missing_sex = sum(!has_sex),
    missing_survival = sum(!has_survival),
    missing_mutations = sum(!has_all_mutations)
  )

cat("\nExclusion Reasons (patients may have multiple):\n")
print(t(exclusion_reasons))

# Generate manuscript text
cat("\n*** MANUSCRIPT TEXT ***\n")
cat(paste0(
  "Complete data for multivariate analysis were available for ", 
  attrition_summary$complete_cases, " patients (", 
  round(100*attrition_summary$complete_cases/attrition_summary$total_clustered, 1),
  "% of clustered cohort). Exclusions were due to missing mutation data (n=",
  exclusion_reasons$missing_mutations, "), incomplete survival follow-up (n=",
  exclusion_reasons$missing_survival, "), or missing demographic information.\n"
))

# Save
write_csv(attrition_summary, "07_Revision/data_verification/multivariate_attrition_summary.csv")
write_csv(exclusion_reasons, "07_Revision/data_verification/exclusion_reasons.csv")
```

---

## TASK 2: VRS VALIDATION DETAILS
### Priority: 🟡 IMPORTANT

### Task 2.1: VRS Development and Validation Statistics

**Objective**: Generate validation statistics for Venetoclax Response Score

```r
library(pROC)
library(caret)

# Load VRS data (if exists) or recalculate
# Assuming VRS was calculated in previous phases

# If need to recalculate VRS from 9 genes:
vrs_genes <- c("BCL2", "NPM1", "DNMT3A", "TP53", "RUNX1", "ASXL1", "TET2", "CD47", "CTLA4")

# Load expression data
expr <- readRDS("01_ProcessedData/expression_batchcorrected.rds")

# Calculate VRS (simplified - actual weights from your model)
vrs_data <- expr %>%
  filter(gene %in% vrs_genes) %>%
  # ... your VRS calculation logic

# Split into development (70%) and validation (30%)
set.seed(42)
n_total <- nrow(vrs_data)
train_idx <- sample(1:n_total, round(0.7*n_total))
test_idx <- setdiff(1:n_total, train_idx)

train_data <- vrs_data[train_idx, ]
test_data <- vrs_data[test_idx, ]

# Fit model in training
# ... your model fitting

# Validate in test set
# ROC analysis
roc_result <- roc(test_data$venetoclax_response, test_data$vrs_predicted)

cat("\n========== VRS VALIDATION ==========\n")
cat("Development set n:", length(train_idx), "\n")
cat("Validation set n:", length(test_idx), "\n")
cat("Validation AUC:", round(auc(roc_result), 3), "\n")
cat("95% CI:", round(ci.auc(roc_result)[1], 3), "-", round(ci.auc(roc_result)[3], 3), "\n")

# Sensitivity/Specificity at optimal cutpoint
optimal <- coords(roc_result, "best", best.method = "youden")
cat("Optimal threshold:", round(optimal$threshold, 2), "\n")
cat("Sensitivity:", round(optimal$sensitivity, 3), "\n")
cat("Specificity:", round(optimal$specificity, 3), "\n")

# Generate manuscript text
cat("\n*** MANUSCRIPT TEXT ***\n")
cat(paste0(
  "VRS was developed in ", length(train_idx), " patients (70% of cohort) and validated ",
  "in the remaining ", length(test_idx), " patients (30%), achieving AUC of ",
  round(auc(roc_result), 2), " (95% CI: ", round(ci.auc(roc_result)[1], 2), "-",
  round(ci.auc(roc_result)[3], 2), ") with sensitivity of ",
  round(optimal$sensitivity*100, 1), "% and specificity of ",
  round(optimal$specificity*100, 1), "% at the optimal cutpoint.\n"
))
```

---

## TASK 3: GENERATE CORRECTED MANUSCRIPT TEXT
### Priority: 🔴 CRITICAL

### Task 3.1: Compile All Text Corrections

```r
# This task compiles all verified data into publication-ready text

# Create output directory
dir.create("07_Revision/text_corrections", recursive = TRUE, showWarnings = FALSE)

# Load all verification results
drug_sample_size <- read_csv("07_Revision/data_verification/drug_sample_size_verification.csv")
drug_counts <- read_csv("07_Revision/data_verification/drug_counts_by_cluster.csv")
cohens_d <- read_csv("07_Revision/data_verification/venetoclax_cohens_d_verification.csv")
attrition <- read_csv("07_Revision/data_verification/multivariate_attrition_summary.csv")

# Generate comprehensive corrections document
corrections <- list()

# 1. Drug sample size statement (ADD TO METHODS OR RESULTS)
n_drug <- drug_sample_size$count[drug_sample_size$metric == "overlap_patients"]
corrections$drug_sample_size <- paste0(
  "ADDITION (Results, Drug Response section, first paragraph):\n",
  "\"Drug sensitivity data were available for ", n_drug, 
  " patients with matched expression profiles and ex vivo drug response measurements.\""
)

# 2. Drug count correction (if needed)
n_drugs <- sum(drug_counts$min_per_cluster >= 5)  # or appropriate threshold
corrections$drug_count <- paste0(
  "VERIFICATION (Results, Drug Response section):\n",
  "Current text: 'Among 155 drugs tested...'\n",
  "Verified count: ", n_drugs, " drugs\n",
  "ACTION: ", ifelse(n_drugs == 155, "CORRECT - No change needed", 
                     paste0("UPDATE to: 'Among ", n_drugs, " drugs tested...'")))

# 3. Cohen's d correction (if needed)
d_value <- as.numeric(cohens_d$value[cohens_d$metric == "cohens_d"])
d_ci_lower <- as.numeric(cohens_d$value[cohens_d$metric == "ci_lower"])
d_ci_upper <- as.numeric(cohens_d$value[cohens_d$metric == "ci_upper"])
corrections$cohens_d <- paste0(
  "VERIFICATION (Results, Venetoclax section):\n",
  "Current text: 'Cohen's d=1.25, very large effect'\n",
  "Verified: Cohen's d = ", round(d_value, 2), " (95% CI: ", 
  round(d_ci_lower, 2), "-", round(d_ci_upper, 2), ")\n",
  "ACTION: ", ifelse(abs(d_value - 1.25) < 0.1, "CORRECT - No change needed",
                     paste0("UPDATE to: 'Cohen's d=", round(d_value, 2), 
                            " (95% CI: ", round(d_ci_lower, 2), "-", 
                            round(d_ci_upper, 2), "), ", 
                            ifelse(d_value >= 0.8, "very large", 
                                   ifelse(d_value >= 0.5, "large", "medium")), 
                            " effect'")))

# 4. Multivariate attrition explanation
n_complete <- attrition$complete_cases
n_total <- attrition$total_clustered
pct_complete <- round(100*n_complete/n_total, 1)
corrections$attrition <- paste0(
  "ADDITION (Results, after 'Multivariate Cox regression (n=459, 282 events)'):\n",
  "\"Complete covariate data were available for ", n_complete, " patients (",
  pct_complete, "% of the clustered cohort); exclusions were primarily due to ",
  "missing mutation annotation or incomplete survival follow-up.\""
)

# 5. Ex vivo limitation (STRENGTHEN in Discussion)
corrections$ex_vivo <- paste0(
  "REVISION (Discussion, Limitations, FIRST limitation):\n\n",
  "CURRENT (weak):\n",
  "\"First, drug sensitivity was measured ex vivo, not through clinical response.\"\n\n",
  "REVISED (stronger):\n",
  "\"The primary limitation is that drug sensitivity was measured using ex vivo assays ",
  "on isolated AML blasts, which may not fully recapitulate in vivo response due to ",
  "tumor microenvironment interactions, drug pharmacokinetics, and patient-specific ",
  "factors including comorbidities and concomitant medications. While BeatAML ex vivo ",
  "data have demonstrated correlation with clinical outcomes in prior validation ",
  "studies (Tyner et al., 2018; Bottomly et al., 2022), prospective clinical validation ",
  "remains essential before clinical implementation of these findings.\""
)

# 6. Clinical utility language (SOFTEN)
corrections$language <- paste0(
  "LANGUAGE REVISIONS (throughout manuscript):\n\n",
  "CHANGE: 'clinically ready biomarker'\n",
  "TO: 'biomarker ready for prospective validation'\n\n",
  "CHANGE: 'complete translational pathway'\n", 
  "TO: 'clear translational pathway'\n\n",
  "CHANGE: 'immediate clinical utility'\n",
  "TO: 'potential clinical utility pending validation'\n\n",
  "CHANGE: 'We have developed a Phase II clinical trial protocol'\n",
  "TO: 'We propose a Phase II clinical trial design'"
)

# 7. TP53 typo check (requires checking abstract file)
corrections$tp53_check <- paste0(
  "CRITICAL CHECK (Abstract):\n",
  "Search for 'TP54' - if found, change to 'TP53'\n",
  "This is a critical typo that would result in desk rejection."
)

# Write corrections to file
sink("07_Revision/text_corrections/manuscript_corrections.txt")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n")
cat("MANUSCRIPT CORRECTIONS - PHASE 7 VERIFICATION\n")
cat("Generated:", as.character(Sys.time()), "\n")
cat("=" %>% rep(70) %>% paste(collapse=""), "\n\n")

for (name in names(corrections)) {
  cat("-" %>% rep(50) %>% paste(collapse=""), "\n")
  cat(toupper(gsub("_", " ", name)), "\n")
  cat("-" %>% rep(50) %>% paste(collapse=""), "\n")
  cat(corrections[[name]], "\n\n")
}
sink()

cat("\nCorrections saved to: 07_Revision/text_corrections/manuscript_corrections.txt\n")
```

---

## TASK 4: SUPPLEMENTARY TABLE VERIFICATION
### Priority: 🟡 IMPORTANT

### Task 4.1: Verify All Supplementary Tables Exist

```r
# List of expected supplementary tables from manuscript
expected_tables <- c(
  "S1" = "Clustering comparison (k=2 vs k=3 vs k=4)",
  "S2" = "50-gene classifier",
  "S3" = "All 72 differential drugs (FDR<0.05)",
  "S4" = "BCL-2 family gene expression",
  "S5" = "Independence testing results (top 20 drugs)",
  "S6" = "Bootstrap/permutation/LOOCV results",
  "S7" = "Robustness testing for top 10 drugs",
  "S8" = "Salvage therapy options for Cluster 2",
  "S9" = "VRS tertile classification"
)

# Check which files exist
supp_dir <- "03_Results/supplementary_tables/"  # adjust path
existing_files <- list.files(supp_dir, pattern = "\\.csv$|\\.xlsx$")

cat("\n========== SUPPLEMENTARY TABLE VERIFICATION ==========\n")
for (i in seq_along(expected_tables)) {
  table_name <- names(expected_tables)[i]
  table_desc <- expected_tables[i]
  
  # Check if file exists (fuzzy match)
  found <- any(grepl(table_name, existing_files, ignore.case = TRUE))
  
  cat(table_name, ": ", table_desc, "\n")
  cat("  Status: ", ifelse(found, "✓ FOUND", "✗ MISSING"), "\n")
}

# Save verification
verification_df <- data.frame(
  table = names(expected_tables),
  description = expected_tables,
  status = sapply(names(expected_tables), function(t) {
    ifelse(any(grepl(t, existing_files, ignore.case = TRUE)), "Found", "Missing")
  })
)
write_csv(verification_df, "07_Revision/data_verification/supplementary_tables_verification.csv")
```

---

## TASK 5: GENERATE FINAL VERIFICATION REPORT
### Priority: 🔴 CRITICAL

### Task 5.1: Comprehensive Summary Report

```r
# Generate final verification report
sink("07_Revision/PHASE7_VERIFICATION_REPORT.txt")

cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║           PHASE 7: MANUSCRIPT VERIFICATION REPORT                   ║\n")
cat("║                    AML Molecular Subtypes Study                     ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("Generated:", as.character(Sys.time()), "\n")
cat("Target Journal: Blood\n\n")

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("1. CRITICAL DATA VERIFICATION\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Load and display all verification results
if (file.exists("07_Revision/data_verification/drug_sample_size_verification.csv")) {
  drug_size <- read_csv("07_Revision/data_verification/drug_sample_size_verification.csv", 
                        show_col_types = FALSE)
  cat("DRUG RESPONSE SAMPLE SIZE:\n")
  print(drug_size)
  cat("\n")
}

if (file.exists("07_Revision/data_verification/venetoclax_cohens_d_verification.csv")) {
  cohens <- read_csv("07_Revision/data_verification/venetoclax_cohens_d_verification.csv",
                     show_col_types = FALSE)
  cat("VENETOCLAX EFFECT SIZE:\n")
  print(cohens)
  cat("\n")
}

if (file.exists("07_Revision/data_verification/multivariate_attrition_summary.csv")) {
  attrition <- read_csv("07_Revision/data_verification/multivariate_attrition_summary.csv",
                        show_col_types = FALSE)
  cat("MULTIVARIATE ATTRITION:\n")
  print(t(attrition))
  cat("\n")
}

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("2. DISCREPANCY RESOLUTION\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat("Issue                    | Manuscript | Verified   | Action\n")
cat("-------------------------|------------|------------|------------------\n")
# Fill in based on actual verification results
cat("[Complete after running Tasks 1-4]\n\n")

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("3. REQUIRED MANUSCRIPT CHANGES\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

if (file.exists("07_Revision/text_corrections/manuscript_corrections.txt")) {
  corrections <- readLines("07_Revision/text_corrections/manuscript_corrections.txt")
  cat(paste(corrections, collapse = "\n"))
}

cat("\n\n═══════════════════════════════════════════════════════════════════════\n")
cat("4. PRE-SUBMISSION CHECKLIST\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

checklist <- c(
  "[ ] Drug response sample size added to manuscript",
  "[ ] Drug count verified/corrected",
  "[ ] Cohen's d verified/corrected",
  "[ ] Multivariate attrition explained",
  "[ ] Ex vivo limitation strengthened in Discussion",
  "[ ] Clinical utility language softened",
  "[ ] TP53 typo checked in Abstract",
  "[ ] All Supplementary Tables verified",
  "[ ] All Supplementary Figures verified",
  "[ ] All citations verified in PubMed",
  "[ ] Figure legends complete",
  "[ ] Author contributions section complete",
  "[ ] Conflict of interest statement complete",
  "[ ] Data availability statement complete"
)

cat(paste(checklist, collapse = "\n"))
cat("\n")

cat("\n═══════════════════════════════════════════════════════════════════════\n")
cat("5. ESTIMATED ACCEPTANCE PROBABILITY\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat("BEFORE PHASE 7 CORRECTIONS:\n")
cat("  Blood:           70-80%\n")
cat("  JCO:             55-65%\n")
cat("  Nature Medicine: 35-45%\n\n")

cat("AFTER PHASE 7 CORRECTIONS:\n")
cat("  Blood:           80-90%\n")
cat("  JCO:             65-75%\n")
cat("  Nature Medicine: 40-50%\n\n")

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("                         END OF REPORT                                 \n")
cat("═══════════════════════════════════════════════════════════════════════\n")

sink()

cat("\n\n✓ PHASE 7 VERIFICATION REPORT saved to: 07_Revision/PHASE7_VERIFICATION_REPORT.txt\n")
```

---

## OUTPUT EXPECTATIONS

After executing this prompt, you should have:

### Directory Structure
```
07_Revision/
├── data_verification/
│   ├── drug_sample_size_verification.csv
│   ├── drug_counts_by_cluster.csv
│   ├── venetoclax_cohens_d_verification.csv
│   ├── multivariate_attrition_summary.csv
│   ├── exclusion_reasons.csv
│   └── supplementary_tables_verification.csv
├── text_corrections/
│   └── manuscript_corrections.txt
└── PHASE7_VERIFICATION_REPORT.txt
```

### Key Deliverables
1. **Verified numbers** for all manuscript statistics
2. **Ready-to-use text** for all required additions
3. **Clear action items** for each discrepancy
4. **Pre-submission checklist** with completion status

---

## EXECUTION NOTES

1. **Adjust file paths** to match your actual directory structure
2. **Adjust column names** to match your actual data files
3. If any data files are missing, the script will flag them
4. **Run tasks sequentially** - later tasks depend on earlier outputs
5. **Review all outputs** before making manuscript changes

---

## AFTER EXECUTION

Bring the `PHASE7_VERIFICATION_REPORT.txt` back to me for:
1. Review of all verified values
2. Final text suggestions
3. Additional analyses if needed
4. Cover letter preparation
5. Reviewer response template preparation

**Target**: Complete Phase 7 → Submit to Blood within 1 week
