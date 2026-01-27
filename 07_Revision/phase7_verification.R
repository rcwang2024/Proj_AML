# Phase 7: Manuscript Revision & Data Verification
# Systematic verification of all manuscript statistics
# Date: 2025-12-10

setwd("D:/Projects/Project_AML")

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

cat("\n")
cat("═══════════════════════════════════════════════════════════════════════\n")
cat("           PHASE 7: MANUSCRIPT DATA VERIFICATION                       \n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# ============================================================================
# TASK 1.1: DRUG RESPONSE SAMPLE SIZE VERIFICATION
# ============================================================================

cat("\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TASK 1.1: DRUG RESPONSE SAMPLE SIZE VERIFICATION                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# Load cluster assignments
clusters <- read_csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv",
                     show_col_types = FALSE)

# Load drug response data
drug_response <- read_tsv("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt",
                          show_col_types = FALSE)

# Identify patients
expr_patients <- unique(clusters$sample_id)
drug_patients <- unique(drug_response$dbgap_rnaseq_sample)

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


# ============================================================================
# TASK 1.2: DRUG COUNT VERIFICATION
# ============================================================================

cat("\n\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TASK 1.2: DRUG COUNT VERIFICATION                                  ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# Count unique drugs in dataset
drugs <- unique(drug_response$inhibitor)
n_drugs_total <- length(drugs)

# Count drugs with sufficient data for analysis
drug_counts <- drug_response %>%
  left_join(clusters %>% select(sample_id, cluster) %>%
              rename(dbgap_rnaseq_sample = sample_id),
            by = "dbgap_rnaseq_sample") %>%
  filter(!is.na(cluster)) %>%
  group_by(inhibitor, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(inhibitor) %>%
  summarise(
    min_per_cluster = min(n),
    total_n = sum(n),
    .groups = "drop"
  )

# Drugs with minimum N per cluster
n_drugs_analyzed_n5 <- sum(drug_counts$min_per_cluster >= 5)
n_drugs_analyzed_n10 <- sum(drug_counts$min_per_cluster >= 10)
n_drugs_analyzed_n30 <- sum(drug_counts$min_per_cluster >= 30)

cat("\n========== DRUG COUNT VERIFICATION ==========\n")
cat("Total unique drugs in database:", n_drugs_total, "\n")
cat("Drugs with n≥5 per cluster:", n_drugs_analyzed_n5, "\n")
cat("Drugs with n≥10 per cluster:", n_drugs_analyzed_n10, "\n")
cat("Drugs with n≥30 per cluster:", n_drugs_analyzed_n30, "\n")

cat("\nDISCREPANCY CHECK:\n")
cat("Manuscript states: 155 drugs\n")
cat("Previous documentation: 122 compounds\n")
cat("Actual count: ", n_drugs_total, " total, ", n_drugs_analyzed_n30, " analyzed (n≥30)\n")

# Save
write_csv(drug_counts, "07_Revision/data_verification/drug_counts_by_cluster.csv")


# ============================================================================
# TASK 1.3: COHEN'S D VERIFICATION FOR VENETOCLAX
# ============================================================================

cat("\n\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TASK 1.3: COHEN'S D VERIFICATION FOR VENETOCLAX                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# Get Venetoclax data
venetoclax_data <- drug_response %>%
  filter(grepl("venetoclax|ABT.199", inhibitor, ignore.case = TRUE)) %>%
  left_join(clusters %>% select(sample_id, cluster) %>%
              rename(dbgap_rnaseq_sample = sample_id),
            by = "dbgap_rnaseq_sample") %>%
  filter(!is.na(cluster), !is.na(auc))

# Sample sizes per cluster
n_per_cluster <- venetoclax_data %>% count(cluster)
cat("\n========== VENETOCLAX SAMPLE SIZES ==========\n")
print(n_per_cluster)

# Calculate Cohen's d
cluster1_auc <- venetoclax_data$auc[venetoclax_data$cluster == 1]
cluster2_auc <- venetoclax_data$auc[venetoclax_data$cluster == 2]

# Manual calculation
mean1 <- mean(cluster1_auc, na.rm = TRUE)
mean2 <- mean(cluster2_auc, na.rm = TRUE)
sd1 <- sd(cluster1_auc, na.rm = TRUE)
sd2 <- sd(cluster2_auc, na.rm = TRUE)
n1 <- length(cluster1_auc)
n2 <- length(cluster2_auc)

# Pooled SD
pooled_sd <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2))
manual_d <- (mean1 - mean2) / pooled_sd

# Magnitude interpretation
magnitude <- ifelse(abs(manual_d) >= 0.8, "very large",
                   ifelse(abs(manual_d) >= 0.5, "large",
                         ifelse(abs(manual_d) >= 0.2, "medium", "small")))

cat("\n========== COHEN'S D VERIFICATION ==========\n")
cat("Manual calculation d:", round(manual_d, 3), "\n")
cat("Magnitude:", magnitude, "\n")

cat("\nDISCREPANCY CHECK:\n")
cat("Manuscript states: d = 1.25\n")
cat("Previous discussions: d > 2.0\n")
cat("Actual calculated d:", round(manual_d, 3), "\n")

# Additional details
cat("\nDetailed Statistics:\n")
cat("Cluster 1 - Mean AUC:", round(mean1, 2), "SD:", round(sd1, 2), "n:", n1, "\n")
cat("Cluster 2 - Mean AUC:", round(mean2, 2), "SD:", round(sd2, 2), "n:", n2, "\n")
cat("Pooled SD:", round(pooled_sd, 2), "\n")

# Calculate 95% CI for Cohen's d
d_se <- sqrt((n1+n2)/(n1*n2) + (manual_d^2)/(2*(n1+n2)))
d_ci_lower <- manual_d - 1.96*d_se
d_ci_upper <- manual_d + 1.96*d_se

cat("95% CI for d: [", round(d_ci_lower, 3), ", ", round(d_ci_upper, 3), "]\n")

# Save verification
verification_df <- data.frame(
  metric = c("cohens_d", "ci_lower", "ci_upper", "magnitude",
             "cluster1_mean", "cluster1_sd", "cluster1_n",
             "cluster2_mean", "cluster2_sd", "cluster2_n"),
  value = c(manual_d, d_ci_lower, d_ci_upper, magnitude,
            mean1, sd1, n1, mean2, sd2, n2)
)
write_csv(verification_df, "07_Revision/data_verification/venetoclax_cohens_d_verification.csv")


# ============================================================================
# TASK 1.4: MULTIVARIATE SAMPLE ATTRITION ANALYSIS
# ============================================================================

cat("\n\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  TASK 1.4: MULTIVARIATE SAMPLE ATTRITION ANALYSIS                   ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# Load clinical data (has some mutations and clinical variables)
clinical <- read_csv("03_Results/05_Analysis_Ready_Data/clinical_gold_standard.csv",
                     show_col_types = FALSE)

# Load mutations data (has TET2 and other mutations)
mutations <- read_csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv",
                      show_col_types = FALSE) %>%
  rename(dbgap_rnaseq_sample = ...1)  # First column is sample ID

# Merge cluster assignments with clinical and mutation data
# Rename sample_id to dbgap_rnaseq_sample for easier joining
clusters_renamed <- clusters %>% rename(dbgap_rnaseq_sample = sample_id)

all_data <- clusters_renamed %>%
  left_join(clinical, by = "dbgap_rnaseq_sample") %>%
  left_join(mutations %>% select(dbgap_rnaseq_sample, TET2), by = "dbgap_rnaseq_sample")

# Check completeness for multivariate model variables
completeness <- all_data %>%
  mutate(
    has_age = !is.na(ageAtDiagnosis),
    has_sex = !is.na(consensus_sex),
    has_survival = !is.na(overallSurvival) & !is.na(vitalStatus),
    has_TP53 = !is.na(TP53),
    has_TET2 = !is.na(TET2),
    has_RUNX1 = !is.na(RUNX1),
    has_ASXL1 = !is.na(ASXL1),
    has_all_mutations = has_TP53 & has_TET2 & has_RUNX1 & has_ASXL1,
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
    has_RUNX1 = sum(has_RUNX1),
    has_ASXL1 = sum(has_ASXL1),
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


# ============================================================================
# SUMMARY OF CRITICAL FINDINGS
# ============================================================================

cat("\n\n╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                  VERIFICATION SUMMARY                                ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

cat("1. Drug response sample size:", n_drug_analysis, "patients\n")
cat("2. Drug count: Total =", n_drugs_total, ", Analyzed (n≥30) =", n_drugs_analyzed_n30, "\n")
cat("3. Venetoclax Cohen's d:", round(manual_d, 3),
    "(", round(d_ci_lower, 2), "-", round(d_ci_upper, 2), ")\n")
cat("4. Multivariate attrition:", attrition_summary$total_clustered, "→",
    attrition_summary$complete_cases, "patients\n")

cat("\n✓ All verification files saved to: 07_Revision/data_verification/\n")
cat("\nNext steps:\n")
cat("1. Review verification results\n")
cat("2. Generate corrected manuscript text\n")
cat("3. Verify supplementary tables\n")
cat("4. Compile final report\n\n")
