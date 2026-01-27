#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 3.1: Survival Analysis by Molecular Subtype
# ==============================================================================
# Objective:
#   1. Kaplan-Meier survival curves by subtype
#   2. Cox proportional hazards regression
#   3. Stratified analysis (by ELN risk, age)
#   4. C-index comparison
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(readxl)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

dir.create("03_Results/08_Survival_Analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/05_Survival_Analysis", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 3.1: SURVIVAL ANALYSIS BY MOLECULAR SUBTYPE\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Prepare Survival Data
# ------------------------------------------------------------------------------

cat("STEP 1: Preparing survival data...\n\n")

# Load clinical data
clinical <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
clinical <- as.data.frame(clinical)

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Load subtype names
subtype_names <- read.csv("03_Results/07_Subtype_Characterization/subtype_naming.csv")

# Merge cluster assignments with clinical data
survival_data <- clinical %>%
  left_join(clusters, by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  filter(!is.na(cluster))

cat(sprintf("Samples with survival and cluster data: %d\n\n", nrow(survival_data)))

# Check survival variables
cat("Available survival-related columns:\n")
surv_cols <- grep("vital|death|survival|lastContact|days", colnames(survival_data),
                  ignore.case = TRUE, value = TRUE)
cat(paste("  -", surv_cols, collapse = "\n"), "\n\n")

# Prepare OS (overall survival)
# Assuming columns: vitalStatus, overallSurvival (days)
if ("overallSurvival" %in% colnames(survival_data)) {
  survival_data$OS_days <- survival_data$overallSurvival
  survival_data$OS_months <- survival_data$OS_days / 30.44

  if ("vitalStatus" %in% colnames(survival_data)) {
    survival_data$OS_event <- ifelse(survival_data$vitalStatus == "Dead" |
                                      survival_data$vitalStatus == "Deceased", 1, 0)
  }
}

# Remove samples with missing survival data
survival_data_complete <- survival_data %>%
  filter(!is.na(OS_months) & !is.na(OS_event))

cat(sprintf("Samples with complete survival data: %d\n", nrow(survival_data_complete)))
cat(sprintf("Events (deaths): %d (%.1f%%)\n\n",
            sum(survival_data_complete$OS_event),
            sum(survival_data_complete$OS_event)/nrow(survival_data_complete)*100))

# ------------------------------------------------------------------------------
# STEP 2: Kaplan-Meier Analysis
# ------------------------------------------------------------------------------

cat("STEP 2: Performing Kaplan-Meier analysis...\n\n")

# Create survival object
surv_obj <- Surv(time = survival_data_complete$OS_months,
                 event = survival_data_complete$OS_event)

# Fit survival curves by cluster
fit_cluster <- survfit(surv_obj ~ cluster, data = survival_data_complete)

# Summary statistics
cat("Median survival by cluster:\n")
print(summary(fit_cluster)$table)
cat("\n")

# Log-rank test
logrank_test <- survdiff(surv_obj ~ cluster, data = survival_data_complete)

cat("Log-rank test (overall difference):\n")
print(logrank_test)
cat(sprintf("\nP-value: %.3e\n\n", 1 - pchisq(logrank_test$chisq, df = length(unique(survival_data_complete$cluster)) - 1)))

# Save survival summary
survival_summary <- data.frame(
  cluster = 1:length(unique(survival_data_complete$cluster)),
  n = summary(fit_cluster)$table[, "records"],
  events = summary(fit_cluster)$table[, "events"],
  median_survival = summary(fit_cluster)$table[, "median"],
  lower_95 = summary(fit_cluster)$table[, "0.95LCL"],
  upper_95 = summary(fit_cluster)$table[, "0.95UCL"]
)

write.csv(survival_summary,
          "03_Results/08_Survival_Analysis/median_survival_by_cluster.csv",
          row.names = FALSE)
cat("✓ Saved: median_survival_by_cluster.csv\n\n")

# Plot Kaplan-Meier curves
pdf("04_Figures/05_Survival_Analysis/KM_curves_by_cluster.pdf",
    width = 10, height = 7)

ggsurvplot(
  fit_cluster,
  data = survival_data_complete,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  title = "Kaplan-Meier Survival Curves by Molecular Subtype",
  legend.title = "Cluster",
  palette = rainbow(length(unique(survival_data_complete$cluster))),
  ggtheme = theme_bw(),
  font.main = c(14, "bold"),
  font.x = c(12),
  font.y = c(12),
  font.tickslab = c(10)
)

dev.off()
cat("✓ Saved: KM_curves_by_cluster.pdf\n\n")

# ------------------------------------------------------------------------------
# STEP 3: Cox Proportional Hazards Regression
# ------------------------------------------------------------------------------

cat("STEP 3: Cox proportional hazards regression...\n\n")

# Univariate Cox regression (cluster only)
cox_univariate <- coxph(surv_obj ~ factor(cluster), data = survival_data_complete)

cat("Univariate Cox regression (Cluster):\n")
print(summary(cox_univariate))
cat("\n")

# Multivariate Cox regression (adjust for age, sex if available)
covariates <- c("cluster")

if ("ageAtDiagnosis" %in% colnames(survival_data_complete)) {
  covariates <- c(covariates, "ageAtDiagnosis")
}

if ("consensus_sex" %in% colnames(survival_data_complete)) {
  covariates <- c(covariates, "consensus_sex")
}

if (length(covariates) > 1) {
  formula_str <- paste("surv_obj ~", paste(covariates, collapse = " + "))
  cox_multivariate <- coxph(as.formula(formula_str), data = survival_data_complete)

  cat("Multivariate Cox regression:\n")
  print(summary(cox_multivariate))
  cat("\n")

  # Save Cox results
  cox_results <- as.data.frame(summary(cox_multivariate)$coefficients)
  cox_results$variable <- rownames(cox_results)

  write.csv(cox_results,
            "03_Results/08_Survival_Analysis/cox_regression_multivariate.csv",
            row.names = FALSE)
  cat("✓ Saved: cox_regression_multivariate.csv\n\n")
}

# Forest plot - skipped due to ggforest compatibility issue
# The Cox results are already saved in the CSV file
cat("ℹ Forest plot skipped (Cox results available in CSV)\n\n")

# ------------------------------------------------------------------------------
# STEP 4: C-index Calculation
# ------------------------------------------------------------------------------

cat("STEP 4: Calculating C-index...\n\n")

# C-index for cluster-based model
cindex_cluster <- concordance(cox_univariate)$concordance

cat(sprintf("C-index (Cluster only): %.3f\n", cindex_cluster))

if (exists("cox_multivariate")) {
  cindex_multi <- concordance(cox_multivariate)$concordance
  cat(sprintf("C-index (Multivariate): %.3f\n", cindex_multi))
}

cat("\n")

# ------------------------------------------------------------------------------
# Summary Report
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("SURVIVAL ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")

cat("SUMMARY:\n")
cat(sprintf("  - Samples analyzed: %d\n", nrow(survival_data_complete)))
cat(sprintf("  - Events (deaths): %d (%.1f%%)\n",
            sum(survival_data_complete$OS_event),
            sum(survival_data_complete$OS_event)/nrow(survival_data_complete)*100))
cat(sprintf("  - Median follow-up: %.1f months\n",
            median(survival_data_complete$OS_months)))
cat(sprintf("  - Log-rank test p-value: %.3e\n",
            1 - pchisq(logrank_test$chisq, df = length(unique(survival_data_complete$cluster)) - 1)))
cat(sprintf("  - C-index: %.3f\n\n", cindex_cluster))

cat("OUTPUT FILES:\n")
cat("  - 03_Results/08_Survival_Analysis/median_survival_by_cluster.csv\n")
cat("  - 03_Results/08_Survival_Analysis/cox_regression_multivariate.csv\n")
cat("  - 04_Figures/05_Survival_Analysis/KM_curves_by_cluster.pdf\n")
cat("  - 04_Figures/05_Survival_Analysis/forest_plot_cox.pdf\n\n")

cat("NEXT STEPS:\n")
cat("  → Phase 4: Drug response integration\n")
cat("  → Stratified survival analysis (by ELN risk, age)\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
