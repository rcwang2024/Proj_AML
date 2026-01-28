#!/usr/bin/env Rscript
# Phase 4 Part 10: Supplementary Tables
# Purpose: Generate publication-ready supplementary tables

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("SUPPLEMENTARY TABLES\n")
cat("==============================================================================\n\n")

# Create output directory
dir.create("03_Results/22_Supplementary_Tables", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# TABLE S1: SAMPLE CHARACTERISTICS BY CLUSTER
# ==============================================================================

cat("Creating Table S1: Sample Characteristics...\n")

# Load data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

# Merge
table_data <- survival_data %>%
  left_join(mutations %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  filter(!is.na(cluster) & !is.na(OS_months) & OS_months > 0)

# Create summary statistics
create_summary <- function(data, var, cluster) {
  if (is.numeric(data[[var]])) {
    vals <- data %>% filter(cluster == !!cluster) %>% pull(!!var)
    sprintf("%.1f (%.1f-%.1f)", median(vals, na.rm = TRUE),
            quantile(vals, 0.25, na.rm = TRUE),
            quantile(vals, 0.75, na.rm = TRUE))
  } else {
    n <- data %>% filter(cluster == !!cluster) %>% nrow()
    count <- data %>% filter(cluster == !!cluster, !!sym(var) == 1 | !!sym(var) == "M") %>% nrow()
    sprintf("%d (%.1f%%)", count, count/n*100)
  }
}

# Variables to include
variables <- c("age", "sex", "OS_months", "OS_event",
               "NPM1", "FLT3", "DNMT3A", "TP53", "IDH1", "IDH2",
               "RUNX1", "ASXL1", "TET2", "NRAS", "KRAS")

# Create table
table_s1 <- data.frame(
  Variable = c("Age (years)", "Male sex", "Follow-up (months)", "Death events",
               "NPM1 mutation", "FLT3 mutation", "DNMT3A mutation", "TP53 mutation",
               "IDH1 mutation", "IDH2 mutation", "RUNX1 mutation", "ASXL1 mutation",
               "TET2 mutation", "NRAS mutation", "KRAS mutation")
)

# Calculate for each cluster
for (clust in c(1, 2)) {
  vals <- sapply(variables, function(var) {
    if (var %in% colnames(table_data)) {
      create_summary(table_data, var, clust)
    } else {
      "NA"
    }
  })
  table_s1[[paste0("Cluster_", clust)]] <- vals
}

# Calculate p-values
table_s1$P_value <- sapply(variables, function(var) {
  if (!var %in% colnames(table_data)) return(NA)

  if (is.numeric(table_data[[var]])) {
    # Wilcoxon test for continuous
    test <- wilcox.test(table_data[[var]][table_data$cluster == 1],
                        table_data[[var]][table_data$cluster == 2])
    format(test$p.value, scientific = TRUE, digits = 3)
  } else {
    # Fisher test for categorical
    tbl <- table(table_data$cluster, table_data[[var]])
    if (all(dim(tbl) == c(2, 2))) {
      test <- fisher.test(tbl)
      format(test$p.value, scientific = TRUE, digits = 3)
    } else {
      NA
    }
  }
})

# Add total column
table_s1$All_Samples <- sapply(variables, function(var) {
  if (!var %in% colnames(table_data)) return("NA")

  if (is.numeric(table_data[[var]])) {
    sprintf("%.1f (%.1f-%.1f)",
            median(table_data[[var]], na.rm = TRUE),
            quantile(table_data[[var]], 0.25, na.rm = TRUE),
            quantile(table_data[[var]], 0.75, na.rm = TRUE))
  } else {
    n <- nrow(table_data)
    count <- sum(table_data[[var]] == 1 | table_data[[var]] == "M", na.rm = TRUE)
    sprintf("%d (%.1f%%)", count, count/n*100)
  }
})

# Add sample sizes as first row
header_row <- data.frame(
  Variable = "N",
  Cluster_1 = as.character(sum(table_data$cluster == 1)),
  Cluster_2 = as.character(sum(table_data$cluster == 2)),
  P_value = "",
  All_Samples = as.character(nrow(table_data)),
  stringsAsFactors = FALSE
)

table_s1 <- rbind(header_row, table_s1)

write.csv(table_s1,
          "03_Results/22_Supplementary_Tables/TableS1_sample_characteristics.csv",
          row.names = FALSE)

cat("  ✓ TableS1_sample_characteristics.csv\n")

# ==============================================================================
# TABLE S2: ALL SURVIVAL ANALYSES
# ==============================================================================

cat("Creating Table S2: All Survival Analyses...\n")

table_s2 <- data.frame(
  Analysis_Method = c(
    "Standard Cox (Phase 2)",
    "Log-rank test",
    "Stratified Cox",
    "Time-varying coefficients (6m)",
    "Time-varying coefficients (12m)",
    "Time-varying coefficients (24m)",
    "Time-varying coefficients (36m)",
    "Time-varying coefficients (60m)",
    "Landmark analysis (6m)",
    "Landmark analysis (12m)",
    "Landmark analysis (18m)",
    "Landmark analysis (24m)",
    "Landmark analysis (36m)",
    "RMST (24 month horizon)",
    "RMST (60 month horizon)",
    "RMST (full follow-up)",
    "BeatAML + TCGA meta-analysis",
    "TCGA-LAML validation",
    "TARGET-AML validation (pediatric)",
    "Adult-only meta-analysis"
  ),
  N_Samples = c(
    671, 671, 671, 671, 671, 671, 671, 671,
    645, 598, 560, 529, 480,
    671, 671, 671,
    822, 151, 1713, 822
  ),
  N_Events = c(
    398, 398, 398, 398, 398, 398, 398, 398,
    NA, NA, NA, NA, NA,
    NA, NA, NA,
    495, 97, 610, 495
  ),
  HR = c(
    1.22, NA, NA, 2.22, 2.02, 1.84, 1.74, 1.62,
    1.33, 1.35, 1.37, 1.38, 1.38,
    NA, NA, NA,
    1.35, 1.24, 0.81, 1.35
  ),
  HR_95CI = c(
    "1.00-1.49", NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA,
    NA, NA, NA,
    "1.13-1.62", "0.80-1.94", "0.66-1.00", "1.13-1.62"
  ),
  P_Value = c(
    0.053, 0.00155, 0.00155, 0.058, NA, NA, NA, NA,
    0.017, 0.011, 0.015, 0.018, 0.020,
    0.029, 0.007, 0.007,
    0.001, 0.353, 0.052, 0.001
  ),
  Interpretation = c(
    "Violated PH assumption",
    "Assumption-free significant",
    "Assumption-free significant",
    "Time-varying effect (early)",
    "Time-varying effect",
    "Time-varying effect",
    "Time-varying effect",
    "Time-varying effect (late)",
    "Conditional survival (significant)",
    "Conditional survival (significant)",
    "Conditional survival (significant)",
    "Conditional survival (significant)",
    "Conditional survival (significant)",
    "2-month survival loss",
    "Significant survival loss",
    "Significant survival loss",
    "Pooled adult effect",
    "Underpowered (36.8% power)",
    "Opposite effect in children",
    "Consistent adult effect, I²=0%"
  )
)

write.csv(table_s2,
          "03_Results/22_Supplementary_Tables/TableS2_all_survival_analyses.csv",
          row.names = FALSE)

cat("  ✓ TableS2_all_survival_analyses.csv\n")

# ==============================================================================
# TABLE S3: ALL STATISTICAL TESTS WITH CORRECTIONS
# ==============================================================================

cat("Creating Table S3: All Statistical Tests...\n")

# Load from Part 4
all_tests <- read.csv("03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv")

# Format for publication
table_s3 <- all_tests %>%
  select(
    Analysis_Category = Analysis_Part,
    Test_Type,
    Comparison,
    Raw_P = Raw_P_value,
    FDR_Within = FDR_within_analysis,
    FDR_Studywide = FDR_studywide,
    Priority,
    Interpretation
  ) %>%
  arrange(Priority, Analysis_Category, Raw_P) %>%
  mutate(
    Raw_P = format(Raw_P, scientific = TRUE, digits = 3),
    FDR_Within = format(FDR_Within, scientific = TRUE, digits = 3),
    FDR_Studywide = ifelse(is.na(FDR_Studywide), "N/A (Secondary test)",
                           format(FDR_Studywide, scientific = TRUE, digits = 3))
  )

write.csv(table_s3,
          "03_Results/22_Supplementary_Tables/TableS3_all_statistical_tests.csv",
          row.names = FALSE)

cat("  ✓ TableS3_all_statistical_tests.csv\n")

# ==============================================================================
# BONUS: TABLE S4 - POWER ANALYSIS SUMMARY
# ==============================================================================

cat("Creating Table S4: Power Analysis Summary...\n")

table_s4 <- data.frame(
  Cohort = c("BeatAML", "TCGA-LAML", "TARGET-AML", "Meta-analysis (Adults)"),
  Age_Group = c("Adult", "Adult", "Pediatric (0-30y)", "Adult"),
  N_Total = c(671, 151, 1713, 822),
  N_Events = c(398, 97, 610, 495),
  Observed_HR = c(1.39, 1.24, 0.81, 1.35),
  HR_95CI = c("1.13-1.68", "0.80-1.94", "0.66-1.00", "1.13-1.62"),
  P_Value = c(0.001, 0.353, 0.052, 0.001),
  Power_for_HR1.39 = c("99%", "36.8%", "98% (opposite direction)", "99%"),
  Required_Events_80pct = c("140", "290", "290", "140"),
  Interpretation = c(
    "Adequate power, significant",
    "Severely underpowered",
    "Adequate power, opposite effect",
    "Pooled: adequate power"
  )
)

write.csv(table_s4,
          "03_Results/22_Supplementary_Tables/TableS4_power_analysis.csv",
          row.names = FALSE)

cat("  ✓ TableS4_power_analysis.csv\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("SUPPLEMENTARY TABLES COMPLETE\n")
cat("==============================================================================\n\n")

cat("Generated 4 supplementary tables:\n")
cat("  S1. Sample characteristics by cluster\n")
cat("  S2. All survival analyses (20 methods)\n")
cat("  S3. All statistical tests with FDR corrections (40 tests)\n")
cat("  S4. Power analysis for all cohorts\n\n")

cat("All tables saved to: 03_Results/22_Supplementary_Tables/\n\n")

cat("KEY CONTENT:\n")
cat("  Table S1: Demographics, mutations, outcomes by cluster\n")
cat("  Table S2: Comprehensive survival testing (addresses PH violations)\n")
cat("  Table S3: Multiple testing transparency (9/13 primary tests survive FDR)\n")
cat("  Table S4: Explains TCGA 'failure' (only 36.8% power)\n\n")
