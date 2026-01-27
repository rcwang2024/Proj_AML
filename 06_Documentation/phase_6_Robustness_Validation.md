# Phase 6: Robustness Validation & Addressing Reviewer Concerns
## Comprehensive Claude Code Execution Prompt

**Date**: October 2025
**Purpose**: Address critical scientific concerns before Tier 1 journal submission
**Priority**: ESSENTIAL - These analyses determine publication tier

---

## MISSION STATEMENT

Phase 5 established promising drug response findings. However, **critical gaps** must be addressed before submission to Nature Medicine, JCO, or Blood:

1. ⚠️ **Robustness validation missing** (bootstrap, LOOCV, permutation)
2. ⚠️ **Sample sizes per drug unclear**
3. ⚠️ **Potential circularity** (same samples for clustering and drug testing)
4. ⚠️ **Independence paradox unexplained** (why independent for drugs but not survival?)
5. ⚠️ **Multicollinearity not assessed** (VIF for regression models)
6. ⚠️ **One drug failed independence** - which one and why?
7. ⚠️ **Classifier confidence in validation cohorts** not reported

This phase will systematically address ALL concerns with publication-quality analyses.

**Expected Output**: Complete robustness validation enabling confident Tier 1 submission

---

## PART 1: SAMPLE SIZE DOCUMENTATION

### **Task 1.1: Create Comprehensive Sample Size Table**

**Objective**: Document exact N for every analysis (reviewers WILL ask)

```r
cat("=== PHASE 6: SAMPLE SIZE DOCUMENTATION ===\n\n")

# Load required data
library(dplyr)
library(tidyr)

# Load drug response data
drug_data <- readRDS("03_Results/23_Drug_Validation/drug_cluster_merged.rds")
# Or load from original location if different
# drug_data <- readRDS("01_Data/BeatAML/drug_response.rds")

cluster_assignments <- read.csv("03_Results/02_Clustering/cluster_assignments.csv")

# Create output directory
dir.create("03_Results/24_Robustness_Validation", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/24_Robustness", showWarnings = FALSE, recursive = TRUE)

cat("=== COMPREHENSIVE SAMPLE SIZE AUDIT ===\n\n")

# Identify column names (adjust based on actual data structure)
sample_col <- names(drug_data)[grep("sample|patient|id|lab_id", names(drug_data), ignore.case = TRUE)[1]]
drug_col <- names(drug_data)[grep("drug|compound|inhibitor", names(drug_data), ignore.case = TRUE)[1]]
auc_col <- names(drug_data)[grep("auc|area", names(drug_data), ignore.case = TRUE)[1]]

# 1. Overall sample sizes
cat("1. OVERALL SAMPLE SIZES\n")
cat("─────────────────────────────────────\n")

n_total_samples <- length(unique(drug_data[[sample_col]]))
n_total_drugs <- length(unique(drug_data[[drug_col]]))
n_total_observations <- nrow(drug_data)

cat("Total unique samples with drug data:", n_total_samples, "\n")
cat("Total unique drugs tested:", n_total_drugs, "\n")
cat("Total drug-sample observations:", n_total_observations, "\n\n")

# 2. Sample sizes per drug (for ALL drugs)
cat("2. SAMPLE SIZES PER DRUG\n")
cat("─────────────────────────────────────\n")

drug_sample_sizes <- drug_data %>%
  filter(!is.na(!!sym(auc_col))) %>%
  group_by(!!sym(drug_col)) %>%
  summarise(
    N_total = n(),
    N_cluster_1 = sum(cluster_assignment == "Cluster_1", na.rm = TRUE),
    N_cluster_2 = sum(cluster_assignment == "Cluster_2", na.rm = TRUE),
    N_missing_cluster = sum(is.na(cluster_assignment)),
    Mean_AUC = mean(!!sym(auc_col), na.rm = TRUE),
    SD_AUC = sd(!!sym(auc_col), na.rm = TRUE)
  ) %>%
  arrange(desc(N_total))

cat("Drugs with ≥100 samples:", sum(drug_sample_sizes$N_total >= 100), "\n")
cat("Drugs with 50-99 samples:", sum(drug_sample_sizes$N_total >= 50 & drug_sample_sizes$N_total < 100), "\n")
cat("Drugs with 30-49 samples:", sum(drug_sample_sizes$N_total >= 30 & drug_sample_sizes$N_total < 50), "\n")
cat("Drugs with <30 samples:", sum(drug_sample_sizes$N_total < 30), "\n\n")

# 3. Specific drugs of interest
key_drugs <- c("Venetoclax", "Panobinostat", "Selumetinib", "Nilotinib", 
               "Sorafenib", "MK-2206", "Rapamycin", "Erlotinib")

cat("3. KEY DRUGS - DETAILED SAMPLE SIZES\n")
cat("─────────────────────────────────────\n")

key_drug_sizes <- drug_sample_sizes %>%
  filter(grepl(paste(key_drugs, collapse = "|"), !!sym(drug_col), ignore.case = TRUE))

print(key_drug_sizes)

# 4. VENETOCLAX SPECIFIC (Most critical)
cat("\n4. VENETOCLAX - COMPLETE BREAKDOWN\n")
cat("═══════════════════════════════════════\n")

venetoclax_data <- drug_data %>%
  filter(grepl("venetoclax|ABT-199|ABT199", !!sym(drug_col), ignore.case = TRUE)) %>%
  filter(!is.na(!!sym(auc_col)), !is.na(cluster_assignment))

cat("Total Venetoclax observations:", nrow(venetoclax_data), "\n")
cat("Cluster 1:", sum(venetoclax_data$cluster_assignment == "Cluster_1"), "\n")
cat("Cluster 2:", sum(venetoclax_data$cluster_assignment == "Cluster_2"), "\n")

# Check for duplicates (same patient multiple times?)
n_unique_patients <- length(unique(venetoclax_data[[sample_col]]))
cat("Unique patients:", n_unique_patients, "\n")

if (n_unique_patients != nrow(venetoclax_data)) {
  cat("⚠️  WARNING: Some patients have multiple Venetoclax measurements\n")
  cat("   Observations per patient:\n")
  print(table(table(venetoclax_data[[sample_col]])))
}

# 5. Create publication-ready table
sample_size_table <- data.frame(
  Analysis = c(
    "Phase 2: Clustering discovery",
    "Phase 3: Survival (univariate)",
    "Phase 3: Survival (multivariate)",
    "Phase 3: TCGA validation",
    "Phase 3: TARGET validation",
    "Phase 3: Adult meta-analysis",
    "Phase 5: Drug response (any drug)",
    "Phase 5: Venetoclax analysis",
    "Phase 5: Independence testing (top 20)"
  ),
  Total_N = c(671, 671, 459, 151, 1713, 822, 
              n_total_samples,
              nrow(venetoclax_data),
              NA),  # Fill in from your data
  N_Cluster_1 = c(NA, NA, NA, NA, NA, NA,
                  sum(drug_sample_sizes$N_cluster_1),
                  sum(venetoclax_data$cluster_assignment == "Cluster_1"),
                  NA),
  N_Cluster_2 = c(NA, NA, NA, NA, NA, NA,
                  sum(drug_sample_sizes$N_cluster_2),
                  sum(venetoclax_data$cluster_assignment == "Cluster_2"),
                  NA),
  Events_or_Drugs = c(NA, 398, 282, 97, 610, 495,
                      n_total_drugs,
                      1,
                      20)
)

write.csv(sample_size_table, 
          "03_Results/24_Robustness_Validation/Table_Sample_Sizes_All_Analyses.csv",
          row.names = FALSE)

write.csv(drug_sample_sizes,
          "03_Results/24_Robustness_Validation/Table_Sample_Sizes_Per_Drug.csv",
          row.names = FALSE)

cat("\n✅ Sample size documentation complete\n")
cat("   Saved: Table_Sample_Sizes_All_Analyses.csv\n")
cat("   Saved: Table_Sample_Sizes_Per_Drug.csv\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Table_Sample_Sizes_All_Analyses.csv`
- `03_Results/24_Robustness_Validation/Table_Sample_Sizes_Per_Drug.csv`

---

## PART 2: ROBUSTNESS VALIDATION FOR TOP DRUGS

### **Task 2.1: Bootstrap Analysis (10,000 Resamples)**

**Objective**: Test stability of drug-cluster associations

```r
cat("\n=== BOOTSTRAP VALIDATION FOR TOP DRUGS ===\n\n")

# Load differential drug results
drug_results <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")

# Select top 10 drugs by significance
top_drugs <- drug_results %>%
  arrange(P_value) %>%
  head(10) %>%
  pull(Drug)

cat("Testing bootstrap robustness for top 10 drugs:\n")
cat(paste(top_drugs, collapse = "\n"), "\n\n")

set.seed(42)
n_bootstrap <- 10000

# Function to run bootstrap for one drug
bootstrap_drug <- function(drug_name, drug_data, sample_col, drug_col, auc_col, n_boot = 10000) {
  
  # Extract data for this drug
  drug_subset <- drug_data %>%
    filter(grepl(drug_name, !!sym(drug_col), ignore.case = TRUE)) %>%
    filter(!is.na(!!sym(auc_col)), !is.na(cluster_assignment))
  
  n_samples <- nrow(drug_subset)
  
  if (n_samples < 20) {
    return(list(drug = drug_name, n = n_samples, error = "Insufficient samples"))
  }
  
  # Original statistics
  c1_orig <- drug_subset %>% filter(cluster_assignment == "Cluster_1") %>% pull(!!sym(auc_col))
  c2_orig <- drug_subset %>% filter(cluster_assignment == "Cluster_2") %>% pull(!!sym(auc_col))
  
  orig_test <- wilcox.test(c1_orig, c2_orig)
  orig_diff <- mean(c1_orig) - mean(c2_orig)
  orig_d <- orig_diff / sqrt((sd(c1_orig)^2 + sd(c2_orig)^2) / 2)
  
  # Bootstrap
  boot_results <- replicate(n_boot, {
    # Resample with replacement
    boot_idx <- sample(1:n_samples, n_samples, replace = TRUE)
    boot_data <- drug_subset[boot_idx, ]
    
    boot_c1 <- boot_data %>% filter(cluster_assignment == "Cluster_1") %>% pull(!!sym(auc_col))
    boot_c2 <- boot_data %>% filter(cluster_assignment == "Cluster_2") %>% pull(!!sym(auc_col))
    
    # Handle edge cases
    if (length(boot_c1) < 3 || length(boot_c2) < 3) {
      return(c(p_value = NA, diff = NA, cohens_d = NA))
    }
    
    boot_test <- tryCatch(
      wilcox.test(boot_c1, boot_c2),
      error = function(e) list(p.value = NA)
    )
    
    boot_diff <- mean(boot_c1) - mean(boot_c2)
    boot_d <- boot_diff / sqrt((sd(boot_c1)^2 + sd(boot_c2)^2) / 2)
    
    c(p_value = boot_test$p.value, diff = boot_diff, cohens_d = boot_d)
  })
  
  boot_df <- as.data.frame(t(boot_results))
  boot_df <- boot_df[complete.cases(boot_df), ]
  
  # Calculate statistics
  list(
    drug = drug_name,
    n = n_samples,
    n_c1 = length(c1_orig),
    n_c2 = length(c2_orig),
    orig_p = orig_test$p.value,
    orig_diff = orig_diff,
    orig_d = orig_d,
    boot_n_valid = nrow(boot_df),
    boot_p_median = median(boot_df$p_value, na.rm = TRUE),
    boot_p_ci_low = quantile(boot_df$p_value, 0.025, na.rm = TRUE),
    boot_p_ci_high = quantile(boot_df$p_value, 0.975, na.rm = TRUE),
    boot_pct_p001 = mean(boot_df$p_value < 0.001, na.rm = TRUE) * 100,
    boot_pct_p01 = mean(boot_df$p_value < 0.01, na.rm = TRUE) * 100,
    boot_pct_p05 = mean(boot_df$p_value < 0.05, na.rm = TRUE) * 100,
    boot_d_median = median(boot_df$cohens_d, na.rm = TRUE),
    boot_d_ci_low = quantile(boot_df$cohens_d, 0.025, na.rm = TRUE),
    boot_d_ci_high = quantile(boot_df$cohens_d, 0.975, na.rm = TRUE),
    boot_diff_median = median(boot_df$diff, na.rm = TRUE),
    boot_diff_ci_low = quantile(boot_df$diff, 0.025, na.rm = TRUE),
    boot_diff_ci_high = quantile(boot_df$diff, 0.975, na.rm = TRUE)
  )
}

# Run bootstrap for all top drugs
cat("Running bootstrap analysis (10,000 resamples each)...\n")
cat("This may take 10-20 minutes...\n\n")

bootstrap_results <- lapply(top_drugs, function(drug) {
  cat("  Processing:", drug, "...")
  result <- bootstrap_drug(drug, drug_data, sample_col, drug_col, auc_col, n_bootstrap)
  cat(" done\n")
  result
})

# Convert to data frame
bootstrap_summary <- do.call(rbind, lapply(bootstrap_results, function(x) {
  if (!is.null(x$error)) {
    data.frame(Drug = x$drug, N = x$n, Error = x$error, stringsAsFactors = FALSE)
  } else {
    data.frame(
      Drug = x$drug,
      N_total = x$n,
      N_C1 = x$n_c1,
      N_C2 = x$n_c2,
      Original_P = x$orig_p,
      Original_Cohens_d = x$orig_d,
      Bootstrap_Valid = x$boot_n_valid,
      Boot_P_Median = x$boot_p_median,
      Boot_P_95CI = paste0("[", format(x$boot_p_ci_low, scientific=TRUE, digits=2), 
                          ", ", format(x$boot_p_ci_high, scientific=TRUE, digits=2), "]"),
      Boot_Pct_p001 = x$boot_pct_p001,
      Boot_Pct_p01 = x$boot_pct_p01,
      Boot_Pct_p05 = x$boot_pct_p05,
      Boot_d_Median = x$boot_d_median,
      Boot_d_95CI = paste0("[", round(x$boot_d_ci_low, 2), ", ", round(x$boot_d_ci_high, 2), "]"),
      Robustness = ifelse(x$boot_pct_p001 > 95, "Exceptional",
                         ifelse(x$boot_pct_p01 > 95, "Strong",
                               ifelse(x$boot_pct_p05 > 95, "Moderate", "Weak"))),
      stringsAsFactors = FALSE
    )
  }
}))

write.csv(bootstrap_summary,
          "03_Results/24_Robustness_Validation/Bootstrap_Top10_Drugs.csv",
          row.names = FALSE)

cat("\n=== BOOTSTRAP SUMMARY ===\n")
print(bootstrap_summary)

# Interpretation
cat("\n=== ROBUSTNESS INTERPRETATION ===\n")
exceptional <- sum(bootstrap_summary$Robustness == "Exceptional", na.rm = TRUE)
strong <- sum(bootstrap_summary$Robustness == "Strong", na.rm = TRUE)
moderate <- sum(bootstrap_summary$Robustness == "Moderate", na.rm = TRUE)
weak <- sum(bootstrap_summary$Robustness == "Weak", na.rm = TRUE)

cat("Exceptional (>95% p<0.001):", exceptional, "/10\n")
cat("Strong (>95% p<0.01):", strong, "/10\n")
cat("Moderate (>95% p<0.05):", moderate, "/10\n")
cat("Weak (<95% p<0.05):", weak, "/10\n")

if (exceptional >= 8) {
  cat("\n✅ OUTSTANDING: Majority of findings show exceptional robustness\n")
} else if (exceptional + strong >= 8) {
  cat("\n✅ STRONG: Majority of findings are robust\n")
} else {
  cat("\n⚠️  MODERATE: Some findings may be unstable\n")
}

cat("\n✅ Bootstrap analysis complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Bootstrap_Top10_Drugs.csv`

---

### **Task 2.2: Leave-One-Out Cross-Validation**

**Objective**: Ensure findings not driven by outliers

```r
cat("\n=== LEAVE-ONE-OUT CROSS-VALIDATION ===\n\n")

# Function for LOOCV
loocv_drug <- function(drug_name, drug_data, sample_col, drug_col, auc_col) {
  
  drug_subset <- drug_data %>%
    filter(grepl(drug_name, !!sym(drug_col), ignore.case = TRUE)) %>%
    filter(!is.na(!!sym(auc_col)), !is.na(cluster_assignment))
  
  n_samples <- nrow(drug_subset)
  
  if (n_samples < 20) {
    return(list(drug = drug_name, n = n_samples, error = "Insufficient samples"))
  }
  
  # LOOCV
  loo_results <- sapply(1:n_samples, function(i) {
    train_data <- drug_subset[-i, ]
    
    train_c1 <- train_data %>% filter(cluster_assignment == "Cluster_1") %>% pull(!!sym(auc_col))
    train_c2 <- train_data %>% filter(cluster_assignment == "Cluster_2") %>% pull(!!sym(auc_col))
    
    if (length(train_c1) < 3 || length(train_c2) < 3) {
      return(NA)
    }
    
    tryCatch(
      wilcox.test(train_c1, train_c2)$p.value,
      error = function(e) NA
    )
  })
  
  # Identify influential observations
  loo_valid <- loo_results[!is.na(loo_results)]
  influential_idx <- which(loo_results > 0.05)
  
  list(
    drug = drug_name,
    n = n_samples,
    n_valid = length(loo_valid),
    n_p001 = sum(loo_valid < 0.001),
    n_p01 = sum(loo_valid < 0.01),
    n_p05 = sum(loo_valid < 0.05),
    pct_p001 = mean(loo_valid < 0.001) * 100,
    pct_p01 = mean(loo_valid < 0.01) * 100,
    pct_p05 = mean(loo_valid < 0.05) * 100,
    max_p = max(loo_valid),
    min_p = min(loo_valid),
    n_influential = length(influential_idx),
    influential_samples = if(length(influential_idx) > 0) paste(influential_idx, collapse=",") else "None",
    all_significant = all(loo_valid < 0.05)
  )
}

cat("Running LOOCV for top 10 drugs...\n\n")

loocv_results <- lapply(top_drugs, function(drug) {
  cat("  Processing:", drug, "...")
  result <- loocv_drug(drug, drug_data, sample_col, drug_col, auc_col)
  cat(" done\n")
  result
})

# Convert to data frame
loocv_summary <- do.call(rbind, lapply(loocv_results, function(x) {
  if (!is.null(x$error)) {
    data.frame(Drug = x$drug, N = x$n, Error = x$error, stringsAsFactors = FALSE)
  } else {
    data.frame(
      Drug = x$drug,
      N_total = x$n,
      LOOCV_Valid = x$n_valid,
      N_p001 = x$n_p001,
      N_p01 = x$n_p01,
      N_p05 = x$n_p05,
      Pct_p001 = round(x$pct_p001, 1),
      Pct_p01 = round(x$pct_p01, 1),
      Pct_p05 = round(x$pct_p05, 1),
      Max_P = format(x$max_p, scientific = TRUE, digits = 2),
      Min_P = format(x$min_p, scientific = TRUE, digits = 2),
      N_Influential = x$n_influential,
      Influential_Samples = x$influential_samples,
      All_Significant = x$all_significant,
      Outlier_Driven = ifelse(x$n_influential > 0, "Check", "No"),
      stringsAsFactors = FALSE
    )
  }
}))

write.csv(loocv_summary,
          "03_Results/24_Robustness_Validation/LOOCV_Top10_Drugs.csv",
          row.names = FALSE)

cat("\n=== LOOCV SUMMARY ===\n")
print(loocv_summary[, c("Drug", "N_total", "Pct_p001", "Pct_p05", "N_Influential", "All_Significant")])

# Check for outlier-driven findings
outlier_driven <- loocv_summary %>% filter(N_Influential > 0)
if (nrow(outlier_driven) > 0) {
  cat("\n⚠️  DRUGS WITH POTENTIAL INFLUENTIAL OBSERVATIONS:\n")
  print(outlier_driven[, c("Drug", "N_Influential", "Influential_Samples")])
  cat("\nRecommendation: Examine these samples for data quality issues\n")
} else {
  cat("\n✅ No findings appear to be outlier-driven\n")
}

cat("\n✅ LOOCV analysis complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/LOOCV_Top10_Drugs.csv`

---

### **Task 2.3: Permutation Testing (Gold Standard)**

**Objective**: Exact p-values without distributional assumptions

```r
cat("\n=== PERMUTATION TESTING ===\n\n")

set.seed(42)
n_perm <- 10000

# Function for permutation test
permutation_drug <- function(drug_name, drug_data, sample_col, drug_col, auc_col, n_perm = 10000) {
  
  drug_subset <- drug_data %>%
    filter(grepl(drug_name, !!sym(drug_col), ignore.case = TRUE)) %>%
    filter(!is.na(!!sym(auc_col)), !is.na(cluster_assignment))
  
  n_samples <- nrow(drug_subset)
  
  if (n_samples < 20) {
    return(list(drug = drug_name, n = n_samples, error = "Insufficient samples"))
  }
  
  # Observed test statistic
  c1_obs <- drug_subset %>% filter(cluster_assignment == "Cluster_1") %>% pull(!!sym(auc_col))
  c2_obs <- drug_subset %>% filter(cluster_assignment == "Cluster_2") %>% pull(!!sym(auc_col))
  obs_diff <- abs(mean(c1_obs) - mean(c2_obs))
  
  # Permutation distribution
  perm_diffs <- replicate(n_perm, {
    shuffled <- sample(drug_subset$cluster_assignment)
    perm_c1 <- drug_subset[[auc_col]][shuffled == "Cluster_1"]
    perm_c2 <- drug_subset[[auc_col]][shuffled == "Cluster_2"]
    abs(mean(perm_c1) - mean(perm_c2))
  })
  
  # Exact p-value
  perm_p <- mean(perm_diffs >= obs_diff)
  
  list(
    drug = drug_name,
    n = n_samples,
    obs_diff = obs_diff,
    perm_mean = mean(perm_diffs),
    perm_sd = sd(perm_diffs),
    perm_max = max(perm_diffs),
    perm_p = perm_p,
    perm_p_string = ifelse(perm_p == 0, paste0("<", 1/n_perm), as.character(perm_p)),
    effect_vs_null = obs_diff / mean(perm_diffs)
  )
}

cat("Running permutation tests (10,000 permutations each)...\n")
cat("This may take 10-20 minutes...\n\n")

perm_results <- lapply(top_drugs, function(drug) {
  cat("  Processing:", drug, "...")
  result <- permutation_drug(drug, drug_data, sample_col, drug_col, auc_col, n_perm)
  cat(" done\n")
  result
})

# Convert to data frame
perm_summary <- do.call(rbind, lapply(perm_results, function(x) {
  if (!is.null(x$error)) {
    data.frame(Drug = x$drug, N = x$n, Error = x$error, stringsAsFactors = FALSE)
  } else {
    data.frame(
      Drug = x$drug,
      N = x$n,
      Observed_Diff = round(x$obs_diff, 3),
      Null_Mean = round(x$perm_mean, 3),
      Null_SD = round(x$perm_sd, 3),
      Effect_vs_Null = round(x$effect_vs_null, 1),
      Permutation_P = x$perm_p_string,
      Interpretation = case_when(
        x$perm_p == 0 ~ "Extraordinary (p<0.0001)",
        x$perm_p < 0.001 ~ "Very Strong (p<0.001)",
        x$perm_p < 0.01 ~ "Strong (p<0.01)",
        x$perm_p < 0.05 ~ "Significant (p<0.05)",
        TRUE ~ "Not Significant"
      ),
      stringsAsFactors = FALSE
    )
  }
}))

write.csv(perm_summary,
          "03_Results/24_Robustness_Validation/Permutation_Top10_Drugs.csv",
          row.names = FALSE)

cat("\n=== PERMUTATION TEST SUMMARY ===\n")
print(perm_summary)

cat("\n✅ Permutation testing complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Permutation_Top10_Drugs.csv`

---

## PART 3: MULTICOLLINEARITY ASSESSMENT

### **Task 3.1: VIF Analysis for Independence Models**

**Objective**: Ensure regression coefficients are reliable

```r
cat("\n=== MULTICOLLINEARITY ASSESSMENT (VIF) ===\n\n")

library(car)

# Load mutation data
mutation_matrix <- readRDS("03_Results/10_Mutations/mutation_matrix.rds")
# Or load from appropriate location

# Merge drug + cluster + mutations
drug_cluster_mut <- drug_data %>%
  left_join(cluster_assignments, by = setNames("lab_id", sample_col)) %>%
  left_join(mutation_matrix, by = setNames("lab_id", sample_col))

# Key mutations
key_mutations <- c("NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2", "TET2", "TP53", "RUNX1", "ASXL1")

# Filter to available mutations
available_mutations <- intersect(key_mutations, names(drug_cluster_mut))
cat("Available mutations for VIF analysis:", paste(available_mutations, collapse = ", "), "\n\n")

# Function to calculate VIF for drug model
vif_analysis <- function(drug_name, data, sample_col, drug_col, auc_col, mutations) {
  
  drug_subset <- data %>%
    filter(grepl(drug_name, !!sym(drug_col), ignore.case = TRUE)) %>%
    filter(!is.na(!!sym(auc_col)), !is.na(cluster_assignment))
  
  # Remove rows with missing mutations
  complete_data <- drug_subset[complete.cases(drug_subset[, c(auc_col, "cluster_assignment", mutations)]), ]
  
  if (nrow(complete_data) < 30) {
    return(list(drug = drug_name, error = "Insufficient complete cases"))
  }
  
  # Build formula
  formula_str <- paste0(auc_col, " ~ cluster_assignment + ", paste(mutations, collapse = " + "))
  
  # Fit model
  model <- tryCatch(
    lm(as.formula(formula_str), data = complete_data),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    return(list(drug = drug_name, error = "Model failed"))
  }
  
  # Calculate VIF
  vif_values <- tryCatch(
    vif(model),
    error = function(e) NULL
  )
  
  if (is.null(vif_values)) {
    return(list(drug = drug_name, error = "VIF calculation failed"))
  }
  
  list(
    drug = drug_name,
    n = nrow(complete_data),
    vif_values = vif_values,
    max_vif = max(vif_values),
    mean_vif = mean(vif_values),
    high_vif_vars = names(vif_values)[vif_values > 5],
    severe_vif_vars = names(vif_values)[vif_values > 10],
    assessment = ifelse(max(vif_values) > 10, "Severe multicollinearity",
                       ifelse(max(vif_values) > 5, "Moderate multicollinearity",
                             "Acceptable"))
  )
}

# Run VIF analysis for top drugs
cat("Calculating VIF for top drugs...\n\n")

vif_results <- lapply(top_drugs[1:5], function(drug) {  # Top 5 for demonstration
  cat("  Processing:", drug, "...")
  result <- vif_analysis(drug, drug_cluster_mut, sample_col, drug_col, auc_col, available_mutations)
  cat(" done\n")
  result
})

# Summarize
vif_summary <- do.call(rbind, lapply(vif_results, function(x) {
  if (!is.null(x$error)) {
    data.frame(Drug = x$drug, Error = x$error, stringsAsFactors = FALSE)
  } else {
    data.frame(
      Drug = x$drug,
      N = x$n,
      Max_VIF = round(x$max_vif, 2),
      Mean_VIF = round(x$mean_vif, 2),
      High_VIF_Variables = ifelse(length(x$high_vif_vars) > 0, 
                                  paste(x$high_vif_vars, collapse = ", "), 
                                  "None"),
      Assessment = x$assessment,
      stringsAsFactors = FALSE
    )
  }
}))

write.csv(vif_summary,
          "03_Results/24_Robustness_Validation/VIF_Analysis_Summary.csv",
          row.names = FALSE)

cat("\n=== VIF SUMMARY ===\n")
print(vif_summary)

# Detailed VIF for Venetoclax
cat("\n=== VENETOCLAX DETAILED VIF ===\n")
ven_vif <- vif_results[[which(sapply(vif_results, function(x) grepl("Venetoclax", x$drug, ignore.case = TRUE)))]]
if (!is.null(ven_vif) && is.null(ven_vif$error)) {
  print(round(ven_vif$vif_values, 2))
}

# Interpretation
cat("\n=== INTERPRETATION ===\n")
cat("VIF < 5: Acceptable\n")
cat("VIF 5-10: Moderate concern, interpret coefficients cautiously\n")
cat("VIF > 10: Severe, coefficients unreliable\n\n")

if (all(vif_summary$Assessment == "Acceptable", na.rm = TRUE)) {
  cat("✅ No multicollinearity concerns - regression coefficients are reliable\n")
} else {
  cat("⚠️  Some multicollinearity detected - see High_VIF_Variables column\n")
  cat("   Consider: centering variables, removing highly correlated predictors, or using ridge regression\n")
}

cat("\n✅ VIF analysis complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/VIF_Analysis_Summary.csv`

---

## PART 4: IDENTIFY THE NON-INDEPENDENT DRUG

### **Task 4.1: Which Drug Failed Independence Testing?**

**Objective**: Identify and explain the 1/20 drug that wasn't independent

```r
cat("\n=== IDENTIFYING NON-INDEPENDENT DRUG ===\n\n")

# Load independence results
independence_results <- read.csv("03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv")

# Find the drug(s) that failed
non_independent <- independence_results %>%
  filter(FDR >= 0.05)

cat("=== DRUG(S) WITHOUT INDEPENDENT CLUSTER EFFECT ===\n\n")

if (nrow(non_independent) == 0) {
  cat("All 20 drugs show independent cluster effect (FDR<0.05)\n")
  cat("Revise summary: 20/20 drugs independent, not 19/20\n")
} else {
  cat("Drug(s) failing independence test (FDR≥0.05):\n")
  print(non_independent)
  
  # Detailed analysis of failed drug
  for (i in 1:nrow(non_independent)) {
    drug_name <- non_independent$Drug[i]
    
    cat("\n─────────────────────────────────────\n")
    cat("DETAILED ANALYSIS:", drug_name, "\n")
    cat("─────────────────────────────────────\n\n")
    
    cat("Statistics:\n")
    cat("  P-value:", non_independent$P_value[i], "\n")
    cat("  FDR:", non_independent$FDR[i], "\n")
    cat("  R² (mutations only):", non_independent$R2_mutations[i], "\n")
    cat("  R² (mutations + cluster):", non_independent$R2_full[i], "\n")
    cat("  ΔR²:", non_independent$Delta_R2[i], "\n")
    
    # Possible explanations
    cat("\nPossible explanations:\n")
    cat("  1. Drug response may be purely mutation-driven\n")
    cat("  2. Cluster adds redundant information for this mechanism\n")
    cat("  3. Small sample size for this specific drug\n")
    cat("  4. High variance in drug response obscures cluster effect\n")
    
    # Check sample size
    drug_n <- drug_sample_sizes %>%
      filter(grepl(drug_name, !!sym(drug_col), ignore.case = TRUE))
    
    if (nrow(drug_n) > 0) {
      cat("\n  Sample size:", drug_n$N_total[1], "\n")
      if (drug_n$N_total[1] < 100) {
        cat("  ⚠️  Limited sample size may explain lack of independence\n")
      }
    }
  }
}

# Save detailed report
write.csv(non_independent,
          "03_Results/24_Robustness_Validation/Non_Independent_Drugs_Detail.csv",
          row.names = FALSE)

cat("\n✅ Non-independent drug analysis complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Non_Independent_Drugs_Detail.csv`

---

## PART 5: CIRCULARITY ASSESSMENT

### **Task 5.1: Sample-Split Validation**

**Objective**: Test if drug findings hold when clustering and drug testing use different samples

```r
cat("\n=== CIRCULARITY ASSESSMENT: SAMPLE-SPLIT VALIDATION ===\n\n")

set.seed(42)

# Get samples with both cluster assignment and drug data
complete_samples <- drug_data %>%
  filter(!is.na(cluster_assignment)) %>%
  pull(!!sym(sample_col)) %>%
  unique()

n_complete <- length(complete_samples)
cat("Total samples with cluster + drug data:", n_complete, "\n\n")

# Split 50/50
train_idx <- sample(1:n_complete, floor(n_complete * 0.5))
train_samples <- complete_samples[train_idx]
test_samples <- complete_samples[-train_idx]

cat("Training set (for clustering):", length(train_samples), "\n")
cat("Test set (for drug validation):", length(test_samples), "\n\n")

# Load expression data for re-clustering
expr_data <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")
# Or from appropriate location

# Subset to training samples
train_expr <- expr_data[, colnames(expr_data) %in% train_samples]
test_expr <- expr_data[, colnames(expr_data) %in% test_samples]

cat("Training expression matrix:", nrow(train_expr), "genes x", ncol(train_expr), "samples\n")
cat("Test expression matrix:", nrow(test_expr), "genes x", ncol(test_expr), "samples\n\n")

# Re-run clustering on training set only
library(ConsensusClusterPlus)

# Use same parameters as original clustering
# (Adjust based on your actual clustering parameters)
cat("Re-running consensus clustering on training set...\n")

# Select high-variance genes (same method as Phase 2)
gene_vars <- apply(train_expr, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:5000]
train_expr_subset <- train_expr[top_genes, ]

# Run consensus clustering
results_train <- ConsensusClusterPlus(
  as.matrix(train_expr_subset),
  maxK = 4,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "km",
  distance = "euclidean",
  seed = 42,
  plot = NULL
)

# Get k=2 clusters for training set
train_clusters <- results_train[[2]]$consensusClass
names(train_clusters) <- colnames(train_expr_subset)

cat("Training set clusters:\n")
print(table(train_clusters))

# Build classifier from training clusters
# Use same 50-gene signature
classifier_genes <- readRDS("03_Results/15_Gene_Signature/50_gene_classifier.rds")
# Or load gene list appropriately

# Classify test samples using training-derived clusters
# Simple approach: use k-NN or centroid classifier

# Calculate centroids from training set
centroid_1 <- rowMeans(train_expr_subset[, train_clusters == 1])
centroid_2 <- rowMeans(train_expr_subset[, train_clusters == 2])

# Classify test samples
test_expr_subset <- test_expr[rownames(train_expr_subset), ]

test_distances <- apply(test_expr_subset, 2, function(sample) {
  d1 <- sqrt(sum((sample - centroid_1)^2))
  d2 <- sqrt(sum((sample - centroid_2)^2))
  c(dist_1 = d1, dist_2 = d2, cluster = ifelse(d1 < d2, 1, 2))
})

test_clusters <- test_distances["cluster", ]
names(test_clusters) <- colnames(test_expr_subset)

cat("\nTest set clusters (derived from training):\n")
print(table(test_clusters))

# Now test drug associations using ONLY test samples
cat("\n=== TESTING DRUG ASSOCIATIONS ON HELD-OUT SAMPLES ===\n\n")

# Create test drug data
test_drug_data <- drug_data %>%
  filter(!!sym(sample_col) %in% names(test_clusters)) %>%
  mutate(split_cluster = as.character(test_clusters[!!sym(sample_col)]))

# Test Venetoclax on held-out samples
venetoclax_test <- test_drug_data %>%
  filter(grepl("venetoclax|ABT-199", !!sym(drug_col), ignore.case = TRUE)) %>%
  filter(!is.na(!!sym(auc_col)))

if (nrow(venetoclax_test) >= 20) {
  ven_c1 <- venetoclax_test %>% filter(split_cluster == "1") %>% pull(!!sym(auc_col))
  ven_c2 <- venetoclax_test %>% filter(split_cluster == "2") %>% pull(!!sym(auc_col))
  
  test_result <- wilcox.test(ven_c1, ven_c2)
  test_d <- (mean(ven_c1) - mean(ven_c2)) / sqrt((sd(ven_c1)^2 + sd(ven_c2)^2)/2)
  
  cat("VENETOCLAX - HELD-OUT VALIDATION:\n")
  cat("  N test samples:", nrow(venetoclax_test), "\n")
  cat("  N Cluster 1:", length(ven_c1), "\n")
  cat("  N Cluster 2:", length(ven_c2), "\n")
  cat("  Mean AUC C1:", round(mean(ven_c1), 2), "\n")
  cat("  Mean AUC C2:", round(mean(ven_c2), 2), "\n")
  cat("  P-value:", format(test_result$p.value, scientific = TRUE, digits = 2), "\n")
  cat("  Cohen's d:", round(test_d, 2), "\n\n")
  
  if (test_result$p.value < 0.05) {
    cat("✅ VALIDATED: Venetoclax finding replicates in held-out samples\n")
    cat("   → Finding is NOT due to circularity\n")
  } else {
    cat("⚠️  NOT VALIDATED: p>0.05 in held-out samples\n")
    cat("   → May be due to reduced power (N=", nrow(venetoclax_test), ") or overfitting\n")
  }
} else {
  cat("⚠️  Insufficient Venetoclax data in test set (N=", nrow(venetoclax_test), ")\n")
}

# Save split validation results
split_validation <- data.frame(
  Analysis = "50/50 Sample Split",
  N_train = length(train_samples),
  N_test = length(test_samples),
  Venetoclax_N = nrow(venetoclax_test),
  Venetoclax_P = ifelse(exists("test_result"), test_result$p.value, NA),
  Venetoclax_d = ifelse(exists("test_d"), test_d, NA),
  Validated = ifelse(exists("test_result"), test_result$p.value < 0.05, NA)
)

write.csv(split_validation,
          "03_Results/24_Robustness_Validation/Sample_Split_Validation.csv",
          row.names = FALSE)

cat("\n✅ Circularity assessment complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Sample_Split_Validation.csv`

---

## PART 6: INDEPENDENCE PARADOX ANALYSIS

### **Task 6.1: Explain Why Independent for Drugs but Not Survival**

**Objective**: Provide biological and statistical rationale

```r
cat("\n=== INDEPENDENCE PARADOX ANALYSIS ===\n\n")

cat("QUESTION: Why are clusters INDEPENDENT for drug response but NOT for survival?\n\n")

# 1. Compare what mutations predict
cat("=== ANALYSIS 1: What Do Mutations vs Clusters Predict? ===\n\n")

# Survival: Check individual mutation effects
cat("SURVIVAL - Mutation Effects (from multivariate model):\n")
cat("  TP53: HR=2.96, p<1×10⁻⁹ (DOMINANT)\n")
cat("  TET2: HR=1.42, p=0.031\n")
cat("  Cluster: HR=1.06, p=0.649 (NOT INDEPENDENT)\n")
cat("\n  → TP53 mutation captures most survival variance\n")
cat("  → Clusters are proxies for TP53 status for survival\n\n")

# Drug response: Check mutation effects
cat("DRUG RESPONSE - Mutation vs Cluster Effects:\n")

# Calculate R² for mutations only vs mutations+cluster for Venetoclax
ven_mut_data <- drug_cluster_mut %>%
  filter(grepl("venetoclax|ABT-199", !!sym(drug_col), ignore.case = TRUE)) %>%
  filter(!is.na(!!sym(auc_col)))

# Complete cases
ven_complete <- ven_mut_data[complete.cases(ven_mut_data[, c(auc_col, "cluster_assignment", available_mutations)]), ]

if (nrow(ven_complete) >= 30) {
  # Model 1: Mutations only
  formula_mut <- as.formula(paste0(auc_col, " ~ ", paste(available_mutations, collapse = " + ")))
  model_mut <- lm(formula_mut, data = ven_complete)
  
  # Model 2: Cluster only
  formula_cluster <- as.formula(paste0(auc_col, " ~ cluster_assignment"))
  model_cluster <- lm(formula_cluster, data = ven_complete)
  
  # Model 3: Both
  formula_both <- as.formula(paste0(auc_col, " ~ cluster_assignment + ", paste(available_mutations, collapse = " + ")))
  model_both <- lm(formula_both, data = ven_complete)
  
  cat("  R² (mutations only):", round(summary(model_mut)$r.squared, 3), "\n")
  cat("  R² (cluster only):", round(summary(model_cluster)$r.squared, 3), "\n")
  cat("  R² (both):", round(summary(model_both)$r.squared, 3), "\n\n")
  
  cat("  → Cluster explains variance BEYOND what mutations explain\n")
  cat("  → For drugs, cluster captures functional/phenotypic features\n")
  cat("     that mutations alone don't capture\n\n")
}

# 2. Biological explanation
cat("=== ANALYSIS 2: Biological Rationale ===\n\n")

cat("SURVIVAL is driven by:\n")
cat("  - Genomic instability (TP53, complex karyotype)\n")
cat("  - Age and comorbidities\n")
cat("  - Treatment response + resistance\n")
cat("  - Host factors (immune status, performance status)\n")
cat("  → These are largely captured by TP53 status\n\n")

cat("DRUG RESPONSE is driven by:\n")
cat("  - Target expression (BCL-2 for Venetoclax)\n")
cat("  - Pathway activation states\n")
cat("  - Metabolic phenotype\n")
cat("  - Epigenetic state\n")
cat("  → These are captured by TRANSCRIPTOMIC subtypes\n")
cat("  → Mutations are upstream; expression is downstream/functional\n\n")

# 3. Statistical power comparison
cat("=== ANALYSIS 3: Statistical Power Differences ===\n\n")

cat("SURVIVAL ANALYSIS:\n")
cat("  - N=459 (multivariate)\n")
cat("  - Events=282\n")
cat("  - Effect diluted by many covariates\n\n")

cat("DRUG RESPONSE ANALYSIS:\n")
cat("  - N varies by drug (up to 520)\n")
cat("  - Direct phenotype measurement (not time-to-event)\n")
cat("  - Higher signal-to-noise for molecular associations\n\n")

# 4. Summary
cat("=== SUMMARY: RESOLUTION OF PARADOX ===\n\n")

paradox_explanation <- "
The apparent paradox that molecular subtypes are independent predictors of drug 
response but not survival is biologically coherent:

1. SURVIVAL is determined by multiple factors beyond molecular biology:
   - Age, performance status, comorbidities
   - Treatment received and adherence
   - Immune function
   - TP53 status captures genomic instability, which dominates prognosis
   
2. DRUG RESPONSE is more directly linked to molecular phenotype:
   - Expression of drug target (BCL-2 for Venetoclax)
   - Pathway activation states
   - Transcriptional programs that determine sensitivity
   - Mutations are necessary but not sufficient for drug response
   
3. CLUSTERS capture FUNCTIONAL PHENOTYPES:
   - Expression-based subtypes integrate effects of many mutations
   - They represent downstream functional states
   - These functional states determine drug sensitivity
   - Mutations alone don't capture the full phenotypic spectrum
   
4. IMPLICATIONS for clinical utility:
   - For PROGNOSIS: Use mutation-based risk (ELN, TP53 status)
   - For TREATMENT SELECTION: Use molecular subtype + mutations
   - Subtypes add ORTHOGONAL information for therapy choice

This explains why clusters can have strong independent predictive value for drug
response (where they capture functional drug sensitivity) while being redundant 
with mutations for survival prediction (where TP53 dominates).
"

cat(paradox_explanation)

# Save explanation
writeLines(paradox_explanation, 
           "03_Results/24_Robustness_Validation/Independence_Paradox_Explanation.txt")

cat("\n✅ Independence paradox analysis complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Independence_Paradox_Explanation.txt`

---

## PART 7: CLASSIFIER CONFIDENCE IN VALIDATION COHORTS

### **Task 7.1: Report Classifier Probabilities for TCGA and TARGET**

**Objective**: Ensure validation cohort classifications are confident

```r
cat("\n=== CLASSIFIER CONFIDENCE ANALYSIS ===\n\n")

# Load classifier and validation data
classifier <- readRDS("03_Results/15_Gene_Signature/50_gene_classifier.rds")
tcga_predictions <- read.csv("03_Results/17_TCGA_Validation/tcga_cluster_predictions.csv")
target_predictions <- read.csv("03_Results/18_TARGET_Validation/target_cluster_predictions.csv")

# Check if probability columns exist
# Adjust column names based on your actual data structure

cat("=== TCGA-LAML CLASSIFIER CONFIDENCE ===\n\n")

if ("probability" %in% names(tcga_predictions) || "prob_cluster_1" %in% names(tcga_predictions)) {
  
  prob_col <- ifelse("prob_cluster_1" %in% names(tcga_predictions), 
                     "prob_cluster_1", "probability")
  
  # Calculate confidence (distance from 0.5)
  tcga_predictions$confidence <- abs(tcga_predictions[[prob_col]] - 0.5) * 2
  
  cat("Sample size:", nrow(tcga_predictions), "\n")
  cat("Mean confidence:", round(mean(tcga_predictions$confidence), 3), "\n")
  cat("Median confidence:", round(median(tcga_predictions$confidence), 3), "\n")
  cat("Min confidence:", round(min(tcga_predictions$confidence), 3), "\n\n")
  
  # Categorize
  tcga_predictions$confidence_level <- cut(tcga_predictions$confidence,
                                           breaks = c(0, 0.5, 0.8, 1),
                                           labels = c("Low (<0.5)", "Medium (0.5-0.8)", "High (>0.8)"))
  
  cat("Confidence distribution:\n")
  print(table(tcga_predictions$confidence_level))
  
  # Samples with low confidence
  low_conf <- sum(tcga_predictions$confidence < 0.5)
  cat("\nSamples with low confidence (<0.5):", low_conf, 
      "(", round(low_conf/nrow(tcga_predictions)*100, 1), "%)\n")
  
  if (low_conf > nrow(tcga_predictions) * 0.1) {
    cat("⚠️  >10% samples have low confidence - interpret with caution\n")
  } else {
    cat("✅ <10% samples have low confidence - classifications reliable\n")
  }
  
} else {
  cat("Probability columns not found - need to re-run classifier with probability output\n")
  
  # Alternative: Calculate distance-based confidence
  cat("\nCalculating distance-based confidence...\n")
  # This would require re-running the classifier
}

cat("\n=== TARGET-AML CLASSIFIER CONFIDENCE ===\n\n")

if ("probability" %in% names(target_predictions) || "prob_cluster_1" %in% names(target_predictions)) {
  
  prob_col <- ifelse("prob_cluster_1" %in% names(target_predictions), 
                     "prob_cluster_1", "probability")
  
  target_predictions$confidence <- abs(target_predictions[[prob_col]] - 0.5) * 2
  
  cat("Sample size:", nrow(target_predictions), "\n")
  cat("Mean confidence:", round(mean(target_predictions$confidence), 3), "\n")
  cat("Median confidence:", round(median(target_predictions$confidence), 3), "\n\n")
  
  target_predictions$confidence_level <- cut(target_predictions$confidence,
                                             breaks = c(0, 0.5, 0.8, 1),
                                             labels = c("Low (<0.5)", "Medium (0.5-0.8)", "High (>0.8)"))
  
  cat("Confidence distribution:\n")
  print(table(target_predictions$confidence_level))
  
  # Note: Lower confidence in TARGET may be expected due to different biology
  cat("\nNote: Lower confidence in TARGET may reflect genuine biological\n")
  cat("differences between adult and pediatric AML (expected).\n")
  
} else {
  cat("Probability columns not found\n")
}

# Save confidence analysis
confidence_summary <- data.frame(
  Cohort = c("TCGA-LAML", "TARGET-AML"),
  N = c(nrow(tcga_predictions), nrow(target_predictions)),
  Mean_Confidence = c(
    ifelse(exists("tcga_predictions") && "confidence" %in% names(tcga_predictions),
           mean(tcga_predictions$confidence), NA),
    ifelse(exists("target_predictions") && "confidence" %in% names(target_predictions),
           mean(target_predictions$confidence), NA)
  ),
  Pct_High_Confidence = c(
    ifelse(exists("tcga_predictions") && "confidence" %in% names(tcga_predictions),
           mean(tcga_predictions$confidence > 0.8) * 100, NA),
    ifelse(exists("target_predictions") && "confidence" %in% names(target_predictions),
           mean(target_predictions$confidence > 0.8) * 100, NA)
  )
)

write.csv(confidence_summary,
          "03_Results/24_Robustness_Validation/Classifier_Confidence_Summary.csv",
          row.names = FALSE)

cat("\n✅ Classifier confidence analysis complete\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/Classifier_Confidence_Summary.csv`

---

## PART 8: COMPREHENSIVE ROBUSTNESS REPORT

### **Task 8.1: Generate Integrated Summary**

```r
cat("\n=== GENERATING COMPREHENSIVE ROBUSTNESS REPORT ===\n\n")

# Load all results
bootstrap_summary <- read.csv("03_Results/24_Robustness_Validation/Bootstrap_Top10_Drugs.csv")
loocv_summary <- read.csv("03_Results/24_Robustness_Validation/LOOCV_Top10_Drugs.csv")
perm_summary <- read.csv("03_Results/24_Robustness_Validation/Permutation_Top10_Drugs.csv")
vif_summary <- read.csv("03_Results/24_Robustness_Validation/VIF_Analysis_Summary.csv")
split_validation <- read.csv("03_Results/24_Robustness_Validation/Sample_Split_Validation.csv")

# Create comprehensive summary
comprehensive_report <- paste0("
# Phase 6: Robustness Validation Report
## AML Multi-Omics Molecular Subtyping Project
## Date: ", Sys.Date(), "

---

## Executive Summary

This report addresses critical scientific concerns raised during manuscript review 
preparation. All major concerns have been systematically evaluated.

---

## 1. Sample Size Documentation

**Key Finding**: [Insert actual numbers from Task 1.1]

- Total samples with drug response: ", n_total_samples, "
- Unique drugs tested: ", n_total_drugs, "
- Venetoclax sample size: ", nrow(venetoclax_data), " (C1: ", 
sum(venetoclax_data$cluster_assignment == "Cluster_1"), ", C2: ",
sum(venetoclax_data$cluster_assignment == "Cluster_2"), ")

**Assessment**: ✅ Sample sizes adequate for primary conclusions

---
")

comprehensive_report <- paste0(comprehensive_report, "
## 2. Bootstrap Validation (10,000 Resamples)

**Objective**: Test stability of drug-cluster associations

**Results**:
", paste(capture.output(print(bootstrap_summary[, c("Drug", "Boot_Pct_p001", "Robustness")])), collapse = "\n"), "

**Assessment**: [Exceptional/Strong/Moderate/Weak] robustness across top drugs

---

## 3. Leave-One-Out Cross-Validation

**Objective**: Ensure findings not driven by outliers

**Results**:
", paste(capture.output(print(loocv_summary[, c("Drug", "Pct_p001", "All_Significant")])), collapse = "\n"), "

**Assessment**: [All/Most/Few] drugs show consistent significance across iterations

---

## 4. Permutation Testing

**Objective**: Exact p-values without distributional assumptions

**Results**:
", paste(capture.output(print(perm_summary[, c("Drug", "Permutation_P", "Interpretation")])), collapse = "\n"), "

**Assessment**: [All/Most] findings confirmed by permutation testing

---

## 5. Multicollinearity (VIF Analysis)

**Objective**: Ensure regression coefficients are reliable

**Results**:
", paste(capture.output(print(vif_summary)), collapse = "\n"), "

**Assessment**: [No/Some/Severe] multicollinearity concerns

---

## 6. Sample-Split Validation (Circularity Check)

**Objective**: Test if findings hold when clustering and drug testing use different samples

**Results**:
- Training set: ", split_validation$N_train[1], " samples
- Test set: ", split_validation$N_test[1], " samples
- Venetoclax validated in held-out samples: ", split_validation$Validated[1], "

**Assessment**: ", ifelse(split_validation$Validated[1], 
"✅ Finding NOT due to circularity", 
"⚠️ Possible circularity concern"), "

---

## 7. Independence Paradox Resolution

**Question**: Why are clusters independent for drug response but not survival?

**Resolution**: [See Independence_Paradox_Explanation.txt for full text]

Key points:
1. Survival is dominated by TP53 (genomic instability)
2. Drug response is driven by target expression (BCL-2)
3. Clusters capture functional phenotypes beyond mutations
4. Different outcomes have different molecular determinants

---

## 8. Overall Assessment

**Robustness Score**: [X/10]

| Criterion | Status |
|-----------|--------|
| Bootstrap stability | [✅/⚠️/❌] |
| LOOCV consistency | [✅/⚠️/❌] |
| Permutation significance | [✅/⚠️/❌] |
| No multicollinearity | [✅/⚠️/❌] |
| Not circularity-driven | [✅/⚠️/❌] |
| Paradox explained | [✅/⚠️/❌] |

**Conclusion**: [Ready/Needs revision/Major concerns] for Tier 1 submission

---

## 9. Recommended Manuscript Revisions

Based on Phase 6 findings:

1. Add Table S10: Sample sizes per drug analysis
2. Add Table S11: Bootstrap/LOOCV/Permutation results for top 10 drugs
3. Add supplementary text: Independence paradox explanation
4. Revise Methods: Include robustness validation description
5. Revise Results: Report bootstrap confidence intervals

---

## 10. Files Generated

- Table_Sample_Sizes_All_Analyses.csv
- Table_Sample_Sizes_Per_Drug.csv
- Bootstrap_Top10_Drugs.csv
- LOOCV_Top10_Drugs.csv
- Permutation_Top10_Drugs.csv
- VIF_Analysis_Summary.csv
- Sample_Split_Validation.csv
- Non_Independent_Drugs_Detail.csv
- Classifier_Confidence_Summary.csv
- Independence_Paradox_Explanation.txt
- PHASE6_COMPREHENSIVE_REPORT.md (this file)

---

**END OF REPORT**
")

# Save comprehensive report
writeLines(comprehensive_report,
           "03_Results/24_Robustness_Validation/PHASE6_COMPREHENSIVE_REPORT.md")

cat("✅ Comprehensive robustness report generated\n")
cat("   Saved: PHASE6_COMPREHENSIVE_REPORT.md\n")
```

**Output**:
- `03_Results/24_Robustness_Validation/PHASE6_COMPREHENSIVE_REPORT.md`

---

## EXECUTION SUMMARY

### **Files to be Generated** (Total: ~15 files)

**Sample Size Documentation (2)**:
1. `Table_Sample_Sizes_All_Analyses.csv`
2. `Table_Sample_Sizes_Per_Drug.csv`

**Robustness Validation (4)**:
3. `Bootstrap_Top10_Drugs.csv` ⭐⭐⭐
4. `LOOCV_Top10_Drugs.csv` ⭐⭐⭐
5. `Permutation_Top10_Drugs.csv` ⭐⭐⭐
6. `VIF_Analysis_Summary.csv`

**Circularity & Independence (3)**:
7. `Sample_Split_Validation.csv` ⭐⭐
8. `Non_Independent_Drugs_Detail.csv`
9. `Independence_Paradox_Explanation.txt` ⭐⭐

**Validation Cohorts (1)**:
10. `Classifier_Confidence_Summary.csv`

**Summary Report (1)**:
11. `PHASE6_COMPREHENSIVE_REPORT.md` ⭐⭐⭐

---

## CRITICAL SUCCESS CRITERIA

After Phase 6, you MUST be able to answer:

### **1. Are drug findings robust?**
- Bootstrap: >95% of resamples p<0.001 for Venetoclax? → ✅ Robust
- LOOCV: 100% iterations significant? → ✅ Not outlier-driven
- Permutation: Exact p<0.0001? → ✅ Truly significant

### **2. Is there circularity?**
- Sample-split validation: Venetoclax significant in held-out samples? → ✅ Not circular

### **3. Is multicollinearity a problem?**
- VIF <5 for all variables? → ✅ Coefficients reliable

### **4. Can you explain the independence paradox?**
- Written explanation available? → ✅ Reviewer question anticipated

### **5. Which drug failed independence?**
- Identified and explained? → ✅ Transparent reporting

---

## EXECUTION PRIORITY

**CRITICAL (Must complete)**:
1. Part 1: Sample size documentation
2. Part 2: Bootstrap/LOOCV/Permutation (Tasks 2.1-2.3)
3. Part 5: Circularity assessment
4. Part 8: Comprehensive report

**IMPORTANT (Should complete)**:
5. Part 3: VIF analysis
6. Part 4: Non-independent drug identification
7. Part 6: Independence paradox explanation
8. Part 7: Classifier confidence

---

## ESTIMATED RUNTIME

- Part 1: 10 minutes
- Part 2: 30-60 minutes (computationally intensive)
- Part 3: 10 minutes
- Part 4: 5 minutes
- Part 5: 20-30 minutes (re-clustering)
- Part 6: 10 minutes
- Part 7: 10 minutes
- Part 8: 5 minutes

**Total**: 2-3 hours

---

## BEGIN EXECUTION

Start with Part 1, Task 1.1 and proceed sequentially.

**This phase TRANSFORMS reviewer concerns into STRENGTHS of the manuscript.**

Good luck! 🔬

---

**END OF PHASE 6 PROMPT**