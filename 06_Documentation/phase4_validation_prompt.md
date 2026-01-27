# Phase 4: Manuscript Preparation & Critical Issues Resolution
## Comprehensive Claude Code Execution Prompt

**Date**: October 12, 2025
**Purpose**: Address remaining critical issues and prepare publication-ready materials
**Priority**: HIGH - Required for manuscript submission

---

## MISSION STATEMENT

Phase 3 revealed critical findings that require manuscript reframing:
- ✅ Molecular subtypes are real and robust
- ✅ Strong prognostic associations in adults
- ❌ **NOT independent of mutations** (p=0.649)
- ❌ **Opposite effect in pediatric AML**

This phase will:
1. Address remaining statistical/technical concerns
2. Perform publication-required sensitivity analyses
3. Generate manuscript-ready figures and tables
4. Document findings transparently
5. Prepare supplementary materials

---

## PART 1: SAMPLE ATTRITION ANALYSIS

### **Background**
Sample sizes vary across analyses:
- Cox PH: 671 samples
- Multivariate: 459 samples (32% reduction!)
- Mutations: 522 samples
- Need to understand why

### **Task 1.1: Create Sample Flow Diagram**

**Objective**: Document every step of sample filtering

**Method**:
```r
library(dplyr)
library(DiagrammeR)

# Load all datasets
expr_data <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")
cluster_assignments <- read.csv("03_Results/02_Clustering/cluster_assignments.csv")
clinical_data <- readRDS("01_Data/BeatAML/cohort1_clinical.rds")
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutation_matrix <- readRDS("03_Results/10_Mutations/mutation_matrix.rds")

# Step-by-step tracking
step1 <- data.frame(
  Step = "1. Expression data",
  N = ncol(expr_data),
  Description = "Raw expression matrix"
)

step2 <- data.frame(
  Step = "2. After clustering",
  N = nrow(cluster_assignments),
  Description = "Samples with cluster assignments"
)

step3 <- data.frame(
  Step = "3. With clinical data",
  N = nrow(clinical_data),
  Description = "Samples with clinical information"
)

step4 <- survival_data %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS)) %>%
  nrow() %>%
  data.frame(
    Step = "4. With survival data",
    N = .,
    Description = "Complete survival outcomes"
  )

step5 <- survival_data %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS) & !is.na(cluster_assignment)) %>%
  nrow() %>%
  data.frame(
    Step = "5. Survival + Clusters",
    N = .,
    Description = "Complete for survival analysis"
  )

# For multivariate analysis
multivar_data <- survival_data %>%
  left_join(mutation_matrix, by="lab_id") %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS), !is.na(cluster_assignment))

# Check completeness for key mutations
multivar_complete <- multivar_data %>%
  filter(!is.na(NPM1), !is.na(TP53), !is.na(RUNX1), !is.na(DNMT3A), 
         !is.na(FLT3), !is.na(TET2), !is.na(ASXL1),
         !is.na(AGE), !is.na(SEX))

step6 <- data.frame(
  Step = "6. With mutation data",
  N = nrow(multivar_data),
  Description = "Expression + survival + mutations"
)

step7 <- data.frame(
  Step = "7. Complete for multivariate",
  N = nrow(multivar_complete),
  Description = "No missing values in key variables"
)

# Combine all steps
attrition_table <- bind_rows(step1, step2, step3, step4, step5, step6, step7)
attrition_table$N_excluded <- c(0, diff(attrition_table$N))
attrition_table$Percent_retained <- round(attrition_table$N / attrition_table$N[1] * 100, 1)

write.csv(attrition_table, "03_Results/21_Manuscript_Prep/sample_attrition_table.csv", row.names=FALSE)

# Print detailed breakdown
cat("\n=== SAMPLE ATTRITION ANALYSIS ===\n")
print(attrition_table)

# Identify specific reasons for exclusion
cat("\n=== REASONS FOR MULTIVARIATE EXCLUSION ===\n")

# Missing mutation data
missing_mutations <- multivar_data %>%
  summarise(
    Missing_NPM1 = sum(is.na(NPM1)),
    Missing_TP53 = sum(is.na(TP53)),
    Missing_RUNX1 = sum(is.na(RUNX1)),
    Missing_TET2 = sum(is.na(TET2)),
    Missing_any = sum(is.na(NPM1) | is.na(TP53) | is.na(RUNX1) | is.na(TET2))
  )

print(missing_mutations)

# Missing clinical data
missing_clinical <- multivar_data %>%
  summarise(
    Missing_AGE = sum(is.na(AGE)),
    Missing_SEX = sum(is.na(SEX))
  )

print(missing_clinical)

# Create visual flow diagram
flow_diagram <- grViz("
digraph sample_flow {
  graph [rankdir=TB, fontname=Arial]
  node [shape=box, style=filled, fillcolor=lightblue, fontname=Arial]
  
  start [label='Expression Data\\nn=707', fillcolor=lightgreen]
  cluster [label='With Cluster Assignment\\nn=707']
  survival [label='With Survival Data\\nn=671\\n(36 excluded)', fillcolor=lightyellow]
  mutation [label='With Mutation Data\\nn=615\\n(56 excluded)', fillcolor=lightyellow]
  complete [label='Complete for Multivariate\\nn=459\\n(156 excluded)', fillcolor=orange]
  
  start -> cluster [label='Clustering']
  cluster -> survival [label='Merge clinical']
  survival -> mutation [label='Merge mutations']
  mutation -> complete [label='Complete case\\nanalysis']
}
")

# Save diagram
flow_diagram %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_pdf("04_Figures/20_Manuscript_Prep/sample_attrition_flowchart.pdf")

cat("\n✅ Sample attrition analysis complete\n")
```

**Output Files**:
- `03_Results/21_Manuscript_Prep/sample_attrition_table.csv`
- `04_Figures/20_Manuscript_Prep/sample_attrition_flowchart.pdf`

---

### **Task 1.2: Sensitivity Analysis - Missing Data Handling**

**Objective**: Test if complete-case analysis biases results

**Method**:
```r
library(mice)

# Option 1: Multiple Imputation
set.seed(42)

# Prepare data for imputation
impute_data <- multivar_data %>%
  select(lab_id, OS_MONTHS, OS_STATUS, cluster_assignment, AGE, SEX,
         NPM1, TP53, RUNX1, DNMT3A, FLT3, TET2, ASXL1)

# Impute missing mutations (binary, so use logistic regression)
imputed <- mice(impute_data, m=5, method='logreg', seed=42, print=FALSE)

# Fit model on each imputed dataset
fits_imputed <- with(imputed, 
  coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + AGE + SEX + 
        TP53 + TET2 + RUNX1 + ASXL1)
)

# Pool results
pooled <- pool(fits_imputed)
summary_pooled <- summary(pooled)

# Compare to complete case
fit_complete <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + AGE + SEX + 
                      TP53 + TET2 + RUNX1 + ASXL1,
                      data = multivar_complete)

# Extract coefficients
comparison <- data.frame(
  Variable = rownames(summary(fit_complete)$coefficients),
  Complete_Case_HR = exp(coef(fit_complete)),
  Complete_Case_p = summary(fit_complete)$coefficients[,"Pr(>|z|)"],
  Imputed_HR = exp(summary_pooled$estimate),
  Imputed_p = summary_pooled$p.value
)

comparison$Difference <- abs(comparison$Complete_Case_HR - comparison$Imputed_HR)
comparison$Percent_Change <- abs(comparison$Difference / comparison$Complete_Case_HR * 100)

write.csv(comparison, "03_Results/21_Manuscript_Prep/imputation_sensitivity_analysis.csv", row.names=FALSE)

cat("\n=== IMPUTATION SENSITIVITY ANALYSIS ===\n")
print(comparison)

# Key finding: Does cluster p-value change?
cluster_row <- comparison[grepl("cluster", comparison$Variable, ignore.case=TRUE), ]
cat("\n📊 Cluster assignment p-value:\n")
cat("  Complete case: p =", cluster_row$Complete_Case_p, "\n")
cat("  With imputation: p =", cluster_row$Imputed_p, "\n")

if(abs(cluster_row$Complete_Case_p - cluster_row$Imputed_p) < 0.05) {
  cat("  ✅ Results are robust to missing data\n")
} else {
  cat("  ⚠️ Results sensitive to missing data handling\n")
}
```

**Output**: 
- `03_Results/21_Manuscript_Prep/imputation_sensitivity_analysis.csv`

---

## PART 2: TARGET-AML SENSITIVITY ANALYSES

### **Background**
- 8 of 50 genes missing in TARGET (16%)
- Imputed using BeatAML means
- Need to test robustness

### **Task 2.1: Classifier Performance with Only Available Genes**

**Objective**: Does classifier still work with only 42 genes?

**Method**:
```r
# Load TARGET data
target_expr <- readRDS("03_Results/18_TARGET_Validation/target_aml_expression_normalized.rds")
target_clinical <- read.csv("03_Results/18_TARGET_Validation/target_aml_clinical.csv")

# Load gene signature
gene_sig <- read.csv("03_Results/15_Gene_Signature/50_gene_signature_symbols.csv")
gene_sig_symbols <- read.csv("03_Results/15_Gene_Signature/50_gene_signature_symbols.csv")

# Identify available vs missing genes
available_genes <- intersect(gene_sig_symbols$Gene_Symbol, rownames(target_expr))
missing_genes <- setdiff(gene_sig_symbols$Gene_Symbol, rownames(target_expr))

cat("\n=== GENE AVAILABILITY IN TARGET ===\n")
cat("Available:", length(available_genes), "/", nrow(gene_sig_symbols), "\n")
cat("Missing:", length(missing_genes), "\n")
cat("Missing genes:", paste(missing_genes, collapse=", "), "\n")

# Load BeatAML classifier
rf_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")

# Prepare TARGET data - 3 scenarios
# Scenario 1: Mean imputation (current)
target_pred_1 <- as.data.frame(t(target_expr))

for(gene in missing_genes) {
  # Use BeatAML mean
  beatAML_expr <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")
  if(gene %in% rownames(beatAML_expr)) {
    target_pred_1[[gene]] <- mean(beatAML_expr[gene, ])
  } else {
    target_pred_1[[gene]] <- 0
  }
}

# Scenario 2: Zero imputation
target_pred_2 <- target_pred_1
for(gene in missing_genes) {
  target_pred_2[[gene]] <- 0
}

# Scenario 3: Median imputation
target_pred_3 <- target_pred_1
for(gene in missing_genes) {
  beatAML_expr <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")
  if(gene %in% rownames(beatAML_expr)) {
    target_pred_3[[gene]] <- median(beatAML_expr[gene, ])
  } else {
    target_pred_3[[gene]] <- 0
  }
}

# Scenario 4: Only 42 available genes (retrain-free prediction)
# Note: Can't use original classifier - would need different model
# Instead, use 42-gene subset with original weights
target_pred_4 <- target_pred_1[, available_genes]

# Predict with each scenario
pred_1 <- predict(rf_classifier, newdata=target_pred_1[, gene_sig_symbols$Gene_Symbol])
pred_2 <- predict(rf_classifier, newdata=target_pred_2[, gene_sig_symbols$Gene_Symbol])
pred_3 <- predict(rf_classifier, newdata=target_pred_3[, gene_sig_symbols$Gene_Symbol])

# Get probabilities
prob_1 <- predict(rf_classifier, newdata=target_pred_1[, gene_sig_symbols$Gene_Symbol], type="prob")
prob_2 <- predict(rf_classifier, newdata=target_pred_2[, gene_sig_symbols$Gene_Symbol], type="prob")
prob_3 <- predict(rf_classifier, newdata=target_pred_3[, gene_sig_symbols$Gene_Symbol], type="prob")

# Compare cluster proportions
scenario_comparison <- data.frame(
  Scenario = c("Mean imputation", "Zero imputation", "Median imputation"),
  Cluster1_n = c(sum(pred_1=="Cluster_1"), sum(pred_2=="Cluster_1"), sum(pred_3=="Cluster_1")),
  Cluster2_n = c(sum(pred_1=="Cluster_2"), sum(pred_2=="Cluster_2"), sum(pred_3=="Cluster_2")),
  Cluster1_pct = c(mean(pred_1=="Cluster_1")*100, mean(pred_2=="Cluster_1")*100, mean(pred_3=="Cluster_1")*100),
  Mean_confidence = c(mean(apply(prob_1, 1, max)), mean(apply(prob_2, 1, max)), mean(apply(prob_3, 1, max))),
  Low_confidence_pct = c(
    mean(apply(prob_1, 1, max) < 0.6)*100,
    mean(apply(prob_2, 1, max) < 0.6)*100,
    mean(apply(prob_3, 1, max) < 0.6)*100
  )
)

write.csv(scenario_comparison, "03_Results/21_Manuscript_Prep/target_imputation_sensitivity.csv", row.names=FALSE)

cat("\n=== IMPUTATION SENSITIVITY IN TARGET ===\n")
print(scenario_comparison)

# Test survival for each scenario
target_surv <- target_clinical %>%
  mutate(
    OS_MONTHS = case_when(
      !is.na(days_to_death) ~ days_to_death / 30.44,
      !is.na(days_to_last_follow_up) ~ days_to_last_follow_up / 30.44,
      TRUE ~ NA_real_
    ),
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))

# Merge with predictions
target_surv$cluster_mean <- pred_1[match(target_surv$case_id, colnames(target_expr))]
target_surv$cluster_zero <- pred_2[match(target_surv$case_id, colnames(target_expr))]
target_surv$cluster_median <- pred_3[match(target_surv$case_id, colnames(target_expr))]

# Cox models
cox_mean <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_mean, data=target_surv)
cox_zero <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_zero, data=target_surv)
cox_median <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_median, data=target_surv)

survival_comparison <- data.frame(
  Scenario = c("Mean imputation", "Zero imputation", "Median imputation"),
  HR = c(exp(coef(cox_mean)), exp(coef(cox_zero)), exp(coef(cox_median))),
  P_value = c(
    summary(cox_mean)$coefficients[,"Pr(>|z|)"],
    summary(cox_zero)$coefficients[,"Pr(>|z|)"],
    summary(cox_median)$coefficients[,"Pr(>|z|)"]
  )
)

write.csv(survival_comparison, "03_Results/21_Manuscript_Prep/target_survival_sensitivity.csv", row.names=FALSE)

cat("\n=== SURVIVAL SENSITIVITY IN TARGET ===\n")
print(survival_comparison)

# Agreement between scenarios
cat("\n=== AGREEMENT BETWEEN SCENARIOS ===\n")
cat("Mean vs Zero:", mean(pred_1 == pred_2) * 100, "% agreement\n")
cat("Mean vs Median:", mean(pred_1 == pred_3) * 100, "% agreement\n")
cat("Zero vs Median:", mean(pred_2 == pred_3) * 100, "% agreement\n")
```

**Output Files**:
- `03_Results/21_Manuscript_Prep/target_imputation_sensitivity.csv`
- `03_Results/21_Manuscript_Prep/target_survival_sensitivity.csv`

---

### **Task 2.2: Identify Which 8 Genes Are Missing**

**Method**:
```r
# Detailed missing gene analysis
gene_sig_full <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")

missing_gene_details <- gene_sig_full %>%
  filter(Gene %in% missing_genes) %>%
  arrange(desc(Importance))

# Add importance rank
missing_gene_details$Importance_Rank <- match(missing_gene_details$Gene, 
                                              gene_sig_full$Gene[order(-gene_sig_full$Importance)])

write.csv(missing_gene_details, "03_Results/21_Manuscript_Prep/missing_genes_in_target.csv", row.names=FALSE)

cat("\n=== MISSING GENES IN TARGET ===\n")
print(missing_gene_details[, c("Gene", "Importance", "Importance_Rank")])

# Are missing genes important?
median_importance_rank <- median(gene_sig_full$Importance_Rank)
cat("\nMedian importance rank of missing genes:", median(missing_gene_details$Importance_Rank), "\n")
cat("Median importance rank of all 50 genes:", median_importance_rank, "\n")

if(median(missing_gene_details$Importance_Rank) > median_importance_rank) {
  cat("✅ Missing genes are less important (higher rank = lower importance)\n")
} else {
  cat("⚠️ Missing genes include important predictors\n")
}
```

**Output**: 
- `03_Results/21_Manuscript_Prep/missing_genes_in_target.csv`

---

## PART 3: TCGA POWER ANALYSIS - IMPROVED VERSION

### **Task 3.1: Calculate Minimum Detectable Effect Size**

**Objective**: What HR could TCGA detect with 80% power?

**Method**:
```r
library(powerSurvEpi)

# TCGA parameters
tcga_n <- 151
tcga_events <- 97
tcga_c1 <- 76
tcga_c2 <- 75
tcga_event_rate <- tcga_events / tcga_n

# Test range of HRs
hr_range <- seq(1.1, 3.0, by=0.05)

power_results <- data.frame(
  HR = hr_range,
  Power = sapply(hr_range, function(hr) {
    tryCatch({
      powerCT.default(
        nE = tcga_c2,
        nC = tcga_c1,
        pE = tcga_event_rate,
        pC = tcga_event_rate,
        RR = hr,
        alpha = 0.05
      )
    }, error = function(e) NA)
  })
)

# Find HR for 80% power
power_80 <- power_results$HR[which(power_results$Power >= 0.80)[1]]

# Find power for observed HRs
beatAML_observed_hr <- 1.39
tcga_observed_hr <- 1.24

power_beatAML_hr <- power_results$Power[which.min(abs(power_results$HR - beatAML_observed_hr))]
power_tcga_hr <- power_results$Power[which.min(abs(power_results$HR - tcga_observed_hr))]

# Calculate required sample size for BeatAML effect
required_events <- function(hr, power, alpha=0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta <- qnorm(power)
  
  n <- (z_alpha + z_beta)^2 * 4 / (log(hr))^2
  return(ceiling(n))
}

required_n_80 <- required_events(beatAML_observed_hr, 0.80)

# Summary table
power_summary <- data.frame(
  Metric = c(
    "TCGA total n",
    "TCGA events",
    "TCGA event rate (%)",
    "BeatAML observed HR",
    "TCGA observed HR",
    "Minimum detectable HR (80% power)",
    "Power to detect BeatAML HR (%)",
    "Power to detect TCGA HR (%)",
    "Required events for 80% power",
    "TCGA has (% of required)"
  ),
  Value = c(
    tcga_n,
    tcga_events,
    round(tcga_event_rate * 100, 1),
    round(beatAML_observed_hr, 2),
    round(tcga_observed_hr, 2),
    round(power_80, 2),
    round(power_beatAML_hr * 100, 1),
    round(power_tcga_hr * 100, 1),
    required_n_80,
    round(tcga_events / required_n_80 * 100, 1)
  )
)

write.csv(power_summary, "03_Results/21_Manuscript_Prep/tcga_power_analysis_detailed.csv", row.names=FALSE)

cat("\n=== TCGA POWER ANALYSIS ===\n")
print(power_summary)

# Plot power curve
pdf("04_Figures/20_Manuscript_Prep/tcga_power_curve.pdf", width=8, height=6)
plot(power_results$HR, power_results$Power, type="l", lwd=2,
     xlab="Hazard Ratio (Cluster 2 vs Cluster 1)",
     ylab="Statistical Power",
     main="TCGA-LAML Power Analysis",
     ylim=c(0, 1))
abline(h=0.80, lty=2, col="red", lwd=1)
abline(v=beatAML_observed_hr, lty=2, col="blue", lwd=1)
abline(v=power_80, lty=2, col="darkgreen", lwd=1)

# Add points
points(beatAML_observed_hr, power_beatAML_hr, pch=19, col="blue", cex=1.5)
points(power_80, 0.80, pch=19, col="darkgreen", cex=1.5)

# Add labels
text(beatAML_observed_hr, power_beatAML_hr + 0.1, 
     sprintf("BeatAML HR=%.2f\nPower=%.1f%%", beatAML_observed_hr, power_beatAML_hr*100),
     col="blue")
text(power_80, 0.85, 
     sprintf("Min detectable HR=%.2f", power_80),
     col="darkgreen")

legend("bottomright", 
       legend=c("Power curve", "80% power threshold", "BeatAML observed HR", "Min detectable HR"),
       lty=c(1, 2, 2, 2), col=c("black", "red", "blue", "darkgreen"), lwd=2)
dev.off()

cat("\n📊 TCGA is underpowered to detect small effects\n")
cat("   Can only detect HR ≥", round(power_80, 2), "with 80% power\n")
cat("   BeatAML HR (", beatAML_observed_hr, ") has only", round(power_beatAML_hr*100), "% power\n")
```

**Output Files**:
- `03_Results/21_Manuscript_Prep/tcga_power_analysis_detailed.csv`
- `04_Figures/20_Manuscript_Prep/tcga_power_curve.pdf`

---

## PART 4: MULTIPLE TESTING CORRECTION FRAMEWORK

### **Background**
Multiple analyses performed without overall correction:
- 11 mutations tested
- 11 immune cells tested  
- 7 Cox models
- 5 clustering solutions
- Need transparent framework

### **Task 4.1: Catalog All Statistical Tests**

**Method**:
```r
# Create comprehensive catalog of all tests
all_tests <- data.frame(
  Analysis_Part = character(),
  Test_Type = character(),
  Comparison = character(),
  Raw_P_value = numeric(),
  Test_Category = character(),
  Priority = character()
)

# Part 1: Mutation enrichment tests
mutation_enrichment <- read.csv("03_Results/10_Mutations/mutation_enrichment_by_cluster.csv")
mutation_tests <- mutation_enrichment %>%
  mutate(
    Analysis_Part = "Mutation Enrichment",
    Test_Type = "Fisher's Exact Test",
    Comparison = paste(Gene, "enrichment by cluster"),
    Raw_P_value = P_Value,
    Test_Category = "Exploratory",
    Priority = "Primary"
  ) %>%
  select(Analysis_Part, Test_Type, Comparison, Raw_P_value, Test_Category, Priority)

all_tests <- rbind(all_tests, mutation_tests)

# Part 2: Immune cell enrichment
immune_enrichment <- read.csv("03_Results/16_Immune_Deconvolution/mcp_counter_cell_enrichment.csv")
immune_tests <- immune_enrichment %>%
  mutate(
    Analysis_Part = "Immune Deconvolution",
    Test_Type = "Wilcoxon Rank-Sum",
    Comparison = paste(Cell_Type, "enrichment by cluster"),
    Raw_P_value = P_Value,
    Test_Category = "Exploratory",
    Priority = "Secondary"
  ) %>%
  select(Analysis_Part, Test_Type, Comparison, Raw_P_value, Test_Category, Priority)

all_tests <- rbind(all_tests, immune_tests)

# Part 3: Survival analyses
survival_tests <- data.frame(
  Analysis_Part = c(
    "Survival - Stratified Cox",
    "Survival - Time-varying",
    "Survival - Landmark (6m)",
    "Survival - Landmark (12m)",
    "Survival - Landmark (18m)",
    "Survival - Landmark (24m)",
    "Survival - RMST (24m)",
    "Survival - RMST (60m)"
  ),
  Test_Type = c(rep("Log-rank / Cox", 6), rep("RMST", 2)),
  Comparison = "Cluster survival difference",
  Raw_P_value = c(0.00155, 0.058, 0.017, 0.011, 0.015, 0.018, 0.029, 0.007),  # From Phase 3 results
  Test_Category = "Confirmatory",
  Priority = "Primary"
)

all_tests <- rbind(all_tests, survival_tests)

# Part 4: Multivariate models
multivar_tests <- data.frame(
  Analysis_Part = "Multivariate Analysis",
  Test_Type = "Likelihood Ratio Test",
  Comparison = c(
    "Cluster vs Clinical only",
    "Cluster vs Mutations only",
    "Cluster in full model"
  ),
  Raw_P_value = c(0.052, NA, 0.649),
  Test_Category = "Confirmatory",
  Priority = "Primary"
)

all_tests <- rbind(all_tests, multivar_tests)

# Part 5: External validation
validation_tests <- data.frame(
  Analysis_Part = c("TCGA Validation", "TARGET Validation", "Meta-analysis"),
  Test_Type = c("Cox regression", "Cox regression", "Random effects"),
  Comparison = "Cluster survival difference",
  Raw_P_value = c(0.353, 0.052, 0.001),
  Test_Category = "Confirmatory",
  Priority = "Primary"
)

all_tests <- rbind(all_tests, validation_tests)

# Add FDR correction WITHIN each analysis part
all_tests <- all_tests %>%
  group_by(Analysis_Part) %>%
  mutate(
    FDR_within_part = p.adjust(Raw_P_value, method="BH"),
    Bonferroni_within_part = p.adjust(Raw_P_value, method="bonferroni")
  ) %>%
  ungroup()

# Add overall study-wide correction for PRIMARY tests only
primary_tests <- all_tests %>% filter(Priority == "Primary")
primary_tests$FDR_studywide <- p.adjust(primary_tests$Raw_P_value, method="BH")
primary_tests$Bonferroni_studywide <- p.adjust(primary_tests$Raw_P_value, method="bonferroni")

# Merge back
all_tests <- all_tests %>%
  left_join(
    primary_tests %>% select(Comparison, FDR_studywide, Bonferroni_studywide),
    by="Comparison"
  )

# Flag significant results
all_tests <- all_tests %>%
  mutate(
    Sig_raw = Raw_P_value < 0.05,
    Sig_FDR_within = FDR_within_part < 0.05,
    Sig_FDR_studywide = FDR_studywide < 0.05 | is.na(FDR_studywide),
    Interpretation = case_when(
      is.na(Raw_P_value) ~ "Not tested",
      Sig_FDR_studywide ~ "Significant (study-wide FDR<0.05)",
      Sig_FDR_within ~ "Significant (within-analysis FDR<0.05)",
      Sig_raw ~ "Nominally significant (p<0.05)",
      TRUE ~ "Not significant"
    )
  )

write.csv(all_tests, "03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv", row.names=FALSE)

# Summary statistics
cat("\n=== MULTIPLE TESTING SUMMARY ===\n")
cat("Total tests performed:", nrow(all_tests), "\n")
cat("Primary confirmatory tests:", sum(all_tests$Priority == "Primary"), "\n")
cat("Exploratory tests:", sum(all_tests$Priority != "Primary"), "\n\n")

cat("Significant at raw p<0.05:", sum(all_tests$Sig_raw, na.rm=TRUE), "\n")
cat("Significant at FDR<0.05 (within analysis):", sum(all_tests$Sig_FDR_within, na.rm=TRUE), "\n")
cat("Significant at study-wide FDR<0.05:", sum(all_tests$Sig_FDR_studywide, na.rm=TRUE), "\n")

# Key findings that survive correction
cat("\n=== KEY FINDINGS (Study-wide FDR < 0.05) ===\n")
key_findings <- all_tests %>%
  filter(Sig_FDR_studywide) %>%
  arrange(FDR_studywide) %>%
  select(Comparison, Raw_P_value, FDR_studywide)

print(key_findings)
```

**Output**:
- `03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv`

---

## PART 5: MUTATION-CLUSTER INTERACTION ANALYSIS

### **Background**
Test if mutations have DIFFERENT prognostic effects in different clusters

### **Task 5.1: Test Interaction Terms**

**Method**:
```r
# Load data
multivar_complete <- read.csv("03_Results/11_Survival_Reanalysis/05_multivariate_complete_data.csv")

# Test key mutation × cluster interactions
mutations_to_test <- c("NPM1", "TP53", "RUNX1", "DNMT3A", "FLT3", "TET2", "ASXL1")

interaction_results <- lapply(mutations_to_test, function(mut) {
  
  formula_main <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + ", 
                                    mut, " + AGE + SEX"))
  formula_interact <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment * ", 
                                        mut, " + AGE + SEX"))
  
  # Fit both models
  fit_main <- coxph(formula_main, data=multivar_complete)
  fit_interact <- coxph(formula_interact, data=multivar_complete)
  
  # Likelihood ratio test
  lrt <- anova(fit_main, fit_interact)
  
  # Extract interaction coefficient
  interact_term <- grep(":", names(coef(fit_interact)), value=TRUE)
  interact_coef <- coef(fit_interact)[interact_term]
  interact_se <- sqrt(vcov(fit_interact)[interact_term, interact_term])
  interact_p <- summary(fit_interact)$coefficients[interact_term, "Pr(>|z|)"]
  
  # HR in each cluster
  # Cluster 1 (reference)
  hr_c1 <- exp(coef(fit_interact)[mut])
  
  # Cluster 2 (reference + interaction)
  hr_c2 <- exp(coef(fit_interact)[mut] + interact_coef)
  
  data.frame(
    Mutation = mut,
    LRT_ChiSq = lrt$Chisq[2],
    LRT_p = lrt$`Pr(>|Chi|)`[2],
    Interaction_HR = exp(interact_coef),
    Interaction_p = interact_p,
    HR_in_Cluster1 = hr_c1,
    HR_in_Cluster2 = hr_c2,
    HR_ratio = hr_c2 / hr_c1
  )
  
}) %>% bind_rows()

# FDR correction
interaction_results$FDR <- p.adjust(interaction_results$LRT_p, method="BH")

# Sort by significance
interaction_results <- interaction_results %>% arrange(LRT_p)

write.csv(interaction_results, "03_Results/21_Manuscript_Prep/mutation_cluster_interactions.csv", row.names=FALSE)

cat("\n=== MUTATION × CLUSTER INTERACTIONS ===\n")
print(interaction_results)

# Identify significant interactions
sig_interactions <- interaction_results %>% filter(FDR < 0.1)

if(nrow(sig_interactions) > 0) {
  cat("\n✅ SIGNIFICANT INTERACTIONS FOUND (FDR < 0.1):\n")
  print(sig_interactions[, c("Mutation", "LRT_p", "FDR", "HR_in_Cluster1", "HR_in_Cluster2")])
  cat("\n🎯 This suggests clusters MODIFY mutation effects - independent value!\n")
} else {
  cat("\n❌ No significant interactions found\n")
  cat("   Mutation effects are similar across clusters\n")
}

# Visualize key interaction
if(nrow(sig_interactions) > 0) {
  top_mutation <- sig_interactions$Mutation[1]
  
  # Create 4-group comparison
  multivar_complete <- multivar_complete %>%
    mutate(
      group_4way = case_when(
        cluster_assignment == "Cluster_1" & get(top_mutation) == 0 ~ paste0("C1/", top_mutation, "-"),
        cluster_assignment == "Cluster_1" & get(top_mutation) == 1 ~ paste0("C1/", top_mutation, "+"),
        cluster_assignment == "Cluster_2" & get(top_mutation) == 0 ~ paste0("C2/", top_mutation, "-"),
        cluster_assignment == "Cluster_2" & get(top_mutation) == 1 ~ paste0("C2/", top_mutation, "+")
      )
    )
  
  fit_4way <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ group_4way, data=multivar_complete)
  
  pdf(paste0("04_Figures/20_Manuscript_Prep/interaction_", top_mutation, ".pdf"), width=10, height=8)
  ggsurvplot(fit_4way, 
             data=multivar_complete,
             pval=TRUE,
             conf.int=FALSE,
             risk.table=TRUE,
             title=paste0("Cluster × ", top_mutation, " Interaction"))
  dev.off()
}
```

**Output Files**:
- `03_Results/21_Manuscript_Prep/mutation_cluster_interactions.csv`
- `04_Figures/20_Manuscript_Prep/interaction_*.pdf` (if significant)

**Interpretation**:
- **If interactions significant**: Clusters have INDEPENDENT value (they modify mutation effects)
- **If no interactions**: Clusters are just proxies for mutations (current finding)

---

## PART 6: TIME-STRATIFIED ANALYSIS

### **Background**
HR decreases from 2.22 → 1.62 over time. Explore why.

### **Task 6.1: Early vs Late Death Analysis**

**Method**:
```r
# Load survival data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Define early vs late deaths
# Early: 0-24 months
# Late: >24 months

survival_stratified <- survival_data %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS)) %>%
  mutate(
    death_timing = case_when(
      OS_STATUS == 0 ~ "Censored",
      OS_MONTHS <= 24 ~ "Early death (≤24m)",
      OS_MONTHS > 24 ~ "Late death (>24m)",
      TRUE ~ NA_character_
    ),
    early_death = ifelse(OS_STATUS == 1 & OS_MONTHS <= 24, 1, 0),
    late_death = ifelse(OS_STATUS == 1 & OS_MONTHS > 24, 1, 0)
  )

# Test cluster association with early vs late
early_by_cluster <- table(survival_stratified$cluster_assignment, 
                          survival_stratified$early_death)
fisher_early <- fisher.test(early_by_cluster)

late_by_cluster <- table(survival_stratified$cluster_assignment,
                         survival_stratified$late_death)
fisher_late <- fisher.test(late_by_cluster)

timing_summary <- data.frame(
  Death_Type = c("Early (≤24m)", "Late (>24m)"),
  Cluster1_n = c(early_by_cluster[1,2], late_by_cluster[1,2]),
  Cluster2_n = c(early_by_cluster[2,2], late_by_cluster[2,2]),
  OR = c(fisher_early$estimate, fisher_late$estimate),
  P_value = c(fisher_early$p.value, fisher_late$p.value)
)

write.csv(timing_summary, "03_Results/21_Manuscript_Prep/early_vs_late_death_analysis.csv", row.names=FALSE)

cat("\n=== EARLY VS LATE DEATH BY CLUSTER ===\n")
print(timing_summary)

# Competing risks analysis
# Treat early death as competing risk for late death
library(cmprsk)

survival_stratified <- survival_stratified %>%
  mutate(
    event_type = case_when(
      OS_STATUS == 0 ~ 0,  # Censored
      OS_MONTHS <= 24 ~ 1,  # Early death
      OS_MONTHS > 24 ~ 2    # Late death
    )
  )

# Cumulative incidence by cluster
cif_result <- cuminc(
  ftime = survival_stratified$OS_MONTHS,
  fstatus = survival_stratified$event_type,
  group = survival_stratified$cluster_assignment
)

# Extract cumulative incidence at key timepoints
cif_summary <- data.frame(
  Time = c(12, 24, 36, 48, 60),
  C1_Early = sapply(c(12, 24, 36, 48, 60), function(t) {
    idx <- max(which(cif_result$`Cluster_1 1`$time <= t))
    cif_result$`Cluster_1 1`$est[idx]
  }),
  C2_Early = sapply(c(12, 24, 36, 48, 60), function(t) {
    idx <- max(which(cif_result$`Cluster_2 1`$time <= t))
    cif_result$`Cluster_2 1`$est[idx]
  })
)

cif_summary$Difference <- cif_summary$C2_Early - cif_summary$C1_Early

write.csv(cif_summary, "03_Results/21_Manuscript_Prep/competing_risks_cumulative_incidence.csv", row.names=FALSE)

# Plot competing risks
pdf("04_Figures/20_Manuscript_Prep/competing_risks_plot.pdf", width=10, height=6)
plot(cif_result, 
     col=rep(c("#E41A1C", "#377EB8"), each=2),
     lty=c(1,2,1,2),
     xlab="Months", ylab="Cumulative Incidence",
     main="Competing Risks: Early vs Late Death by Cluster")
legend("topleft",
       legend=c("C1 Early death", "C1 Late death", "C2 Early death", "C2 Late death"),
       col=rep(c("#E41A1C", "#377EB8"), each=2),
       lty=c(1,2,1,2), lwd=2)
dev.off()

cat("\n📊 Interpretation:\n")
if(fisher_early$p.value < 0.05 & fisher_late$p.value > 0.05) {
  cat("  Cluster effect is EARLY-SPECIFIC (first 24 months)\n")
  cat("  Likely reflects treatment response differences\n")
} else if(fisher_early$p.value > 0.05 & fisher_late$p.value < 0.05) {
  cat("  Cluster effect is LATE-SPECIFIC (after 24 months)\n")
  cat("  Likely reflects relapse/progression differences\n")
} else {
  cat("  Cluster effect affects both early and late mortality\n")
}
```

**Output Files**:
- `03_Results/21_Manuscript_Prep/early_vs_late_death_analysis.csv`
- `03_Results/21_Manuscript_Prep/competing_risks_cumulative_incidence.csv`
- `04_Figures/20_Manuscript_Prep/competing_risks_plot.pdf`

---

## PART 7: TREATMENT RESPONSE ANALYSIS (IF DATA AVAILABLE)

### **Task 7.1: Check for Treatment Data in BeatAML**

**Method**:
```r
# Check clinical data for treatment information
clinical_data <- readRDS("01_Data/BeatAML/cohort1_clinical.rds")

cat("\n=== CHECKING FOR TREATMENT DATA ===\n")
cat("Available clinical variables:\n")
print(names(clinical_data))

# Look for treatment-related columns
treatment_cols <- grep("treat|therapy|chemo|drug|response|remission|CR|PR", 
                       names(clinical_data), 
                       ignore.case=TRUE, 
                       value=TRUE)

if(length(treatment_cols) > 0) {
  cat("\n✅ Found treatment-related variables:\n")
  print(treatment_cols)
  
  # Analyze response by cluster
  for(col in treatment_cols) {
    cat("\n--- ", col, " ---\n")
    print(table(clinical_data$cluster_assignment, clinical_data[[col]], useNA="ifany"))
  }
  
  # If complete remission (CR) data available
  if("response_CR" %in% names(clinical_data) | "CR" %in% names(clinical_data)) {
    response_col <- ifelse("response_CR" %in% names(clinical_data), "response_CR", "CR")
    
    # Merge with clusters
    response_data <- clinical_data %>%
      select(lab_id, all_of(response_col)) %>%
      left_join(cluster_assignments, by="lab_id") %>%
      filter(!is.na(cluster_assignment), !is.na(get(response_col)))
    
    # Test association
    response_table <- table(response_data$cluster_assignment, response_data[[response_col]])
    fisher_response <- fisher.test(response_table)
    
    response_summary <- data.frame(
      Cluster = c("Cluster 1", "Cluster 2"),
      CR_n = response_table[,2],
      Total_n = rowSums(response_table),
      CR_rate = response_table[,2] / rowSums(response_table) * 100,
      OR = fisher_response$estimate,
      P_value = fisher_response$p.value
    )
    
    write.csv(response_summary, "03_Results/21_Manuscript_Prep/treatment_response_by_cluster.csv", row.names=FALSE)
    
    cat("\n=== TREATMENT RESPONSE BY CLUSTER ===\n")
    print(response_summary)
  }
  
} else {
  cat("\n❌ No treatment response data found in clinical file\n")
  cat("   Treatment response analysis cannot be performed\n")
}

# Check drug sensitivity data
if(file.exists("01_Data/BeatAML/drug_response.rds")) {
  drug_data <- readRDS("01_Data/BeatAML/drug_response.rds")
  
  cat("\n✅ Found drug sensitivity data\n")
  cat("Drug response dimensions:", dim(drug_data), "\n")
  cat("Number of drugs:", ncol(drug_data) - 1, "\n")  # Assuming first col is sample ID
  
  # Merge with clusters
  drug_cluster <- drug_data %>%
    left_join(cluster_assignments, by="lab_id") %>%
    filter(!is.na(cluster_assignment))
  
  cat("Samples with both drug and cluster data:", nrow(drug_cluster), "\n")
  
  if(nrow(drug_cluster) >= 20) {
    # Test each drug
    drug_cols <- setdiff(names(drug_data), "lab_id")
    
    drug_results <- lapply(drug_cols, function(drug) {
      wilcox_result <- wilcox.test(
        drug_cluster[[drug]][drug_cluster$cluster_assignment == "Cluster_1"],
        drug_cluster[[drug]][drug_cluster$cluster_assignment == "Cluster_2"]
      )
      
      data.frame(
        Drug = drug,
        C1_mean = mean(drug_cluster[[drug]][drug_cluster$cluster_assignment == "Cluster_1"], na.rm=TRUE),
        C2_mean = mean(drug_cluster[[drug]][drug_cluster$cluster_assignment == "Cluster_2"], na.rm=TRUE),
        Difference = mean(drug_cluster[[drug]][drug_cluster$cluster_assignment == "Cluster_2"], na.rm=TRUE) -
                     mean(drug_cluster[[drug]][drug_cluster$cluster_assignment == "Cluster_1"], na.rm=TRUE),
        P_value = wilcox_result$p.value
      )
    }) %>% bind_rows()
    
    drug_results$FDR <- p.adjust(drug_results$P_value, method="BH")
    drug_results <- drug_results %>% arrange(P_value)
    
    write.csv(drug_results, "03_Results/21_Manuscript_Prep/drug_sensitivity_by_cluster.csv", row.names=FALSE)
    
    cat("\n=== DRUG SENSITIVITY DIFFERENCES ===\n")
    cat("Significant at FDR<0.05:", sum(drug_results$FDR < 0.05), "drugs\n")
    print(head(drug_results, 10))
  }
}
```

**Output** (if data available):
- `03_Results/21_Manuscript_Prep/treatment_response_by_cluster.csv`
- `03_Results/21_Manuscript_Prep/drug_sensitivity_by_cluster.csv`

---

## PART 8: IMMUNE CHECKPOINT EXPRESSION ANALYSIS

### **Task 8.1: Test Immune Checkpoint Gene Expression**

**Method**:
```r
# Load expression data
expr_data <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")
cluster_assignments <- read.csv("03_Results/02_Clustering/cluster_assignments.csv")

# Immune checkpoint genes to test
checkpoint_genes <- c(
  "CD274",    # PD-L1
  "PDCD1",    # PD-1
  "CTLA4",    # CTLA-4
  "LAG3",     # LAG-3
  "HAVCR2",   # TIM-3
  "TIGIT",    # TIGIT
  "BTLA",     # BTLA
  "CD276",    # B7-H3
  "VTCN1",    # B7-H4
  "IDO1",     # IDO1
  "CD47"      # CD47
)

# Check availability
available_checkpoints <- intersect(checkpoint_genes, rownames(expr_data))
cat("\n=== IMMUNE CHECKPOINT GENES AVAILABLE ===\n")
cat("Found", length(available_checkpoints), "of", length(checkpoint_genes), "genes\n")
cat("Available:", paste(available_checkpoints, collapse=", "), "\n")

if(length(available_checkpoints) == 0) {
  cat("❌ No checkpoint genes found - check gene identifiers\n")
  
  # Try Ensembl IDs if symbols not found
  checkpoint_ensembl <- c(
    "ENSG00000120217",  # CD274 (PD-L1)
    "ENSG00000188389",  # PDCD1 (PD-1)
    "ENSG00000163599",  # CTLA4
    "ENSG00000089692",  # LAG3
    "ENSG00000135077",  # HAVCR2 (TIM-3)
    "ENSG00000171813",  # TIGIT
    "ENSG00000186265",  # BTLA
    "ENSG00000103855",  # CD276 (B7-H3)
    "ENSG00000134256",  # VTCN1 (B7-H4)
    "ENSG00000131203",  # IDO1
    "ENSG00000091972"   # CD47
  )
  
  available_checkpoints <- intersect(checkpoint_ensembl, rownames(expr_data))
  cat("Trying Ensembl IDs... Found", length(available_checkpoints), "genes\n")
}

if(length(available_checkpoints) > 0) {
  # Extract checkpoint expression
  checkpoint_expr <- expr_data[available_checkpoints, , drop=FALSE]
  checkpoint_df <- as.data.frame(t(checkpoint_expr))
  checkpoint_df$lab_id <- rownames(checkpoint_df)
  
  # Merge with clusters
  checkpoint_cluster <- checkpoint_df %>%
    left_join(cluster_assignments, by="lab_id") %>%
    filter(!is.na(cluster_assignment))
  
  # Test each checkpoint
  checkpoint_results <- lapply(available_checkpoints, function(gene) {
    wilcox_result <- wilcox.test(
      checkpoint_cluster[[gene]][checkpoint_cluster$cluster_assignment == "Cluster_1"],
      checkpoint_cluster[[gene]][checkpoint_cluster$cluster_assignment == "Cluster_2"]
    )
    
    c1_vals <- checkpoint_cluster[[gene]][checkpoint_cluster$cluster_assignment == "Cluster_1"]
    c2_vals <- checkpoint_cluster[[gene]][checkpoint_cluster$cluster_assignment == "Cluster_2"]
    
    data.frame(
      Gene = gene,
      C1_mean = mean(c1_vals, na.rm=TRUE),
      C1_median = median(c1_vals, na.rm=TRUE),
      C2_mean = mean(c2_vals, na.rm=TRUE),
      C2_median = median(c2_vals, na.rm=TRUE),
      Fold_change = mean(c2_vals, na.rm=TRUE) / mean(c1_vals, na.rm=TRUE),
      P_value = wilcox_result$p.value,
      Higher_in = ifelse(mean(c2_vals) > mean(c1_vals), "Cluster 2", "Cluster 1")
    )
  }) %>% bind_rows()
  
  checkpoint_results$FDR <- p.adjust(checkpoint_results$P_value, method="BH")
  checkpoint_results <- checkpoint_results %>% arrange(P_value)
  
  write.csv(checkpoint_results, "03_Results/21_Manuscript_Prep/immune_checkpoint_expression.csv", row.names=FALSE)
  
  cat("\n=== IMMUNE CHECKPOINT EXPRESSION ===\n")
  print(checkpoint_results)
  
  # Plot top checkpoints
  sig_checkpoints <- checkpoint_results %>% filter(FDR < 0.05)
  
  if(nrow(sig_checkpoints) > 0) {
    cat("\n✅ Significant checkpoint differences (FDR<0.05):", nrow(sig_checkpoints), "\n")
    
    # Create boxplots
    pdf("04_Figures/20_Manuscript_Prep/immune_checkpoints_boxplots.pdf", width=12, height=8)
    par(mfrow=c(2, ceiling(nrow(sig_checkpoints)/2)))
    
    for(gene in sig_checkpoints$Gene) {
      boxplot(checkpoint_cluster[[gene]] ~ checkpoint_cluster$cluster_assignment,
              main=gene,
              xlab="Cluster",
              ylab="Expression",
              col=c("#E41A1C", "#377EB8"))
    }
    dev.off()
    
    # Interpretation
    cat("\n📊 Clinical implications:\n")
    if(sum(sig_checkpoints$Higher_in == "Cluster 2") > 0) {
      cat("  Cluster 2 has higher checkpoint expression\n")
      cat("  → May respond better to immune checkpoint inhibitors\n")
      cat("  → Suggests immune exhaustion phenotype\n")
    }
  } else {
    cat("\n❌ No significant checkpoint differences\n")
  }
}
```

**Output Files**:
- `03_Results/21_Manuscript_Prep/immune_checkpoint_expression.csv`
- `04_Figures/20_Manuscript_Prep/immune_checkpoints_boxplots.pdf`

---

## PART 9: GENERATE MANUSCRIPT-READY FIGURES

### **Task 9.1: Main Figure 1 - Mutation Landscape**

**Method**:
```r
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Load data
mutation_enrichment <- read.csv("03_Results/10_Mutations/mutation_enrichment_by_cluster.csv")
mutation_matrix <- readRDS("03_Results/10_Mutations/mutation_matrix.rds")
cluster_assignments <- read.csv("03_Results/02_Clustering/cluster_assignments.csv")

# Merge
mut_cluster <- mutation_matrix %>%
  left_join(cluster_assignments, by="lab_id") %>%
  filter(!is.na(cluster_assignment))

# Panel A: Mutation frequencies by cluster
top_mutations <- mutation_enrichment %>%
  filter(P_Value < 0.05) %>%
  arrange(P_Value) %>%
  head(10) %>%
  pull(Gene)

freq_data <- mut_cluster %>%
  select(cluster_assignment, all_of(top_mutations)) %>%
  pivot_longer(cols = all_of(top_mutations), names_to="Gene", values_to="Mutated") %>%
  group_by(cluster_assignment, Gene) %>%
  summarise(Frequency = mean(Mutated) * 100, .groups="drop")

pdf("04_Figures/21_Main_Figures/Figure1_mutation_landscape.pdf", width=14, height=10)

layout(matrix(c(1,1,2,2,3,3), nrow=3, byrow=TRUE))

# Panel A: Frequency barplot
p1 <- ggplot(freq_data, aes(x=Gene, y=Frequency, fill=cluster_assignment)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("#E41A1C", "#377EB8"),
                    labels=c("Proliferative (C1)", "Immune-Inflammatory (C2)")) +
  labs(x="", y="Mutation Frequency (%)", 
       title="A. Mutation Frequencies by Molecular Subtype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="top",
        legend.title=element_blank(),
        text=element_text(size=12))

print(p1)

# Panel B: Forest plot of enrichment
enrichment_plot_data <- mutation_enrichment %>%
  filter(Gene %in% top_mutations) %>%
  mutate(
    Enriched_in = ifelse(Cluster_1_Freq > Cluster_2_Freq, "Cluster 1", "Cluster 2"),
    log_OR = log(Odds_Ratio),
    log_OR_lower = log(Odds_Ratio) - 1.96*sqrt(1/Cluster_1_Mutated + 1/Cluster_1_WT + 
                                                 1/Cluster_2_Mutated + 1/Cluster_2_WT),
    log_OR_upper = log(Odds_Ratio) + 1.96*sqrt(1/Cluster_1_Mutated + 1/Cluster_1_WT + 
                                                 1/Cluster_2_Mutated + 1/Cluster_2_WT)
  ) %>%
  arrange(desc(log_OR))

p2 <- ggplot(enrichment_plot_data, aes(x=reorder(Gene, log_OR), y=log_OR)) +
  geom_point(aes(color=Enriched_in), size=4) +
  geom_errorbar(aes(ymin=log_OR_lower, ymax=log_OR_upper, color=Enriched_in), width=0.2) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_color_manual(values=c("Cluster 1"="#E41A1C", "Cluster 2"="#377EB8")) +
  coord_flip() +
  labs(y="Log Odds Ratio (C1 vs C2)", x="",
       title="B. Mutation Enrichment by Subtype") +
  theme_bw() +
  theme(legend.position="top",
        text=element_text(size=12))

print(p2)

# Panel C: Oncoprint
mat <- as.matrix(mut_cluster[, top_mutations])
rownames(mat) <- mut_cluster$lab_id
mat[is.na(mat)] <- 0

# Order by cluster
cluster_order <- order(mut_cluster$cluster_assignment)
mat_ordered <- mat[cluster_order, ]

# Annotation
ha <- HeatmapAnnotation(
  Cluster = mut_cluster$cluster_assignment[cluster_order],
  col = list(Cluster = c("Cluster_1"="#E41A1C", "Cluster_2"="#377EB8")),
  show_legend=TRUE
)

ht <- Heatmap(t(mat_ordered),
              name = "Mutation",
              col = c("0"="white", "1"="black"),
              top_annotation = ha,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              row_names_side = "left",
              column_title = "C. Mutation Matrix by Subtype (n=522 samples)",
              heatmap_legend_param = list(
                at = c(0, 1),
                labels = c("Wild-type", "Mutated")
              ))

print(ht)

dev.off()

cat("\n✅ Figure 1 created: Mutation landscape\n")
```

---

### **Task 9.2: Main Figure 2 - Survival Meta-Analysis**

**Method**:
```r
library(meta)
library(grid)
library(gridExtra)

pdf("04_Figures/21_Main_Figures/Figure2_survival_meta_analysis.pdf", width=14, height=10)

layout(matrix(c(1,1,2,2,3,3), nrow=3, byrow=TRUE))

# Panel A: BeatAML KM curve
beatAML_surv <- survival_data %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS), !is.na(cluster_assignment))

fit_beatAML <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment, data=beatAML_surv)

p1 <- ggsurvplot(fit_beatAML,
                 data = beatAML_surv,
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 legend.labs = c("Proliferative (C1)", "Immune-Inflammatory (C2)"),
                 palette = c("#E41A1C", "#377EB8"),
                 title = "A. BeatAML Discovery Cohort (n=671)",
                 xlab = "Months",
                 ylab = "Overall Survival",
                 break.time.by = 12,
                 ggtheme = theme_bw())

print(p1)

# Panel B: TCGA KM curve
tcga_surv <- read.csv("03_Results/17_TCGA_Validation/tcga_survival_validation.csv")
fit_tcga <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster, data=tcga_surv)

p2 <- ggsurvplot(fit_tcga,
                 data = tcga_surv,
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 legend.labs = c("Proliferative (C1)", "Immune-Inflammatory (C2)"),
                 palette = c("#E41A1C", "#377EB8"),
                 title = "B. TCGA-LAML Validation (n=151)",
                 xlab = "Months",
                 ylab = "Overall Survival",
                 break.time.by = 12,
                 ggtheme = theme_bw())

print(p2)

# Panel C: Forest plot meta-analysis
meta_data <- data.frame(
  Cohort = c("BeatAML", "TCGA-LAML"),
  HR = c(1.38, 1.24),
  Lower = c(1.13, 0.80),
  Upper = c(1.68, 1.94),
  N = c(671, 151),
  Events = c(398, 97)
)

# Meta-analysis
meta_result <- metagen(
  TE = log(meta_data$HR),
  seTE = (log(meta_data$Upper) - log(meta_data$Lower)) / (2*1.96),
  studlab = meta_data$Cohort,
  sm = "HR",
  fixed = TRUE,
  random = TRUE
)

forest(meta_result,
       leftcols = c("studlab", "N", "Events"),
       leftlabs = c("Cohort", "N", "Events"),
       xlab = "Hazard Ratio (C2 vs C1)",
       xlim = c(0.5, 3),
       at = c(0.5, 1, 1.5, 2, 2.5, 3),
       digits = 2,
       col.square = "navy",
       col.square.lines = "navy",
       col.diamond = "darkred",
       col.diamond.lines = "darkred",
       print.tau2 = FALSE,
       print.I2 = TRUE,
       print.pval.Q = TRUE)

title("C. Meta-Analysis of Adult Cohorts", line=2.5)

dev.off()

cat("\n✅ Figure 2 created: Survival meta-analysis\n")
```

---

### **Task 9.3: Main Figure 3 - Age-Specific Heterogeneity**

**Method**:
```r
pdf("04_Figures/21_Main_Figures/Figure3_age_heterogeneity.pdf", width=14, height=12)

# Three KM curves stacked
par(mfrow=c(3,1))

# Panel A: BeatAML (adults)
fit_beatAML <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment, data=beatAML_surv)
plot(fit_beatAML, col=c("#E41A1C", "#377EB8"), lwd=2,
     xlab="Months", ylab="Survival Probability",
     main="A. BeatAML (Adult, median age 62y) - HR=1.38, p=0.001")
legend("topright", legend=c("C1 (better)", "C2 (worse)"),
       col=c("#E41A1C", "#377EB8"), lwd=2)

# Panel B: TCGA (adults)
plot(fit_tcga, col=c("#E41A1C", "#377EB8"), lwd=2,
     xlab="Months", ylab="Survival Probability",
     main="B. TCGA-LAML (Adult, median age 58y) - HR=1.24, p=0.35")
legend("topright", legend=c("C1 (better)", "C2 (worse)"),
       col=c("#E41A1C", "#377EB8"), lwd=2)

# Panel C: TARGET (pediatric)
target_surv <- read.csv("03_Results/18_TARGET_Validation/target_survival_validation.csv")
fit_target <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster, data=target_surv)
plot(fit_target, col=c("#E41A1C", "#377EB8"), lwd=2,
     xlab="Months", ylab="Survival Probability",
     main="C. TARGET-AML (Pediatric, median age 10y) - HR=0.81, p=0.052 (OPPOSITE!)")
legend("topright", legend=c("C1 (WORSE in children)", "C2 (BETTER in children)"),
       col=c("#E41A1C", "#377EB8"), lwd=2)

# Panel D: Forest plot showing heterogeneity
par(mfrow=c(1,1))

meta_all <- data.frame(
  Cohort = c("BeatAML", "TCGA", "TARGET"),
  Age_Group = c("Adult", "Adult", "Pediatric"),
  HR = c(1.38, 1.24, 0.81),
  Lower = c(1.13, 0.80, 0.66),
  Upper = c(1.68, 1.94, 1.00)
)

meta_all3 <- metagen(
  TE = log(meta_all$HR),
  seTE = (log(meta_all$Upper) - log(meta_all$Lower)) / (2*1.96),
  studlab = paste0(meta_all$Cohort, " (", meta_all$Age_Group, ")"),
  sm = "HR"
)

forest(meta_all3,
       xlab = "Hazard Ratio (C2 vs C1)",
       leftcols = c("studlab"),
       leftlabs = c("Cohort"),
       rightcols = c("effect", "ci"),
       rightlabs = c("HR", "95% CI"),
       col.square = c("navy", "navy", "darkred"),
       print.I2 = TRUE,
       print.tau2 = TRUE,
       print.pval.Q = TRUE,
       xlim = c(0.5, 2))

title("D. Cross-Cohort Meta-Analysis Showing Age Heterogeneity", line=2)
text(1, 0.5, paste0("Heterogeneity: I² = 84.8%, p = 0.001"), pos=4)

dev.off()

cat("\n✅ Figure 3 created: Age-specific heterogeneity\n")
```

---

### **Task 9.4: Main Figure 4 - Multivariate Analysis**

**Method**:
```r
pdf("04_Figures/21_Main_Figures/Figure4_multivariate_analysis.pdf", width=12, height=8)

# Load multivariate results
multivar_full <- read.csv("03_Results/11_Survival_Reanalysis/05_full_model_coefficients.csv")

# Forest plot
multivar_full <- multivar_full %>%
  mutate(
    Significant = P_value < 0.05,
    Variable_clean = case_when(
      Variable == "cluster_assignmentCluster_2" ~ "Cluster (C2 vs C1)",
      Variable == "AGE" ~ "Age (per year)",
      Variable == "SEXMale" ~ "Sex (Male vs Female)",
      Variable == "TP53" ~ "TP53 mutation",
      Variable == "TET2" ~ "TET2 mutation",
      Variable == "RUNX1" ~ "RUNX1 mutation",
      Variable == "ASXL1" ~ "ASXL1 mutation",
      TRUE ~ Variable
    )
  ) %>%
  arrange(desc(HR))

ggplot(multivar_full, aes(x=HR, y=reorder(Variable_clean, HR))) +
  geom_vline(xintercept=1, linetype="dashed", color="gray50") +
  geom_point(aes(color=Significant), size=4) +
  geom_errorbarh(aes(xmin=HR_lower, xmax=HR_upper, color=Significant), height=0.2, linewidth=1) +
  scale_color_manual(values=c("FALSE"="gray70", "TRUE"="darkred"),
                     labels=c("FALSE"="p≥0.05", "TRUE"="p<0.05")) +
  scale_x_log10(breaks=c(0.5, 1, 2, 3, 4)) +
  labs(x="Hazard Ratio (log scale)", y="",
       title="Multivariate Cox Regression - Full Model",
       subtitle=paste0("n=459 samples, 282 events\n",
                      "Cluster HR=1.06, p=0.649 (NOT significant)\n",
                      "TP53 HR=2.96, p<1e-9 (highly significant)")) +
  theme_bw() +
  theme(text=element_text(size=12),
        legend.position="top",
        legend.title=element_blank())

dev.off()

cat("\n✅ Figure 4 created: Multivariate analysis\n")
```

---

## PART 10: SUPPLEMENTARY TABLES

### **Task 10.1: Table S1 - Sample Characteristics**

**Method**:
```r
library(tableone)

# Load data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Create Table 1
vars <- c("AGE", "SEX", "OS_MONTHS", "OS_STATUS", 
          "NPM1", "TP53", "RUNX1", "DNMT3A", "FLT3", "TET2", "ASXL1")

table1 <- CreateTableOne(
  vars = vars,
  strata = "cluster_assignment",
  data = survival_data,
  factorVars = c("SEX", "OS_STATUS", "NPM1", "TP53", "RUNX1", "DNMT3A", "FLT3", "TET2", "ASXL1"),
  test = TRUE
)

table1_df <- print(table1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

write.csv(table1_df, "03_Results/22_Supplementary_Tables/TableS1_sample_characteristics.csv")

cat("\n✅ Table S1 created: Sample characteristics\n")
```

---

### **Task 10.2: Table S2 - All Survival Models Comparison**

**Method**:
```r
# Compile all survival analyses
survival_summary <- data.frame(
  Method = c(
    "Standard Cox (Phase 2)",
    "Stratified Cox",
    "Time-varying coefficients",
    "Landmark 6m",
    "Landmark 12m",
    "Landmark 18m",
    "Landmark 24m",
    "RMST (24m)",
    "RMST (60m)",
    "BeatAML meta-analysis",
    "TCGA validation",
    "TARGET validation (pediatric)"
  ),
  HR = c(1.22, NA, 2.22, 1.33, 1.35, 1.37, 1.38, NA, NA, 1.35, 1.24, 0.81),
  Lower_CI = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 1.13, 0.80, 0.66),
  Upper_CI = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 1.62, 1.94, 1.00),
  P_value = c(0.053, 0.00155, 0.058, 0.017, 0.011, 0.015, 0.018, 0.029, 0.007, 0.001, 0.353, 0.052),
  N_samples = c(671, 671, 671, 645, 598, 560, 529, 671, 671, 822, 151, 1713),
  Interpretation = c(
    "Violated PH assumption",
    "Assumption-free, significant",
    "Time-varying effect",
    "Conditional survival significant",
    "Conditional survival significant",
    "Conditional survival significant",
    "Conditional survival significant",
    "Survival time difference",
    "Survival time difference",
    "Pooled adult effect",
    "Underpowered validation",
    "Opposite effect in children"
  )
)

write.csv(survival_summary, "03_Results/22_Supplementary_Tables/TableS2_all_survival_analyses.csv", row.names=FALSE)

cat("\n✅ Table S2 created: All survival models\n")
```

---

### **Task 10.3: Table S3 - All Statistical Tests with Corrections**

**Method**:
```r
# Load from Part 4
all_tests <- read.csv("03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv")

# Format for publication
table_s3 <- all_tests %>%
  select(
    Analysis = Analysis_Part,
    Test = Test_Type,
    Comparison,
    Raw_p = Raw_P_value,
    FDR_within = FDR_within_part,
    FDR_studywide,
    Interpretation
  ) %>%
  arrange(Analysis, Raw_p)

write.csv(table_s3, "03_Results/22_Supplementary_Tables/TableS3_all_statistical_tests.csv", row.names=FALSE)

cat("\n✅ Table S3 created: All statistical tests\n")
```

---

## EXECUTION SUMMARY

### **Files to be Generated** (Total: ~30 files)

**Results Files** (15):
- `03_Results/21_Manuscript_Prep/sample_attrition_table.csv`
- `03_Results/21_Manuscript_Prep/imputation_sensitivity_analysis.csv`
- `03_Results/21_Manuscript_Prep/target_imputation_sensitivity.csv`
- `03_Results/21_Manuscript_Prep/target_survival_sensitivity.csv`
- `03_Results/21_Manuscript_Prep/missing_genes_in_target.csv`
- `03_Results/21_Manuscript_Prep/tcga_power_analysis_detailed.csv`
- `03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv`
- `03_Results/21_Manuscript_Prep/mutation_cluster_interactions.csv`
- `03_Results/21_Manuscript_Prep/early_vs_late_death_analysis.csv`
- `03_Results/21_Manuscript_Prep/competing_risks_cumulative_incidence.csv`
- `03_Results/21_Manuscript_Prep/treatment_response_by_cluster.csv` (if available)
- `03_Results/21_Manuscript_Prep/drug_sensitivity_by_cluster.csv` (if available)
- `03_Results/21_Manuscript_Prep/immune_checkpoint_expression.csv`
- `03_Results/22_Supplementary_Tables/TableS1_sample_characteristics.csv`
- `03_Results/22_Supplementary_Tables/TableS2_all_survival_analyses.csv`
- `03_Results/22_Supplementary_Tables/TableS3_all_statistical_tests.csv`

**Figure Files** (12):
- `04_Figures/20_Manuscript_Prep/sample_attrition_flowchart.pdf`
- `04_Figures/20_Manuscript_Prep/tcga_power_curve.pdf`
- `04_Figures/20_Manuscript_Prep/interaction_*.pdf` (if significant)
- `04_Figures/20_Manuscript_Prep/competing_risks_plot.pdf`
- `04_Figures/20_Manuscript_Prep/immune_checkpoints_boxplots.pdf`
- `04_Figures/21_Main_Figures/Figure1_mutation_landscape.pdf`
- `04_Figures/21_Main_Figures/Figure2_survival_meta_analysis.pdf`
- `04_Figures/21_Main_Figures/Figure3_age_heterogeneity.pdf`
- `04_Figures/21_Main_Figures/Figure4_multivariate_analysis.pdf`

---

## EXECUTION PRIORITY

**MUST DO** (Critical for manuscript):
1. ✅ Part 1: Sample attrition (explains n=459)
2. ✅ Part 3: TCGA power analysis (defends validation)
3. ✅ Part 4: Multiple testing catalog (transparency)
4. ✅ Part 9: Main figures (publication-ready)
5. ✅ Part 10: Supplementary tables (required)

**SHOULD DO** (Strengthens manuscript):
6. ✅ Part 2: TARGET sensitivity (validates robustness)
7. ✅ Part 5: Mutation interactions (tests independence)
8. ✅ Part 6: Time-stratified analysis (explains HR decay)

**NICE TO HAVE** (If data available):
9. ⭕ Part 7: Treatment response (adds clinical value)
10. ⭕ Part 8: Immune checkpoints (therapeutic implications)

---

## EXPECTED RUNTIME

- Part 1-6: ~4 hours
- Part 7-8: ~2 hours (if data available)
- Part 9-10: ~3 hours
- **Total: 8-10 hours**

---

## BEGIN EXECUTION

Start with Part 1 and proceed sequentially. Report completion of each part before moving to next.

**Priority order**: 1 → 3 → 4 → 9 → 10 → 2 → 5 → 6 → (7, 8 if applicable)

Good luck! 🚀