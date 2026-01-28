# PHASE 3 - PART 1.1: STRATIFIED COX REGRESSION
# Fix proportional hazards violations by treating cluster as stratification variable
# This addresses the PH violation (global test p=0.0002) found in Phase 2

library(tidyverse)
library(survival)
library(survminer)

cat("=== PHASE 3: PART 1.1 - STRATIFIED COX REGRESSION ===\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("=== STEP 1: LOADING DATA ===\n\n")

# Load survival data with cluster assignments
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

cat("Loaded survival data:", nrow(survival_data), "samples\n")
cat("Variables:", paste(colnames(survival_data), collapse=", "), "\n\n")

# Check for required columns
required_cols <- c("OS_months", "OS_event", "cluster")
if (!all(required_cols %in% colnames(survival_data))) {
  cat("ERROR: Missing required columns\n")
  cat("Required:", paste(required_cols, collapse=", "), "\n")
  cat("Found:", paste(colnames(survival_data), collapse=", "), "\n")
  quit(save = "no", status = 1)
}

# Rename columns for consistency with script
survival_data <- survival_data %>%
  rename(
    cluster_assignment = cluster,
    OS_MONTHS = OS_months,
    OS_STATUS = OS_event
  )

# Summary
cat("Cluster distribution:\n")
print(table(survival_data$cluster_assignment))
cat("\nEvents by cluster:\n")
print(table(survival_data$cluster_assignment, survival_data$OS_STATUS))
cat("\n")

# ============================================================================
# 2. ASSUMPTION-FREE LOG-RANK TEST
# ============================================================================

cat("=== STEP 2: ASSUMPTION-FREE LOG-RANK TEST ===\n\n")

# Create survival object
surv_obj <- Surv(time = survival_data$OS_MONTHS, event = survival_data$OS_STATUS)

# Log-rank test (makes no assumptions about hazard ratios)
logrank_test <- survdiff(surv_obj ~ cluster_assignment, data = survival_data)

# Extract p-value
logrank_pvalue <- 1 - pchisq(logrank_test$chisq, df = 1)

cat("Log-Rank Test Results:\n")
cat("  Chi-square statistic:", round(logrank_test$chisq, 3), "\n")
cat("  Degrees of freedom:", 1, "\n")
cat("  P-value:", format(logrank_pvalue, scientific = TRUE, digits = 3), "\n")

if (logrank_pvalue < 0.05) {
  cat("  ✓ SIGNIFICANT: Cluster assignment is associated with survival\n")
} else {
  cat("  ✗ NOT SIGNIFICANT: No survival difference between clusters\n")
}
cat("\n")

# ============================================================================
# 3. KAPLAN-MEIER ANALYSIS (Stratified by Cluster)
# ============================================================================

cat("=== STEP 3: KAPLAN-MEIER ANALYSIS ===\n\n")

# Fit Kaplan-Meier curves
km_fit <- survfit(surv_obj ~ cluster_assignment, data = survival_data)

# Get median survival times
surv_summary <- summary(km_fit)$table

cat("Median Survival Times:\n")
print(surv_summary[, c("median", "0.95LCL", "0.95UCL")])
cat("\n")

# Get survival probabilities at key timepoints
timepoints <- c(6, 12, 24, 36)
surv_probs <- summary(km_fit, times = timepoints)

cat("Survival Probabilities at Key Timepoints:\n")
for (i in 1:length(timepoints)) {
  cat("\n", timepoints[i], "months:\n", sep="")
  for (j in 1:length(unique(survival_data$cluster_assignment))) {
    idx <- (i-1) * length(unique(survival_data$cluster_assignment)) + j
    if (idx <= length(surv_probs$surv)) {
      cat("  Cluster", j, ":", round(surv_probs$surv[idx], 3),
          "(95% CI:", round(surv_probs$lower[idx], 3), "-",
          round(surv_probs$upper[idx], 3), ")\n")
    }
  }
}
cat("\n")

# ============================================================================
# 4. STRATIFIED COX REGRESSION
# ============================================================================

cat("=== STEP 4: STRATIFIED COX REGRESSION ===\n\n")

# Treat cluster as stratification variable (not covariate)
# This allows baseline hazard to differ between clusters
# No PH assumption needed for stratification variable

cat("Fitting stratified Cox model: cluster as strata\n")
cat("Formula: Surv(OS_MONTHS, OS_STATUS) ~ strata(cluster_assignment)\n\n")

stratified_cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ strata(cluster_assignment),
                        data = survival_data)

cat("Stratified Cox Model Results:\n")
print(summary(stratified_cox))
cat("\n")

# Note: Stratified model has no coefficients for strata variable
# This is correct - we're allowing different baseline hazards

# ============================================================================
# 5. COMPARE WITH STANDARD COX (for reference)
# ============================================================================

cat("=== STEP 5: COMPARISON WITH STANDARD COX ===\n\n")

# Standard Cox (cluster as covariate) - for comparison only
cat("Standard Cox model (for comparison):\n")
standard_cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                      data = survival_data)

# Test PH assumption
cox_zph <- cox.zph(standard_cox)

cat("\nProportional Hazards Test (Standard Cox):\n")
print(cox_zph)

cat("\nPH Assumption Test Results:\n")
cat("  Global test p-value:", format(cox_zph$table["GLOBAL", "p"], scientific = TRUE, digits = 3), "\n")

if (cox_zph$table["GLOBAL", "p"] < 0.05) {
  cat("  ✗ PH ASSUMPTION VIOLATED - Standard Cox is inappropriate\n")
  cat("  ✓ Stratified Cox (without PH assumption) is the correct approach\n")
} else {
  cat("  ✓ PH assumption holds - Standard Cox is appropriate\n")
}
cat("\n")

# Extract hazard ratio from standard Cox (for reference)
hr <- exp(coef(standard_cox))
hr_ci <- exp(confint(standard_cox))

cat("\nHazard Ratio from Standard Cox (Cluster 2 vs Cluster 1):\n")
cat("  HR:", round(hr, 3), "\n")
cat("  95% CI:", round(hr_ci[1], 3), "-", round(hr_ci[2], 3), "\n")
cat("  P-value:", format(summary(standard_cox)$coefficients[5], scientific = TRUE, digits = 3), "\n\n")

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

cat("=== STEP 6: SAVING RESULTS ===\n\n")

# Create output directory
dir.create("03_Results/11_Survival_Reanalysis", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/11_Survival_Reanalysis", showWarnings = FALSE, recursive = TRUE)

# Save stratified survival results
results_df <- data.frame(
  analysis = "Stratified Cox Regression",
  method = "Log-Rank Test",
  test_statistic = logrank_test$chisq,
  df = 1,
  pvalue = logrank_pvalue,
  median_surv_cluster1 = surv_summary["cluster_assignment=1", "median"],
  median_surv_cluster2 = surv_summary["cluster_assignment=2", "median"],
  ph_assumption_pvalue = cox_zph$table["GLOBAL", "p"],
  ph_assumption_violated = cox_zph$table["GLOBAL", "p"] < 0.05,
  standard_cox_hr = hr,
  standard_cox_hr_lower = hr_ci[1],
  standard_cox_hr_upper = hr_ci[2],
  note = "Stratified Cox treats cluster as strata (no PH assumption needed)"
)

write.csv(results_df,
          "03_Results/11_Survival_Reanalysis/01_stratified_cox_results.csv",
          row.names = FALSE)
cat("✓ Saved: 01_stratified_cox_results.csv\n")

# Save detailed survival summary
detailed_summary <- as.data.frame(surv_summary)
write.csv(detailed_summary,
          "03_Results/11_Survival_Reanalysis/01_km_survival_summary.csv",
          row.names = TRUE)
cat("✓ Saved: 01_km_survival_summary.csv\n\n")

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

cat("=== STEP 7: CREATING VISUALIZATIONS ===\n\n")

# Create cluster labels for plotting
survival_data$cluster_label <- factor(
  survival_data$cluster_assignment,
  levels = c(1, 2),
  labels = c("Cluster 1 (Proliferative)", "Cluster 2 (Immune-Inflammatory)")
)

# Kaplan-Meier plot with risk table
p_km <- ggsurvplot(
  km_fit,
  data = survival_data,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  ncensor.plot = FALSE,
  legend.title = "Molecular Subtype",
  legend.labs = c("Cluster 1 (Proliferative)", "Cluster 2 (Immune-Inflammatory)"),
  palette = c("#2E9FDF", "#E7B800"),
  title = "Stratified Survival Analysis (Assumption-Free Log-Rank Test)",
  subtitle = paste0("Log-rank p = ", format(logrank_pvalue, scientific = TRUE, digits = 3)),
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  ggtheme = theme_bw(base_size = 12),
  font.legend = 11,
  font.x = 11,
  font.y = 11,
  font.tickslab = 10
)

pdf("04_Figures/11_Survival_Reanalysis/01_stratified_kaplan_meier.pdf",
    width = 10, height = 8)
print(p_km)
dev.off()

cat("✓ Saved: 01_stratified_kaplan_meier.pdf\n")

# Schoenfeld residuals plot (visualize PH violation)
pdf("04_Figures/11_Survival_Reanalysis/01_schoenfeld_residuals.pdf",
    width = 8, height = 6)
plot(cox_zph, main = "Schoenfeld Residuals Test\n(Standard Cox Model)")
abline(h = 0, col = "red", lty = 2)
dev.off()

cat("✓ Saved: 01_schoenfeld_residuals.pdf\n\n")

# ============================================================================
# 8. SUMMARY AND INTERPRETATION
# ============================================================================

cat("=== SUMMARY AND INTERPRETATION ===\n\n")

cat("PART 1.1: STRATIFIED COX REGRESSION COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. PROPORTIONAL HAZARDS ASSUMPTION:\n")
cat("   Global test p-value:", format(cox_zph$table["GLOBAL", "p"], scientific = TRUE, digits = 3), "\n")
if (cox_zph$table["GLOBAL", "p"] < 0.05) {
  cat("   ✗ PH assumption is VIOLATED\n")
  cat("   → Standard Cox regression is INAPPROPRIATE\n")
  cat("   → Stratified Cox or alternative methods required\n\n")
} else {
  cat("   ✓ PH assumption holds\n\n")
}

cat("2. ASSUMPTION-FREE LOG-RANK TEST:\n")
cat("   P-value:", format(logrank_pvalue, scientific = TRUE, digits = 3), "\n")
if (logrank_pvalue < 0.05) {
  cat("   ✓ SIGNIFICANT survival difference between clusters\n")
  cat("   → Molecular subtypes are prognostic\n\n")
} else {
  cat("   ✗ No significant survival difference\n\n")
}

cat("3. MEDIAN SURVIVAL TIMES:\n")
cat("   Cluster 1:", round(surv_summary["cluster_assignment=1", "median"], 1), "months\n")
cat("   Cluster 2:", round(surv_summary["cluster_assignment=2", "median"], 1), "months\n")
cat("   Difference:", round(abs(surv_summary["cluster_assignment=1", "median"] -
                                 surv_summary["cluster_assignment=2", "median"]), 1), "months\n\n")

cat("4. STATISTICAL APPROACH:\n")
cat("   ✓ Log-rank test: Makes NO assumptions about hazard ratios\n")
cat("   ✓ Stratified Cox: Allows different baseline hazards per cluster\n")
cat("   ✓ Both methods are appropriate when PH assumption is violated\n\n")

cat("INTERPRETATION:\n")
if (logrank_pvalue < 0.05) {
  cat("The molecular subtypes show SIGNIFICANT prognostic value for overall survival.\n")
  cat("The log-rank test confirms this association WITHOUT requiring the proportional\n")
  cat("hazards assumption that was violated in the standard Cox model.\n\n")
} else {
  cat("No significant survival difference was detected between molecular subtypes.\n")
  cat("This may indicate limited prognostic value or insufficient power.\n\n")
}

cat("NEXT STEPS:\n")
cat("  → Part 1.2: Time-varying coefficient model\n")
cat("  → Part 1.3: Landmark analysis\n")
cat("  → Part 1.4: RMST analysis\n\n")

cat("### PART 1.1 COMPLETE ###\n")
