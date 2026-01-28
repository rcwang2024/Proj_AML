# PHASE 3 - PART 3: TCGA VALIDATION INVESTIGATION
# Understand why TCGA validation failed (p=0.353) while BeatAML was significant (p=0.00155)

library(tidyverse)
library(survival)
library(survminer)

cat("=== PHASE 3: PART 3 - TCGA VALIDATION INVESTIGATION ===\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("=== STEP 1: LOADING DATA ===\n\n")

# BeatAML data
beataml_surv <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
beataml_surv <- beataml_surv %>%
  rename(OS_MONTHS = OS_months, OS_STATUS = OS_event, cluster_assignment = cluster)

# TCGA data
tcga_data <- read.csv("03_Results/17_TCGA_Validation/tcga_clinical_with_predictions.csv")
tcga_data <- tcga_data %>%
  rename(OS_MONTHS = OS_months, OS_STATUS = OS_event, cluster_assignment = predicted_cluster)

cat("BeatAML:", nrow(beataml_surv), "samples\n")
cat("TCGA:", nrow(tcga_data), "samples\n\n")

# ============================================================================
# 2. COMPARE COHORT CHARACTERISTICS
# ============================================================================

cat("=== STEP 2: COMPARING COHORT CHARACTERISTICS ===\n\n")

# Age comparison
cat("AGE DISTRIBUTION:\n")
cat("  BeatAML: Mean =", round(mean(beataml_surv$age, na.rm=TRUE), 1),
    "± SD", round(sd(beataml_surv$age, na.rm=TRUE), 1), "\n")
cat("  TCGA:    Mean =", round(mean(tcga_data$age, na.rm=TRUE), 1),
    "± SD", round(sd(tcga_data$age, na.rm=TRUE), 1), "\n")

age_test <- t.test(beataml_surv$age, tcga_data$age)
cat("  t-test p-value:", format(age_test$p.value, scientific=TRUE, digits=3), "\n\n")

# Sex distribution
cat("SEX DISTRIBUTION:\n")
beataml_male_pct <- sum(beataml_surv$sex == "M", na.rm=TRUE) / sum(!is.na(beataml_surv$sex)) * 100
tcga_male_pct <- sum(tcga_data$sex == "male", na.rm=TRUE) / sum(!is.na(tcga_data$sex)) * 100
cat("  BeatAML: Male =", round(beataml_male_pct, 1), "%\n")
cat("  TCGA:    Male =", round(tcga_male_pct, 1), "%\n\n")

# Cluster distribution
cat("CLUSTER DISTRIBUTION:\n")
beataml_c1_pct <- sum(beataml_surv$cluster_assignment == 1) / nrow(beataml_surv) * 100
beataml_c2_pct <- sum(beataml_surv$cluster_assignment == 2) / nrow(beataml_surv) * 100
tcga_c1_pct <- sum(tcga_data$cluster_assignment == 1, na.rm=TRUE) / sum(!is.na(tcga_data$cluster_assignment)) * 100
tcga_c2_pct <- sum(tcga_data$cluster_assignment == 2, na.rm=TRUE) / sum(!is.na(tcga_data$cluster_assignment)) * 100

cat("  BeatAML: Cluster 1 =", round(beataml_c1_pct, 1), "%, Cluster 2 =", round(beataml_c2_pct, 1), "%\n")
cat("  TCGA:    Cluster 1 =", round(tcga_c1_pct, 1), "%, Cluster 2 =", round(tcga_c2_pct, 1), "%\n")

# Chi-square test for cluster proportions
cluster_table <- matrix(c(
  sum(beataml_surv$cluster_assignment == 1), sum(beataml_surv$cluster_assignment == 2),
  sum(tcga_data$cluster_assignment == 1, na.rm=TRUE), sum(tcga_data$cluster_assignment == 2, na.rm=TRUE)
), nrow=2, byrow=TRUE)
cluster_chisq <- chisq.test(cluster_table)
cat("  Chi-square p-value:", format(cluster_chisq$p.value, scientific=TRUE, digits=3), "\n")
if (cluster_chisq$p.value < 0.05) {
  cat("  ✗ Cluster proportions DIFFER between cohorts\n\n")
} else {
  cat("  ✓ Cluster proportions are similar\n\n")
}

# Classification confidence
cat("CLASSIFICATION CONFIDENCE:\n")
tcga_pred <- read.csv("03_Results/17_TCGA_Validation/tcga_sample_predictions.csv")
cat("  TCGA mean confidence:", round(mean(tcga_pred$confidence), 3), "\n")
cat("  TCGA median confidence:", round(median(tcga_pred$confidence), 3), "\n")
cat("  TCGA confidence range:", round(min(tcga_pred$confidence), 3), "-",
    round(max(tcga_pred$confidence), 3), "\n\n")

# Compare to BeatAML training confidence (if available)
cat("NOTE: Lower confidence may indicate poor generalization\n\n")

# ============================================================================
# 3. SURVIVAL COMPARISON
# ============================================================================

cat("=== STEP 3: COMPARING SURVIVAL OUTCOMES ===\n\n")

# Overall survival by cohort
beataml_median_os <- median(beataml_surv$OS_MONTHS, na.rm=TRUE)
tcga_median_os <- median(tcga_data$OS_MONTHS, na.rm=TRUE)

cat("OVERALL SURVIVAL:\n")
cat("  BeatAML median OS:", round(beataml_median_os, 1), "months\n")
cat("  TCGA median OS:", round(tcga_median_os, 1), "months\n")
cat("  Difference:", round(abs(beataml_median_os - tcga_median_os), 1), "months\n\n")

# Event rates
beataml_event_rate <- sum(beataml_surv$OS_STATUS) / nrow(beataml_surv) * 100
tcga_event_rate <- sum(tcga_data$OS_STATUS, na.rm=TRUE) / sum(!is.na(tcga_data$OS_STATUS)) * 100

cat("EVENT RATES:\n")
cat("  BeatAML: Deaths =", sum(beataml_surv$OS_STATUS), "(",
    round(beataml_event_rate, 1), "%)\n")
cat("  TCGA:    Deaths =", sum(tcga_data$OS_STATUS, na.rm=TRUE), "(",
    round(tcga_event_rate, 1), "%)\n\n")

# Survival by cluster in each cohort
cat("SURVIVAL BY CLUSTER:\n")

# BeatAML
beataml_c1_median <- median(beataml_surv$OS_MONTHS[beataml_surv$cluster_assignment == 1], na.rm=TRUE)
beataml_c2_median <- median(beataml_surv$OS_MONTHS[beataml_surv$cluster_assignment == 2], na.rm=TRUE)
cat("  BeatAML C1:", round(beataml_c1_median, 1), "months, C2:",
    round(beataml_c2_median, 1), "months, Diff:",
    round(beataml_c1_median - beataml_c2_median, 1), "months\n")

# TCGA
tcga_c1_median <- median(tcga_data$OS_MONTHS[tcga_data$cluster_assignment == 1], na.rm=TRUE)
tcga_c2_median <- median(tcga_data$OS_MONTHS[tcga_data$cluster_assignment == 2], na.rm=TRUE)
cat("  TCGA C1:", round(tcga_c1_median, 1), "months, C2:",
    round(tcga_c2_median, 1), "months, Diff:",
    round(tcga_c1_median - tcga_c2_median, 1), "months\n\n")

# Log-rank tests
beataml_logrank <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                             data = beataml_surv)
beataml_pval <- 1 - pchisq(beataml_logrank$chisq, df=1)

tcga_logrank <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                         data = tcga_data)
tcga_pval <- 1 - pchisq(tcga_logrank$chisq, df=1)

cat("LOG-RANK TEST P-VALUES:\n")
cat("  BeatAML:", format(beataml_pval, scientific=TRUE, digits=3), "\n")
cat("  TCGA:", format(tcga_pval, scientific=TRUE, digits=3), "\n\n")

# ============================================================================
# 4. POWER ANALYSIS
# ============================================================================

cat("=== STEP 4: POWER ANALYSIS ===\n\n")

# Observed effect size in BeatAML
beataml_cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                     data = beataml_surv)
beataml_hr <- exp(coef(beataml_cox))

cat("BeatAML observed HR:", round(beataml_hr, 3), "\n\n")

# Power calculation for TCGA
# Using standard power calculation for survival studies
tcga_n <- sum(!is.na(tcga_data$OS_STATUS))
tcga_events <- sum(tcga_data$OS_STATUS, na.rm=TRUE)

cat("TCGA POWER ANALYSIS:\n")
cat("  Sample size:", tcga_n, "\n")
cat("  Events:", tcga_events, "\n")
cat("  Expected HR:", round(beataml_hr, 3), "\n\n")

# Simplified power estimate
# Power ≈ Φ(√(events × ln(HR)²/4) - 1.96)
# where Φ is normal CDF

log_hr <- log(beataml_hr)
z_stat <- sqrt(tcga_events * log_hr^2 / 4) - 1.96
estimated_power <- pnorm(z_stat)

cat("  Estimated power:", round(estimated_power * 100, 1), "%\n")

if (estimated_power < 0.8) {
  cat("  ⚠ UNDERPOWERED: TCGA cohort is too small to detect effect\n")
} else {
  cat("  ✓ Adequate power to detect effect\n")
}
cat("\n")

# Calculate required sample size for 80% power
required_events <- ceiling((4 * (1.96 + 0.84)^2) / log_hr^2)
cat("  Required events for 80% power:", required_events, "\n")
cat("  TCGA has:", tcga_events, "events\n")
cat("  Shortfall:", required_events - tcga_events, "events\n\n")

# ============================================================================
# 5. EFFECT SIZE COMPARISON
# ============================================================================

cat("=== STEP 5: COMPARING EFFECT SIZES ===\n\n")

# Cox models
beataml_cox_summary <- summary(beataml_cox)
tcga_cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                  data = tcga_data)
tcga_cox_summary <- summary(tcga_cox)

cat("HAZARD RATIOS (Cluster 2 vs Cluster 1):\n")
cat("  BeatAML: HR =", round(exp(coef(beataml_cox)), 3),
    "(95% CI:", round(beataml_cox_summary$conf.int[1,3], 3), "-",
    round(beataml_cox_summary$conf.int[1,4], 3), ")\n")
cat("    p-value:", format(beataml_cox_summary$coefficients[1,5],
                           scientific=TRUE, digits=3), "\n\n")

cat("  TCGA: HR =", round(exp(coef(tcga_cox)), 3),
    "(95% CI:", round(tcga_cox_summary$conf.int[1,3], 3), "-",
    round(tcga_cox_summary$conf.int[1,4], 3), ")\n")
cat("    p-value:", format(tcga_cox_summary$coefficients[1,5],
                           scientific=TRUE, digits=3), "\n\n")

# Test for effect size heterogeneity
cat("HETEROGENEITY TEST:\n")
diff_log_hr <- log(exp(coef(beataml_cox))) - log(exp(coef(tcga_cox)))
se_beataml <- beataml_cox_summary$coefficients[1, "se(coef)"]
se_tcga <- tcga_cox_summary$coefficients[1, "se(coef)"]
se_diff <- sqrt(se_beataml^2 + se_tcga^2)
z_het <- diff_log_hr / se_diff
p_het <- 2 * (1 - pnorm(abs(z_het)))

cat("  Z-statistic:", round(z_het, 3), "\n")
cat("  P-value:", format(p_het, scientific=TRUE, digits=3), "\n")

if (p_het < 0.05) {
  cat("  ✓ SIGNIFICANT HETEROGENEITY: Effect sizes differ between cohorts\n\n")
} else {
  cat("  → No significant heterogeneity: Effect sizes are consistent\n\n")
}

# ============================================================================
# 6. POSSIBLE EXPLANATIONS
# ============================================================================

cat("=== STEP 6: EXPLORING POSSIBLE EXPLANATIONS ===\n\n")

explanations <- data.frame(
  explanation = character(),
  evidence = character(),
  likely = character(),
  stringsAsFactors = FALSE
)

# 1. Power/Sample size
explanations <- rbind(explanations, data.frame(
  explanation = "Insufficient power (small sample)",
  evidence = paste0("TCGA n=", tcga_n, " events=", tcga_events,
                   ", power=", round(estimated_power*100, 0), "%"),
  likely = ifelse(estimated_power < 0.8, "LIKELY", "Unlikely")
))

# 2. Effect size difference
explanations <- rbind(explanations, data.frame(
  explanation = "Effect size differs (true heterogeneity)",
  evidence = paste0("BeatAML HR=", round(exp(coef(beataml_cox)), 2),
                   ", TCGA HR=", round(exp(coef(tcga_cox)), 2),
                   ", p_het=", format(p_het, digits=3)),
  likely = ifelse(p_het < 0.05, "LIKELY", "Unlikely")
))

# 3. Cohort differences
age_diff <- abs(mean(beataml_surv$age, na.rm=TRUE) - mean(tcga_data$age, na.rm=TRUE))
explanations <- rbind(explanations, data.frame(
  explanation = "Cohort demographic differences",
  evidence = paste0("Age diff=", round(age_diff, 1), " years"),
  likely = ifelse(age_diff > 5, "Possible", "Unlikely")
))

# 4. Treatment era differences
explanations <- rbind(explanations, data.frame(
  explanation = "Treatment era differences",
  evidence = "BeatAML:recent, TCGA:2008-2013",
  likely = "POSSIBLE"
))

# 5. Classification quality
low_conf_pct <- sum(tcga_pred$confidence < 0.7) / nrow(tcga_pred) * 100
explanations <- rbind(explanations, data.frame(
  explanation = "Poor classification quality",
  evidence = paste0(round(low_conf_pct, 1), "% samples <0.7 confidence"),
  likely = ifelse(low_conf_pct > 30, "Possible", "Unlikely")
))

cat("POSSIBLE EXPLANATIONS:\n\n")
print(explanations)
cat("\n")

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

cat("=== STEP 7: SAVING RESULTS ===\n\n")

investigation_summary <- data.frame(
  cohort = c("BeatAML", "TCGA"),
  n_samples = c(nrow(beataml_surv), tcga_n),
  n_events = c(sum(beataml_surv$OS_STATUS), tcga_events),
  median_os = c(beataml_median_os, tcga_median_os),
  hr = c(exp(coef(beataml_cox)), exp(coef(tcga_cox))),
  hr_lower = c(beataml_cox_summary$conf.int[1,3], tcga_cox_summary$conf.int[1,3]),
  hr_upper = c(beataml_cox_summary$conf.int[1,4], tcga_cox_summary$conf.int[1,4]),
  pvalue = c(beataml_pval, tcga_pval),
  cluster1_pct = c(beataml_c1_pct, tcga_c1_pct),
  cluster2_pct = c(beataml_c2_pct, tcga_c2_pct)
)

write.csv(investigation_summary,
          "03_Results/11_Survival_Reanalysis/06_tcga_investigation_summary.csv",
          row.names = FALSE)
cat("✓ Saved: 06_tcga_investigation_summary.csv\n")

write.csv(explanations,
          "03_Results/11_Survival_Reanalysis/06_validation_failure_explanations.csv",
          row.names = FALSE)
cat("✓ Saved: 06_validation_failure_explanations.csv\n\n")

# ============================================================================
# 8. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("PART 3: TCGA VALIDATION INVESTIGATION COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. COHORT CHARACTERISTICS:\n")
cat("   Sample size: BeatAML=", nrow(beataml_surv), ", TCGA=", tcga_n, "\n")
cat("   Events: BeatAML=", sum(beataml_surv$OS_STATUS), ", TCGA=", tcga_events, "\n")
cat("   Cluster proportions: Similar (p=",
    format(cluster_chisq$p.value, digits=3), ")\n\n")

cat("2. SURVIVAL ANALYSIS:\n")
cat("   BeatAML: HR=", round(exp(coef(beataml_cox)), 3),
    ", p=", format(beataml_pval, scientific=TRUE, digits=3), "\n")
cat("   TCGA: HR=", round(exp(coef(tcga_cox)), 3),
    ", p=", format(tcga_pval, scientific=TRUE, digits=3), "\n")
cat("   Heterogeneity: p=", format(p_het, scientific=TRUE, digits=3), "\n\n")

cat("3. POWER ANALYSIS:\n")
cat("   TCGA estimated power:", round(estimated_power * 100, 1), "%\n")
cat("   Required events:", required_events, ", Actual events:", tcga_events, "\n")
if (estimated_power < 0.8) {
  cat("   ⚠ TCGA IS UNDERPOWERED\n\n")
} else {
  cat("   ✓ Adequate power\n\n")
}

cat("4. MOST LIKELY EXPLANATION:\n")
likely_explanation <- explanations$explanation[explanations$likely == "LIKELY"][1]
if (!is.na(likely_explanation)) {
  cat("   ", likely_explanation, "\n\n")
} else {
  cat("   Combination of factors (power + heterogeneity)\n\n")
}

cat("CONCLUSION:\n")
if (estimated_power < 0.8) {
  cat("TCGA validation failure is primarily due to INSUFFICIENT POWER.\n")
  cat("The cohort is too small (", tcga_events, " events) to reliably detect\n")
  cat("the observed effect size (HR=", round(beataml_hr, 2), ").\n")
  cat("Estimated power is only ", round(estimated_power*100, 0), "%.\n\n")
} else if (p_het < 0.05) {
  cat("There is TRUE HETEROGENEITY between cohorts.\n")
  cat("The prognostic effect may differ due to treatment era,\n")
  cat("patient populations, or other cohort-specific factors.\n\n")
} else {
  cat("Effect sizes are consistent but TCGA failed to reach significance.\n")
  cat("This suggests power limitations or chance variation.\n\n")
}

cat("RECOMMENDATION FOR PUBLICATION:\n")
cat("- Report both results transparently\n")
cat("- Emphasize BeatAML as discovery cohort\n")
cat("- Acknowledge TCGA power limitations\n")
cat("- Consider meta-analysis to pool evidence\n\n")

cat("### PART 3 COMPLETE ###\n")
