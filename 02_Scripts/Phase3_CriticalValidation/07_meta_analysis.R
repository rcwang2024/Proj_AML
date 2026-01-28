# PHASE 3 - PART 5: META-ANALYSIS
# Pool evidence from BeatAML and TCGA cohorts using fixed and random effects models

library(tidyverse)
library(survival)

cat("=== PHASE 3: PART 5 - META-ANALYSIS ===\n\n")

# ============================================================================
# 1. LOAD COHORT RESULTS
# ============================================================================

cat("=== STEP 1: LOADING COHORT RESULTS ===\n\n")

# Load investigation summary from Part 3
cohort_summary <- read.csv("03_Results/11_Survival_Reanalysis/06_tcga_investigation_summary.csv")

cat("Cohort results loaded:\n")
print(cohort_summary[, c("cohort", "n_samples", "n_events", "hr", "hr_lower", "hr_upper", "pvalue")])
cat("\n")

# ============================================================================
# 2. PREPARE META-ANALYSIS DATA
# ============================================================================

cat("=== STEP 2: PREPARING META-ANALYSIS DATA ===\n\n")

# Extract log-HRs and standard errors
log_hr_beataml <- log(cohort_summary$hr[cohort_summary$cohort == "BeatAML"])
log_hr_tcga <- log(cohort_summary$hr[cohort_summary$cohort == "TCGA"])

# Calculate SE from confidence intervals
# SE = (log(upper) - log(lower)) / (2 * 1.96)
se_beataml <- (log(cohort_summary$hr_upper[cohort_summary$cohort == "BeatAML"]) -
               log(cohort_summary$hr_lower[cohort_summary$cohort == "BeatAML"])) / (2 * 1.96)
se_tcga <- (log(cohort_summary$hr_upper[cohort_summary$cohort == "TCGA"]) -
            log(cohort_summary$hr_lower[cohort_summary$cohort == "TCGA"])) / (2 * 1.96)

meta_data <- data.frame(
  study = c("BeatAML", "TCGA"),
  log_hr = c(log_hr_beataml, log_hr_tcga),
  se = c(se_beataml, se_tcga),
  n = cohort_summary$n_samples,
  events = cohort_summary$n_events
)

cat("Meta-analysis data:\n")
print(meta_data)
cat("\n")

# ============================================================================
# 3. FIXED EFFECTS META-ANALYSIS
# ============================================================================

cat("=== STEP 3: FIXED EFFECTS META-ANALYSIS ===\n\n")

# Manual calculation: inverse-variance weighted
weights_fixed <- 1 / meta_data$se^2
pooled_log_hr_fixed <- sum(meta_data$log_hr * weights_fixed) / sum(weights_fixed)
pooled_se_fixed <- sqrt(1 / sum(weights_fixed))
pooled_hr_fixed <- exp(pooled_log_hr_fixed)
pooled_lower_fixed <- exp(pooled_log_hr_fixed - 1.96 * pooled_se_fixed)
pooled_upper_fixed <- exp(pooled_log_hr_fixed + 1.96 * pooled_se_fixed)
z_fixed <- pooled_log_hr_fixed / pooled_se_fixed
p_fixed <- 2 * (1 - pnorm(abs(z_fixed)))

cat("FIXED EFFECTS MODEL (Inverse Variance):\n")
cat("  Pooled HR:", round(pooled_hr_fixed, 3), "\n")
cat("  95% CI:", round(pooled_lower_fixed, 3), "-", round(pooled_upper_fixed, 3), "\n")
cat("  Z-statistic:", round(z_fixed, 3), "\n")
cat("  P-value:", format(p_fixed, scientific=TRUE, digits=3), "\n")

if (p_fixed < 0.05) {
  cat("  ✓ SIGNIFICANT pooled effect\n\n")
} else {
  cat("  ✗ Not significant\n\n")
}

# ============================================================================
# 4. RANDOM EFFECTS META-ANALYSIS
# ============================================================================

cat("=== STEP 4: RANDOM EFFECTS META-ANALYSIS ===\n\n")

# Calculate Q statistic for heterogeneity
Q <- sum(weights_fixed * (meta_data$log_hr - pooled_log_hr_fixed)^2)
df <- length(meta_data$log_hr) - 1
Q_pvalue <- 1 - pchisq(Q, df)

cat("HETEROGENEITY TEST:\n")
cat("  Q statistic:", round(Q, 3), "\n")
cat("  df:", df, "\n")
cat("  P-value:", format(Q_pvalue, scientific=TRUE, digits=3), "\n")

# Calculate I² (proportion of variation due to heterogeneity)
I2 <- max(0, (Q - df) / Q * 100)
cat("  I² statistic:", round(I2, 1), "%\n")

if (I2 < 25) {
  cat("  → Low heterogeneity\n")
} else if (I2 < 50) {
  cat("  → Moderate heterogeneity\n")
} else {
  cat("  → High heterogeneity\n")
}
cat("\n")

# DerSimonian-Laird random effects
# τ² = max(0, (Q - df) / C)
# where C = sum(weights) - sum(weights²)/sum(weights)
C <- sum(weights_fixed) - sum(weights_fixed^2) / sum(weights_fixed)
tau2 <- max(0, (Q - df) / C)

cat("Between-study variance (τ²):", round(tau2, 4), "\n\n")

# Random effects pooled estimate
weights_random <- 1 / (meta_data$se^2 + tau2)
pooled_log_hr_random <- sum(meta_data$log_hr * weights_random) / sum(weights_random)
pooled_se_random <- sqrt(1 / sum(weights_random))
pooled_hr_random <- exp(pooled_log_hr_random)
pooled_lower_random <- exp(pooled_log_hr_random - 1.96 * pooled_se_random)
pooled_upper_random <- exp(pooled_log_hr_random + 1.96 * pooled_se_random)
z_random <- pooled_log_hr_random / pooled_se_random
p_random <- 2 * (1 - pnorm(abs(z_random)))

cat("RANDOM EFFECTS MODEL (DerSimonian-Laird):\n")
cat("  Pooled HR:", round(pooled_hr_random, 3), "\n")
cat("  95% CI:", round(pooled_lower_random, 3), "-", round(pooled_upper_random, 3), "\n")
cat("  Z-statistic:", round(z_random, 3), "\n")
cat("  P-value:", format(p_random, scientific=TRUE, digits=3), "\n")

if (p_random < 0.05) {
  cat("  ✓ SIGNIFICANT pooled effect\n\n")
} else {
  cat("  ✗ Not significant\n\n")
}

# ============================================================================
# 5. INDIVIDUAL STUDY CONTRIBUTIONS
# ============================================================================

cat("=== STEP 5: INDIVIDUAL STUDY CONTRIBUTIONS ===\n\n")

# Fixed effects weights
meta_data$weight_fixed_pct <- weights_fixed / sum(weights_fixed) * 100

# Random effects weights
meta_data$weight_random_pct <- weights_random / sum(weights_random) * 100

cat("STUDY CONTRIBUTIONS:\n")
cat("\nFixed Effects:\n")
for (i in 1:nrow(meta_data)) {
  cat("  ", meta_data$study[i], ": ", round(meta_data$weight_fixed_pct[i], 1), "%\n", sep="")
}

cat("\nRandom Effects:\n")
for (i in 1:nrow(meta_data)) {
  cat("  ", meta_data$study[i], ": ", round(meta_data$weight_random_pct[i], 1), "%\n", sep="")
}
cat("\n")

# ============================================================================
# 6. POOLED INDIVIDUAL PATIENT DATA ANALYSIS
# ============================================================================

cat("=== STEP 6: POOLED INDIVIDUAL PATIENT DATA ANALYSIS ===\n\n")

# Load raw data
beataml_surv <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
tcga_data <- read.csv("03_Results/17_TCGA_Validation/tcga_clinical_with_predictions.csv")

# Standardize column names
beataml_surv <- beataml_surv %>%
  mutate(cohort = "BeatAML") %>%
  rename(OS_MONTHS = OS_months, OS_STATUS = OS_event, cluster_assignment = cluster) %>%
  select(cohort, cluster_assignment, OS_MONTHS, OS_STATUS, age, sex)

tcga_data <- tcga_data %>%
  mutate(cohort = "TCGA") %>%
  rename(OS_MONTHS = OS_months, OS_STATUS = OS_event, cluster_assignment = predicted_cluster) %>%
  select(cohort, cluster_assignment, OS_MONTHS, OS_STATUS, age, sex)

# Pool datasets
pooled_data <- rbind(beataml_surv, tcga_data)

cat("Pooled dataset:\n")
cat("  Total samples:", nrow(pooled_data), "\n")
cat("  BeatAML:", sum(pooled_data$cohort == "BeatAML"), "\n")
cat("  TCGA:", sum(pooled_data$cohort == "TCGA"), "\n")
cat("  Total events:", sum(pooled_data$OS_STATUS), "\n\n")

# IPD meta-analysis (stratified by cohort)
cat("IPD Meta-Analysis (stratified by cohort):\n")
ipd_cox <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + strata(cohort),
                 data = pooled_data)

print(summary(ipd_cox))
cat("\n")

ipd_hr <- exp(coef(ipd_cox))
ipd_ci <- exp(confint(ipd_cox))
ipd_pval <- summary(ipd_cox)$coefficients[1, "Pr(>|z|)"]

cat("IPD POOLED ESTIMATE (Stratified Cox):\n")
cat("  Pooled HR:", round(ipd_hr, 3), "\n")
cat("  95% CI:", round(ipd_ci[1], 3), "-", round(ipd_ci[2], 3), "\n")
cat("  P-value:", format(ipd_pval, scientific=TRUE, digits=3), "\n\n")

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

cat("=== STEP 7: SAVING RESULTS ===\n\n")

meta_results <- data.frame(
  model = c("Fixed Effects", "Random Effects", "IPD Stratified Cox"),
  pooled_hr = c(pooled_hr_fixed, pooled_hr_random, ipd_hr),
  lower_95 = c(pooled_lower_fixed, pooled_lower_random, ipd_ci[1]),
  upper_95 = c(pooled_upper_fixed, pooled_upper_random, ipd_ci[2]),
  pvalue = c(p_fixed, p_random, ipd_pval),
  heterogeneity_Q = c(Q, Q, NA),
  heterogeneity_pval = c(Q_pvalue, Q_pvalue, NA),
  I2_percent = c(I2, I2, NA)
)

write.csv(meta_results,
          "03_Results/11_Survival_Reanalysis/07_meta_analysis_results.csv",
          row.names = FALSE)
cat("✓ Saved: 07_meta_analysis_results.csv\n")

study_contributions <- data.frame(
  study = meta_data$study,
  hr = exp(meta_data$log_hr),
  se = meta_data$se,
  n = meta_data$n,
  events = meta_data$events,
  weight_fixed_pct = meta_data$weight_fixed_pct,
  weight_random_pct = meta_data$weight_random_pct
)

write.csv(study_contributions,
          "03_Results/11_Survival_Reanalysis/07_study_contributions.csv",
          row.names = FALSE)
cat("✓ Saved: 07_study_contributions.csv\n\n")

# ============================================================================
# 8. VISUALIZATION
# ============================================================================

cat("=== STEP 8: CREATING FOREST PLOT ===\n\n")

# Create forest plot data
forest_data <- data.frame(
  study = c("BeatAML", "TCGA", "Pooled (Fixed)", "Pooled (Random)", "IPD Stratified"),
  hr = c(exp(meta_data$log_hr), pooled_hr_fixed, pooled_hr_random, ipd_hr),
  lower = c(exp(meta_data$log_hr - 1.96 * meta_data$se),
            pooled_lower_fixed, pooled_lower_random, ipd_ci[1]),
  upper = c(exp(meta_data$log_hr + 1.96 * meta_data$se),
            pooled_upper_fixed, pooled_upper_random, ipd_ci[2]),
  type = c(rep("Individual", 2), rep("Meta", 3)),
  weight = c(meta_data$weight_fixed_pct, NA, NA, NA)
)

forest_data$study <- factor(forest_data$study,
                             levels = rev(c("BeatAML", "TCGA", "Pooled (Fixed)",
                                            "Pooled (Random)", "IPD Stratified")))

library(ggplot2)

p_forest <- ggplot(forest_data, aes(x = hr, y = study)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_point(aes(color = type, size = type)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = type), height = 0.2) +
  scale_color_manual(values = c("Individual" = "steelblue", "Meta" = "#E7B800")) +
  scale_size_manual(values = c("Individual" = 3, "Meta" = 4)) +
  labs(
    title = "Meta-Analysis: Molecular Subtypes and Survival",
    subtitle = paste0("Pooled HR = ", round(pooled_hr_random, 2),
                     " (95% CI: ", round(pooled_lower_random, 2), "-",
                     round(pooled_upper_random, 2), "), p = ",
                     format(p_random, scientific=TRUE, digits=2)),
    x = "Hazard Ratio (Cluster 2 vs Cluster 1)",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("04_Figures/11_Survival_Reanalysis/07_meta_analysis_forest_plot.pdf",
       p_forest, width = 10, height = 6)
cat("✓ Saved: 07_meta_analysis_forest_plot.pdf\n\n")

# ============================================================================
# 9. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("PART 5: META-ANALYSIS COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. FIXED EFFECTS META-ANALYSIS:\n")
cat("   Pooled HR:", round(pooled_hr_fixed, 3), "(95% CI:", round(pooled_lower_fixed, 3),
    "-", round(pooled_upper_fixed, 3), ")\n")
cat("   P-value:", format(p_fixed, scientific=TRUE, digits=3), "\n")
if (p_fixed < 0.05) {
  cat("   ✓ Significant pooled effect\n\n")
} else {
  cat("   ✗ Not significant\n\n")
}

cat("2. RANDOM EFFECTS META-ANALYSIS:\n")
cat("   Pooled HR:", round(pooled_hr_random, 3), "(95% CI:", round(pooled_lower_random, 3),
    "-", round(pooled_upper_random, 3), ")\n")
cat("   P-value:", format(p_random, scientific=TRUE, digits=3), "\n")
if (p_random < 0.05) {
  cat("   ✓ Significant pooled effect\n\n")
} else {
  cat("   ✗ Not significant\n\n")
}

cat("3. HETEROGENEITY:\n")
cat("   Q statistic:", round(Q, 3), ", p =", format(Q_pvalue, scientific=TRUE, digits=3), "\n")
cat("   I² statistic:", round(I2, 1), "%\n")
if (I2 < 25) {
  cat("   → Low heterogeneity between studies\n\n")
} else {
  cat("   → Moderate/High heterogeneity\n\n")
}

cat("4. IPD POOLED ANALYSIS:\n")
cat("   Stratified Cox HR:", round(ipd_hr, 3), "(95% CI:", round(ipd_ci[1], 3),
    "-", round(ipd_ci[2], 3), ")\n")
cat("   P-value:", format(ipd_pval, scientific=TRUE, digits=3), "\n\n")

cat("5. STUDY CONTRIBUTIONS:\n")
cat("   BeatAML weight:", round(meta_data$weight_fixed_pct[1], 1), "%\n")
cat("   TCGA weight:", round(meta_data$weight_fixed_pct[2], 1), "%\n\n")

cat("INTERPRETATION:\n")
if (p_random < 0.05) {
  cat("Meta-analysis provides STRONGER EVIDENCE for prognostic effect by pooling\n")
  cat("both cohorts. The random effects model accounts for between-study variability\n")
  cat("and provides a more conservative estimate. The effect is significant when\n")
  cat("combining evidence across studies.\n\n")
} else {
  cat("Despite pooling evidence, the effect remains non-significant in meta-analysis.\n")
  cat("This suggests either true lack of effect or substantial heterogeneity.\n\n")
}

cat("RECOMMENDATION:\n")
cat("- Report random effects estimate (more conservative)\n")
cat("- Acknowledge BeatAML drives most of the evidence (larger sample)\n")
cat("- Note that TCGA was underpowered (35% power)\n")
cat("- Meta-analysis provides best overall estimate of effect\n\n")

cat("### PART 5 COMPLETE ###\n")
