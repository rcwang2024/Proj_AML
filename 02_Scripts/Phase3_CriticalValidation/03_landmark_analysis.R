# PHASE 3 - PART 1.3: LANDMARK ANALYSIS
# Perform landmark analysis at multiple timepoints to avoid guaranteeist time bias
# and address PH violations. Analyze survival FROM specific landmarks among patients
# who survived TO those landmarks.

library(tidyverse)
library(survival)
library(survminer)

cat("=== PHASE 3: PART 1.3 - LANDMARK ANALYSIS ===\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("=== STEP 1: LOADING DATA ===\n\n")

survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Rename columns for consistency
survival_data <- survival_data %>%
  rename(
    cluster_assignment = cluster,
    OS_MONTHS = OS_months,
    OS_STATUS = OS_event
  )

cat("Loaded survival data:", nrow(survival_data), "samples\n")
cat("Total events:", sum(survival_data$OS_STATUS), "\n")
cat("Median follow-up:", round(median(survival_data$OS_MONTHS), 1), "months\n\n")

# ============================================================================
# 2. DEFINE LANDMARK TIMEPOINTS
# ============================================================================

cat("=== STEP 2: DEFINING LANDMARK TIMEPOINTS ===\n\n")

# Define landmark timepoints (in months)
landmarks <- c(0, 6, 12, 18, 24, 36)

cat("Landmark timepoints:", paste(landmarks, collapse=", "), "months\n\n")

# ============================================================================
# 3. PERFORM LANDMARK ANALYSES
# ============================================================================

cat("=== STEP 3: PERFORMING LANDMARK ANALYSES ===\n\n")

# Initialize results storage
landmark_results <- data.frame()

# Loop through each landmark
for (lm in landmarks) {

  cat("### LANDMARK:", lm, "MONTHS ###\n\n")

  # Select patients who survived past landmark
  lm_data <- survival_data %>%
    filter(OS_MONTHS > lm)

  cat("Patients alive at", lm, "months:", nrow(lm_data), "\n")

  if (nrow(lm_data) < 50) {
    cat("  → Too few patients, skipping this landmark\n\n")
    next
  }

  # Adjust survival times (time since landmark)
  lm_data$OS_MONTHS_ADJUSTED <- lm_data$OS_MONTHS - lm

  # Check cluster distribution
  cluster_dist <- table(lm_data$cluster_assignment)
  cat("  Cluster 1:", cluster_dist[1], "\n")
  cat("  Cluster 2:", cluster_dist[2], "\n")

  # Events after landmark
  events_after <- sum(lm_data$OS_STATUS)
  cat("  Events after landmark:", events_after, "\n")

  if (events_after < 20) {
    cat("  → Too few events, skipping this landmark\n\n")
    next
  }

  # Create survival object (from landmark forward)
  surv_obj_lm <- Surv(time = lm_data$OS_MONTHS_ADJUSTED, event = lm_data$OS_STATUS)

  # Log-rank test
  logrank_lm <- survdiff(surv_obj_lm ~ cluster_assignment, data = lm_data)
  logrank_pval_lm <- 1 - pchisq(logrank_lm$chisq, df = 1)

  cat("  Log-rank p-value:", format(logrank_pval_lm, scientific = TRUE, digits = 3), "\n")

  # Cox model from landmark
  cox_lm <- coxph(surv_obj_lm ~ cluster_assignment, data = lm_data)

  # Extract HR and CI
  hr_lm <- exp(coef(cox_lm))
  hr_ci_lm <- exp(confint(cox_lm))
  cox_pval_lm <- summary(cox_lm)$coefficients[1, "Pr(>|z|)"]

  cat("  Hazard Ratio:", round(hr_lm, 3),
      "(95% CI:", round(hr_ci_lm[1], 3), "-", round(hr_ci_lm[2], 3), ")\n")
  cat("  Cox p-value:", format(cox_pval_lm, scientific = TRUE, digits = 3), "\n")

  # Test PH assumption at this landmark
  cox_zph_lm <- cox.zph(cox_lm)
  ph_pval_lm <- cox_zph_lm$table["GLOBAL", "p"]

  cat("  PH assumption p-value:", format(ph_pval_lm, scientific = TRUE, digits = 3), "\n")

  if (ph_pval_lm < 0.05) {
    cat("  ✗ PH violated even at this landmark\n")
  } else {
    cat("  ✓ PH assumption holds at this landmark\n")
  }

  # Median survival from landmark
  km_lm <- survfit(surv_obj_lm ~ cluster_assignment, data = lm_data)
  median_surv_lm <- summary(km_lm)$table[, "median"]

  cat("  Median survival from landmark:\n")
  cat("    Cluster 1:", round(median_surv_lm[1], 1), "months\n")
  cat("    Cluster 2:", round(median_surv_lm[2], 1), "months\n")

  # Store results
  landmark_results <- rbind(landmark_results, data.frame(
    landmark_months = lm,
    n_patients = nrow(lm_data),
    n_cluster1 = cluster_dist[1],
    n_cluster2 = cluster_dist[2],
    n_events = events_after,
    logrank_pvalue = logrank_pval_lm,
    hazard_ratio = hr_lm,
    hr_95_lower = hr_ci_lm[1],
    hr_95_upper = hr_ci_lm[2],
    cox_pvalue = cox_pval_lm,
    ph_assumption_pvalue = ph_pval_lm,
    ph_holds = ph_pval_lm >= 0.05,
    median_surv_cluster1 = median_surv_lm[1],
    median_surv_cluster2 = median_surv_lm[2]
  ))

  cat("\n")
}

# ============================================================================
# 4. ANALYZE CONSISTENCY ACROSS LANDMARKS
# ============================================================================

cat("=== STEP 4: ANALYZING CONSISTENCY ACROSS LANDMARKS ===\n\n")

cat("Landmark Analysis Summary:\n")
print(landmark_results[, c("landmark_months", "n_patients", "n_events",
                            "hazard_ratio", "logrank_pvalue")])
cat("\n")

# Check if results are consistent
significant_count <- sum(landmark_results$logrank_pvalue < 0.05, na.rm = TRUE)
total_landmarks <- nrow(landmark_results)

cat("Statistical Consistency:\n")
cat("  Significant landmarks:", significant_count, "of", total_landmarks, "\n")
cat("  Proportion:", round(significant_count / total_landmarks, 2), "\n\n")

if (significant_count / total_landmarks >= 0.5) {
  cat("  ✓ CONSISTENT: Prognostic effect observed across most landmarks\n")
} else {
  cat("  ✗ INCONSISTENT: Prognostic effect varies across landmarks\n")
}
cat("\n")

# Check HR stability
hr_range <- range(landmark_results$hazard_ratio, na.rm = TRUE)
hr_cv <- sd(landmark_results$hazard_ratio, na.rm = TRUE) /
         mean(landmark_results$hazard_ratio, na.rm = TRUE)

cat("Hazard Ratio Stability:\n")
cat("  Range:", round(hr_range[1], 2), "-", round(hr_range[2], 2), "\n")
cat("  Coefficient of variation:", round(hr_cv, 3), "\n")

if (hr_cv < 0.3) {
  cat("  ✓ STABLE: HRs are consistent across landmarks\n")
} else {
  cat("  ⚠ VARIABLE: HRs show substantial variation\n")
}
cat("\n")

# ============================================================================
# 5. CONDITIONAL SURVIVAL PROBABILITIES
# ============================================================================

cat("=== STEP 5: CONDITIONAL SURVIVAL PROBABILITIES ===\n\n")

# Calculate 12-month conditional survival from each landmark
conditional_surv <- data.frame()

for (lm in landmarks) {

  lm_data <- survival_data %>%
    filter(OS_MONTHS > lm)

  if (nrow(lm_data) < 50) next

  lm_data$OS_MONTHS_ADJUSTED <- lm_data$OS_MONTHS - lm

  surv_obj_lm <- Surv(time = lm_data$OS_MONTHS_ADJUSTED, event = lm_data$OS_STATUS)
  km_lm <- survfit(surv_obj_lm ~ cluster_assignment, data = lm_data)

  # Get 12-month survival probability from landmark
  surv_12m <- summary(km_lm, times = 12)

  if (length(surv_12m$surv) >= 2) {
    conditional_surv <- rbind(conditional_surv, data.frame(
      landmark_months = lm,
      cluster = c(1, 2),
      surv_12m_conditional = surv_12m$surv,
      lower_95 = surv_12m$lower,
      upper_95 = surv_12m$upper
    ))
  }
}

cat("12-Month Conditional Survival Probabilities:\n")
print(conditional_surv)
cat("\n")

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

cat("=== STEP 6: SAVING RESULTS ===\n\n")

write.csv(landmark_results,
          "03_Results/11_Survival_Reanalysis/03_landmark_analysis_results.csv",
          row.names = FALSE)
cat("✓ Saved: 03_landmark_analysis_results.csv\n")

write.csv(conditional_surv,
          "03_Results/11_Survival_Reanalysis/03_conditional_survival.csv",
          row.names = FALSE)
cat("✓ Saved: 03_conditional_survival.csv\n\n")

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

cat("=== STEP 7: CREATING VISUALIZATIONS ===\n\n")

# Plot 1: HR across landmarks
p_hr_landmark <- ggplot(landmark_results,
                        aes(x = landmark_months, y = hazard_ratio)) +
  geom_line(size = 1.2, color = "#E7B800") +
  geom_point(size = 3, color = "#E7B800") +
  geom_errorbar(aes(ymin = hr_95_lower, ymax = hr_95_upper),
                width = 2, color = "#E7B800") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Hazard Ratios Across Landmark Timepoints",
    subtitle = "Cluster 2 vs Cluster 1 (from each landmark forward)",
    x = "Landmark Timepoint (months)",
    y = "Hazard Ratio (95% CI)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("04_Figures/11_Survival_Reanalysis/03_hr_by_landmark.pdf",
       p_hr_landmark, width = 10, height = 6)
cat("✓ Saved: 03_hr_by_landmark.pdf\n")

# Plot 2: P-values across landmarks
landmark_results$significance <- ifelse(landmark_results$logrank_pvalue < 0.05,
                                       "Significant", "Not Significant")

p_pval_landmark <- ggplot(landmark_results,
                          aes(x = landmark_months, y = -log10(logrank_pvalue),
                              color = significance)) +
  geom_line(size = 1.2, color = "gray50") +
  geom_point(aes(color = significance), size = 4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Significant" = "#00BA38",
                                 "Not Significant" = "#F8766D")) +
  labs(
    title = "Statistical Significance Across Landmarks",
    subtitle = "Log-rank test p-values",
    x = "Landmark Timepoint (months)",
    y = "-log10(P-value)",
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("04_Figures/11_Survival_Reanalysis/03_pvalue_by_landmark.pdf",
       p_pval_landmark, width = 10, height = 6)
cat("✓ Saved: 03_pvalue_by_landmark.pdf\n")

# Plot 3: Conditional 12-month survival
if (nrow(conditional_surv) > 0) {
  p_conditional <- ggplot(conditional_surv,
                          aes(x = landmark_months, y = surv_12m_conditional,
                              color = factor(cluster), group = factor(cluster))) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 2) +
    scale_color_manual(values = c("1" = "#2E9FDF", "2" = "#E7B800"),
                       labels = c("Cluster 1 (Proliferative)",
                                  "Cluster 2 (Immune-Inflammatory)")) +
    labs(
      title = "12-Month Conditional Survival Probabilities",
      subtitle = "Probability of surviving 12 more months given survival to landmark",
      x = "Landmark Timepoint (months)",
      y = "12-Month Conditional Survival Probability",
      color = "Molecular Subtype"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave("04_Figures/11_Survival_Reanalysis/03_conditional_survival.pdf",
         p_conditional, width = 10, height = 6)
  cat("✓ Saved: 03_conditional_survival.pdf\n")
}

cat("\n")

# ============================================================================
# 8. KAPLAN-MEIER PLOTS FOR KEY LANDMARKS
# ============================================================================

cat("=== STEP 8: CREATING KM PLOTS FOR KEY LANDMARKS ===\n\n")

# Generate KM plots for 0, 12, and 24 month landmarks
key_landmarks <- c(0, 12, 24)

pdf("04_Figures/11_Survival_Reanalysis/03_km_plots_by_landmark.pdf",
    width = 12, height = 10)

par(mfrow = c(2, 2))

for (lm in key_landmarks) {

  lm_data <- survival_data %>%
    filter(OS_MONTHS > lm)

  if (nrow(lm_data) < 50) next

  lm_data$OS_MONTHS_ADJUSTED <- lm_data$OS_MONTHS - lm
  lm_data$cluster_label <- factor(lm_data$cluster_assignment,
                                   labels = c("Cluster 1", "Cluster 2"))

  surv_obj_lm <- Surv(time = lm_data$OS_MONTHS_ADJUSTED, event = lm_data$OS_STATUS)
  km_lm <- survfit(surv_obj_lm ~ cluster_label, data = lm_data)

  # Get p-value
  logrank_lm <- survdiff(surv_obj_lm ~ cluster_assignment, data = lm_data)
  pval_lm <- 1 - pchisq(logrank_lm$chisq, df = 1)

  p_km_lm <- ggsurvplot(
    km_lm,
    data = lm_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = FALSE,
    legend.title = "",
    legend.labs = c("Cluster 1 (Proliferative)", "Cluster 2 (Immune-Inflammatory)"),
    palette = c("#2E9FDF", "#E7B800"),
    title = paste0("Landmark: ", lm, " months (n=", nrow(lm_data), ")"),
    xlab = paste0("Time since ", lm, "-month landmark (months)"),
    ylab = "Survival Probability",
    ggtheme = theme_bw(base_size = 11)
  )

  print(p_km_lm)
}

dev.off()

cat("✓ Saved: 03_km_plots_by_landmark.pdf\n\n")

# ============================================================================
# 9. SUMMARY
# ============================================================================

cat("=== SUMMARY AND INTERPRETATION ===\n\n")

cat("PART 1.3: LANDMARK ANALYSIS COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. LANDMARK ANALYSES PERFORMED:\n")
cat("   Total landmarks:", total_landmarks, "\n")
cat("   Significant landmarks:", significant_count,
    "(", round(100 * significant_count / total_landmarks, 1), "%)\\n\n")

cat("2. HAZARD RATIO TRENDS:\n")
cat("   HR range:", round(hr_range[1], 2), "-", round(hr_range[2], 2), "\n")
cat("   Mean HR:", round(mean(landmark_results$hazard_ratio, na.rm = TRUE), 2), "\n")
cat("   Stability (CV):", round(hr_cv, 3), "\n\n")

cat("3. PROPORTIONAL HAZARDS AT LANDMARKS:\n")
ph_ok_count <- sum(landmark_results$ph_holds, na.rm = TRUE)
cat("   Landmarks where PH holds:", ph_ok_count, "of", total_landmarks, "\n")
if (ph_ok_count / total_landmarks >= 0.7) {
  cat("   ✓ PH assumption generally satisfied from landmarks\n")
} else {
  cat("   ⚠ PH issues persist even at landmarks\n")
}
cat("\n")

cat("INTERPRETATION:\n")
cat("Landmark analysis eliminates guaranteeist time bias by only including\n")
cat("patients who survived to each landmark timepoint. This approach is robust\n")
cat("to PH violations and provides time-conditional assessments of prognosis.\n\n")

if (significant_count / total_landmarks >= 0.5) {
  cat("The molecular subtypes show CONSISTENT prognostic value across multiple\n")
  cat("landmark timepoints, supporting their clinical utility for risk stratification.\n\n")
} else {
  cat("Prognostic value varies across landmarks, suggesting the effect may be\n")
  cat("time-dependent or limited to specific periods of follow-up.\n\n")
}

cat("NEXT STEPS:\n")
cat("  → Part 1.4: RMST analysis (assumption-free survival measure)\n\n")

cat("### PART 1.3 COMPLETE ###\n")
