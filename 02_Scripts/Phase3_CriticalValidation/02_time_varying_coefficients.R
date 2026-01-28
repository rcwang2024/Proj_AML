# PHASE 3 - PART 1.2: TIME-VARYING COEFFICIENT MODEL
# Model how hazard ratio changes over time when PH assumption is violated
# Uses time-dependent coefficients to capture temporal variation

library(tidyverse)
library(survival)
library(survminer)

cat("=== PHASE 3: PART 1.2 - TIME-VARYING COEFFICIENT MODEL ===\n\n")

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
cat("Cluster 1:", sum(survival_data$cluster_assignment == 1), "\n")
cat("Cluster 2:", sum(survival_data$cluster_assignment == 2), "\n")
cat("Events:", sum(survival_data$OS_STATUS), "\n\n")

# ============================================================================
# 2. TIME-VARYING COEFFICIENT MODEL
# ============================================================================

cat("=== STEP 2: FITTING TIME-VARYING COEFFICIENT MODEL ===\n\n")

# Method 1: Time-dependent coefficient using tt() function
# This allows the coefficient to vary as a function of time

cat("Fitting Cox model with time-varying coefficient...\n")
cat("Formula: coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster + tt(cluster))\n")
cat("where tt(cluster) = cluster * log(time)\n\n")

# Create time-dependent Cox model
# The coefficient β(t) = β0 + β1*f(t)
# Here we use f(t) = log(t) as a common choice
time_varying_cox <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + tt(cluster_assignment),
  data = survival_data,
  tt = function(x, t, ...) x * log(t)
)

cat("Time-Varying Cox Model Results:\n")
print(summary(time_varying_cox))
cat("\n")

# Extract coefficients
beta0 <- coef(time_varying_cox)[1]  # Main effect
beta1 <- coef(time_varying_cox)[2]  # Time interaction

cat("Time-Varying Hazard Ratio Model:\n")
cat("  HR(t) = exp(β₀ + β₁ × log(t))\n")
cat("  β₀ (main effect):", round(beta0, 4), "\n")
cat("  β₁ (time interaction):", round(beta1, 4), "\n\n")

# ============================================================================
# 3. CALCULATE HAZARD RATIOS OVER TIME
# ============================================================================

cat("=== STEP 3: CALCULATING TIME-SPECIFIC HAZARD RATIOS ===\n\n")

# Calculate HR at different timepoints
timepoints <- c(6, 12, 24, 36, 48, 60)

hr_over_time <- data.frame(
  time_months = timepoints,
  log_hr = beta0 + beta1 * log(timepoints),
  hr = exp(beta0 + beta1 * log(timepoints))
)

# Calculate confidence intervals
se_beta0 <- summary(time_varying_cox)$coefficients[1, "se(coef)"]
se_beta1 <- summary(time_varying_cox)$coefficients[2, "se(coef)"]
cov_beta <- vcov(time_varying_cox)[1, 2]

hr_over_time$se_log_hr <- sqrt(
  se_beta0^2 +
  (log(timepoints))^2 * se_beta1^2 +
  2 * log(timepoints) * cov_beta
)

hr_over_time$hr_lower <- exp(hr_over_time$log_hr - 1.96 * hr_over_time$se_log_hr)
hr_over_time$hr_upper <- exp(hr_over_time$log_hr + 1.96 * hr_over_time$se_log_hr)

cat("Hazard Ratios at Key Timepoints (Cluster 2 vs Cluster 1):\n")
print(hr_over_time[, c("time_months", "hr", "hr_lower", "hr_upper")])
cat("\n")

# Interpret trend
if (beta1 > 0) {
  cat("TREND: Hazard ratio INCREASES over time\n")
  cat("  → Cluster 2's higher risk becomes MORE pronounced with longer follow-up\n\n")
} else if (beta1 < 0) {
  cat("TREND: Hazard ratio DECREASES over time\n")
  cat("  → Cluster 2's higher risk becomes LESS pronounced with longer follow-up\n\n")
} else {
  cat("TREND: Hazard ratio is CONSTANT over time\n")
  cat("  → PH assumption may actually hold\n\n")
}

# Test if time interaction is significant
pval_interaction <- summary(time_varying_cox)$coefficients[2, "Pr(>|z|)"]
cat("Statistical Significance of Time Interaction:\n")
cat("  P-value:", format(pval_interaction, scientific = TRUE, digits = 3), "\n")

if (pval_interaction < 0.05) {
  cat("  ✓ SIGNIFICANT: Hazard ratio varies significantly over time\n")
  cat("  → Confirms PH violation\n")
} else {
  cat("  ✗ NOT SIGNIFICANT: No strong evidence of time-varying effect\n")
  cat("  → PH violation may be minor\n")
}
cat("\n")

# ============================================================================
# 4. ALTERNATIVE TIME FUNCTIONS
# ============================================================================

cat("=== STEP 4: TESTING ALTERNATIVE TIME FUNCTIONS ===\n\n")

# Test different time transformation functions

# Linear time
cat("Model 1: Linear time (tt = cluster * t)\n")
model_linear <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + tt(cluster_assignment),
  data = survival_data,
  tt = function(x, t, ...) x * t
)
cat("  Time interaction p-value:", format(summary(model_linear)$coefficients[2, "Pr(>|z|)"], scientific = TRUE, digits = 3), "\n")
cat("  AIC:", AIC(model_linear), "\n\n")

# Square root time
cat("Model 2: Square root time (tt = cluster * sqrt(t))\n")
model_sqrt <- coxph(
  Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment + tt(cluster_assignment),
  data = survival_data,
  tt = function(x, t, ...) x * sqrt(t)
)
cat("  Time interaction p-value:", format(summary(model_sqrt)$coefficients[2, "Pr(>|z|)"], scientific = TRUE, digits = 3), "\n")
cat("  AIC:", AIC(model_sqrt), "\n\n")

# Log time (original)
cat("Model 3: Log time (tt = cluster * log(t)) [Original]\n")
cat("  Time interaction p-value:", format(pval_interaction, scientific = TRUE, digits = 3), "\n")
cat("  AIC:", AIC(time_varying_cox), "\n\n")

# Select best model by AIC
aic_values <- c(
  linear = AIC(model_linear),
  sqrt = AIC(model_sqrt),
  log = AIC(time_varying_cox)
)

best_model <- names(which.min(aic_values))
cat("Best-fitting time function:", best_model, "(lowest AIC)\n\n")

# ============================================================================
# 5. PIECEWISE CONSTANT HAZARDS
# ============================================================================

cat("=== STEP 5: PIECEWISE CONSTANT HAZARDS MODEL ===\n\n")

# Divide follow-up into time intervals
# Test if HR differs between early and late periods

cat("Splitting follow-up at 12 months (early vs late)\n")

# Create time-split dataset
# Early period: 0-12 months
# Late period: >12 months

# Filter out any patients with 0 survival time (if any)
survival_data_filtered <- survival_data %>%
  filter(OS_MONTHS > 0)

cat("Samples after filtering zero times:", nrow(survival_data_filtered), "\n")

# For each patient, create records for each period they experienced
# Add small epsilon to avoid zero time issues
survival_split <- survSplit(
  Surv(OS_MONTHS, OS_STATUS) ~ .,
  data = survival_data_filtered,
  cut = c(12),
  episode = "time_period",
  start = "tstart",
  end = "tstop"
)

survival_split$period <- factor(survival_split$time_period,
                                labels = c("Early (0-12m)", "Late (>12m)"))

# Fit model with period interaction
piecewise_cox <- coxph(
  Surv(tstart, tstop, OS_STATUS) ~ cluster_assignment * period,
  data = survival_split
)

cat("Piecewise Cox Model Results:\n")
print(summary(piecewise_cox))
cat("\n")

# Extract HRs for each period
# Early period: exp(β_cluster)
# Late period: exp(β_cluster + β_interaction)
beta_cluster <- coef(piecewise_cox)["cluster_assignment"]
beta_interaction <- coef(piecewise_cox)["cluster_assignment:periodLate (>12m)"]

hr_early <- exp(beta_cluster)
hr_late <- exp(beta_cluster + beta_interaction)

cat("Hazard Ratios by Time Period:\n")
cat("  Early period (0-12 months):", round(hr_early, 3), "\n")
cat("  Late period (>12 months):", round(hr_late, 3), "\n")
cat("  Ratio (Late/Early):", round(hr_late / hr_early, 3), "\n\n")

# Test interaction significance
pval_period_interaction <- summary(piecewise_cox)$coefficients["cluster_assignment:periodLate (>12m)", "Pr(>|z|)"]

cat("Period Interaction Test:\n")
cat("  P-value:", format(pval_period_interaction, scientific = TRUE, digits = 3), "\n")

if (pval_period_interaction < 0.05) {
  cat("  ✓ SIGNIFICANT: Hazard ratio differs between early and late periods\n")
} else {
  cat("  ✗ NOT SIGNIFICANT: Hazard ratio is similar across time periods\n")
}
cat("\n")

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

cat("=== STEP 6: SAVING RESULTS ===\n\n")

# Save time-varying HR results
write.csv(hr_over_time,
          "03_Results/11_Survival_Reanalysis/02_time_varying_hr.csv",
          row.names = FALSE)
cat("✓ Saved: 02_time_varying_hr.csv\n")

# Save model comparison
model_comparison <- data.frame(
  model = c("Linear time", "Square root time", "Log time"),
  aic = aic_values,
  time_interaction_pvalue = c(
    summary(model_linear)$coefficients[2, "Pr(>|z|)"],
    summary(model_sqrt)$coefficients[2, "Pr(>|z|)"],
    pval_interaction
  ),
  best_model = names(aic_values) == best_model
)

write.csv(model_comparison,
          "03_Results/11_Survival_Reanalysis/02_time_function_comparison.csv",
          row.names = FALSE)
cat("✓ Saved: 02_time_function_comparison.csv\n")

# Save piecewise results
piecewise_results <- data.frame(
  period = c("Early (0-12m)", "Late (>12m)"),
  hr = c(hr_early, hr_late),
  interaction_pvalue = pval_period_interaction
)

write.csv(piecewise_results,
          "03_Results/11_Survival_Reanalysis/02_piecewise_hr.csv",
          row.names = FALSE)
cat("✓ Saved: 02_piecewise_hr.csv\n\n")

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

cat("=== STEP 7: CREATING VISUALIZATIONS ===\n\n")

# Plot 1: Hazard ratio over time
p_hr_time <- ggplot(hr_over_time, aes(x = time_months, y = hr)) +
  geom_line(size = 1.2, color = "#E7B800") +
  geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper), alpha = 0.2, fill = "#E7B800") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Time-Varying Hazard Ratio (Cluster 2 vs Cluster 1)",
    subtitle = paste0("HR(t) = exp(", round(beta0, 3), " + ", round(beta1, 3), " × log(t))"),
    x = "Time (months)",
    y = "Hazard Ratio"
  ) +
  scale_x_continuous(breaks = timepoints) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("04_Figures/11_Survival_Reanalysis/02_hr_over_time.pdf",
       p_hr_time, width = 10, height = 6)
cat("✓ Saved: 02_hr_over_time.pdf\n")

# Plot 2: Piecewise HRs
piecewise_plot_data <- data.frame(
  period = c("Early\n(0-12m)", "Late\n(>12m)"),
  hr = c(hr_early, hr_late)
)

p_piecewise <- ggplot(piecewise_plot_data, aes(x = period, y = hr, fill = period)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = round(hr, 2)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#2E9FDF", "#E7B800")) +
  labs(
    title = "Piecewise Hazard Ratios by Time Period",
    x = "Time Period",
    y = "Hazard Ratio (Cluster 2 vs Cluster 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave("04_Figures/11_Survival_Reanalysis/02_piecewise_hr.pdf",
       p_piecewise, width = 8, height = 6)
cat("✓ Saved: 02_piecewise_hr.pdf\n\n")

# ============================================================================
# 8. SUMMARY
# ============================================================================

cat("=== SUMMARY AND INTERPRETATION ===\n\n")

cat("PART 1.2: TIME-VARYING COEFFICIENT MODEL COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. TIME-VARYING HAZARD RATIO:\n")
cat("   HR at 6 months:", round(hr_over_time$hr[1], 3), "\n")
cat("   HR at 12 months:", round(hr_over_time$hr[2], 3), "\n")
cat("   HR at 24 months:", round(hr_over_time$hr[3], 3), "\n")
cat("   HR at 36 months:", round(hr_over_time$hr[4], 3), "\n\n")

cat("2. TEMPORAL TREND:\n")
if (beta1 > 0) {
  cat("   ↑ Hazard ratio INCREASES over time\n")
  cat("   → Risk difference between clusters grows with longer follow-up\n\n")
} else if (beta1 < 0) {
  cat("   ↓ Hazard ratio DECREASES over time\n")
  cat("   → Risk difference between clusters shrinks with longer follow-up\n\n")
} else {
  cat("   → Hazard ratio is relatively constant\n\n")
}

cat("3. STATISTICAL SIGNIFICANCE:\n")
cat("   Time interaction p-value:", format(pval_interaction, scientific = TRUE, digits = 3), "\n")
if (pval_interaction < 0.05) {
  cat("   ✓ Significant time-varying effect confirms PH violation\n\n")
} else {
  cat("   → Time-varying effect not statistically significant\n\n")
}

cat("4. PIECEWISE ANALYSIS:\n")
cat("   Early HR (0-12m):", round(hr_early, 3), "\n")
cat("   Late HR (>12m):", round(hr_late, 3), "\n")
cat("   Period interaction p-value:", format(pval_period_interaction, scientific = TRUE, digits = 3), "\n\n")

cat("INTERPRETATION:\n")
cat("The time-varying coefficient model reveals how the hazard ratio changes\n")
cat("throughout the follow-up period. This analysis is appropriate when the\n")
cat("proportional hazards assumption is violated and provides a more nuanced\n")
cat("understanding of the temporal dynamics of risk.\n\n")

cat("NEXT STEPS:\n")
cat("  → Part 1.3: Landmark analysis\n")
cat("  → Part 1.4: RMST analysis\n\n")

cat("### PART 1.2 COMPLETE ###\n")
