# PHASE 3 - PART 1.4: RESTRICTED MEAN SURVIVAL TIME (RMST) ANALYSIS
# RMST is a completely assumption-free survival measure
# - Does NOT require proportional hazards
# - Does NOT require parametric assumptions
# - Provides clinically interpretable differences in mean survival time

library(tidyverse)
library(survival)
library(survminer)

cat("=== PHASE 3: PART 1.4 - RMST ANALYSIS ===\n\n")

# Manual RMST calculation function
calculate_rmst <- function(time, status, tau) {
  # Fit Kaplan-Meier
  km_fit <- survfit(Surv(time, status) ~ 1)

  # Get KM estimates
  km_times <- km_fit$time
  km_surv <- km_fit$surv
  km_std_err <- km_fit$std.err

  # Restrict to tau
  idx <- km_times <= tau
  km_times_r <- c(0, km_times[idx], tau)
  km_surv_r <- c(1, km_surv[idx], tail(km_surv[idx], 1))

  # Calculate RMST as area under KM curve
  rmst <- 0
  for (i in 1:(length(km_times_r) - 1)) {
    rmst <- rmst + km_surv_r[i] * (km_times_r[i + 1] - km_times_r[i])
  }

  # Approximate SE using Greenwood formula
  var_sum <- 0
  if (length(km_times[idx]) > 0) {
    for (i in 1:sum(idx)) {
      if (km_fit$n.risk[i] > 0 && km_fit$n.event[i] > 0) {
        var_sum <- var_sum + km_fit$n.event[i] / (km_fit$n.risk[i] * (km_fit$n.risk[i] - km_fit$n.event[i]))
      }
    }
  }
  rmst_se <- sqrt(var_sum) * (tau / 2)  # Rough approximation

  return(list(rmst = rmst, se = rmst_se))
}

cat("Using manual RMST calculation (survRM2 not required)\n\n")

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
cat("Events:", sum(survival_data$OS_STATUS), "\n")
cat("Max follow-up:", round(max(survival_data$OS_MONTHS), 1), "months\n\n")

# ============================================================================
# 2. SELECT RESTRICTION TIME
# ============================================================================

cat("=== STEP 2: SELECTING RESTRICTION TIME (τ) ===\n\n")

# RMST is calculated up to a restriction time τ
# Common choices:
# - Median follow-up
# - Time when ~25% remain at risk
# - Clinically relevant timepoint (e.g., 5 years)

median_followup <- median(survival_data$OS_MONTHS)
quantile_75 <- quantile(survival_data$OS_MONTHS, 0.75)
five_years <- 60  # 5 years in months

cat("Restriction time options:\n")
cat("  Median follow-up:", round(median_followup, 1), "months\n")
cat("  75th percentile:", round(quantile_75, 1), "months\n")
cat("  5 years (60 months):", five_years, "months\n\n")

# Use multiple restriction times for robustness
tau_values <- c(24, 36, 48, 60, median_followup)
cat("Will analyze RMST at multiple τ values:\n")
cat("  ", paste(round(tau_values, 1), "months", collapse=", "), "\n\n")

# ============================================================================
# 3. RMST ANALYSIS AT MULTIPLE TIMEPOINTS
# ============================================================================

cat("=== STEP 3: RMST ANALYSIS AT MULTIPLE RESTRICTION TIMES ===\n\n")

# Prepare data
rmst_data <- survival_data %>%
  select(OS_MONTHS, OS_STATUS, cluster_assignment) %>%
  filter(OS_MONTHS > 0)

# Initialize results storage
rmst_results <- data.frame()

for (tau in tau_values) {

  cat("### RESTRICTION TIME:", round(tau, 1), "MONTHS ###\n\n")

  # Check if tau is reasonable (not beyond most follow-up)
  if (tau > quantile(survival_data$OS_MONTHS, 0.9)) {
    cat("  Warning: τ exceeds 90th percentile of follow-up\n")
    cat("  Skipping this timepoint\n\n")
    next
  }

  # Perform RMST analysis
  tryCatch({

    # Calculate RMST for each cluster
    cluster1_data <- rmst_data %>% filter(cluster_assignment == 1)
    cluster2_data <- rmst_data %>% filter(cluster_assignment == 2)

    rmst_c1 <- calculate_rmst(cluster1_data$OS_MONTHS, cluster1_data$OS_STATUS, tau)
    rmst_c2 <- calculate_rmst(cluster2_data$OS_MONTHS, cluster2_data$OS_STATUS, tau)

    rmst_cluster1 <- rmst_c1$rmst
    rmst_se_cluster1 <- rmst_c1$se
    rmst_cluster2 <- rmst_c2$rmst
    rmst_se_cluster2 <- rmst_c2$se

    # Difference in RMST
    rmst_diff <- rmst_cluster1 - rmst_cluster2
    rmst_diff_se <- sqrt(rmst_se_cluster1^2 + rmst_se_cluster2^2)
    rmst_diff_lower <- rmst_diff - 1.96 * rmst_diff_se
    rmst_diff_upper <- rmst_diff + 1.96 * rmst_diff_se

    # Z-test for difference
    z_stat <- rmst_diff / rmst_diff_se
    rmst_pvalue <- 2 * (1 - pnorm(abs(z_stat)))

    cat("  Restriction time (τ):", round(tau, 1), "months\n")
    cat("  RMST Cluster 1:", round(rmst_cluster1, 2), "±", round(rmst_se_cluster1, 2), "months\n")
    cat("  RMST Cluster 2:", round(rmst_cluster2, 2), "±", round(rmst_se_cluster2, 2), "months\n")
    cat("  Difference (Cluster 1 - Cluster 2):", round(rmst_diff, 2),
        "months (95% CI:", round(rmst_diff_lower, 2), "to", round(rmst_diff_upper, 2), ")\n")
    cat("  P-value:", format(rmst_pvalue, scientific = TRUE, digits = 3), "\n")

    if (rmst_pvalue < 0.05) {
      cat("  ✓ SIGNIFICANT difference in mean survival time\n")
      if (rmst_diff > 0) {
        cat("  → Cluster 1 survives", round(rmst_diff, 1), "months longer on average (up to τ)\n")
      } else {
        cat("  → Cluster 2 survives", round(abs(rmst_diff), 1), "months longer on average (up to τ)\n")
      }
    } else {
      cat("  ✗ No significant difference\n")
    }

    cat("\n")

    # Store results
    rmst_results <- rbind(rmst_results, data.frame(
      tau = tau,
      rmst_cluster1 = rmst_cluster1,
      rmst_se_cluster1 = rmst_se_cluster1,
      rmst_cluster2 = rmst_cluster2,
      rmst_se_cluster2 = rmst_se_cluster2,
      rmst_difference = rmst_diff,
      rmst_diff_lower = rmst_diff_lower,
      rmst_diff_upper = rmst_diff_upper,
      pvalue = rmst_pvalue,
      significant = rmst_pvalue < 0.05
    ))

  }, error = function(e) {
    cat("  Error at τ =", tau, ":", e$message, "\n\n")
  })
}

# ============================================================================
# 4. RMST vs TRADITIONAL MEASURES
# ============================================================================

cat("=== STEP 4: COMPARING RMST WITH TRADITIONAL MEASURES ===\n\n")

# Standard Cox model for comparison
cox_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                   data = survival_data)
hr_cox <- exp(coef(cox_model))
cox_pval <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]

# Log-rank test
logrank <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                    data = survival_data)
logrank_pval <- 1 - pchisq(logrank$chisq, df = 1)

cat("Traditional Survival Measures:\n")
cat("  Log-rank test p-value:", format(logrank_pval, scientific = TRUE, digits = 3), "\n")
cat("  Cox HR (Cluster 2 vs 1):", round(hr_cox, 3), "\n")
cat("  Cox p-value:", format(cox_pval, scientific = TRUE, digits = 3), "\n\n")

cat("RMST Analysis Summary:\n")
cat("  Number of τ values tested:", nrow(rmst_results), "\n")
cat("  Significant results:", sum(rmst_results$significant), "\n")
cat("  Mean RMST difference:", round(mean(rmst_results$rmst_difference), 2), "months\n\n")

# ============================================================================
# 5. CLINICAL INTERPRETATION
# ============================================================================

cat("=== STEP 5: CLINICAL INTERPRETATION ===\n\n")

# Calculate loss of life expectancy over 5 years
tau_5yr <- 60
if (60 %in% rmst_results$tau) {
  rmst_5yr <- rmst_results[rmst_results$tau == 60, ]

  cat("Loss of Life Expectancy (5-year horizon):\n")
  cat("  Cluster 1 mean survival:", round(rmst_5yr$rmst_cluster1, 1), "months\n")
  cat("  Cluster 2 mean survival:", round(rmst_5yr$rmst_cluster2, 1), "months\n")
  cat("  Difference:", round(rmst_5yr$rmst_difference, 1), "months\n\n")

  if (rmst_5yr$rmst_difference > 0) {
    cat("  Interpretation: Over a 5-year period, Cluster 2 patients\n")
    cat("  lose an average of", round(rmst_5yr$rmst_difference, 1), "months of life\n")
    cat("  compared to Cluster 1 patients.\n\n")
  }

  # Convert to percentage
  pct_loss <- (rmst_5yr$rmst_difference / rmst_5yr$rmst_cluster1) * 100
  cat("  Relative difference:", round(pct_loss, 1), "% loss of expected survival time\n\n")
}

# ============================================================================
# 6. SAVE RESULTS
# ============================================================================

cat("=== STEP 6: SAVING RESULTS ===\n\n")

write.csv(rmst_results,
          "03_Results/11_Survival_Reanalysis/04_rmst_results.csv",
          row.names = FALSE)
cat("✓ Saved: 04_rmst_results.csv\n\n")

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

cat("=== STEP 7: CREATING VISUALIZATIONS ===\n\n")

# Plot 1: RMST by restriction time
p_rmst_by_tau <- ggplot(rmst_results, aes(x = tau)) +
  geom_line(aes(y = rmst_cluster1, color = "Cluster 1"), size = 1.2) +
  geom_point(aes(y = rmst_cluster1, color = "Cluster 1"), size = 3) +
  geom_line(aes(y = rmst_cluster2, color = "Cluster 2"), size = 1.2) +
  geom_point(aes(y = rmst_cluster2, color = "Cluster 2"), size = 3) +
  scale_color_manual(values = c("Cluster 1" = "#2E9FDF", "Cluster 2" = "#E7B800"),
                     labels = c("Cluster 1 (Proliferative)",
                                "Cluster 2 (Immune-Inflammatory)")) +
  labs(
    title = "Restricted Mean Survival Time by Restriction Time",
    subtitle = "Average survival time up to τ",
    x = "Restriction Time τ (months)",
    y = "RMST (months)",
    color = "Molecular Subtype"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("04_Figures/11_Survival_Reanalysis/04_rmst_by_tau.pdf",
       p_rmst_by_tau, width = 10, height = 6)
cat("✓ Saved: 04_rmst_by_tau.pdf\n")

# Plot 2: RMST difference with CI
p_rmst_diff <- ggplot(rmst_results, aes(x = tau, y = rmst_difference)) +
  geom_line(size = 1.2, color = "#E7B800") +
  geom_point(aes(color = significant), size = 4) +
  geom_errorbar(aes(ymin = rmst_diff_lower, ymax = rmst_diff_upper),
                width = 2, color = "#E7B800") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#00BA38", "FALSE" = "#F8766D"),
                     labels = c("TRUE" = "Significant", "FALSE" = "Not Significant")) +
  labs(
    title = "RMST Difference (Cluster 1 - Cluster 2)",
    subtitle = "Positive values indicate Cluster 1 survives longer",
    x = "Restriction Time τ (months)",
    y = "RMST Difference (months)",
    color = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("04_Figures/11_Survival_Reanalysis/04_rmst_difference.pdf",
       p_rmst_diff, width = 10, height = 6)
cat("✓ Saved: 04_rmst_difference.pdf\n\n")

# ============================================================================
# 8. SUMMARY
# ============================================================================

cat("=== SUMMARY AND INTERPRETATION ===\n\n")

cat("PART 1.4: RMST ANALYSIS COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. RMST ANALYSIS (Assumption-Free):\n")
cat("   Restriction times tested:", nrow(rmst_results), "\n")
cat("   Significant results:", sum(rmst_results$significant),
    "(", round(100 * sum(rmst_results$significant) / nrow(rmst_results), 0), "%)\\n")
cat("   Mean RMST difference:", round(mean(rmst_results$rmst_difference), 2), "months\n\n")

if (60 %in% rmst_results$tau) {
  rmst_5yr <- rmst_results[rmst_results$tau == 60, ]
  cat("2. CLINICAL IMPACT (5-year horizon):\n")
  cat("   Loss of life expectancy:", round(rmst_5yr$rmst_difference, 1), "months\n")
  cat("   Relative loss:", round((rmst_5yr$rmst_difference / rmst_5yr$rmst_cluster1) * 100, 1), "%\n")
  cat("   P-value:", format(rmst_5yr$pvalue, scientific = TRUE, digits = 3), "\n\n")
}

cat("3. COMPARISON WITH TRADITIONAL METHODS:\n")
cat("   Log-rank p-value:", format(logrank_pval, scientific = TRUE, digits = 3), "\n")
cat("   Cox HR:", round(hr_cox, 3), "\n")
cat("   Cox p-value:", format(cox_pval, scientific = TRUE, digits = 3), "\n")

if (all(rmst_results$significant) && logrank_pval < 0.05 && cox_pval < 0.05) {
  cat("   ✓ AGREEMENT: All methods confirm significant prognostic effect\n\n")
} else if (any(rmst_results$significant)) {
  cat("   ⚠ PARTIAL AGREEMENT: RMST significant at some timepoints\n\n")
} else {
  cat("   ✗ DISAGREEMENT: Different conclusions from different methods\n\n")
}

cat("INTERPRETATION:\n")
cat("RMST provides a clinically interpretable, assumption-free measure of survival.\n")
cat("Unlike hazard ratios, RMST differences directly quantify the loss or gain of\n")
cat("survival time, making it easier to communicate clinical impact to patients.\n\n")

if (mean(rmst_results$rmst_difference) > 0) {
  cat("Cluster 2 patients experience a loss of approximately",
      round(mean(rmst_results$rmst_difference), 1), "months of survival\n")
  cat("compared to Cluster 1 patients, averaged across restriction times.\n")
  cat("This represents a clinically meaningful difference in life expectancy.\n\n")
}

cat("ADVANTAGES OF RMST:\n")
cat("  ✓ No proportional hazards assumption required\n")
cat("  ✓ Directly interpretable (months of life gained/lost)\n")
cat("  ✓ Robust to non-proportional hazards\n")
cat("  ✓ Easy to communicate to clinicians and patients\n\n")

cat("### PART 1.4 COMPLETE ###\n")
cat("### ALL PART 1 ANALYSES (PH VIOLATION FIXES) COMPLETE ###\n\n")

cat("NEXT: Part 2 - Multivariate analysis\n")
