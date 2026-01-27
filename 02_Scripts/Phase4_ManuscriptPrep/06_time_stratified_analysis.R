#!/usr/bin/env Rscript
# Phase 4 Part 6: Time-Stratified Analysis
# Purpose: Explore why HR decreases over time (2.22 → 1.62)

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("TIME-STRATIFIED ANALYSIS\n")
cat("==============================================================================\n\n")

cat("Phase 3 showed HR decreases over time:\n")
cat("  - 6 months: HR=2.22\n")
cat("  - 12 months: HR=2.02\n")
cat("  - 24 months: HR=1.84\n")
cat("  - 60 months: HR=1.62\n\n")

cat("Question: Is cluster effect EARLY-specific or LATE-specific?\n\n")

# Load data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

analysis_data <- survival_data %>%
  filter(!is.na(OS_months) & OS_months > 0 & !is.na(OS_event) & !is.na(cluster))

cat(sprintf("Analysis dataset: %d samples, %d events\n\n", nrow(analysis_data), sum(analysis_data$OS_event)))

# ==============================================================================
# EARLY VS LATE DEATH ANALYSIS
# ==============================================================================

cat("=== EARLY VS LATE DEATH PATTERNS ===\n\n")

# Define cutoff (24 months = 2 years)
cutoff <- 24

analysis_data <- analysis_data %>%
  mutate(
    death_timing = case_when(
      OS_event == 0 ~ "Censored",
      OS_months <= cutoff ~ "Early death (≤24m)",
      OS_months > cutoff ~ "Late death (>24m)"
    ),
    early_death = ifelse(OS_event == 1 & OS_months <= cutoff, 1, 0),
    late_death = ifelse(OS_event == 1 & OS_months > cutoff, 1, 0)
  )

# Counts by cluster
counts <- analysis_data %>%
  group_by(cluster, death_timing) %>%
  summarise(N = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = death_timing, values_from = N, values_fill = 0)

cat("Death patterns by cluster:\n")
print(counts)
cat("\n")

# Test cluster association with early vs late death
early_by_cluster <- table(analysis_data$cluster, analysis_data$early_death)
late_by_cluster <- table(analysis_data$cluster, analysis_data$late_death)

fisher_early <- fisher.test(early_by_cluster)
fisher_late <- fisher.test(late_by_cluster)

timing_summary <- data.frame(
  Death_Type = c("Early (≤24m)", "Late (>24m)"),
  Cluster1_Deaths = c(early_by_cluster[1, 2], late_by_cluster[1, 2]),
  Cluster2_Deaths = c(early_by_cluster[2, 2], late_by_cluster[2, 2]),
  Cluster1_Pct = c(
    early_by_cluster[1, 2] / sum(early_by_cluster[1, ]) * 100,
    late_by_cluster[1, 2] / sum(late_by_cluster[1, ]) * 100
  ),
  Cluster2_Pct = c(
    early_by_cluster[2, 2] / sum(early_by_cluster[2, ]) * 100,
    late_by_cluster[2, 2] / sum(late_by_cluster[2, ]) * 100
  ),
  OR = c(fisher_early$estimate, fisher_late$estimate),
  P_value = c(fisher_early$p.value, fisher_late$p.value)
)

cat("Early vs Late Death Analysis:\n")
print(timing_summary, row.names = FALSE)

write.csv(timing_summary,
          "03_Results/21_Manuscript_Prep/early_vs_late_death_analysis.csv",
          row.names = FALSE)

cat("\n")

# ==============================================================================
# CUMULATIVE INCIDENCE CURVES
# ==============================================================================

cat("=== CUMULATIVE INCIDENCE BY TIME PERIOD ===\n\n")

# Calculate cumulative incidence at key timepoints
timepoints <- c(6, 12, 18, 24, 36, 48, 60)

km_fit <- survfit(Surv(OS_months, OS_event) ~ cluster, data = analysis_data)

# Get summary at specific timepoints
km_summary <- summary(km_fit, times = timepoints)

# Extract cumulative incidence by cluster (1 - survival)
# The strata are ordered: first all cluster=1 timepoints, then all cluster=2 timepoints
n_times <- length(timepoints)
n_strata <- 2

cum_inc <- data.frame(
  Time = timepoints,
  C1_CumInc = 1 - km_summary$surv[km_summary$strata == "cluster=1"],
  C2_CumInc = 1 - km_summary$surv[km_summary$strata == "cluster=2"]
)

cum_inc$Difference <- cum_inc$C2_CumInc - cum_inc$C1_CumInc
cum_inc$Ratio <- cum_inc$C2_CumInc / cum_inc$C1_CumInc

cat("Cumulative Incidence (proportion with events):\n")
print(cum_inc, row.names = FALSE)

write.csv(cum_inc,
          "03_Results/21_Manuscript_Prep/cumulative_incidence_by_time.csv",
          row.names = FALSE)

cat("\n")

# ==============================================================================
# CONDITIONAL SURVIVAL ANALYSIS
# ==============================================================================

cat("=== CONDITIONAL SURVIVAL (FOR THOSE WHO SURVIVE TO LANDMARK) ===\n\n")

landmarks <- c(6, 12, 24)

conditional_results <- lapply(landmarks, function(landmark) {

  # Subset to those alive at landmark
  cond_data <- analysis_data %>% filter(OS_months >= landmark)

  n_c1 <- sum(cond_data$cluster == 1)
  n_c2 <- sum(cond_data$cluster == 2)

  if (n_c1 < 10 || n_c2 < 10) {
    cat(sprintf("Landmark %dm: Insufficient data\n", landmark))
    return(NULL)
  }

  # Cox from landmark
  cox_fit <- coxph(Surv(OS_months - landmark, OS_event) ~ cluster, data = cond_data)
  hr <- exp(coef(cox_fit))
  hr_ci <- exp(confint(cox_fit))
  p_val <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]

  cat(sprintf("Landmark %dm:\n", landmark))
  cat(sprintf("  N at risk: C1=%d, C2=%d\n", n_c1, n_c2))
  cat(sprintf("  HR: %.2f (%.2f-%.2f), p=%.4f\n", hr, hr_ci[1], hr_ci[2], p_val))
  cat("\n")

  data.frame(
    Landmark_Months = landmark,
    N_Cluster1 = n_c1,
    N_Cluster2 = n_c2,
    HR = hr,
    HR_Lower = hr_ci[1],
    HR_Upper = hr_ci[2],
    P_Value = p_val
  )

}) %>% bind_rows()

write.csv(conditional_results,
          "03_Results/21_Manuscript_Prep/conditional_survival_by_landmark.csv",
          row.names = FALSE)

# ==============================================================================
# VISUALIZATION
# ==============================================================================

cat("Creating visualization...\n")

pdf("04_Figures/20_Manuscript_Prep/time_stratified_analysis.pdf", width = 12, height = 10)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 1))

# Panel A: Early vs Late death rates
par(mar = c(5, 5, 4, 2))
barplot_data <- as.matrix(timing_summary[, c("Cluster1_Pct", "Cluster2_Pct")])
rownames(barplot_data) <- timing_summary$Death_Type

barplot(barplot_data, beside = TRUE,
        col = c("lightblue", "lightcoral"),
        ylab = "Mortality Rate (%)",
        main = "A. Early vs Late Mortality by Cluster",
        legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n"),
        ylim = c(0, 70),
        cex.lab = 1.2, cex.main = 1.3)

# Add p-values
text(c(2, 5), c(65, 65),
     labels = sprintf("p=%.3f", timing_summary$P_value),
     cex = 1.1)

# Panel B: Cumulative incidence over time
par(mar = c(5, 5, 4, 2))
plot(cum_inc$Time, cum_inc$C1_CumInc * 100,
     type = "l", lwd = 3, col = "#E41A1C",
     ylim = c(0, 100),
     xlab = "Months",
     ylab = "Cumulative Incidence (%)",
     main = "B. Cumulative Incidence Over Time",
     cex.lab = 1.2, cex.main = 1.3)

lines(cum_inc$Time, cum_inc$C2_CumInc * 100, lwd = 3, col = "#377EB8")

legend("topleft",
       legend = c("Cluster 1", "Cluster 2"),
       col = c("#E41A1C", "#377EB8"),
       lwd = 3, bty = "n", cex = 1.2)

# Add shaded regions for early/late
abline(v = cutoff, lty = 2, col = "gray50", lwd = 2)
text(12, 90, "Early\n(0-24m)", cex = 1.1, col = "gray40")
text(42, 90, "Late\n(>24m)", cex = 1.1, col = "gray40")

# Panel C: HR over time
par(mar = c(5, 5, 4, 2))

# HR decreases over time (from Phase 3)
hr_time <- data.frame(
  Time = c(6, 12, 24, 36, 60),
  HR = c(2.22, 2.02, 1.84, 1.74, 1.62)
)

plot(hr_time$Time, hr_time$HR,
     type = "b", lwd = 3, col = "darkred", pch = 19, cex = 1.5,
     ylim = c(1, 2.5),
     xlab = "Months from Diagnosis",
     ylab = "Hazard Ratio (C2 vs C1)",
     main = "C. Time-Varying Hazard Ratio",
     cex.lab = 1.2, cex.main = 1.3)

abline(h = 1, lty = 2, col = "gray50", lwd = 2)
abline(v = cutoff, lty = 2, col = "gray50", lwd = 1)

# Add annotation
text(40, 2.3, "HR decreases\nover time", cex = 1.1, col = "darkred", font = 2)

dev.off()

cat("  ✓ time_stratified_analysis.pdf\n\n")

# ==============================================================================
# INTERPRETATION
# ==============================================================================

cat("==============================================================================\n")
cat("INTERPRETATION\n")
cat("==============================================================================\n\n")

# Determine pattern
early_sig <- timing_summary$P_value[1] < 0.05
late_sig <- timing_summary$P_value[2] < 0.05

early_diff <- timing_summary$Cluster2_Pct[1] - timing_summary$Cluster1_Pct[1]
late_diff <- timing_summary$Cluster2_Pct[2] - timing_summary$Cluster1_Pct[2]

if (early_sig && !late_sig) {
  cat("✓ EARLY-SPECIFIC EFFECT\n\n")
  cat("Cluster effect is strongest in the first 24 months:\n")
  cat(sprintf("  - Early death rate difference: %.1f percentage points (p=%.4f)\n",
              early_diff, timing_summary$P_value[1]))
  cat(sprintf("  - Late death rate difference: %.1f percentage points (p=%.4f)\n\n",
              late_diff, timing_summary$P_value[2]))
  cat("LIKELY MECHANISMS:\n")
  cat("  1. Treatment response differences (induction/consolidation)\n")
  cat("  2. Cluster 2 enriched for aggressive mutations (TP53)\n")
  cat("  3. Early refractory disease in high-risk cluster\n\n")

} else if (!early_sig && late_sig) {
  cat("✓ LATE-SPECIFIC EFFECT\n\n")
  cat("Cluster effect emerges after 24 months:\n")
  cat(sprintf("  - Early death rate similar between clusters (p=%.4f)\n",
              timing_summary$P_value[1]))
  cat(sprintf("  - Late death rate difference: %.1f percentage points (p=%.4f)\n\n",
              late_diff, timing_summary$P_value[2]))
  cat("LIKELY MECHANISMS:\n")
  cat("  1. Relapse patterns differ by cluster\n")
  cat("  2. Long-term remission duration differences\n")
  cat("  3. Different clonal evolution trajectories\n\n")

} else if (early_sig && late_sig) {
  cat("✓ SUSTAINED EFFECT\n\n")
  cat("Cluster effect present throughout follow-up:\n")
  cat(sprintf("  - Early deaths higher in Cluster 2 (p=%.4f)\n",
              timing_summary$P_value[1]))
  cat(sprintf("  - Late deaths also higher in Cluster 2 (p=%.4f)\n\n",
              timing_summary$P_value[2]))
  cat("Interpretation: Persistent biological difference\n\n")

} else {
  cat("✗ NO CLEAR TEMPORAL PATTERN\n\n")
  cat("Neither early nor late deaths significantly different\n")
  cat("HR decrease may reflect survivor bias or small effects\n\n")
}

cat("WHY DOES HR DECREASE OVER TIME?\n\n")
cat("Possible explanations:\n")
cat("1. Survivor selection bias:\n")
cat("   - High-risk Cluster 2 patients die early\n")
cat("   - Survivors in both clusters are more similar\n")
cat("   - HR attenuates among long-term survivors\n\n")

cat("2. Treatment phase effects:\n")
cat("   - Induction/consolidation: Large effect (0-24m)\n")
cat("   - Maintenance/surveillance: Smaller effect (>24m)\n\n")

cat("3. Competing risks:\n")
cat("   - Early deaths dominated by disease aggression\n")
cat("   - Late deaths include treatment complications, relapse\n\n")

cat("✅ Time-stratified analysis complete\n\n")
cat("Files generated:\n")
cat("  - 03_Results/21_Manuscript_Prep/early_vs_late_death_analysis.csv\n")
cat("  - 03_Results/21_Manuscript_Prep/cumulative_incidence_by_time.csv\n")
cat("  - 03_Results/21_Manuscript_Prep/conditional_survival_by_landmark.csv\n")
cat("  - 04_Figures/20_Manuscript_Prep/time_stratified_analysis.pdf\n\n")
