#!/usr/bin/env Rscript
# Phase 4 Part 3: TCGA Power Analysis
# Purpose: Calculate minimum detectable effect size and power curves

suppressPackageStartupMessages({
  library(dplyr)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("TCGA POWER ANALYSIS\n")
cat("==============================================================================\n\n")

# Create output directories
dir.create("03_Results/21_Manuscript_Prep", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/20_Manuscript_Prep", recursive = TRUE, showWarnings = FALSE)

# TCGA parameters
tcga_n <- 151
tcga_events <- 97
tcga_c1 <- 76
tcga_c2 <- 75
tcga_event_rate <- tcga_events / tcga_n

cat("TCGA-LAML Cohort:\n")
cat(sprintf("  Total n: %d\n", tcga_n))
cat(sprintf("  Events: %d (%.1f%%)\n", tcga_events, tcga_event_rate * 100))
cat(sprintf("  Cluster 1: %d\n", tcga_c1))
cat(sprintf("  Cluster 2: %d\n\n", tcga_c2))

# BeatAML observed effect
beatAML_observed_hr <- 1.39
tcga_observed_hr <- 1.24

# Function to calculate power for survival studies
# Schoenfeld formula: N_events = (Z_alpha + Z_beta)^2 * 4 / (log(HR))^2
calculate_power <- function(hr, n_events, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha/2)
  required_events <- (z_alpha + qnorm(0.80))^2 * 4 / (log(hr))^2

  # Observed vs required ratio gives approximate power
  # Using non-central chi-square approximation
  ncp <- (log(hr)^2 * n_events) / 4
  power <- 1 - pchisq(qchisq(1 - alpha, df = 1), df = 1, ncp = ncp)

  return(list(
    power = power,
    required_events = required_events
  ))
}

# Calculate power for range of HRs
hr_range <- seq(1.05, 3.0, by = 0.05)
power_results <- data.frame(
  HR = hr_range,
  Power = sapply(hr_range, function(hr) {
    result <- calculate_power(hr, tcga_events)
    return(result$power)
  }),
  Required_events_80pct = sapply(hr_range, function(hr) {
    (qnorm(0.975) + qnorm(0.80))^2 * 4 / (log(hr))^2
  })
)

# Find HR for 80% power
power_80_idx <- which(power_results$Power >= 0.80)[1]
if (is.na(power_80_idx)) {
  power_80_hr <- NA
  cat("⚠ TCGA cannot detect any HR up to 3.0 with 80% power\n\n")
} else {
  power_80_hr <- power_results$HR[power_80_idx]
  cat(sprintf("Minimum detectable HR (80%% power): %.2f\n\n", power_80_hr))
}

# Power for observed HRs
power_beatAML_hr <- calculate_power(beatAML_observed_hr, tcga_events)$power
power_tcga_hr <- calculate_power(tcga_observed_hr, tcga_events)$power

cat("Power to detect observed effects:\n")
cat(sprintf("  BeatAML HR=%.2f: %.1f%% power\n", beatAML_observed_hr, power_beatAML_hr * 100))
cat(sprintf("  TCGA HR=%.2f: %.1f%% power\n\n", tcga_observed_hr, power_tcga_hr * 100))

# Required sample size for 80% power at BeatAML effect
required_events_beatAML <- (qnorm(0.975) + qnorm(0.80))^2 * 4 / (log(beatAML_observed_hr))^2
required_n_beatAML <- required_events_beatAML / tcga_event_rate

cat("Required sample size to detect BeatAML effect (HR=1.39):\n")
cat(sprintf("  Required events: %d\n", ceiling(required_events_beatAML)))
cat(sprintf("  Required total n: %d (at %.1f%% event rate)\n",
            ceiling(required_n_beatAML), tcga_event_rate * 100))
cat(sprintf("  TCGA has: %d events (%.1f%% of required)\n\n",
            tcga_events, tcga_events / required_events_beatAML * 100))

# Create power summary table
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
    "Required events for 80% power (HR=1.39)",
    "TCGA has (% of required)"
  ),
  Value = c(
    tcga_n,
    tcga_events,
    round(tcga_event_rate * 100, 1),
    round(beatAML_observed_hr, 2),
    round(tcga_observed_hr, 2),
    ifelse(is.na(power_80_hr), ">3.0", round(power_80_hr, 2)),
    round(power_beatAML_hr * 100, 1),
    round(power_tcga_hr * 100, 1),
    ceiling(required_events_beatAML),
    round(tcga_events / required_events_beatAML * 100, 1)
  )
)

write.csv(power_summary,
          "03_Results/21_Manuscript_Prep/tcga_power_analysis_detailed.csv",
          row.names = FALSE)

cat("=== POWER SUMMARY ===\n\n")
print(power_summary, row.names = FALSE)

# Create power curve plot
pdf("04_Figures/20_Manuscript_Prep/tcga_power_curve.pdf", width = 10, height = 7)
par(mar = c(5, 5, 4, 2))

plot(power_results$HR, power_results$Power,
     type = "l", lwd = 3, col = "navy",
     xlab = "Hazard Ratio (Cluster 2 vs Cluster 1)",
     ylab = "Statistical Power",
     main = "TCGA-LAML Power Analysis (n=151, 97 events)",
     ylim = c(0, 1),
     xlim = c(1, 3),
     cex.lab = 1.3,
     cex.axis = 1.2,
     cex.main = 1.4)

# Add reference lines
abline(h = 0.80, lty = 2, col = "red", lwd = 2)
abline(h = 0.05, lty = 3, col = "gray50", lwd = 1)
abline(v = beatAML_observed_hr, lty = 2, col = "blue", lwd = 2)
if (!is.na(power_80_hr) && power_80_hr <= 3.0) {
  abline(v = power_80_hr, lty = 2, col = "darkgreen", lwd = 2)
}

# Add points for observed HRs
points(beatAML_observed_hr, power_beatAML_hr, pch = 19, col = "blue", cex = 2)
points(tcga_observed_hr, power_tcga_hr, pch = 19, col = "purple", cex = 2)
if (!is.na(power_80_hr) && power_80_hr <= 3.0) {
  points(power_80_hr, 0.80, pch = 19, col = "darkgreen", cex = 2)
}

# Add labels
text(beatAML_observed_hr, power_beatAML_hr + 0.12,
     sprintf("BeatAML HR=%.2f\nPower=%.1f%%", beatAML_observed_hr, power_beatAML_hr * 100),
     col = "blue", cex = 1.1)

text(tcga_observed_hr, power_tcga_hr - 0.12,
     sprintf("TCGA HR=%.2f\nPower=%.1f%%", tcga_observed_hr, power_tcga_hr * 100),
     col = "purple", cex = 1.1)

if (!is.na(power_80_hr) && power_80_hr <= 3.0) {
  text(power_80_hr + 0.3, 0.85,
       sprintf("Min detectable\nHR=%.2f", power_80_hr),
       col = "darkgreen", cex = 1.1)
}

# Add legend
legend("bottomright",
       legend = c("Power curve", "80% power", "5% power (type I error)",
                  "BeatAML observed", "TCGA observed",
                  if (!is.na(power_80_hr) && power_80_hr <= 3.0) "Min detectable" else NULL),
       lty = c(1, 2, 3, 2, 2, if (!is.na(power_80_hr) && power_80_hr <= 3.0) 2 else NULL),
       col = c("navy", "red", "gray50", "blue", "purple",
               if (!is.na(power_80_hr) && power_80_hr <= 3.0) "darkgreen" else NULL),
       lwd = c(3, 2, 1, 2, 2, if (!is.na(power_80_hr) && power_80_hr <= 3.0) 2 else NULL),
       cex = 1.1,
       bty = "n")

# Add text box with interpretation
text(2.5, 0.3,
     sprintf("TCGA is severely underpowered\nto detect moderate effects"),
     cex = 1.2, col = "red", font = 2)

dev.off()

cat("\n✅ Power analysis complete\n")
cat("\nFiles generated:\n")
cat("  - 03_Results/21_Manuscript_Prep/tcga_power_analysis_detailed.csv\n")
cat("  - 04_Figures/20_Manuscript_Prep/tcga_power_curve.pdf\n\n")

cat("==============================================================================\n")
cat("INTERPRETATION\n")
cat("==============================================================================\n\n")

if (power_beatAML_hr < 0.50) {
  cat("TCGA is SEVERELY UNDERPOWERED to detect the BeatAML effect\n")
  cat(sprintf("  - Only %.1f%% power to detect HR=%.2f\n", power_beatAML_hr * 100, beatAML_observed_hr))
  cat(sprintf("  - Would need %d events (%.0fx more) for 80%% power\n",
              ceiling(required_events_beatAML),
              required_events_beatAML / tcga_events))
  cat("\n✓ TCGA 'failure' is expected due to insufficient sample size\n")
  cat("✓ Effect size is consistent (HR=1.24 vs 1.39), just underpowered\n")
} else {
  cat("TCGA has moderate power to detect the BeatAML effect\n")
}

cat("\n")
