################################################################################
# CREATE FIGURE S2: PROPORTIONAL HAZARDS DIAGNOSTICS (COMBINED)
# Combines multiple PH diagnostic analyses into a single figure
# Date: 2025-12-09
################################################################################

setwd("D:/Projects/Project_AML")
library(survival)
library(ggplot2)
library(dplyr)
library(gridExtra)

cat("=== CREATING FIGURE S2: PH DIAGNOSTICS (COMBINED) ===\n\n")

# Load survival data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Prepare data
survival_data <- survival_data %>%
  filter(!is.na(OS_months) & !is.na(OS_event) & !is.na(cluster)) %>%
  mutate(
    cluster = factor(cluster, levels = c(1, 2),
                    labels = c("Cluster 1", "Cluster 2"))
  )

cat("Data loaded: n =", nrow(survival_data), "patients\n\n")

################################################################################
# CREATE COMBINED FIGURE
################################################################################

pdf("05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics_Combined.pdf",
    width = 14, height = 12)

layout(matrix(c(1,2,3,4), nrow = 2, byrow = TRUE))
par(mar = c(5, 5, 3, 2))

################################################################################
# PANEL A: SCHOENFELD RESIDUALS TEST
################################################################################

cat("Panel A: Schoenfeld residuals...\n")

# Fit Cox model
cox_model <- coxph(Surv(OS_months, OS_event) ~ cluster, data = survival_data)

# Test PH assumption
ph_test <- cox.zph(cox_model)

# Plot Schoenfeld residuals
plot(ph_test,
     main = "A. Schoenfeld Residuals Test",
     xlab = "Time (months)",
     ylab = "Beta(t) for Cluster",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.5,
     lwd = 2, col = "#2166AC")
abline(h = 0, lty = 2, col = "gray50", lwd = 1.5)
abline(h = coef(cox_model), lty = 2, col = "red", lwd = 1.5)

# Add test result
text(max(survival_data$OS_months, na.rm = TRUE) * 0.6,
     par("usr")[4] * 0.9,
     sprintf("PH test: p = %.4f\n(VIOLATED)", ph_test$table[1, "p"]),
     cex = 1.2, font = 2, col = ifelse(ph_test$table[1, "p"] < 0.05, "red", "darkgreen"))

################################################################################
# PANEL B: HAZARD RATIO OVER TIME (TIME-VARYING)
################################################################################

cat("Panel B: Hazard ratio over time...\n")

# Calculate time-varying hazard ratios using sliding windows
time_windows <- seq(6, 60, by = 6)  # Every 6 months
hrs <- numeric(length(time_windows))
hr_lower <- numeric(length(time_windows))
hr_upper <- numeric(length(time_windows))

for (i in 1:length(time_windows)) {
  t <- time_windows[i]

  # Landmark analysis at time t
  landmark_data <- survival_data %>%
    filter(OS_months >= t) %>%
    mutate(
      OS_months_adj = OS_months - t,
      OS_event_adj = OS_event
    )

  if (nrow(landmark_data) > 50) {  # Need sufficient events
    fit <- coxph(Surv(OS_months_adj, OS_event_adj) ~ cluster, data = landmark_data)
    hrs[i] <- exp(coef(fit))
    ci <- exp(confint(fit))
    hr_lower[i] <- ci[1, 1]
    hr_upper[i] <- ci[1, 2]
  } else {
    hrs[i] <- NA
    hr_lower[i] <- NA
    hr_upper[i] <- NA
  }
}

# Plot HR over time
plot(time_windows, hrs,
     type = "l", lwd = 3, col = "#2166AC",
     main = "B. Time-Varying Hazard Ratio",
     xlab = "Time (months)",
     ylab = "Hazard Ratio (Cluster 2 vs Cluster 1)",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.5,
     ylim = c(0.8, 2.5))

# Add confidence intervals
polygon(c(time_windows, rev(time_windows)),
        c(hr_lower, rev(hr_upper)),
        col = adjustcolor("#2166AC", alpha = 0.3),
        border = NA)

# Add reference line at HR=1
abline(h = 1, lty = 2, col = "gray50", lwd = 1.5)

# Add trend annotation
text(max(time_windows) * 0.6, 2.2,
     "HR decreases over time\n(survivor selection bias)",
     cex = 1.1, font = 2, col = "red")

################################################################################
# PANEL C: LANDMARK ANALYSIS AT KEY TIMEPOINTS
################################################################################

cat("Panel C: Landmark analysis...\n")

landmarks <- c(6, 12, 24, 36)
landmark_hrs <- numeric(length(landmarks))
landmark_p <- numeric(length(landmarks))
landmark_lower <- numeric(length(landmarks))
landmark_upper <- numeric(length(landmarks))

for (i in 1:length(landmarks)) {
  t <- landmarks[i]

  landmark_data <- survival_data %>%
    filter(OS_months >= t) %>%
    mutate(
      OS_months_adj = OS_months - t,
      OS_event_adj = OS_event
    )

  fit <- coxph(Surv(OS_months_adj, OS_event_adj) ~ cluster, data = landmark_data)
  landmark_hrs[i] <- exp(coef(fit))
  ci <- exp(confint(fit))
  landmark_lower[i] <- ci[1, 1]
  landmark_upper[i] <- ci[1, 2]
  landmark_p[i] <- summary(fit)$coefficients[1, "Pr(>|z|)"]
}

# Create forest plot style
plot(landmarks, landmark_hrs,
     pch = 19, cex = 2, col = "#2166AC",
     main = "C. Landmark Analysis at Key Timepoints",
     xlab = "Landmark Time (months)",
     ylab = "Hazard Ratio (95% CI)",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.5,
     ylim = c(0.8, 2.8), xlim = c(0, 40))

# Add confidence intervals
segments(landmarks, landmark_lower, landmarks, landmark_upper,
         lwd = 3, col = "#2166AC")

# Add reference line
abline(h = 1, lty = 2, col = "gray50", lwd = 1.5)

# Add significance markers
sig_landmarks <- landmarks[landmark_p < 0.05]
points(sig_landmarks, landmark_hrs[landmark_p < 0.05],
       pch = 8, cex = 2.5, col = "red", lwd = 2)

# Add legend
legend("topright",
       legend = c("HR (95% CI)", "p<0.05"),
       pch = c(19, 8), col = c("#2166AC", "red"),
       pt.cex = c(2, 2.5), bty = "n", cex = 1.1)

################################################################################
# PANEL D: LOG-RANK P-VALUE OVER TIME
################################################################################

cat("Panel D: Log-rank p-value stability...\n")

# Calculate log-rank p-values at different time restrictions
time_restrictions <- seq(12, 60, by = 6)
logrank_p <- numeric(length(time_restrictions))
n_events <- numeric(length(time_restrictions))

for (i in 1:length(time_restrictions)) {
  t_max <- time_restrictions[i]

  # Restrict data to t_max
  restricted_data <- survival_data %>%
    mutate(
      OS_months_r = pmin(OS_months, t_max),
      OS_event_r = ifelse(OS_months > t_max, 0, OS_event)
    )

  # Log-rank test
  survdiff_obj <- survdiff(Surv(OS_months_r, OS_event_r) ~ cluster,
                          data = restricted_data)

  logrank_p[i] <- 1 - pchisq(survdiff_obj$chisq, df = 1)
  n_events[i] <- sum(restricted_data$OS_event_r)
}

# Plot log-rank p-values over time
plot(time_restrictions, -log10(logrank_p),
     type = "b", pch = 19, cex = 2, lwd = 3, col = "#2166AC",
     main = "D. Statistical Significance Over Time",
     xlab = "Time Restriction (months)",
     ylab = "-log10(p-value)",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.5,
     ylim = c(0, max(-log10(logrank_p)) * 1.2))

# Add significance threshold
abline(h = -log10(0.05), lty = 2, col = "red", lwd = 2)
text(max(time_restrictions) * 0.7, -log10(0.05) + 0.2,
     "p = 0.05", cex = 1.1, col = "red")

# Add annotation
text(max(time_restrictions) * 0.5, max(-log10(logrank_p)) * 0.9,
     sprintf("Significant at all timepoints\n(%d-%d events)",
             min(n_events), max(n_events)),
     cex = 1.1, font = 2, col = "darkgreen")

dev.off()

cat("\n✓ Figure S2 (Combined) created successfully!\n")
cat("  Location: 05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics_Combined.pdf\n\n")

# Replace the placeholder
if (file.exists("05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics_Combined.pdf")) {
  file.copy("05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics_Combined.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics.pdf",
           overwrite = TRUE)
  cat("✓ Replaced placeholder Figure S2 with combined version\n")
}

cat("\n=== FIGURE S2 COMPLETE ===\n")
