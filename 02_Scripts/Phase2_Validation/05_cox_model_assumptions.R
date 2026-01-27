# TASK 5: Check Cox Model Assumptions
# Verify proportional hazards assumption and check for influential observations

library(survival)
library(tidyverse)

cat("=== TASK 5: COX MODEL ASSUMPTION CHECKS ===\n\n")

# Load survival data with clusters
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

cat("Loaded survival data:", nrow(survival_data), "samples\n")
cat("Variables:", paste(colnames(survival_data), collapse = ", "), "\n\n")

# Fit Cox model
cat("Fitting Cox proportional hazards model...\n")
cox_model <- coxph(Surv(OS_months, OS_event) ~ factor(cluster) + age + sex,
                   data = survival_data)

cat("\n=== COX MODEL SUMMARY ===\n")
print(summary(cox_model))

cat("\n\n=== 1. TESTING PROPORTIONAL HAZARDS ASSUMPTION ===\n")
cat("(p > 0.05 indicates assumption is met)\n\n")

ph_test <- cox.zph(cox_model)
print(ph_test)

if (any(ph_test$table[, "p"] < 0.05)) {
  cat("\n⚠️ WARNING: Proportional hazards assumption violated for some variables\n")
  cat("Consider stratified Cox model or time-varying effects\n\n")

  # Plot Schoenfeld residuals
  cat("Saving Schoenfeld residual plots...\n")
  pdf("04_Figures/10_Model_Diagnostics/schoenfeld_residuals.pdf", width = 12, height = 8)
  plot(ph_test)
  dev.off()

  cat("✓ Saved: 04_Figures/10_Model_Diagnostics/schoenfeld_residuals.pdf\n")
} else {
  cat("\n✓ Proportional hazards assumption MET for all variables\n")
}

cat("\n\n=== 2. CHECKING FOR INFLUENTIAL OBSERVATIONS ===\n\n")

dfbeta_vals <- residuals(cox_model, type = "dfbeta")

# Plot dfbeta values
pdf("04_Figures/10_Model_Diagnostics/dfbeta_plots.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))
for (i in 1:min(ncol(dfbeta_vals), 4)) {
  plot(dfbeta_vals[, i],
       ylab = paste("dfbeta -", colnames(dfbeta_vals)[i]),
       main = colnames(dfbeta_vals)[i])
  abline(h = 0, lty = 2)
}
dev.off()

cat("✓ Saved: 04_Figures/10_Model_Diagnostics/dfbeta_plots.pdf\n\n")

# Identify influential points (>3 SD)
influential <- which(abs(dfbeta_vals[, 1]) > 3 * sd(dfbeta_vals[, 1]))
if (length(influential) > 0) {
  cat("Found", length(influential), "potentially influential observations\n")
  cat("Indices:", head(influential, 10), "...\n")
} else {
  cat("✓ No highly influential observations detected\n")
}

cat("\n\n=== 3. CHECKING LINEARITY OF AGE EFFECT ===\n\n")

# Compare continuous vs categorical age
survival_data$age_cat <- cut(survival_data$age,
                             breaks = c(0, 40, 60, 80, 100),
                             labels = c("<40", "40-60", "60-80", ">80"))

cox_age_cat <- coxph(Surv(OS_months, OS_event) ~ factor(cluster) + age_cat + sex,
                     data = survival_data)

aic_continuous <- AIC(cox_model)
aic_categorical <- AIC(cox_age_cat)

cat("AIC with continuous age:", round(aic_continuous, 1), "\n")
cat("AIC with categorical age:", round(aic_categorical, 1), "\n\n")

if (aic_continuous < aic_categorical) {
  cat("✓ Continuous age is appropriate (lower AIC)\n")
} else {
  cat("⚠️ Consider using categorical age (lower AIC)\n")
}

cat("\n\n=== 4. TESTING FOR INTERACTIONS ===\n\n")

# Test cluster:age interaction
cox_interaction <- coxph(Surv(OS_months, OS_event) ~
                           factor(cluster) * age + sex,
                         data = survival_data)

anova_int <- anova(cox_model, cox_interaction)
cat("Cluster:Age Interaction Test:\n")
print(anova_int)

if (!is.null(anova_int$`P(>|Chi|)`) && length(anova_int$`P(>|Chi|)`) >= 2) {
  if (anova_int$`P(>|Chi|)`[2] < 0.05) {
    cat("\n⚠️ Significant interaction detected (p < 0.05)\n")
    cat("Effect of cluster varies by age\n")
  } else {
    cat("\n✓ No significant interaction (p >= 0.05)\n")
  }
}

# Summary
cat("\n\n=== SUMMARY OF DIAGNOSTICS ===\n")
cat("1. Proportional Hazards: ", ifelse(any(ph_test$table[, "p"] < 0.05), "⚠️ Some violations", "✓ Passed"), "\n")
cat("2. Influential Observations: ", ifelse(length(influential) > 0, paste("⚠️", length(influential), "found"), "✓ None"), "\n")
cat("3. Age Linearity: ", ifelse(aic_continuous < aic_categorical, "✓ Linear OK", "⚠️ Consider categorical"), "\n")

cat("\n### Task 5 COMPLETE ###\n")
