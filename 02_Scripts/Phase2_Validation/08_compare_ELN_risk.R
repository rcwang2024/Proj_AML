# TASK 3: Compare Molecular Subtypes to ELN 2017 Risk Classification
# Determine if molecular subtypes add prognostic value beyond ELN risk

library(tidyverse)
library(readxl)
library(survival)
library(ggplot2)

cat("=== TASK 3: COMPARE SUBTYPES TO ELN RISK CLASSIFICATION ===\n\n")

# Load clinical data
cat("Loading clinical data...\n")
clinical <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
cat("Clinical data loaded:", nrow(clinical), "samples\n")
cat("ELN2017 column present:", "ELN2017" %in% colnames(clinical), "\n\n")

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
cat("Cluster assignments:", nrow(clusters), "samples\n\n")

# === 1. MERGE DATA ===
cat("=== 1. MERGING CLINICAL AND CLUSTER DATA ===\n\n")

# Merge on RNA sample ID
merged_data <- clinical %>%
  filter(!is.na(dbgap_rnaseq_sample)) %>%
  inner_join(clusters, by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  select(sample_id = dbgap_rnaseq_sample,
         cluster,
         ELN2017,
         age = ageAtDiagnosis,
         sex = consensus_sex,
         vitalStatus,
         overallSurvival) %>%
  filter(!is.na(ELN2017)) %>%  # Remove samples without ELN classification
  filter(!is.na(cluster))       # Remove samples without cluster assignment

cat("Matched samples with both ELN risk and cluster:", nrow(merged_data), "\n")
cat("Samples by ELN2017 category:\n")
print(table(merged_data$ELN2017))
cat("\n")

# === 2. CONTINGENCY TABLE ===
cat("=== 2. ELN RISK DISTRIBUTION BY MOLECULAR SUBTYPE ===\n\n")

# Create contingency table
contingency <- table(Cluster = merged_data$cluster, ELN_Risk = merged_data$ELN2017)
cat("Contingency table:\n")
print(contingency)
cat("\n")

# Proportions within each cluster
cat("Proportions within each cluster:\n")
prop_table <- prop.table(contingency, margin = 1) * 100
print(round(prop_table, 1))
cat("\n")

# Chi-square test for association
chi_test <- chisq.test(contingency)
cat("Chi-square test for independence:\n")
cat("  Chi-square statistic:", round(chi_test$statistic, 3), "\n")
cat("  p-value:", format(chi_test$p.value, scientific = TRUE), "\n")

if (chi_test$p.value < 0.001) {
  cat("  ✓ HIGHLY SIGNIFICANT association (p < 0.001)\n")
} else if (chi_test$p.value < 0.05) {
  cat("  ✓ Significant association (p < 0.05)\n")
} else {
  cat("  ~ No significant association (p >= 0.05)\n")
}
cat("\n")

# === 3. SURVIVAL ANALYSIS ===
cat("=== 3. SURVIVAL STRATIFICATION COMPARISON ===\n\n")

# Prepare survival data
survival_data <- merged_data %>%
  mutate(
    OS_event = ifelse(vitalStatus == "Dead", 1, 0),
    OS_months = overallSurvival,
    ELN_favorable = ifelse(ELN2017 == "Favorable", "Favorable", "Intermediate/Adverse")
  ) %>%
  filter(!is.na(OS_months) & !is.na(OS_event))

cat("Samples with complete survival data:", nrow(survival_data), "\n\n")

# Create survival object
surv_obj <- Surv(time = survival_data$OS_months, event = survival_data$OS_event)

# --- A. ELN Risk Alone ---
cat("A. ELN 2017 Risk Stratification:\n")
fit_eln <- survfit(surv_obj ~ ELN2017, data = survival_data)
survdiff_eln <- survdiff(surv_obj ~ ELN2017, data = survival_data)
pval_eln <- 1 - pchisq(survdiff_eln$chisq, df = length(unique(survival_data$ELN2017)) - 1)

cat("  p-value (log-rank):", format(pval_eln, scientific = TRUE), "\n")
cat("  Median survival by ELN risk:\n")
print(summary(fit_eln)$table[, "median"])
cat("\n")

# --- B. Molecular Subtype Alone ---
cat("B. Molecular Subtype Stratification:\n")
fit_cluster <- survfit(surv_obj ~ cluster, data = survival_data)
survdiff_cluster <- survdiff(surv_obj ~ cluster, data = survival_data)
pval_cluster <- 1 - pchisq(survdiff_cluster$chisq, df = 1)

cat("  p-value (log-rank):", format(pval_cluster, scientific = TRUE), "\n")
cat("  Median survival by cluster:\n")
print(summary(fit_cluster)$table[, "median"])
cat("\n")

# --- C. Combined ELN + Molecular Subtype ---
cat("C. Combined ELN Risk + Molecular Subtype:\n")
survival_data$combined <- paste0("ELN_", survival_data$ELN2017, "_Cluster_", survival_data$cluster)
fit_combined <- survfit(surv_obj ~ combined, data = survival_data)
survdiff_combined <- survdiff(surv_obj ~ combined, data = survival_data)
pval_combined <- 1 - pchisq(survdiff_combined$chisq,
                            df = length(unique(survival_data$combined)) - 1)

cat("  p-value (log-rank):", format(pval_combined, scientific = TRUE), "\n")
cat("  Median survival by combined groups:\n")
medians_combined <- summary(fit_combined)$table[, "median"]
print(sort(medians_combined, decreasing = TRUE))
cat("\n")

# === 4. COX MODELS FOR INDEPENDENT VALUE ===
cat("=== 4. MULTIVARIATE ANALYSIS ===\n\n")

# Model 1: ELN alone
cox_eln <- coxph(surv_obj ~ ELN2017, data = survival_data)
cat("Model 1: ELN2017 alone\n")
print(summary(cox_eln)$coefficients)
cat("  C-index:", summary(cox_eln)$concordance[1], "\n\n")

# Model 2: Cluster alone
cox_cluster <- coxph(surv_obj ~ factor(cluster), data = survival_data)
cat("Model 2: Molecular subtype alone\n")
print(summary(cox_cluster)$coefficients)
cat("  C-index:", summary(cox_cluster)$concordance[1], "\n\n")

# Model 3: ELN + Cluster (test if cluster adds value)
cox_both <- coxph(surv_obj ~ ELN2017 + factor(cluster), data = survival_data)
cat("Model 3: ELN2017 + Molecular subtype\n")
print(summary(cox_both)$coefficients)
cat("  C-index:", summary(cox_both)$concordance[1], "\n\n")

# Likelihood ratio test
cat("Likelihood ratio test: Does molecular subtype add value beyond ELN?\n")
lrt <- anova(cox_eln, cox_both, test = "Chisq")
print(lrt)
cat("\n")

lrt_pval <- lrt$`Pr(>|Chi|)`[2]
if (!is.na(lrt_pval) && lrt_pval < 0.05) {
  cat("✓ Molecular subtypes ADD INDEPENDENT PROGNOSTIC VALUE (p < 0.05)\n")
} else {
  cat("~ Molecular subtypes do not add significant independent value (p >= 0.05)\n")
}
cat("\n")

# Model 4: Full model with age and sex
cox_full <- coxph(surv_obj ~ ELN2017 + factor(cluster) + age + sex, data = survival_data)
cat("Model 4: Full model (ELN + Cluster + Age + Sex)\n")
print(summary(cox_full)$coefficients)
cat("  C-index:", summary(cox_full)$concordance[1], "\n\n")

# === 5. SUBGROUP ANALYSIS ===
cat("=== 5. MOLECULAR SUBTYPE VALUE WITHIN ELN CATEGORIES ===\n\n")

# For each ELN category, test if molecular subtype stratifies survival
eln_categories <- unique(survival_data$ELN2017)

subgroup_results <- data.frame()

for (eln_cat in eln_categories) {
  cat("--- ELN", eln_cat, "patients ---\n")

  subset_data <- survival_data %>% filter(ELN2017 == eln_cat)

  if (length(unique(subset_data$cluster)) < 2) {
    cat("  Only one cluster present, skipping\n\n")
    next
  }

  n_samples <- nrow(subset_data)
  cat("  N =", n_samples, "\n")

  # Test survival difference by cluster
  surv_subset <- Surv(time = subset_data$OS_months, event = subset_data$OS_event)
  fit_subset <- survfit(surv_subset ~ cluster, data = subset_data)
  survdiff_subset <- survdiff(surv_subset ~ cluster, data = subset_data)
  pval_subset <- 1 - pchisq(survdiff_subset$chisq, df = 1)

  cat("  p-value:", format(pval_subset, scientific = TRUE), "\n")

  medians_subset <- summary(fit_subset)$table[, "median"]
  cat("  Median survival by cluster:\n")
  print(medians_subset)
  cat("\n")

  subgroup_results <- rbind(subgroup_results, data.frame(
    ELN_category = eln_cat,
    n_samples = n_samples,
    p_value = pval_subset,
    median_cluster1 = medians_subset[1],
    median_cluster2 = medians_subset[2],
    stringsAsFactors = FALSE
  ))
}

cat("Subgroup analysis summary:\n")
print(subgroup_results, row.names = FALSE)
cat("\n")

# === 6. CREATE VISUALIZATIONS ===
cat("=== 6. CREATING VISUALIZATIONS ===\n\n")

# Create output directory
dir.create("04_Figures/12_ELN_Comparison", showWarnings = FALSE, recursive = TRUE)

# Figure 1: Stacked bar plot of ELN distribution by cluster
p1 <- ggplot(merged_data, aes(x = factor(cluster), fill = ELN2017)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Favorable" = "#00BA38",
                               "Intermediate" = "#F8766D",
                               "Adverse" = "#619CFF"),
                   name = "ELN 2017 Risk") +
  labs(x = "Molecular Subtype",
       y = "Proportion",
       title = "ELN 2017 Risk Distribution by Molecular Subtype") +
  scale_x_discrete(labels = c("1" = "Proliferative", "2" = "Immune-Inflammatory")) +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("04_Figures/12_ELN_Comparison/eln_distribution_by_cluster.pdf",
       p1, width = 8, height = 6)
cat("✓ Saved: 04_Figures/12_ELN_Comparison/eln_distribution_by_cluster.pdf\n")

# Figure 2: Survival curves - Combined stratification
survival_data_plot <- survival_data %>%
  mutate(group = paste0(ifelse(cluster == 1, "Proliferative", "Immune-Inflammatory"),
                        "\n", ELN2017))

fit_plot <- survfit(Surv(OS_months, OS_event) ~ group, data = survival_data_plot)

# Extract survival data for plotting
surv_df <- data.frame(
  time = fit_plot$time,
  surv = fit_plot$surv,
  group = rep(names(fit_plot$strata), fit_plot$strata)
)

p2 <- ggplot(surv_df, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 1) +
  labs(x = "Months",
       y = "Overall Survival Probability",
       title = "Survival by Molecular Subtype and ELN Risk",
       subtitle = paste0("Combined stratification p = ",
                        format(pval_combined, scientific = TRUE, digits = 2))) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank())

ggsave("04_Figures/12_ELN_Comparison/survival_combined_stratification.pdf",
       p2, width = 10, height = 6)
cat("✓ Saved: 04_Figures/12_ELN_Comparison/survival_combined_stratification.pdf\n")

# === 7. SAVE RESULTS ===
cat("\n=== 7. SAVING RESULTS ===\n\n")

dir.create("03_Results/12_ELN_Comparison", showWarnings = FALSE, recursive = TRUE)

# Save contingency table
write.csv(contingency,
          "03_Results/12_ELN_Comparison/eln_cluster_contingency.csv")

# Save merged data
write.csv(merged_data,
          "03_Results/12_ELN_Comparison/samples_with_eln_and_cluster.csv",
          row.names = FALSE)

# Save survival comparison summary
comparison_summary <- data.frame(
  model = c("ELN alone", "Cluster alone", "ELN + Cluster", "Full (ELN+Cluster+Age+Sex)"),
  p_value = c(pval_eln, pval_cluster, pval_combined, NA),
  c_index = c(summary(cox_eln)$concordance[1],
              summary(cox_cluster)$concordance[1],
              summary(cox_both)$concordance[1],
              summary(cox_full)$concordance[1]),
  stringsAsFactors = FALSE
)

write.csv(comparison_summary,
          "03_Results/12_ELN_Comparison/survival_model_comparison.csv",
          row.names = FALSE)

# Save subgroup analysis
write.csv(subgroup_results,
          "03_Results/12_ELN_Comparison/subgroup_analysis_by_eln.csv",
          row.names = FALSE)

cat("✓ Saved: 03_Results/12_ELN_Comparison/eln_cluster_contingency.csv\n")
cat("✓ Saved: 03_Results/12_ELN_Comparison/samples_with_eln_and_cluster.csv\n")
cat("✓ Saved: 03_Results/12_ELN_Comparison/survival_model_comparison.csv\n")
cat("✓ Saved: 03_Results/12_ELN_Comparison/subgroup_analysis_by_eln.csv\n")

# === 8. SUMMARY ===
cat("\n=== SUMMARY ===\n\n")

cat("Samples analyzed:", nrow(merged_data), "\n")
cat("Samples with survival data:", nrow(survival_data), "\n\n")

cat("Key Findings:\n")
cat("1. Association between ELN risk and molecular subtypes:\n")
cat("   Chi-square p-value:", format(chi_test$p.value, scientific = TRUE), "\n\n")

cat("2. Survival stratification:\n")
cat("   ELN alone p-value:", format(pval_eln, scientific = TRUE), "\n")
cat("   Cluster alone p-value:", format(pval_cluster, scientific = TRUE), "\n")
cat("   Combined p-value:", format(pval_combined, scientific = TRUE), "\n\n")

cat("3. Independent prognostic value:\n")
cat("   Likelihood ratio test p-value:", format(lrt_pval, scientific = TRUE), "\n")

if (!is.na(lrt_pval) && lrt_pval < 0.05) {
  cat("   ✓ Molecular subtypes ADD INDEPENDENT VALUE beyond ELN risk\n")
} else {
  cat("   ~ Molecular subtypes captured by ELN risk classification\n")
}

cat("\n### Task 3 COMPLETE ###\n")
