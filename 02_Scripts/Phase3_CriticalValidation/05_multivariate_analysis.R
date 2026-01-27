# PHASE 3 - PART 2: MULTIVARIATE ANALYSIS
# Test if molecular subtypes have INDEPENDENT prognostic value
# beyond clinical variables (age, sex) and key mutations

library(tidyverse)
library(survival)
library(survminer)
library(MASS)  # For stepAIC

cat("=== PHASE 3: PART 2 - MULTIVARIATE ANALYSIS ===\n\n")

# ============================================================================
# 1. LOAD AND MERGE DATA
# ============================================================================

cat("=== STEP 1: LOADING AND MERGING DATA ===\n\n")

# Load survival + cluster data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Load mutation data
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

cat("Survival data:", nrow(survival_data), "samples\n")
cat("Mutation data:", nrow(mutations), "samples x", ncol(mutations), "genes\n\n")

# Merge datasets
merged_data <- survival_data %>%
  left_join(
    mutations %>% rownames_to_column("sample_id"),
    by = "sample_id"
  )

cat("Merged data:", nrow(merged_data), "samples\n")
cat("Variables:", ncol(merged_data), "\n")
cat("Samples with complete survival:", sum(!is.na(merged_data$OS_months) & merged_data$OS_months > 0), "\n\n")

# Clean data
analysis_data <- merged_data %>%
  filter(!is.na(OS_months) & OS_months > 0) %>%
  rename(
    OS_MONTHS = OS_months,
    OS_STATUS = OS_event,
    cluster_assignment = cluster,
    AGE = age,
    SEX = sex
  ) %>%
  mutate(
    SEX = factor(SEX, levels = c("F", "M")),
    cluster_assignment = factor(cluster_assignment, levels = c(1, 2),
                                labels = c("Cluster1", "Cluster2"))
  ) %>%
  # Filter to complete clinical data
  filter(!is.na(AGE) & !is.na(SEX) & !is.na(cluster_assignment) &
         SEX %in% c("F", "M"))

cat("Analysis dataset: ", nrow(analysis_data), "samples\n")
cat("Events:", sum(analysis_data$OS_STATUS), "\n")
cat("Missing AGE:", sum(is.na(analysis_data$AGE)), "\n")
cat("Missing SEX:", sum(is.na(analysis_data$SEX)), "\n")
cat("Missing cluster:", sum(is.na(analysis_data$cluster_assignment)), "\n\n")

# ============================================================================
# 2. KEY MUTATIONS FOR AML PROGNOSIS
# ============================================================================

cat("=== STEP 2: SELECTING KEY PROGNOSTIC MUTATIONS ===\n\n")

# Key mutations known to be prognostic in AML
key_mutations <- c("FLT3", "NPM1", "DNMT3A", "TP53", "IDH1", "IDH2",
                   "RUNX1", "ASXL1", "TET2", "NRAS", "KRAS")

# Check which are available
available_mutations <- key_mutations[key_mutations %in% colnames(analysis_data)]
cat("Key mutations available:", length(available_mutations), "of", length(key_mutations), "\n")
cat("  ", paste(available_mutations, collapse=", "), "\n\n")

# Mutation frequencies
mut_freq <- colSums(analysis_data[, available_mutations], na.rm = TRUE) / nrow(analysis_data) * 100

cat("Mutation frequencies:\n")
for (gene in available_mutations) {
  cat("  ", gene, ": ", round(mut_freq[gene], 1), "%\n", sep="")
}
cat("\n")

# ============================================================================
# 3. UNIVARIATE ANALYSES
# ============================================================================

cat("=== STEP 3: UNIVARIATE SURVIVAL ANALYSES ===\n\n")

# Test each variable individually
univariate_results <- data.frame()

# Clinical variables
for (var in c("cluster_assignment", "AGE", "SEX")) {

  if (var == "AGE") {
    # Continuous variable
    cox_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ AGE, data = analysis_data)
  } else {
    # Categorical variable
    formula_str <- paste("Surv(OS_MONTHS, OS_STATUS) ~", var)
    cox_model <- coxph(as.formula(formula_str), data = analysis_data)
  }

  hr <- exp(coef(cox_model)[1])
  ci <- exp(confint(cox_model)[1, ])
  pval <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]

  univariate_results <- rbind(univariate_results, data.frame(
    variable = var,
    hazard_ratio = hr,
    hr_lower = ci[1],
    hr_upper = ci[2],
    pvalue = pval,
    significant = pval < 0.05
  ))
}

# Mutation variables
for (gene in available_mutations) {

  # Check for sufficient non-NA data
  non_na_count <- sum(!is.na(analysis_data[[gene]]))
  mut_count <- sum(analysis_data[[gene]], na.rm = TRUE)

  if (non_na_count < 100) {
    cat("  Skipping", gene, "(too many NAs:", (nrow(analysis_data) - non_na_count), ")\n")
    next
  }

  if (mut_count < 10) {
    cat("  Skipping", gene, "(too few mutations:", mut_count, ")\n")
    next
  }

  # Filter to complete cases for this gene
  gene_data <- analysis_data %>%
    filter(!is.na(!!sym(gene)))

  formula_str <- paste("Surv(OS_MONTHS, OS_STATUS) ~", gene)
  cox_model <- coxph(as.formula(formula_str), data = gene_data)

  hr <- exp(coef(cox_model)[1])
  ci <- exp(confint(cox_model)[1, ])
  pval <- summary(cox_model)$coefficients[1, "Pr(>|z|)"]

  univariate_results <- rbind(univariate_results, data.frame(
    variable = gene,
    hazard_ratio = hr,
    hr_lower = ci[1],
    hr_upper = ci[2],
    pvalue = pval,
    significant = pval < 0.05
  ))
}

cat("\nUnivariate Analysis Results:\n")
print(univariate_results %>% arrange(pvalue))
cat("\n")

cat("Significant variables (p<0.05):", sum(univariate_results$significant), "\n\n")

# ============================================================================
# 4. MULTIVARIATE MODEL: CLINICAL + CLUSTER
# ============================================================================

cat("=== STEP 4: MULTIVARIATE MODEL - CLINICAL + CLUSTER ===\n\n")

# Model 1: Clinical variables only
cat("Model 1: Clinical variables only (age, sex)\n")
model_clinical <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ AGE + SEX,
                        data = analysis_data)
print(summary(model_clinical))
cat("\n")

# Model 2: Clinical + Cluster
cat("Model 2: Clinical + Molecular Subtype\n")
model_clinical_cluster <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ AGE + SEX + cluster_assignment,
                                data = analysis_data)
print(summary(model_clinical_cluster))
cat("\n")

# Likelihood ratio test
lrt <- anova(model_clinical, model_clinical_cluster)
cat("Likelihood Ratio Test (adding cluster to clinical model):\n")
print(lrt)
cat("\n")

lrt_pvalue <- lrt$`Pr(>|Chi|)`[2]
cat("LRT p-value:", format(lrt_pvalue, scientific = TRUE, digits = 3), "\n")

if (lrt_pvalue < 0.05) {
  cat("✓ Cluster adds SIGNIFICANT prognostic information beyond clinical variables\n\n")
} else {
  cat("✗ Cluster does NOT significantly improve model beyond clinical variables\n\n")
}

# ============================================================================
# 5. MULTIVARIATE MODEL: MUTATIONS ONLY
# ============================================================================

cat("=== STEP 5: MULTIVARIATE MODEL - KEY MUTATIONS ===\n\n")

# Select mutations with sufficient frequency (>=5%) and prognostic in univariate
significant_mutations <- univariate_results %>%
  filter(variable %in% available_mutations, pvalue < 0.1) %>%
  pull(variable)

cat("Mutations with p<0.1 in univariate:", paste(significant_mutations, collapse=", "), "\n\n")

# Only include mutations present in >=5% of samples
sufficient_mutations <- available_mutations[mut_freq >= 5]
cat("Mutations with >=5% frequency:", paste(sufficient_mutations, collapse=", "), "\n\n")

# Use intersection
model_mutations <- intersect(significant_mutations, sufficient_mutations)

if (length(model_mutations) == 0) {
  cat("No mutations meet both criteria, using top 5 most frequent mutations\n")
  model_mutations <- names(sort(mut_freq, decreasing = TRUE))[1:min(5, length(mut_freq))]
}

cat("Mutations included in model:", paste(model_mutations, collapse=", "), "\n\n")

# Fit mutation-only model
formula_mut <- paste("Surv(OS_MONTHS, OS_STATUS) ~",
                     paste(model_mutations, collapse = " + "))
model_mutations_only <- coxph(as.formula(formula_mut), data = analysis_data)

cat("Model 3: Key mutations only\n")
print(summary(model_mutations_only))
cat("\n")

# ============================================================================
# 6. FULL MULTIVARIATE MODEL: CLINICAL + MUTATIONS + CLUSTER
# ============================================================================

cat("=== STEP 6: FULL MULTIVARIATE MODEL ===\n\n")

# Model 4: Clinical + Mutations
formula_clin_mut <- paste("Surv(OS_MONTHS, OS_STATUS) ~ AGE + SEX +",
                          paste(model_mutations, collapse = " + "))
model_clinical_mutations <- coxph(as.formula(formula_clin_mut), data = analysis_data)

cat("Model 4: Clinical + Mutations\n")
print(summary(model_clinical_mutations))
cat("\n")

# Model 5: Clinical + Mutations + Cluster
formula_full <- paste(formula_clin_mut, "+ cluster_assignment")
model_full <- coxph(as.formula(formula_full), data = analysis_data)

cat("Model 5: Clinical + Mutations + Molecular Subtype (FULL MODEL)\n")
print(summary(model_full))
cat("\n")

# Test if cluster adds value beyond clinical + mutations
lrt_full <- anova(model_clinical_mutations, model_full)
cat("Likelihood Ratio Test (adding cluster to clinical + mutations):\n")
print(lrt_full)
cat("\n")

lrt_full_pvalue <- lrt_full$`Pr(>|Chi|)`[2]
cat("LRT p-value:", format(lrt_full_pvalue, scientific = TRUE, digits = 3), "\n")

if (lrt_full_pvalue < 0.05) {
  cat("✓✓✓ CRITICAL FINDING: Cluster adds INDEPENDENT prognostic value\n")
  cat("    beyond both clinical variables AND key mutations\n")
  cat("    → Molecular subtypes capture biology not explained by standard markers\n\n")
} else {
  cat("✗ Cluster does NOT add significant value beyond clinical + mutations\n")
  cat("  → Subtypes may be proxies for known mutations\n\n")
}

# ============================================================================
# 7. MODEL COMPARISON & SELECTION
# ============================================================================

cat("=== STEP 7: MODEL COMPARISON ===\n\n")

# Compare models by AIC, concordance
model_comparison <- data.frame(
  model = c("Clinical only", "Clinical + Cluster", "Mutations only",
            "Clinical + Mutations", "Full (Clinical + Mutations + Cluster)"),
  variables = c(
    "Age, Sex",
    "Age, Sex, Cluster",
    paste(model_mutations, collapse=", "),
    paste("Age, Sex,", paste(model_mutations, collapse=", ")),
    paste("Age, Sex,", paste(model_mutations, collapse=", "), ", Cluster")
  ),
  aic = c(
    AIC(model_clinical),
    AIC(model_clinical_cluster),
    AIC(model_mutations_only),
    AIC(model_clinical_mutations),
    AIC(model_full)
  ),
  concordance = c(
    summary(model_clinical)$concordance[1],
    summary(model_clinical_cluster)$concordance[1],
    summary(model_mutations_only)$concordance[1],
    summary(model_clinical_mutations)$concordance[1],
    summary(model_full)$concordance[1]
  ),
  loglik = c(
    model_clinical$loglik[2],
    model_clinical_cluster$loglik[2],
    model_mutations_only$loglik[2],
    model_clinical_mutations$loglik[2],
    model_full$loglik[2]
  )
)

model_comparison$delta_aic <- model_comparison$aic - min(model_comparison$aic)
model_comparison$best_model <- model_comparison$delta_aic < 2

cat("Model Comparison Table:\n")
print(model_comparison)
cat("\n")

best_model_name <- model_comparison$model[which.min(model_comparison$aic)]
cat("Best model by AIC:", best_model_name, "\n\n")

# ============================================================================
# 8. MUTATION ENRICHMENT BY CLUSTER
# ============================================================================

cat("=== STEP 8: MUTATION ENRICHMENT BY CLUSTER ===\n\n")

enrichment_results <- data.frame()

for (gene in available_mutations) {

  # Create contingency table
  cluster1_mut <- sum(analysis_data$cluster_assignment == "Cluster1" &
                      analysis_data[[gene]] == 1, na.rm = TRUE)
  cluster1_wt <- sum(analysis_data$cluster_assignment == "Cluster1" &
                     analysis_data[[gene]] == 0, na.rm = TRUE)
  cluster2_mut <- sum(analysis_data$cluster_assignment == "Cluster2" &
                      analysis_data[[gene]] == 1, na.rm = TRUE)
  cluster2_wt <- sum(analysis_data$cluster_assignment == "Cluster2" &
                     analysis_data[[gene]] == 0, na.rm = TRUE)

  # Fisher's exact test
  if (cluster1_mut + cluster2_mut >= 5) {
    fisher_test <- fisher.test(matrix(c(cluster1_mut, cluster1_wt,
                                         cluster2_mut, cluster2_wt), nrow = 2))

    enrichment_results <- rbind(enrichment_results, data.frame(
      gene = gene,
      cluster1_freq = round(cluster1_mut / (cluster1_mut + cluster1_wt) * 100, 1),
      cluster2_freq = round(cluster2_mut / (cluster2_mut + cluster2_wt) * 100, 1),
      pvalue = fisher_test$p.value,
      or = fisher_test$estimate,
      significant = fisher_test$p.value < 0.05
    ))
  }
}

enrichment_results <- enrichment_results %>% arrange(pvalue)

cat("Mutation Enrichment by Cluster:\n")
print(enrichment_results)
cat("\n")

cat("Significant enrichments (p<0.05):", sum(enrichment_results$significant, na.rm = TRUE), "\n\n")

# ============================================================================
# 9. SAVE RESULTS
# ============================================================================

cat("=== STEP 9: SAVING RESULTS ===\n\n")

write.csv(univariate_results,
          "03_Results/11_Survival_Reanalysis/05_univariate_results.csv",
          row.names = FALSE)
cat("✓ Saved: 05_univariate_results.csv\n")

write.csv(model_comparison,
          "03_Results/11_Survival_Reanalysis/05_model_comparison.csv",
          row.names = FALSE)
cat("✓ Saved: 05_model_comparison.csv\n")

write.csv(enrichment_results,
          "03_Results/11_Survival_Reanalysis/05_mutation_enrichment.csv",
          row.names = FALSE)
cat("✓ Saved: 05_mutation_enrichment.csv\n")

# Save full model coefficients
full_model_coef <- as.data.frame(summary(model_full)$coefficients)
full_model_coef$variable <- rownames(full_model_coef)
full_model_coef <- full_model_coef %>%
  mutate(
    HR = exp(coef),
    HR_lower = exp(coef - 1.96 * `se(coef)`),
    HR_upper = exp(coef + 1.96 * `se(coef)`)
  ) %>%
  dplyr::select(variable, coef, HR, HR_lower, HR_upper, `Pr(>|z|)`) %>%
  rename(pvalue = `Pr(>|z|)`)

write.csv(full_model_coef,
          "03_Results/11_Survival_Reanalysis/05_full_model_coefficients.csv",
          row.names = FALSE)
cat("✓ Saved: 05_full_model_coefficients.csv\n\n")

# ============================================================================
# 10. VISUALIZATION
# ============================================================================

cat("=== STEP 10: CREATING VISUALIZATIONS ===\n\n")

# Forest plot of full model
full_model_plot_data <- full_model_coef %>%
  filter(variable != "AGE") %>%
  mutate(
    variable = gsub("cluster_assignment", "Subtype: ", variable),
    variable = gsub("SEX", "Sex: ", variable)
  )

p_forest <- ggplot(full_model_plot_data, aes(x = HR, y = reorder(variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.2) +
  labs(
    title = "Full Multivariate Model: Hazard Ratios",
    subtitle = "Adjusted for age, sex, key mutations, and molecular subtype",
    x = "Hazard Ratio (95% CI)",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("04_Figures/11_Survival_Reanalysis/05_multivariate_forest_plot.pdf",
       p_forest, width = 10, height = 6)
cat("✓ Saved: 05_multivariate_forest_plot.pdf\n\n")

# ============================================================================
# 11. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("PART 2: MULTIVARIATE ANALYSIS COMPLETE\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. UNIVARIATE ANALYSIS:\n")
cat("   Variables tested:", nrow(univariate_results), "\n")
cat("   Significant (p<0.05):", sum(univariate_results$significant), "\n")
cat("   Cluster univariate p-value:",
    format(univariate_results$pvalue[univariate_results$variable == "cluster_assignment"], scientific = TRUE, digits = 3), "\n\n")

cat("2. CLINICAL + CLUSTER MODEL:\n")
cat("   LRT p-value:", format(lrt_pvalue, scientific = TRUE, digits = 3), "\n")
if (lrt_pvalue < 0.05) {
  cat("   ✓ Cluster adds value beyond clinical variables\n\n")
} else {
  cat("   ✗ Cluster does not add value beyond clinical variables\n\n")
}

cat("3. FULL MODEL (Clinical + Mutations + Cluster):\n")
cat("   LRT p-value (cluster):", format(lrt_full_pvalue, scientific = TRUE, digits = 3), "\n")
if (lrt_full_pvalue < 0.05) {
  cat("   ✓✓✓ CLUSTER HAS INDEPENDENT PROGNOSTIC VALUE\n")
  cat("       → Molecular subtypes are NOT just proxies for mutations\n")
  cat("       → They capture additional biological information\n\n")
} else {
  cat("   ✗ Cluster is explained by known mutations\n\n")
}

cat("4. MODEL SELECTION:\n")
cat("   Best model:", best_model_name, "\n")
cat("   Concordance:", round(model_comparison$concordance[which.min(model_comparison$aic)], 3), "\n\n")

cat("5. MUTATION ENRICHMENT:\n")
if (nrow(enrichment_results) > 0) {
  sig_enrich <- enrichment_results %>% filter(significant)
  if (nrow(sig_enrich) > 0) {
    cat("   Significant enrichments:", nrow(sig_enrich), "\n")
    for (i in 1:min(3, nrow(sig_enrich))) {
      cat("     ", sig_enrich$gene[i], ": Cluster1=", sig_enrich$cluster1_freq[i],
          "%, Cluster2=", sig_enrich$cluster2_freq[i], "%, p=",
          format(sig_enrich$pvalue[i], scientific = TRUE, digits = 3), "\n", sep="")
    }
  } else {
    cat("   No significant enrichments\n")
  }
}

cat("\n### PART 2 COMPLETE ###\n")
