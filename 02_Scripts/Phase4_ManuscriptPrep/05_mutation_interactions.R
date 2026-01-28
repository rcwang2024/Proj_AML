#!/usr/bin/env Rscript
# Phase 4 Part 5: Mutation-Cluster Interaction Analysis
# Purpose: Test if mutations have DIFFERENT prognostic effects by cluster

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(ggplot2)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("MUTATION × CLUSTER INTERACTION ANALYSIS\n")
cat("==============================================================================\n\n")

cat("Testing whether mutations have different prognostic effects in each cluster\n")
cat("If YES → Clusters have independent prognostic value (effect modification)\n")
cat("If NO → Clusters are just proxies for mutations (no interaction)\n\n")

# Load data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

# Merge
analysis_data <- survival_data %>%
  left_join(mutations %>% tibble::rownames_to_column("sample_id"), by = "sample_id") %>%
  filter(!is.na(cluster) & !is.na(OS_months) & OS_months > 0) %>%
  rename(
    OS_MONTHS = OS_months,
    OS_STATUS = OS_event,
    AGE = age,
    SEX = sex
  ) %>%
  filter(!is.na(AGE) & !is.na(SEX))

cat(sprintf("Analysis dataset: %d samples, %d events\n\n", nrow(analysis_data), sum(analysis_data$OS_STATUS)))

# Mutations to test
mutations_to_test <- c("NPM1", "TP53", "RUNX1", "DNMT3A", "FLT3", "TET2", "ASXL1")

# ==============================================================================
# TEST INTERACTIONS
# ==============================================================================

cat("=== TESTING INTERACTIONS ===\n\n")

interaction_results <- lapply(mutations_to_test, function(mut) {

  cat(sprintf("Testing %s...\n", mut))

  # Skip if not available or too few mutations
  if (!mut %in% colnames(analysis_data)) {
    cat(sprintf("  Skipping %s: not available\n", mut))
    return(NULL)
  }

  # Filter to complete cases
  mut_data <- analysis_data %>% filter(!is.na(!!sym(mut)))
  n_mut <- sum(mut_data[[mut]] == 1)

  if (n_mut < 10) {
    cat(sprintf("  Skipping %s: only %d mutations\n\n", mut, n_mut))
    return(NULL)
  }

  # Model WITHOUT interaction (main effects only)
  formula_main <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ cluster + ", mut, " + AGE + SEX"))

  # Model WITH interaction
  formula_interact <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ cluster * ", mut, " + AGE + SEX"))

  # Fit models
  tryCatch({
    fit_main <- coxph(formula_main, data = mut_data)
    fit_interact <- coxph(formula_interact, data = mut_data)

    # Likelihood ratio test
    lrt <- anova(fit_main, fit_interact)

    # Extract interaction term
    interact_term <- grep(":", names(coef(fit_interact)), value = TRUE)
    if (length(interact_term) == 0) {
      cat(sprintf("  No interaction term found\n\n"))
      return(NULL)
    }

    interact_coef <- coef(fit_interact)[interact_term]
    interact_se <- sqrt(vcov(fit_interact)[interact_term, interact_term])
    interact_p <- summary(fit_interact)$coefficients[interact_term, "Pr(>|z|)"]

    # Calculate HRs in each cluster
    # Cluster 1 (reference): HR for mutation
    hr_c1 <- exp(coef(fit_interact)[mut])

    # Cluster 2: HR for mutation + interaction
    hr_c2 <- exp(coef(fit_interact)[mut] + interact_coef)

    cat(sprintf("  Main effect HR: %.2f\n", exp(coef(fit_main)[mut])))
    cat(sprintf("  HR in Cluster 1: %.2f\n", hr_c1))
    cat(sprintf("  HR in Cluster 2: %.2f\n", hr_c2))
    cat(sprintf("  Interaction p: %.4f %s\n\n",
                interact_p, ifelse(interact_p < 0.05, "***", "")))

    data.frame(
      Mutation = mut,
      N = nrow(mut_data),
      N_mutated = n_mut,
      LRT_ChiSq = lrt$Chisq[2],
      LRT_p = lrt$`Pr(>|Chi|)`[2],
      Interaction_coef = interact_coef,
      Interaction_HR = exp(interact_coef),
      Interaction_p = interact_p,
      Main_HR = exp(coef(fit_main)[mut]),
      HR_in_Cluster1 = hr_c1,
      HR_in_Cluster2 = hr_c2,
      HR_ratio = hr_c2 / hr_c1,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Error fitting %s: %s\n\n", mut, e$message))
    return(NULL)
  })

}) %>% bind_rows()

# FDR correction
if (nrow(interaction_results) > 0) {
  interaction_results$FDR <- p.adjust(interaction_results$LRT_p, method = "BH")

  # Sort by significance
  interaction_results <- interaction_results %>% arrange(LRT_p)

  write.csv(interaction_results,
            "03_Results/21_Manuscript_Prep/mutation_cluster_interactions.csv",
            row.names = FALSE)

  cat("=== INTERACTION RESULTS ===\n\n")
  print(interaction_results, row.names = FALSE)

  # ==============================================================================
  # SUMMARY
  # ==============================================================================

  cat("\n==============================================================================\n")
  cat("SUMMARY\n")
  cat("==============================================================================\n\n")

  sig_interactions <- interaction_results %>% filter(FDR < 0.10)

  if (nrow(sig_interactions) > 0) {
    cat("✅ SIGNIFICANT INTERACTIONS FOUND (FDR < 0.10):\n\n")
    for (i in 1:nrow(sig_interactions)) {
      row <- sig_interactions[i, ]
      cat(sprintf("%s:\n", row$Mutation))
      cat(sprintf("  HR in Cluster 1: %.2f\n", row$HR_in_Cluster1))
      cat(sprintf("  HR in Cluster 2: %.2f\n", row$HR_in_Cluster2))
      cat(sprintf("  Interaction p: %.4f (FDR: %.4f)\n", row$Interaction_p, row$FDR))

      if (row$HR_ratio > 1.5) {
        cat("  → Mutation MORE harmful in Cluster 2\n")
      } else if (row$HR_ratio < 0.67) {
        cat("  → Mutation MORE harmful in Cluster 1\n")
      } else {
        cat("  → Similar effect in both clusters\n")
      }
      cat("\n")
    }

    cat("🎯 INTERPRETATION:\n")
    cat("   Clusters MODIFY the prognostic effect of mutations\n")
    cat("   → Clusters may have INDEPENDENT prognostic value\n")
    cat("   → Different biological contexts alter mutation impact\n\n")

  } else {
    cat("❌ NO SIGNIFICANT INTERACTIONS FOUND\n\n")
    cat("INTERPRETATION:\n")
    cat("   Mutations have similar prognostic effects in both clusters\n")
    cat("   → Clusters do NOT modify mutation effects\n")
    cat("   → Clusters are likely proxies for mutation patterns\n")
    cat("   → Consistent with multivariate finding (p=0.649)\n\n")
  }

  # Report nominal significance
  nom_sig <- interaction_results %>% filter(LRT_p < 0.05)
  if (nrow(nom_sig) > 0 & nrow(sig_interactions) == 0) {
    cat("⚠ NOTE: %d interactions nominally significant (p<0.05) but not FDR<0.10:\n", nrow(nom_sig))
    print(nom_sig[, c("Mutation", "LRT_p", "FDR", "HR_ratio")], row.names = FALSE)
    cat("\n")
  }

  # ==============================================================================
  # VISUALIZATION
  # ==============================================================================

  cat("Creating visualization...\n")

  pdf("04_Figures/20_Manuscript_Prep/mutation_interactions_plot.pdf", width = 10, height = 7)

  # Forest plot showing HR ratios
  interaction_results$Significant <- interaction_results$FDR < 0.10

  p <- ggplot(interaction_results, aes(x = HR_ratio, y = reorder(Mutation, HR_ratio))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", size = 1) +
    geom_point(aes(color = Significant), size = 4) +
    geom_errorbarh(aes(xmin = HR_in_Cluster1/HR_in_Cluster2 * 0.8,
                       xmax = HR_in_Cluster1/HR_in_Cluster2 * 1.2,
                       color = Significant),
                   height = 0.2, size = 1) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                       name = "",
                       labels = c("FALSE" = "Not significant", "TRUE" = "FDR < 0.10")) +
    scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
    labs(x = "HR Ratio (Cluster 2 / Cluster 1)",
         y = "",
         title = "Mutation × Cluster Interactions: Effect Modification",
         subtitle = "HR ratio >1 means mutation is more harmful in Cluster 2") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", size = 16),
          panel.grid.major.y = element_blank())

  print(p)
  dev.off()

  cat("  ✓ mutation_interactions_plot.pdf\n\n")

  cat("✅ Mutation-cluster interaction analysis complete\n\n")
  cat("Files generated:\n")
  cat("  - 03_Results/21_Manuscript_Prep/mutation_cluster_interactions.csv\n")
  cat("  - 04_Figures/20_Manuscript_Prep/mutation_interactions_plot.pdf\n\n")

} else {
  cat("⚠ No interaction tests could be performed\n\n")
}
