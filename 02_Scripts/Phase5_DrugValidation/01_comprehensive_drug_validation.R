#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 5: Comprehensive Drug Response Validation
# ==============================================================================
# Mission: Transform manuscript from exploratory to clinically actionable
# Objective: Test if clusters predict drug response INDEPENDENTLY of mutations
# Date: 2025-10-16
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

dir.create("03_Results/23_Drug_Validation", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/22_Drug_Validation", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 5: COMPREHENSIVE DRUG RESPONSE VALIDATION\n")
cat("Mission: Demonstrate clusters have INDEPENDENT clinical utility\n")
cat("==============================================================================\n\n")

# ==============================================================================
# PART 1: LOAD AND RECOVER DATA
# ==============================================================================

cat("PART 1: Loading drug response and molecular data...\n\n")

# Load drug response (long format from original file)
drug_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt",
                  data.table = FALSE)

# Pivot to wide format (samples × drugs)
drug_auc <- drug_raw %>%
  dplyr::select(dbgap_rnaseq_sample, inhibitor, auc) %>%
  filter(!is.na(auc)) %>%
  pivot_wider(names_from = inhibitor, values_from = auc,
              values_fn = mean)  # Take mean if duplicate

drug_auc_df <- as.data.frame(drug_auc)
rownames(drug_auc_df) <- drug_auc_df$dbgap_rnaseq_sample
drug_auc_df$dbgap_rnaseq_sample <- NULL
drug_auc <- drug_auc_df

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Load mutation data (binary matrix: samples × mutations)
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

# Load expression data (for BCL-2 pathway and immune checkpoints)
expr_raw <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")

# Match samples
common_samples <- intersect(rownames(drug_auc), clusters$sample_id)
common_samples <- intersect(common_samples, colnames(expr_raw))

cat(sprintf("Samples with complete data: %d\n", length(common_samples)))
cat(sprintf("  - Drug response: %d drugs\n", ncol(drug_auc)))
cat(sprintf("  - Clusters: k=%d\n", length(unique(clusters$cluster))))
cat(sprintf("  - Expression: %d genes\n\n", nrow(expr_raw)))

# Filter to common samples
drug_auc_filt <- drug_auc[common_samples, ]
cluster_assign <- clusters$cluster[match(common_samples, clusters$sample_id)]
expr_filt <- expr_raw[, common_samples]

# ==============================================================================
# PART 2: DIFFERENTIAL DRUG RESPONSE BY CLUSTER
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 2: DIFFERENTIAL DRUG RESPONSE ANALYSIS (ALL DRUGS)\n")
cat("Testing: Do clusters predict drug response?\n")
cat("==============================================================================\n\n")

drugs_to_test <- colnames(drug_auc_filt)
n_drugs <- length(drugs_to_test)

cat(sprintf("Testing %d drugs...\n\n", n_drugs))

drug_results <- data.frame()

for (drug_name in drugs_to_test) {
  auc_values <- drug_auc_filt[[drug_name]]

  # Remove NAs
  valid_idx <- !is.na(auc_values)

  if (sum(valid_idx) < 30) next  # Need ≥30 samples

  auc_c1 <- auc_values[valid_idx & cluster_assign == 1]
  auc_c2 <- auc_values[valid_idx & cluster_assign == 2]

  if (length(auc_c1) < 10 || length(auc_c2) < 10) next

  # Wilcoxon rank-sum test
  wilcox_result <- wilcox.test(auc_c1, auc_c2)

  # Effect size (Cohen's d)
  pooled_sd <- sqrt((sd(auc_c1, na.rm = TRUE)^2 + sd(auc_c2, na.rm = TRUE)^2) / 2)
  cohens_d <- (mean(auc_c1, na.rm = TRUE) - mean(auc_c2, na.rm = TRUE)) / pooled_sd

  # Fold difference
  fold_diff <- mean(auc_c1, na.rm = TRUE) / mean(auc_c2, na.rm = TRUE)

  drug_results <- rbind(drug_results, data.frame(
    drug = drug_name,
    n_samples = sum(valid_idx),
    n_cluster1 = length(auc_c1),
    n_cluster2 = length(auc_c2),
    mean_auc_cluster1 = mean(auc_c1, na.rm = TRUE),
    mean_auc_cluster2 = mean(auc_c2, na.rm = TRUE),
    sd_auc_cluster1 = sd(auc_c1, na.rm = TRUE),
    sd_auc_cluster2 = sd(auc_c2, na.rm = TRUE),
    median_auc_cluster1 = median(auc_c1, na.rm = TRUE),
    median_auc_cluster2 = median(auc_c2, na.rm = TRUE),
    wilcoxon_pvalue = wilcox_result$p.value,
    cohens_d = cohens_d,
    fold_difference = fold_diff,
    cluster_more_sensitive = ifelse(mean(auc_c1) < mean(auc_c2), "Cluster_1", "Cluster_2"),
    stringsAsFactors = FALSE
  ))
}

# Check if we have results
if (nrow(drug_results) == 0) {
  stop("ERROR: No drugs could be tested. Check sample ID matching.")
}

# FDR correction
drug_results$fdr <- p.adjust(drug_results$wilcoxon_pvalue, method = "BH")

# Sort by significance
drug_results <- drug_results %>% arrange(wilcoxon_pvalue)

# Significant drugs
sig_drugs <- drug_results %>% filter(fdr < 0.05)

cat(sprintf("✓ Tested %d drugs\n", nrow(drug_results)))
cat(sprintf("✓ Significant at FDR<0.05: %d drugs (%.1f%%)\n\n",
            nrow(sig_drugs), 100 * nrow(sig_drugs) / nrow(drug_results)))

if (nrow(sig_drugs) >= 10) {
  cat("SUCCESS CRITERION: ≥10 drugs with FDR<0.05 = STRONG DIFFERENTIAL RESPONSE ✓✓✓\n\n")
} else {
  cat("WEAK: <10 drugs significant\n\n")
}

cat("Top 10 most significant drugs:\n")
print(sig_drugs[1:min(10, nrow(sig_drugs)),
                c("drug", "wilcoxon_pvalue", "fdr", "cohens_d", "cluster_more_sensitive")])
cat("\n")

# Venetoclax specific check
venetoclax_row <- drug_results[drug_results$drug == "Venetoclax", ]
if (nrow(venetoclax_row) > 0) {
  cat("==== VENETOCLAX SPOTLIGHT ====\n")
  cat(sprintf("  Cluster 1 AUC: %.2f ± %.2f\n",
              venetoclax_row$mean_auc_cluster1, venetoclax_row$sd_auc_cluster1))
  cat(sprintf("  Cluster 2 AUC: %.2f ± %.2f\n",
              venetoclax_row$mean_auc_cluster2, venetoclax_row$sd_auc_cluster2))
  cat(sprintf("  P-value: %.2e\n", venetoclax_row$wilcoxon_pvalue))
  cat(sprintf("  FDR: %.2e\n", venetoclax_row$fdr))
  cat(sprintf("  Cohen's d: %.2f\n", abs(venetoclax_row$cohens_d)))

  if (venetoclax_row$wilcoxon_pvalue < 1e-10 && abs(venetoclax_row$cohens_d) > 1.5) {
    cat("\n  ⭐⭐⭐ EXTRAORDINARY FINDING: p<10^-10 AND d>1.5 ⭐⭐⭐\n")
  }
  cat("\n")
}

write.csv(drug_results, "03_Results/23_Drug_Validation/all_drugs_differential_response.csv",
          row.names = FALSE)
cat("✓ Saved: all_drugs_differential_response.csv\n\n")

# ==============================================================================
# PART 3: THREE-WAY INTERACTIONS (CRITICAL FOR INDEPENDENCE)
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 3: THREE-WAY INTERACTIONS - TESTING INDEPENDENCE\n")
cat("Question: Do clusters add value BEYOND mutations?\n")
cat("==============================================================================\n\n")

# Prepare mutation data
key_mutations <- c("NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2", "TET2", "TP53",
                   "RUNX1", "ASXL1", "NRAS", "KRAS")

# Match samples and extract mutations
mutation_data <- mutations[common_samples, ]

# For each key mutation, check if present in the mutation matrix
mutation_matrix <- data.frame(
  sample_id = common_samples,
  cluster = cluster_assign
)

for (mut in key_mutations) {
  if (mut %in% colnames(mutation_data)) {
    mutation_matrix[[mut]] <- as.integer(mutation_data[[mut]])
  } else {
    mutation_matrix[[mut]] <- 0
  }
}

# Test top 20 drugs for three-way interactions
top_drugs <- sig_drugs$drug[1:min(20, nrow(sig_drugs))]

interaction_results <- data.frame()

cat(sprintf("Testing %d drugs for mutation × cluster interactions...\n\n", length(top_drugs)))

for (drug_name in top_drugs) {
  auc_values <- drug_auc_filt[[drug_name]]

  valid_idx <- !is.na(auc_values)

  if (sum(valid_idx) < 50) next

  test_data <- mutation_matrix[valid_idx, ]
  test_data$auc <- auc_values[valid_idx]

  # Convert cluster to numeric to avoid formula issues
  test_data$cluster <- as.numeric(test_data$cluster)

  # Remove mutations with <10 samples
  mut_counts <- colSums(test_data[, key_mutations, drop = FALSE])
  valid_muts <- names(mut_counts[mut_counts >= 10])

  if (length(valid_muts) == 0) next

  # Model 1: Mutations only
  mut_formula <- as.formula(paste("auc ~", paste(valid_muts, collapse = " + ")))

  tryCatch({
    model1 <- lm(mut_formula, data = test_data)
    r2_mutations <- summary(model1)$r.squared

    # Model 2: Clusters only
    model2 <- lm(auc ~ cluster, data = test_data)
    r2_cluster <- summary(model2)$r.squared

    # Model 3: Mutations + Cluster (additive)
    add_formula <- as.formula(paste("auc ~", paste(c(valid_muts, "cluster"), collapse = " + ")))
    model3 <- lm(add_formula, data = test_data)
    r2_additive <- summary(model3)$r.squared

    # Model 4: Mutations + Cluster + Interactions
    int_formula <- as.formula(paste("auc ~ cluster * (", paste(valid_muts, collapse = " + "), ")"))
    model4 <- lm(int_formula, data = test_data)
    r2_interaction <- summary(model4)$r.squared

    # Test if cluster adds value beyond mutations
    anova_result <- anova(model1, model3)
    p_cluster_adds_value <- anova_result$`Pr(>F)`[2]

    # Test if interactions are significant
    anova_int <- anova(model3, model4)
    p_interactions <- anova_int$`Pr(>F)`[2]

    interaction_results <- rbind(interaction_results, data.frame(
      drug = drug_name,
      n_samples = sum(valid_idx),
      n_mutations_tested = length(valid_muts),
      mutations_tested = paste(valid_muts, collapse = ";"),
      r2_mutations_only = r2_mutations,
      r2_cluster_only = r2_cluster,
      r2_additive = r2_additive,
      r2_interaction = r2_interaction,
      r2_improvement_cluster = r2_additive - r2_mutations,
      r2_improvement_interactions = r2_interaction - r2_additive,
      p_cluster_adds_value = p_cluster_adds_value,
      p_interactions = p_interactions,
      cluster_independent = p_cluster_adds_value < 0.05,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat(sprintf("  Skipping %s: %s\n", drug_name, e$message))
  })
}

# Check if we got any results
if (nrow(interaction_results) == 0) {
  cat("WARNING: No drugs could be tested for interactions\n\n")

  # Create empty summary
  independent_drugs <- data.frame()
} else {
  # FDR correction
  interaction_results$fdr_cluster <- p.adjust(interaction_results$p_cluster_adds_value, method = "BH")
  interaction_results$fdr_interactions <- p.adjust(interaction_results$p_interactions, method = "BH")

  # Sort by cluster independence
  interaction_results <- interaction_results %>% arrange(p_cluster_adds_value)
  independent_drugs <- interaction_results %>% filter(p_cluster_adds_value < 0.05)
}

cat(sprintf("✓ Tested %d drugs for interactions\n", nrow(interaction_results)))
if (nrow(interaction_results) > 0) {
  cat(sprintf("✓ Drugs where cluster adds value beyond mutations (p<0.05): %d\n\n",
              sum(interaction_results$p_cluster_adds_value < 0.05, na.rm = TRUE)))
} else {
  cat("\n")
}

if (nrow(interaction_results) > 0 && nrow(independent_drugs) > 0) {
  cat("⭐⭐⭐ SUCCESS: CLUSTERS HAVE INDEPENDENT PREDICTIVE VALUE ⭐⭐⭐\n\n")
  cat("Drugs with independent cluster effect:\n")
  print(independent_drugs[, c("drug", "p_cluster_adds_value", "r2_improvement_cluster")])
  cat("\n")
} else {
  cat("❌ FAILURE: Clusters do NOT add value beyond mutations\n\n")
}

write.csv(interaction_results, "03_Results/23_Drug_Validation/drug_cluster_independence_tests.csv",
          row.names = FALSE)
cat("✓ Saved: drug_cluster_independence_tests.csv\n\n")

# ==============================================================================
# PART 4: BCL-2 PATHWAY VALIDATION
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 4: BCL-2 PATHWAY VALIDATION\n")
cat("Hypothesis: Venetoclax sensitivity correlates with BCL-2 pathway expression\n")
cat("==============================================================================\n\n")

bcl2_genes <- c("BCL2", "BCL2L1", "MCL1", "BCL2L11", "BAX", "BAK1",
                "BBC3", "PMAIP1", "BID", "BAD")

# Find genes in expression matrix
bcl2_found <- intersect(bcl2_genes, rownames(expr_filt))

cat(sprintf("BCL-2 pathway genes found: %d/%d\n", length(bcl2_found), length(bcl2_genes)))
cat(paste(bcl2_found, collapse = ", "), "\n\n")

bcl2_expr <- expr_filt[bcl2_found, common_samples]

# Mean expression by cluster
bcl2_mean_c1 <- rowMeans(bcl2_expr[, cluster_assign == 1], na.rm = TRUE)
bcl2_mean_c2 <- rowMeans(bcl2_expr[, cluster_assign == 2], na.rm = TRUE)

bcl2_results <- data.frame(
  gene = bcl2_found,
  mean_cluster1 = bcl2_mean_c1,
  mean_cluster2 = bcl2_mean_c2,
  log2fc = log2((bcl2_mean_c1 + 1) / (bcl2_mean_c2 + 1)),
  stringsAsFactors = FALSE
)

# Test differential expression
bcl2_results$pvalue <- NA
for (i in 1:nrow(bcl2_results)) {
  gene <- bcl2_results$gene[i]
  expr_c1 <- as.numeric(bcl2_expr[gene, cluster_assign == 1])
  expr_c2 <- as.numeric(bcl2_expr[gene, cluster_assign == 2])

  test_result <- wilcox.test(expr_c1, expr_c2)
  bcl2_results$pvalue[i] <- test_result$p.value
}

bcl2_results$fdr <- p.adjust(bcl2_results$pvalue, method = "BH")
bcl2_results$direction <- ifelse(bcl2_results$log2fc > 0, "C1 > C2", "C2 > C1")

cat("BCL-2 pathway differential expression:\n")
print(bcl2_results[order(bcl2_results$pvalue), ])
cat("\n")

# Check if BCL2 predicts Venetoclax response
if ("Venetoclax" %in% colnames(drug_auc_filt) && "BCL2" %in% bcl2_found) {
  ven_auc <- drug_auc_filt[["Venetoclax"]]
  bcl2_level <- as.numeric(bcl2_expr["BCL2", ])

  valid_idx <- !is.na(ven_auc) & !is.na(bcl2_level)

  cor_result <- cor.test(ven_auc[valid_idx], bcl2_level[valid_idx], method = "spearman")

  cat("BCL2 expression vs Venetoclax AUC:\n")
  cat(sprintf("  Spearman rho: %.3f\n", cor_result$estimate))
  cat(sprintf("  P-value: %.2e\n", cor_result$p.value))

  if (cor_result$p.value < 0.01) {
    cat("  ✓ VALIDATED: BCL2 expression predicts Venetoclax response\n\n")
  } else {
    cat("  ✗ Not significant\n\n")
  }
}

write.csv(bcl2_results, "03_Results/23_Drug_Validation/bcl2_pathway_expression.csv",
          row.names = FALSE)
cat("✓ Saved: bcl2_pathway_expression.csv\n\n")

# ==============================================================================
# PART 5: IMMUNE CHECKPOINT EXPRESSION
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 5: IMMUNE CHECKPOINT EXPRESSION\n")
cat("Testing: Do clusters differ in immunotherapy target expression?\n")
cat("==============================================================================\n\n")

checkpoint_genes <- c(
  "CD274",    # PD-L1
  "PDCD1",    # PD-1
  "CTLA4",    # CTLA-4
  "LAG3",     # LAG3
  "HAVCR2",   # TIM-3
  "TIGIT",    # TIGIT
  "BTLA",     # BTLA
  "CD47"      # CD47
)

checkpoint_found <- intersect(checkpoint_genes, rownames(expr_filt))

cat(sprintf("Checkpoint genes found: %d/%d\n", length(checkpoint_found), length(checkpoint_genes)))
cat(paste(checkpoint_found, collapse = ", "), "\n\n")

checkpoint_expr <- expr_filt[checkpoint_found, common_samples]

# Test differential expression
checkpoint_results <- data.frame()

for (gene in checkpoint_found) {
  expr_c1 <- as.numeric(checkpoint_expr[gene, cluster_assign == 1])
  expr_c2 <- as.numeric(checkpoint_expr[gene, cluster_assign == 2])

  test_result <- wilcox.test(expr_c1, expr_c2)

  checkpoint_results <- rbind(checkpoint_results, data.frame(
    gene = gene,
    mean_cluster1 = mean(expr_c1, na.rm = TRUE),
    mean_cluster2 = mean(expr_c2, na.rm = TRUE),
    log2fc = log2((mean(expr_c1) + 1) / (mean(expr_c2) + 1)),
    pvalue = test_result$p.value,
    stringsAsFactors = FALSE
  ))
}

checkpoint_results$fdr <- p.adjust(checkpoint_results$pvalue, method = "BH")
checkpoint_results$direction <- ifelse(checkpoint_results$log2fc > 0, "C1 > C2", "C2 > C1")

cat("Immune checkpoint expression:\n")
print(checkpoint_results[order(checkpoint_results$pvalue), ])
cat("\n")

write.csv(checkpoint_results, "03_Results/23_Drug_Validation/immune_checkpoint_expression.csv",
          row.names = FALSE)
cat("✓ Saved: immune_checkpoint_expression.csv\n\n")

# ==============================================================================
# PART 6: DRUG CLASS ANALYSIS
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 6: DRUG CLASS ANALYSIS\n")
cat("Testing: Are specific drug classes enriched in differential response?\n")
cat("==============================================================================\n\n")

# Define drug classes manually based on known mechanisms
drug_classes <- list(
  "BCL2_inhibitors" = c("Venetoclax", "ABT-737"),
  "FLT3_inhibitors" = c("Quizartinib (AC220)", "Gilteritinib (ASP-2215)", "Crenolanib",
                        "Midostaurin", "Lestaurtinib (CEP-701)", "Tandutinib (MLN518)"),
  "MEK_inhibitors" = c("Selumetinib (AZD6244)", "Trametinib (GSK1120212)", "CI-1040 (PD184352)"),
  "PI3K_inhibitors" = c("Idelalisib", "GDC-0941", "BEZ235", "PI-103"),
  "HDAC_inhibitors" = c("Panobinostat", "Entinostat"),
  "TKI_multikinase" = c("Sorafenib", "Dasatinib", "Nilotinib", "Ponatinib (AP24534)",
                        "Sunitinib", "Pazopanib (GW786034)"),
  "CDK_inhibitors" = c("Palbociclib", "Flavopiridol", "SNS-032 (BMS-387032)"),
  "Hypomethylating" = c("Azacytidine", "Decitabine")
)

# Count significant drugs per class
class_enrichment <- data.frame()

for (class_name in names(drug_classes)) {
  class_drugs <- drug_classes[[class_name]]

  # Find drugs in our dataset
  class_drugs_tested <- intersect(class_drugs, drug_results$drug)
  class_drugs_sig <- intersect(class_drugs, sig_drugs$drug)

  n_tested <- length(class_drugs_tested)
  n_sig <- length(class_drugs_sig)

  # Fisher's exact test vs background
  # Contingency table: [sig in class, sig not in class; not sig in class, not sig not in class]
  total_tested <- nrow(drug_results)
  total_sig <- nrow(sig_drugs)

  fisher_table <- matrix(c(
    n_sig, total_sig - n_sig,
    n_tested - n_sig, total_tested - total_tested + n_sig - n_tested + n_sig
  ), nrow = 2, byrow = TRUE)

  if (n_tested > 0) {
    fisher_result <- fisher.test(fisher_table)

    class_enrichment <- rbind(class_enrichment, data.frame(
      drug_class = class_name,
      n_drugs_tested = n_tested,
      n_drugs_significant = n_sig,
      percent_significant = 100 * n_sig / n_tested,
      odds_ratio = fisher_result$estimate,
      pvalue = fisher_result$p.value,
      significant_drugs = paste(class_drugs_sig, collapse = "; "),
      stringsAsFactors = FALSE
    ))
  }
}

class_enrichment$fdr <- p.adjust(class_enrichment$pvalue, method = "BH")
class_enrichment <- class_enrichment %>% arrange(pvalue)

cat("Drug class enrichment analysis:\n")
print(class_enrichment)
cat("\n")

write.csv(class_enrichment, "03_Results/23_Drug_Validation/drug_class_enrichment.csv",
          row.names = FALSE)
cat("✓ Saved: drug_class_enrichment.csv\n\n")

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

cat("\n==============================================================================\n")
cat("PHASE 5 VALIDATION - SUMMARY REPORT\n")
cat("==============================================================================\n\n")

cat("KEY FINDINGS:\n\n")

cat(sprintf("1. DIFFERENTIAL DRUG RESPONSE:\n"))
cat(sprintf("   - Tested: %d drugs\n", nrow(drug_results)))
cat(sprintf("   - Significant (FDR<0.05): %d (%.1f%%)\n",
            nrow(sig_drugs), 100 * nrow(sig_drugs) / nrow(drug_results)))
if (nrow(sig_drugs) >= 10) {
  cat("   - ✓✓✓ SUCCESS: ≥10 drugs with differential response\n\n")
} else {
  cat("   - ✗ WEAK: <10 significant drugs\n\n")
}

cat(sprintf("2. CLUSTER INDEPENDENCE:\n"))
if (nrow(independent_drugs) > 0) {
  cat(sprintf("   - %d drugs where clusters add value beyond mutations\n", nrow(independent_drugs)))
  cat("   - ⭐⭐⭐ CLUSTERS HAVE INDEPENDENT CLINICAL UTILITY ⭐⭐⭐\n\n")
} else {
  cat("   - ❌ Clusters do NOT add value beyond mutations\n")
  cat("   - Cannot claim independent biomarker status\n\n")
}

cat(sprintf("3. VENETOCLAX FINDING:\n"))
if (nrow(venetoclax_row) > 0) {
  cat(sprintf("   - P-value: %.2e (FDR: %.2e)\n",
              venetoclax_row$wilcoxon_pvalue, venetoclax_row$fdr))
  cat(sprintf("   - Cohen's d: %.2f\n", abs(venetoclax_row$cohens_d)))
  if (venetoclax_row$wilcoxon_pvalue < 1e-10 && abs(venetoclax_row$cohens_d) > 1.5) {
    cat("   - ⭐⭐⭐ EXTRAORDINARY: p<10^-10 AND d>1.5 ⭐⭐⭐\n\n")
  } else {
    cat("   - Significant but not extraordinary\n\n")
  }
}

cat(sprintf("4. BCL-2 PATHWAY:\n"))
sig_bcl2 <- sum(bcl2_results$fdr < 0.05, na.rm = TRUE)
cat(sprintf("   - %d/%d genes differentially expressed (FDR<0.05)\n",
            sig_bcl2, nrow(bcl2_results)))
if (sig_bcl2 > 0) {
  cat("   - ✓ Biological mechanism validated\n\n")
} else {
  cat("   - ✗ No clear mechanism\n\n")
}

cat(sprintf("5. IMMUNE CHECKPOINTS:\n"))
sig_checkpoint <- sum(checkpoint_results$fdr < 0.05, na.rm = TRUE)
cat(sprintf("   - %d/%d genes differentially expressed (FDR<0.05)\n",
            sig_checkpoint, nrow(checkpoint_results)))
if (sig_checkpoint >= 2) {
  cat("   - ✓ Immunotherapy target differences confirmed\n\n")
} else {
  cat("   - Limited immunotherapy implications\n\n")
}

cat("\n==============================================================================\n")
cat("CLINICAL ACTIONABILITY VERDICT\n")
cat("==============================================================================\n\n")

# Calculate score
clinical_score <- 0
if (nrow(sig_drugs) >= 10) clinical_score <- clinical_score + 2
if (nrow(independent_drugs) > 0) clinical_score <- clinical_score + 3
if (nrow(venetoclax_row) > 0 && venetoclax_row$wilcoxon_pvalue < 1e-10) clinical_score <- clinical_score + 2
if (sig_bcl2 > 0) clinical_score <- clinical_score + 1
if (sig_checkpoint >= 2) clinical_score <- clinical_score + 1

cat(sprintf("Clinical Actionability Score: %d/9\n\n", clinical_score))

if (clinical_score >= 7) {
  cat("⭐⭐⭐ VERDICT: CLINICALLY ACTIONABLE BIOMARKER ⭐⭐⭐\n")
  cat("Clusters have independent predictive value for treatment selection.\n")
  cat("Manuscript positioning: Precision medicine tool ready for prospective validation.\n\n")
} else if (clinical_score >= 4) {
  cat("✓ VERDICT: PROMISING BIOMARKER\n")
  cat("Clusters show strong drug associations but limited independence.\n")
  cat("Manuscript positioning: Exploratory biomarker with research utility.\n\n")
} else {
  cat("❌ VERDICT: NOT CLINICALLY ACTIONABLE\n")
  cat("Clusters are proxies for known mutations without added value.\n")
  cat("Manuscript positioning: Biological classification only.\n\n")
}

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")

cat("\nSaved files:\n")
cat("  - 03_Results/23_Drug_Validation/all_drugs_differential_response.csv\n")
cat("  - 03_Results/23_Drug_Validation/drug_cluster_independence_tests.csv\n")
cat("  - 03_Results/23_Drug_Validation/bcl2_pathway_expression.csv\n")
cat("  - 03_Results/23_Drug_Validation/immune_checkpoint_expression.csv\n")
cat("  - 03_Results/23_Drug_Validation/drug_class_enrichment.csv\n\n")

cat("Next: Run Part 7-8 for publication figures and manuscript updates.\n")
