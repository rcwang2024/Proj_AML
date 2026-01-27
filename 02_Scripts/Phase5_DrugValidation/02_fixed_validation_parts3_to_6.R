#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 5: FIXED Parts 3-6
# ==============================================================================
# Fixes:
#   - Part 3: Sanitize drug names for formula compatibility
#   - Parts 4-6: Proper gene symbol mapping from Ensembl IDs
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("PHASE 5: FIXED ANALYSIS - PARTS 3-6\n")
cat("==============================================================================\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n\n")

# Load drug results from Part 2
drug_results <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")

# Load drug response
drug_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt", data.table = FALSE)
drug_auc <- drug_raw %>%
  dplyr::select(dbgap_rnaseq_sample, inhibitor, auc) %>%
  filter(!is.na(auc)) %>%
  pivot_wider(names_from = inhibitor, values_from = auc, values_fn = mean)
drug_auc_df <- as.data.frame(drug_auc)
rownames(drug_auc_df) <- drug_auc_df$dbgap_rnaseq_sample
drug_auc_df$dbgap_rnaseq_sample <- NULL
drug_auc <- drug_auc_df

# Load clusters
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Load mutations
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

# Load expression data (Ensembl IDs)
expr_raw <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")

# Load gene annotations (Ensembl ID to gene symbol mapping)
gene_annot <- read.csv("03_Results/01_Processed_Data/gene_annotations.csv")

# Match samples
common_samples <- intersect(rownames(drug_auc), clusters$sample_id)
common_samples <- intersect(common_samples, colnames(expr_raw))

cluster_assign <- clusters$cluster[match(common_samples, clusters$sample_id)]
drug_auc_filt <- drug_auc[common_samples, ]
expr_filt <- expr_raw[, common_samples]

cat(sprintf("Samples with complete data: %d\n\n", length(common_samples)))

# ==============================================================================
# PART 3: THREE-WAY INTERACTIONS (FIXED)
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 3: THREE-WAY INTERACTIONS - TESTING INDEPENDENCE (FIXED)\n")
cat("Question: Do clusters add value BEYOND mutations?\n")
cat("==============================================================================\n\n")

# Prepare mutation data
key_mutations <- c("NPM1", "FLT3", "DNMT3A", "IDH1", "IDH2", "TET2", "TP53",
                   "RUNX1", "ASXL1", "NRAS", "KRAS")

mutation_data <- mutations[common_samples, ]
mutation_matrix <- data.frame(
  sample_id = common_samples,
  cluster = as.numeric(cluster_assign)
)

for (mut in key_mutations) {
  if (mut %in% colnames(mutation_data)) {
    mutation_matrix[[mut]] <- as.integer(mutation_data[[mut]])
  } else {
    mutation_matrix[[mut]] <- 0
  }
}

# Get top 20 significant drugs
sig_drugs <- drug_results %>% filter(fdr < 0.05) %>% arrange(wilcoxon_pvalue)
top_drugs <- sig_drugs$drug[1:min(20, nrow(sig_drugs))]

cat(sprintf("Testing %d drugs for mutation × cluster interactions...\n\n", length(top_drugs)))

interaction_results <- data.frame()

for (drug_name in top_drugs) {
  # Create a clean variable name for the drug
  drug_var <- paste0("drug_", gsub("[^A-Za-z0-9]", "_", drug_name))

  auc_values <- drug_auc_filt[[drug_name]]
  valid_idx <- !is.na(auc_values)

  if (sum(valid_idx) < 50) next

  test_data <- mutation_matrix[valid_idx, ]
  test_data[[drug_var]] <- auc_values[valid_idx]

  # Remove mutations with <10 samples
  mut_counts <- colSums(test_data[, key_mutations, drop = FALSE])
  valid_muts <- names(mut_counts[mut_counts >= 10])

  if (length(valid_muts) == 0) next

  tryCatch({
    # Model 1: Mutations only
    mut_formula <- as.formula(paste(drug_var, "~", paste(valid_muts, collapse = " + ")))
    model1 <- lm(mut_formula, data = test_data)
    r2_mutations <- summary(model1)$r.squared

    # Model 2: Clusters only
    cluster_formula <- as.formula(paste(drug_var, "~ cluster"))
    model2 <- lm(cluster_formula, data = test_data)
    r2_cluster <- summary(model2)$r.squared

    # Model 3: Mutations + Cluster (additive)
    add_formula <- as.formula(paste(drug_var, "~", paste(c(valid_muts, "cluster"), collapse = " + ")))
    model3 <- lm(add_formula, data = test_data)
    r2_additive <- summary(model3)$r.squared

    # Model 4: Mutations + Cluster + Interactions
    int_formula <- as.formula(paste(drug_var, "~ cluster * (", paste(valid_muts, collapse = " + "), ")"))
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
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  ✓ %s: R² mutations=%.3f, +cluster=%.3f (p=%.3f)\n",
                drug_name, r2_mutations, r2_additive, p_cluster_adds_value))

  }, error = function(e) {
    cat(sprintf("  ✗ Skipping %s: %s\n", drug_name, e$message))
  })
}

cat("\n")

if (nrow(interaction_results) > 0) {
  # FDR correction
  interaction_results$fdr_cluster <- p.adjust(interaction_results$p_cluster_adds_value, method = "BH")
  interaction_results$fdr_interactions <- p.adjust(interaction_results$p_interactions, method = "BH")

  # Sort by cluster independence
  interaction_results <- interaction_results %>% arrange(p_cluster_adds_value)

  independent_drugs <- interaction_results %>% filter(p_cluster_adds_value < 0.05)

  cat(sprintf("✓ Tested %d drugs for interactions\n", nrow(interaction_results)))
  cat(sprintf("✓ Drugs where cluster adds value beyond mutations (p<0.05): %d\n\n",
              nrow(independent_drugs)))

  if (nrow(independent_drugs) > 0) {
    cat("⭐⭐⭐ SUCCESS: CLUSTERS HAVE INDEPENDENT PREDICTIVE VALUE ⭐⭐⭐\n\n")
    cat("Drugs with independent cluster effect:\n")
    print(independent_drugs[, c("drug", "p_cluster_adds_value", "fdr_cluster",
                                "r2_improvement_cluster", "r2_mutations_only", "r2_additive")])
    cat("\n")
  } else {
    cat("❌ RESULT: Clusters do NOT add significant value beyond mutations\n")
    cat("   (All p-values > 0.05 for cluster contribution)\n\n")
  }

} else {
  cat("WARNING: No drugs could be tested for interactions\n\n")
  interaction_results <- data.frame()
}

write.csv(interaction_results, "03_Results/23_Drug_Validation/drug_cluster_independence_tests_FIXED.csv",
          row.names = FALSE)
cat("✓ Saved: drug_cluster_independence_tests_FIXED.csv\n\n")

# ==============================================================================
# PART 4: BCL-2 PATHWAY VALIDATION (FIXED)
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 4: BCL-2 PATHWAY VALIDATION (FIXED)\n")
cat("Hypothesis: Venetoclax sensitivity correlates with BCL-2 pathway expression\n")
cat("==============================================================================\n\n")

bcl2_genes <- c("BCL2", "BCL2L1", "MCL1", "BCL2L11", "BAX", "BAK1",
                "BBC3", "PMAIP1", "BID", "BAD")

# Map gene symbols to Ensembl IDs
bcl2_ensembl <- gene_annot %>%
  filter(gene_symbol %in% bcl2_genes) %>%
  pull(ensembl_id)

bcl2_found <- intersect(bcl2_ensembl, rownames(expr_filt))

# Get corresponding gene symbols
bcl2_genes_found <- gene_annot %>%
  filter(ensembl_id %in% bcl2_found) %>%
  pull(gene_symbol)

cat(sprintf("BCL-2 pathway genes found: %d/%d\n", length(bcl2_found), length(bcl2_genes)))
cat(paste(bcl2_genes_found, collapse = ", "), "\n\n")

if (length(bcl2_found) > 0) {
  bcl2_expr <- expr_filt[bcl2_found, common_samples]

  # Create results with gene symbols
  bcl2_results <- data.frame()

  for (i in 1:length(bcl2_found)) {
    ensembl_id <- bcl2_found[i]
    gene_symbol <- gene_annot$gene_symbol[gene_annot$ensembl_id == ensembl_id]

    expr_c1 <- as.numeric(bcl2_expr[ensembl_id, cluster_assign == 1])
    expr_c2 <- as.numeric(bcl2_expr[ensembl_id, cluster_assign == 2])

    test_result <- wilcox.test(expr_c1, expr_c2)

    bcl2_results <- rbind(bcl2_results, data.frame(
      gene = gene_symbol,
      ensembl_id = ensembl_id,
      mean_cluster1 = mean(expr_c1, na.rm = TRUE),
      mean_cluster2 = mean(expr_c2, na.rm = TRUE),
      log2fc = log2((mean(expr_c1) + 1) / (mean(expr_c2) + 1)),
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }

  bcl2_results$fdr <- p.adjust(bcl2_results$pvalue, method = "BH")
  bcl2_results$direction <- ifelse(bcl2_results$log2fc > 0, "C1 > C2", "C2 > C1")
  bcl2_results <- bcl2_results %>% arrange(pvalue)

  cat("BCL-2 pathway differential expression:\n")
  print(bcl2_results)
  cat("\n")

  # Check if BCL2 predicts Venetoclax response
  if ("Venetoclax" %in% colnames(drug_auc_filt) && "BCL2" %in% bcl2_results$gene) {
    bcl2_ensembl_id <- bcl2_results$ensembl_id[bcl2_results$gene == "BCL2"]
    ven_auc <- drug_auc_filt[["Venetoclax"]]
    bcl2_level <- as.numeric(bcl2_expr[bcl2_ensembl_id, ])

    valid_idx <- !is.na(ven_auc) & !is.na(bcl2_level)

    cor_result <- cor.test(ven_auc[valid_idx], bcl2_level[valid_idx], method = "spearman")

    cat("BCL2 expression vs Venetoclax AUC:\n")
    cat(sprintf("  Spearman rho: %.3f\n", cor_result$estimate))
    cat(sprintf("  P-value: %.2e\n", cor_result$p.value))

    if (cor_result$p.value < 0.01) {
      if (cor_result$estimate < 0) {
        cat("  ✓ VALIDATED: Higher BCL2 expression → Lower AUC (more sensitive) ⭐\n\n")
      } else {
        cat("  ✓ Significant correlation (opposite direction expected)\n\n")
      }
    } else {
      cat("  ✗ Not significant\n\n")
    }
  }

  write.csv(bcl2_results, "03_Results/23_Drug_Validation/bcl2_pathway_expression_FIXED.csv",
            row.names = FALSE)
  cat("✓ Saved: bcl2_pathway_expression_FIXED.csv\n\n")

} else {
  cat("ERROR: No BCL-2 pathway genes found in expression data\n\n")
}

# ==============================================================================
# PART 5: IMMUNE CHECKPOINT EXPRESSION (FIXED)
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 5: IMMUNE CHECKPOINT EXPRESSION (FIXED)\n")
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

# Map to Ensembl IDs
checkpoint_ensembl <- gene_annot %>%
  filter(gene_symbol %in% checkpoint_genes) %>%
  pull(ensembl_id)

checkpoint_found <- intersect(checkpoint_ensembl, rownames(expr_filt))

checkpoint_genes_found <- gene_annot %>%
  filter(ensembl_id %in% checkpoint_found) %>%
  pull(gene_symbol)

cat(sprintf("Checkpoint genes found: %d/%d\n", length(checkpoint_found), length(checkpoint_genes)))
cat(paste(checkpoint_genes_found, collapse = ", "), "\n\n")

if (length(checkpoint_found) > 0) {
  checkpoint_expr <- expr_filt[checkpoint_found, common_samples]

  checkpoint_results <- data.frame()

  for (i in 1:length(checkpoint_found)) {
    ensembl_id <- checkpoint_found[i]
    gene_symbol <- gene_annot$gene_symbol[gene_annot$ensembl_id == ensembl_id]

    expr_c1 <- as.numeric(checkpoint_expr[ensembl_id, cluster_assign == 1])
    expr_c2 <- as.numeric(checkpoint_expr[ensembl_id, cluster_assign == 2])

    test_result <- wilcox.test(expr_c1, expr_c2)

    checkpoint_results <- rbind(checkpoint_results, data.frame(
      gene = gene_symbol,
      ensembl_id = ensembl_id,
      mean_cluster1 = mean(expr_c1, na.rm = TRUE),
      mean_cluster2 = mean(expr_c2, na.rm = TRUE),
      log2fc = log2((mean(expr_c1) + 1) / (mean(expr_c2) + 1)),
      pvalue = test_result$p.value,
      stringsAsFactors = FALSE
    ))
  }

  checkpoint_results$fdr <- p.adjust(checkpoint_results$pvalue, method = "BH")
  checkpoint_results$direction <- ifelse(checkpoint_results$log2fc > 0, "C1 > C2", "C2 > C1")
  checkpoint_results <- checkpoint_results %>% arrange(pvalue)

  cat("Immune checkpoint expression:\n")
  print(checkpoint_results)
  cat("\n")

  write.csv(checkpoint_results, "03_Results/23_Drug_Validation/immune_checkpoint_expression_FIXED.csv",
            row.names = FALSE)
  cat("✓ Saved: immune_checkpoint_expression_FIXED.csv\n\n")

} else {
  cat("ERROR: No immune checkpoint genes found in expression data\n\n")
}

# ==============================================================================
# PART 6: DRUG CLASS ANALYSIS (FIXED)
# ==============================================================================

cat("\n==============================================================================\n")
cat("PART 6: DRUG CLASS ANALYSIS\n")
cat("Testing: Are specific drug classes enriched in differential response?\n")
cat("==============================================================================\n\n")

drug_classes <- list(
  "BCL2_inhibitors" = c("Venetoclax", "ABT-737"),
  "FLT3_inhibitors" = c("Quizartinib (AC220)", "Gilteritinib", "Crenolanib",
                        "Midostaurin", "Lestaurtinib (CEP-701)", "Tandutinib (MLN518)"),
  "MEK_inhibitors" = c("Selumetinib (AZD6244)", "Trametinib (GSK1120212)", "CI-1040 (PD184352)"),
  "PI3K_inhibitors" = c("Idelalisib", "GDC-0941", "BEZ235", "PI-103"),
  "HDAC_inhibitors" = c("Panobinostat"),
  "TKI_multikinase" = c("Sorafenib", "Dasatinib", "Nilotinib", "Ponatinib (AP24534)",
                        "Sunitinib", "Pazopanib (GW786034)"),
  "CDK_inhibitors" = c("Palbociclib", "Flavopiridol", "SNS-032 (BMS-387032)"),
  "AKT_inhibitors" = c("MK-2206", "GSK690693"),
  "mTOR_inhibitors" = c("Rapamycin", "INK-128")
)

sig_drugs <- drug_results %>% filter(fdr < 0.05)
class_enrichment <- data.frame()

for (class_name in names(drug_classes)) {
  class_drugs <- drug_classes[[class_name]]

  class_drugs_tested <- intersect(class_drugs, drug_results$drug)
  class_drugs_sig <- intersect(class_drugs, sig_drugs$drug)

  n_tested <- length(class_drugs_tested)
  n_sig <- length(class_drugs_sig)

  total_tested <- nrow(drug_results)
  total_sig <- nrow(sig_drugs)

  if (n_tested > 0) {
    # Fisher's exact test
    fisher_table <- matrix(c(
      n_sig, total_sig - n_sig,
      n_tested - n_sig, total_tested - total_tested + n_sig - n_tested + n_sig
    ), nrow = 2, byrow = TRUE)

    fisher_result <- tryCatch({
      fisher.test(fisher_table)
    }, error = function(e) {
      list(estimate = NA, p.value = 1)
    })

    class_enrichment <- rbind(class_enrichment, data.frame(
      drug_class = class_name,
      n_drugs_tested = n_tested,
      n_drugs_significant = n_sig,
      percent_significant = 100 * n_sig / n_tested,
      odds_ratio = ifelse(is.na(fisher_result$estimate), NA, as.numeric(fisher_result$estimate)),
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

write.csv(class_enrichment, "03_Results/23_Drug_Validation/drug_class_enrichment_FIXED.csv",
          row.names = FALSE)
cat("✓ Saved: drug_class_enrichment_FIXED.csv\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("FIXED ANALYSIS COMPLETE - PARTS 3-6\n")
cat("==============================================================================\n\n")

cat("RESULTS SUMMARY:\n\n")

cat(sprintf("PART 3 - Three-way interactions:\n"))
if (nrow(interaction_results) > 0) {
  n_independent <- sum(interaction_results$p_cluster_adds_value < 0.05, na.rm = TRUE)
  cat(sprintf("  - Tested: %d drugs\n", nrow(interaction_results)))
  cat(sprintf("  - Independent cluster effect (p<0.05): %d drugs\n", n_independent))
  if (n_independent > 0) {
    cat("  - ⭐ CLUSTERS HAVE INDEPENDENT PREDICTIVE VALUE\n\n")
  } else {
    cat("  - ❌ Clusters do NOT add value beyond mutations\n\n")
  }
} else {
  cat("  - ✗ No drugs tested successfully\n\n")
}

if (exists("bcl2_results") && nrow(bcl2_results) > 0) {
  n_sig_bcl2 <- sum(bcl2_results$fdr < 0.05, na.rm = TRUE)
  cat(sprintf("PART 4 - BCL-2 pathway:\n"))
  cat(sprintf("  - Genes found: %d/%d\n", nrow(bcl2_results), length(bcl2_genes)))
  cat(sprintf("  - Differentially expressed (FDR<0.05): %d\n", n_sig_bcl2))
  if (n_sig_bcl2 > 0) {
    cat("  - ✓ Biological mechanism validated\n\n")
  } else {
    cat("  - ✗ No significant differences\n\n")
  }
}

if (exists("checkpoint_results") && nrow(checkpoint_results) > 0) {
  n_sig_checkpoint <- sum(checkpoint_results$fdr < 0.05, na.rm = TRUE)
  cat(sprintf("PART 5 - Immune checkpoints:\n"))
  cat(sprintf("  - Genes found: %d/%d\n", nrow(checkpoint_results), length(checkpoint_genes)))
  cat(sprintf("  - Differentially expressed (FDR<0.05): %d\n", n_sig_checkpoint))
  if (n_sig_checkpoint >= 2) {
    cat("  - ✓ Immunotherapy target differences confirmed\n\n")
  } else {
    cat("  - ✗ Limited immunotherapy implications\n\n")
  }
}

cat(sprintf("PART 6 - Drug class enrichment:\n"))
cat(sprintf("  - Classes tested: %d\n", nrow(class_enrichment)))
n_enriched <- sum(class_enrichment$pvalue < 0.05, na.rm = TRUE)
cat(sprintf("  - Enriched classes (p<0.05): %d\n", n_enriched))
if (n_enriched > 0) {
  cat("  - Top enriched classes:\n")
  top_classes <- class_enrichment %>% filter(pvalue < 0.05) %>% head(5)
  for (i in 1:nrow(top_classes)) {
    cat(sprintf("    • %s: %d/%d drugs (%.1f%%, p=%.3f)\n",
                top_classes$drug_class[i],
                top_classes$n_drugs_significant[i],
                top_classes$n_drugs_tested[i],
                top_classes$percent_significant[i],
                top_classes$pvalue[i]))
  }
  cat("\n")
}

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
