#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 2.3: Molecular Subtype Characterization
# ==============================================================================
# Objective:
#   1. Differential expression analysis per cluster
#   2. Mutation enrichment analysis
#   3. Clinical variable associations
#   4. Biological naming of subtypes
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(data.table)
  library(readxl)
  library(pheatmap)
  library(ggplot2)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

dir.create("03_Results/07_Subtype_Characterization", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/04_Subtype_Characterization", recursive = TRUE, showWarnings = FALSE)

cat("==============================================================================\n")
cat("PHASE 2.3: MOLECULAR SUBTYPE CHARACTERIZATION\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load Data
# ------------------------------------------------------------------------------

cat("STEP 1: Loading data...\n\n")

# Expression data
expr_data <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")

# Cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
n_clusters <- length(unique(clusters$cluster))

# Gold standard cohort for integrated analyses
gold_data <- readRDS("03_Results/05_Analysis_Ready_Data/gold_standard_cohort.rds")

# Mutation data
mutations <- gold_data$mutations

# Clinical data
clinical <- gold_data$clinical

cat(sprintf("Expression: %d genes × %d samples\n", nrow(expr_data), ncol(expr_data)))
cat(sprintf("Clusters: %d samples in %d clusters\n", nrow(clusters), n_clusters))
cat(sprintf("Mutations: %d samples × %d genes\n", nrow(mutations), ncol(mutations)))
cat(sprintf("Clinical: %d samples × %d variables\n\n", nrow(clinical), ncol(clinical)))

# ------------------------------------------------------------------------------
# STEP 2: Differential Expression Analysis
# ------------------------------------------------------------------------------

cat("STEP 2: Performing differential expression analysis...\n\n")

# Align expression data with cluster assignments
expr_clustered <- expr_data[, clusters$sample_id]

# Design matrix for limma (each cluster vs all others)
all_deg_results <- list()

for (k in 1:n_clusters) {
  cat(sprintf("Analyzing Cluster %d vs all others...\n", k))

  # Create design: 1 = in cluster k, 0 = not in cluster k
  group <- ifelse(clusters$cluster == k, "InCluster", "Others")
  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- c("InCluster", "Others")

  # Fit model
  fit <- lmFit(expr_clustered, design)

  # Contrast: InCluster vs Others
  contrast.matrix <- makeContrasts(InCluster - Others, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)

  # Extract results
  results <- topTable(fit2, number = Inf, adjust.method = "BH")
  results$gene_id <- rownames(results)

  # Add significance flags
  results$significant <- results$adj.P.Val < 0.05 & abs(results$logFC) > 1

  all_deg_results[[paste0("Cluster_", k)]] <- results

  # Summary
  n_up <- sum(results$logFC > 1 & results$adj.P.Val < 0.05)
  n_down <- sum(results$logFC < -1 & results$adj.P.Val < 0.05)

  cat(sprintf("  Upregulated genes (logFC>1, FDR<0.05): %d\n", n_up))
  cat(sprintf("  Downregulated genes (logFC<-1, FDR<0.05): %d\n\n", n_down))
}

# Save all DEG results
saveRDS(all_deg_results,
        "03_Results/07_Subtype_Characterization/differential_expression_results.rds")
cat("✓ Saved: differential_expression_results.rds\n\n")

# Extract top genes per cluster
top_genes_per_cluster <- list()

for (k in 1:n_clusters) {
  deg <- all_deg_results[[paste0("Cluster_", k)]]

  top_up <- deg %>%
    filter(logFC > 1, adj.P.Val < 0.05) %>%
    arrange(desc(logFC)) %>%
    head(50)

  top_down <- deg %>%
    filter(logFC < -1, adj.P.Val < 0.05) %>%
    arrange(logFC) %>%
    head(50)

  top_genes_per_cluster[[paste0("Cluster_", k)]] <- list(
    upregulated = top_up,
    downregulated = top_down
  )
}

saveRDS(top_genes_per_cluster,
        "03_Results/07_Subtype_Characterization/top_genes_per_cluster.rds")
cat("✓ Saved: top_genes_per_cluster.rds\n\n")

# ------------------------------------------------------------------------------
# STEP 3: Mutation Enrichment Analysis
# ------------------------------------------------------------------------------

cat("STEP 3: Analyzing mutation enrichment by cluster...\n\n")

# Match mutation data to cluster assignments
mutation_cluster_match <- data.frame(
  sample_id = rownames(mutations),
  stringsAsFactors = FALSE
)

# Map mutation sample IDs to expression IDs via gold standard mapping
sample_mapping <- read.csv("03_Results/01_Processed_Data/master_sample_id_mapping.csv")

mutation_cluster_match <- mutation_cluster_match %>%
  left_join(sample_mapping %>%
              dplyr::select(mutation_id, expression_id),
            by = c("sample_id" = "mutation_id")) %>%
  left_join(clusters %>% dplyr::select(sample_id, cluster),
            by = c("expression_id" = "sample_id"))

# Remove samples without cluster assignment
mutation_cluster_match <- mutation_cluster_match %>%
  filter(!is.na(cluster))

cat(sprintf("Samples with both mutation and cluster data: %d\n\n",
            nrow(mutation_cluster_match)))

# Fisher's exact test for each mutation in each cluster
mutation_enrichment_results <- list()

mutation_genes <- colnames(mutations)

for (k in 1:n_clusters) {
  enrichment_df <- data.frame(
    gene = character(),
    cluster = integer(),
    freq_in_cluster = numeric(),
    freq_other = numeric(),
    OR = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
  )

  samples_in_cluster <- mutation_cluster_match$sample_id[mutation_cluster_match$cluster == k]
  samples_other <- mutation_cluster_match$sample_id[mutation_cluster_match$cluster != k]

  for (gene in mutation_genes) {
    # Contingency table
    mut_in_cluster <- sum(mutations[samples_in_cluster, gene])
    no_mut_in_cluster <- length(samples_in_cluster) - mut_in_cluster

    mut_other <- sum(mutations[samples_other, gene])
    no_mut_other <- length(samples_other) - mut_other

    if (mut_in_cluster + mut_other >= 5) {  # At least 5 total mutations
      cont_table <- matrix(c(mut_in_cluster, no_mut_in_cluster,
                             mut_other, no_mut_other), nrow = 2)

      fisher_result <- fisher.test(cont_table)

      enrichment_df <- rbind(enrichment_df, data.frame(
        gene = gene,
        cluster = k,
        n_mut_in_cluster = mut_in_cluster,
        n_in_cluster = length(samples_in_cluster),
        freq_in_cluster = mut_in_cluster / length(samples_in_cluster) * 100,
        freq_other = mut_other / length(samples_other) * 100,
        OR = as.numeric(fisher_result$estimate),
        pvalue = fisher_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Adjust p-values
  enrichment_df$adj_pvalue <- p.adjust(enrichment_df$pvalue, method = "BH")

  mutation_enrichment_results[[paste0("Cluster_", k)]] <- enrichment_df

  # Display significant enrichments
  sig_enriched <- enrichment_df %>%
    filter(adj_pvalue < 0.05, OR > 1) %>%
    arrange(adj_pvalue)

  if (nrow(sig_enriched) > 0) {
    cat(sprintf("Cluster %d - Enriched mutations (FDR<0.05):\n", k))
    print(head(sig_enriched[, c("gene", "freq_in_cluster", "freq_other", "OR", "adj_pvalue")], 5))
    cat("\n")
  }
}

saveRDS(mutation_enrichment_results,
        "03_Results/07_Subtype_Characterization/mutation_enrichment_results.rds")
cat("✓ Saved: mutation_enrichment_results.rds\n\n")

# ------------------------------------------------------------------------------
# STEP 4: Clinical Variable Associations
# ------------------------------------------------------------------------------

cat("STEP 4: Analyzing clinical variable associations...\n\n")

# Match clinical data to clusters
clinical_cluster_match <- clinical %>%
  dplyr::select(dbgap_rnaseq_sample, ageAtDiagnosis, consensus_sex) %>%
  left_join(clusters, by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  filter(!is.na(cluster))

cat(sprintf("Samples with both clinical and cluster data: %d\n\n",
            nrow(clinical_cluster_match)))

# Age by cluster
cat("Age distribution by cluster:\n")
age_by_cluster <- tapply(clinical_cluster_match$ageAtDiagnosis,
                         clinical_cluster_match$cluster, summary)
print(age_by_cluster)
cat("\n")

# ANOVA for age
age_anova <- aov(ageAtDiagnosis ~ factor(cluster), data = clinical_cluster_match)
cat("ANOVA: Age ~ Cluster\n")
print(summary(age_anova))
cat("\n")

# Sex by cluster
cat("Sex distribution by cluster:\n")
sex_by_cluster <- table(clinical_cluster_match$cluster, clinical_cluster_match$consensus_sex)
print(sex_by_cluster)
cat("\n")

# Chi-square test for sex
if (all(sex_by_cluster >= 5)) {
  sex_chisq <- chisq.test(sex_by_cluster)
  cat("Chi-square test: Sex ~ Cluster\n")
  cat(sprintf("  X-squared = %.2f, p-value = %.3f\n\n", sex_chisq$statistic, sex_chisq$p.value))
}

# Save clinical associations
clinical_summary <- data.frame(
  cluster = 1:n_clusters,
  n_samples = as.numeric(table(clinical_cluster_match$cluster)),
  mean_age = tapply(clinical_cluster_match$ageAtDiagnosis,
                    clinical_cluster_match$cluster, mean, na.rm = TRUE),
  median_age = tapply(clinical_cluster_match$ageAtDiagnosis,
                      clinical_cluster_match$cluster, median, na.rm = TRUE),
  pct_female = tapply(clinical_cluster_match$consensus_sex == "Female",
                      clinical_cluster_match$cluster, mean, na.rm = TRUE) * 100
)

write.csv(clinical_summary,
          "03_Results/07_Subtype_Characterization/clinical_summary_by_cluster.csv",
          row.names = FALSE)
cat("✓ Saved: clinical_summary_by_cluster.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 5: Biological Naming
# ------------------------------------------------------------------------------

cat("STEP 5: Proposing biological names for subtypes...\n\n")

# Load pathway profiles
pathway_profiles <- read.csv("03_Results/06_Molecular_Subtypes/pathway_profiles_by_cluster.csv",
                             row.names = 1)

subtype_names <- character(n_clusters)

for (k in 1:n_clusters) {
  # Get top pathways and mutations for this cluster
  top_pathways <- rownames(pathway_profiles)[order(pathway_profiles[, k], decreasing = TRUE)[1:3]]

  sig_mutations <- mutation_enrichment_results[[paste0("Cluster_", k)]] %>%
    filter(adj_pvalue < 0.05, OR > 2) %>%
    arrange(desc(OR)) %>%
    head(3) %>%
    pull(gene)

  # Propose name based on characteristics
  if (length(sig_mutations) > 0) {
    # Name by dominant mutation
    subtype_names[k] <- paste0("Cluster", k, "_", sig_mutations[1], "-enriched")
  } else if (grepl("IMMUNE|INFLAMMATORY", top_pathways[1])) {
    subtype_names[k] <- paste0("Cluster", k, "_Immune-enriched")
  } else if (grepl("PROLIFERATION|E2F|MYC", top_pathways[1])) {
    subtype_names[k] <- paste0("Cluster", k, "_Proliferative")
  } else if (grepl("METABOLISM|GLYCOLYSIS|OXPHOS", top_pathways[1])) {
    subtype_names[k] <- paste0("Cluster", k, "_Metabolic")
  } else {
    subtype_names[k] <- paste0("Cluster", k, "_Mixed")
  }

  cat(sprintf("Cluster %d → %s\n", k, subtype_names[k]))
  cat(sprintf("  Top pathways: %s\n", paste(gsub("HALLMARK_", "", top_pathways), collapse = ", ")))
  if (length(sig_mutations) > 0) {
    cat(sprintf("  Enriched mutations: %s\n", paste(sig_mutations, collapse = ", ")))
  }
  cat("\n")
}

# Save subtype names
subtype_naming <- data.frame(
  cluster_number = 1:n_clusters,
  proposed_name = subtype_names,
  stringsAsFactors = FALSE
)

write.csv(subtype_naming,
          "03_Results/07_Subtype_Characterization/subtype_naming.csv",
          row.names = FALSE)
cat("✓ Saved: subtype_naming.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 6: Create Summary Figures
# ------------------------------------------------------------------------------

cat("STEP 6: Creating summary figures...\n\n")

# Figure 1: Cluster size bar plot
pdf("04_Figures/04_Subtype_Characterization/cluster_sizes.pdf",
    width = 8, height = 5)

cluster_sizes <- table(clusters$cluster)
barplot(cluster_sizes,
        main = "Molecular Subtype Sample Sizes",
        xlab = "Cluster",
        ylab = "Number of Samples",
        col = rainbow(n_clusters),
        las = 1)

dev.off()
cat("✓ Saved: cluster_sizes.pdf\n\n")

# ------------------------------------------------------------------------------
# Summary Report
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("SUBTYPE CHARACTERIZATION COMPLETE\n")
cat("==============================================================================\n\n")

cat("SUMMARY:\n")
cat(sprintf("  - Number of molecular subtypes: %d\n", n_clusters))
cat(sprintf("  - Samples characterized: %d\n", nrow(clusters)))
cat(sprintf("  - DEGs identified per cluster: Range %d-%d\n",
            min(sapply(all_deg_results, function(x) sum(x$significant))),
            max(sapply(all_deg_results, function(x) sum(x$significant)))))
cat("\n")

cat("SUBTYPE NAMES:\n")
for (i in 1:n_clusters) {
  cat(sprintf("  %s\n", subtype_names[i]))
}
cat("\n")

cat("OUTPUT FILES:\n")
cat("  - 03_Results/07_Subtype_Characterization/differential_expression_results.rds\n")
cat("  - 03_Results/07_Subtype_Characterization/mutation_enrichment_results.rds\n")
cat("  - 03_Results/07_Subtype_Characterization/clinical_summary_by_cluster.csv\n")
cat("  - 03_Results/07_Subtype_Characterization/subtype_naming.csv\n\n")

cat("NEXT STEPS:\n")
cat("  → Phase 3: Survival analysis by molecular subtype\n")
cat("  → Phase 4: Drug response integration\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
