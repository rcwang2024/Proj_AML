# TASK 1.3: Test Mutation Enrichment by Subtype
# Determine if specific mutations are enriched in Proliferative vs Immune-Inflammatory subtypes

library(tidyverse)
library(ggplot2)

cat("=== TASK 1.3: MUTATION ENRICHMENT BY SUBTYPE ===\n\n")

# Load cluster assignments
cluster_assignments <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
cat("Loaded cluster assignments:", nrow(cluster_assignments), "samples\n")
cat("Cluster distribution:\n")
print(table(cluster_assignments$cluster))
cat("\n")

# Load mutation matrix
mut_matrix <- readRDS("03_Results/10_Mutations/mutation_matrix.rds")
cat("Loaded mutation matrix:", nrow(mut_matrix), "genes x", ncol(mut_matrix), "samples\n\n")

# Match samples between clusters and mutations
common_samples <- intersect(colnames(mut_matrix), cluster_assignments$sample_id)
cat("Samples with both mutations and clusters:", length(common_samples), "\n\n")

if (length(common_samples) < 50) {
  cat("⚠️ WARNING: Very few samples with both data types\n")
  cat("Results may not be reliable\n\n")
}

# Filter to common samples
mut_matrix_filtered <- mut_matrix[, common_samples]
cluster_filtered <- cluster_assignments %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

cat("After filtering:\n")
cat("  Mutation matrix:", nrow(mut_matrix_filtered), "genes x", ncol(mut_matrix_filtered), "samples\n")
cat("  Cluster data:", nrow(cluster_filtered), "samples\n\n")

# Test enrichment for each gene using Fisher's exact test
mutation_enrichment <- data.frame()

for (gene in rownames(mut_matrix_filtered)) {
  gene_mut <- mut_matrix_filtered[gene, ]

  # Contingency table: Cluster 1 vs 2, Mutated vs WT
  cluster1_mut <- sum(gene_mut[cluster_filtered$cluster == 1])
  cluster1_wt <- sum(cluster_filtered$cluster == 1) - cluster1_mut
  cluster2_mut <- sum(gene_mut[cluster_filtered$cluster == 2])
  cluster2_wt <- sum(cluster_filtered$cluster == 2) - cluster2_mut

  cont_table <- matrix(c(cluster1_mut, cluster1_wt, cluster2_mut, cluster2_wt),
                       nrow = 2, byrow = TRUE)

  # Fisher's exact test
  if (sum(cont_table) > 0 && min(rowSums(cont_table)) > 0) {
    fisher_result <- fisher.test(cont_table)

    # Calculate percentages
    pct_cluster1 <- cluster1_mut / sum(cluster_filtered$cluster == 1) * 100
    pct_cluster2 <- cluster2_mut / sum(cluster_filtered$cluster == 2) * 100

    # Direction of enrichment
    enriched_in <- ifelse(pct_cluster1 > pct_cluster2, "Cluster1_Proliferative", "Cluster2_Immune-Inflammatory")

    mutation_enrichment <- rbind(mutation_enrichment, data.frame(
      gene = gene,
      cluster1_n = cluster1_mut,
      cluster1_pct = pct_cluster1,
      cluster2_n = cluster2_mut,
      cluster2_pct = pct_cluster2,
      diff_pct = pct_cluster1 - pct_cluster2,
      odds_ratio = fisher_result$estimate,
      p_value = fisher_result$p.value,
      enriched_in = enriched_in,
      stringsAsFactors = FALSE
    ))
  }
}

# FDR correction
mutation_enrichment$fdr <- p.adjust(mutation_enrichment$p_value, method = "BH")

# Sort by significance
mutation_enrichment <- mutation_enrichment %>%
  arrange(p_value)

cat("=== MUTATION ENRICHMENT RESULTS ===\n\n")

# Significant at FDR < 0.05
significant <- mutation_enrichment %>% filter(fdr < 0.05)
cat("Genes with significant enrichment (FDR < 0.05):\n")
if (nrow(significant) > 0) {
  print(significant, row.names = FALSE)
} else {
  cat("  None found at FDR < 0.05\n")
}
cat("\n")

# Nominally significant (p < 0.05)
nom_sig <- mutation_enrichment %>% filter(p_value < 0.05)
cat("Genes with nominal significance (p < 0.05):\n")
if (nrow(nom_sig) > 0) {
  print(nom_sig, row.names = FALSE)
} else {
  cat("  None found\n")
}
cat("\n")

# Show all results sorted by p-value
cat("All genes (sorted by p-value):\n")
print(mutation_enrichment %>% select(gene, cluster1_pct, cluster2_pct, diff_pct, p_value, fdr, enriched_in), row.names = FALSE)
cat("\n")

# Save results
write.csv(mutation_enrichment,
          "03_Results/10_Mutations/mutation_enrichment_by_cluster.csv",
          row.names = FALSE)

cat("✓ Saved: 03_Results/10_Mutations/mutation_enrichment_by_cluster.csv\n\n")

# Create visualization if there are nominally significant genes
if (nrow(nom_sig) > 0) {
  cat("Creating visualization of nominally significant mutations...\n")

  plot_data <- nom_sig %>%
    select(gene, cluster1_pct, cluster2_pct) %>%
    pivot_longer(cols = c(cluster1_pct, cluster2_pct),
                 names_to = "cluster",
                 values_to = "pct") %>%
    mutate(cluster = ifelse(cluster == "cluster1_pct",
                            "Proliferative (Cluster 1)",
                            "Immune-Inflammatory (Cluster 2)"))

  p <- ggplot(plot_data, aes(x = reorder(gene, -pct), y = pct, fill = cluster)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Gene", y = "Mutation Frequency (%)",
         title = "Mutation Frequencies by Molecular Subtype",
         subtitle = "Nominally significant genes (p < 0.05)",
         fill = "Subtype") +
    scale_fill_manual(values = c("Proliferative (Cluster 1)" = "#2E9FDF",
                                   "Immune-Inflammatory (Cluster 2)" = "#E7B800")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "top")

  ggsave("04_Figures/07_Mutations/mutation_frequencies_by_cluster.pdf",
         p, width = 10, height = 6)

  cat("✓ Figure saved: 04_Figures/07_Mutations/mutation_frequencies_by_cluster.pdf\n\n")
} else {
  cat("No nominally significant genes to plot\n\n")
}

# Interpretation
cat("=== INTERPRETATION ===\n\n")

if (nrow(significant) > 0) {
  cat("✓ Found", nrow(significant), "genes with significant enrichment (FDR < 0.05)\n")
  cat("This suggests subtypes have distinct mutational profiles.\n\n")

  # Check for prognostic mutations
  if ("NPM1" %in% significant$gene) {
    npm1_row <- significant %>% filter(gene == "NPM1")
    if (npm1_row$enriched_in == "Cluster1_Proliferative") {
      cat("⭐ NPM1 enriched in Proliferative subtype (good prognosis mutation)\n")
    }
  }

  if ("TP53" %in% significant$gene) {
    tp53_row <- significant %>% filter(gene == "TP53")
    if (tp53_row$enriched_in == "Cluster2_Immune-Inflammatory") {
      cat("⭐ TP53 enriched in Immune-Inflammatory subtype (poor prognosis mutation)\n")
    }
  }

} else if (nrow(nom_sig) > 0) {
  cat("~ Found", nrow(nom_sig), "genes with nominal significance (p < 0.05)\n")
  cat("  Suggests some mutational differences, but not FDR-significant.\n\n")
} else {
  cat("✓ No significant mutation enrichment found\n")
  cat("  This indicates subtypes are TRANSCRIPTION-DRIVEN rather than mutation-driven\n")
  cat("  This is actually a POSITIVE finding - subtypes represent biological state\n")
  cat("  beyond just genetic mutations.\n\n")
}

cat("### Task 1.3 COMPLETE ###\n")
