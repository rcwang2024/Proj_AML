# TASK: Immune Cell Deconvolution Analysis (Simplified)
# Validate "immune-inflammatory" label using MCP-counter method
# Note: Other methods (EPIC, quanTIseq, xCell) failed due to gene name mismatches

library(tidyverse)
library(immunedeconv)
library(ggplot2)
library(pheatmap)

cat("=== IMMUNE CELL DECONVOLUTION ANALYSIS (SIMPLIFIED) ===\n\n")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

cat("Loading expression data...\n")

# Load batch-corrected expression data
expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
cat("Expression data:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")

# Convert Ensembl IDs to gene symbols
cat("Loading gene annotations...\n")
original_expr <- read.delim("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
                            stringsAsFactors = FALSE, sep = "\t")
gene_annotations <- original_expr %>%
  select(Gene_ID = stable_id, Gene_Symbol = display_label) %>%
  distinct()

# Match and convert
ensembl_ids <- rownames(expr_data)
matched_symbols <- gene_annotations$Gene_Symbol[match(ensembl_ids, gene_annotations$Gene_ID)]

has_symbol <- !is.na(matched_symbols) & matched_symbols != ""
expr_data_symbol <- expr_data[has_symbol, ]
rownames(expr_data_symbol) <- matched_symbols[has_symbol]

# Aggregate duplicates
expr_aggregated <- aggregate(expr_data_symbol,
                             by = list(gene = rownames(expr_data_symbol)),
                             FUN = mean)
rownames(expr_aggregated) <- expr_aggregated$gene
expr_aggregated$gene <- NULL
expr_data <- as.matrix(expr_aggregated)

cat("After gene symbol conversion:", nrow(expr_data), "genes\n")

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Match samples
common_samples <- intersect(colnames(expr_data), clusters$sample_id)
expr_matched <- expr_data[, common_samples]
cluster_matched <- clusters %>% filter(sample_id %in% common_samples)

cat("Matched samples:", length(common_samples), "\n")
cat("Cluster 1 (Proliferative):", sum(cluster_matched$cluster == 1), "\n")
cat("Cluster 2 (Immune-Inflammatory):", sum(cluster_matched$cluster == 2), "\n\n")

# ============================================================================
# 2. PREPARE DATA FOR DECONVOLUTION
# ============================================================================

cat("Preparing data for deconvolution...\n")

# Convert to linear scale if needed
if (max(expr_matched, na.rm = TRUE) < 30) {
  cat("Data appears log-transformed, converting to linear scale...\n")
  expr_linear <- 2^expr_matched - 1
  expr_linear[expr_linear < 0] <- 0
} else {
  cat("Data appears to be in linear scale\n")
  expr_linear <- expr_matched
}

cat("Expression range after conversion:", round(min(expr_linear), 2),
    "to", round(max(expr_linear), 2), "\n\n")

# ============================================================================
# 3. RUN MCP-COUNTER DECONVOLUTION
# ============================================================================

cat("=== RUNNING MCP-COUNTER DECONVOLUTION ===\n\n")

dir.create("03_Results/16_Immune_Deconvolution", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/15_Immune_Deconvolution", showWarnings = FALSE, recursive = TRUE)

cat("Running MCP-counter...\n")
mcp_res <- deconvolute(expr_linear, "mcp_counter")
cat("✓ MCP-counter completed\n")
cat("Result dimensions:", nrow(mcp_res), "cell types x", ncol(mcp_res), "columns\n\n")

# Inspect structure
cat("MCP-counter result structure:\n")
print(str(mcp_res))
cat("\n")
cat("First few rows:\n")
print(head(mcp_res))
cat("\n")

# ============================================================================
# 4. PROCESS MCP-COUNTER RESULTS
# ============================================================================

cat("=== PROCESSING MCP-COUNTER RESULTS ===\n\n")

# MCP-counter returns a data frame with cell_type column
# Convert to long format
result_long <- mcp_res %>%
  pivot_longer(cols = -cell_type, names_to = "sample_id", values_to = "score") %>%
  left_join(cluster_matched %>% select(sample_id, cluster), by = "sample_id")

cat("Long format data:", nrow(result_long), "rows\n")
cat("Unique cell types:", length(unique(result_long$cell_type)), "\n\n")

# Calculate mean scores by cluster
cluster_summary <- result_long %>%
  group_by(cell_type, cluster) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  )

cat("Mean immune cell scores by cluster:\n")
print(cluster_summary, n = 20)
cat("\n")

# ============================================================================
# 5. STATISTICAL TESTING
# ============================================================================

cat("=== TESTING FOR DIFFERENTIAL IMMUNE CELL ENRICHMENT ===\n\n")

diff_test_results <- data.frame()

for (cell in unique(result_long$cell_type)) {
  cluster1_scores <- result_long %>%
    filter(cell_type == cell, cluster == 1) %>%
    pull(score)

  cluster2_scores <- result_long %>%
    filter(cell_type == cell, cluster == 2) %>%
    pull(score)

  if (length(cluster1_scores) > 0 && length(cluster2_scores) > 0) {
    wilcox_test <- wilcox.test(cluster1_scores, cluster2_scores)

    diff_test_results <- rbind(diff_test_results, data.frame(
      cell_type = cell,
      mean_cluster1 = mean(cluster1_scores, na.rm = TRUE),
      mean_cluster2 = mean(cluster2_scores, na.rm = TRUE),
      diff = mean(cluster2_scores, na.rm = TRUE) - mean(cluster1_scores, na.rm = TRUE),
      p_value = wilcox_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# FDR correction
diff_test_results$fdr <- p.adjust(diff_test_results$p_value, method = "BH")

# Sort by FDR
diff_test_results <- diff_test_results %>%
  arrange(fdr) %>%
  mutate(
    enriched_in = ifelse(diff > 0, "Cluster2_Immune-Inflammatory", "Cluster1_Proliferative")
  )

# Save results
write.csv(diff_test_results,
          "03_Results/16_Immune_Deconvolution/mcp_counter_cell_enrichment.csv",
          row.names = FALSE)

cat("Significant cell types (FDR < 0.05):\n")
sig_cells <- diff_test_results %>% filter(fdr < 0.05)
if (nrow(sig_cells) > 0) {
  print(sig_cells %>% select(cell_type, diff, p_value, fdr, enriched_in), row.names = FALSE)
} else {
  cat("  None at FDR < 0.05\n")
}
cat("\n")

cat("All cell types (sorted by FDR):\n")
print(diff_test_results, row.names = FALSE)
cat("\n")

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

cat("=== CREATING VISUALIZATIONS ===\n\n")

# Heatmap
result_matrix <- mcp_res %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Match samples
result_matrix_matched <- result_matrix[, colnames(result_matrix) %in% common_samples]

# Order by cluster
cluster_vec <- cluster_matched$cluster[match(colnames(result_matrix_matched), cluster_matched$sample_id)]
sample_order <- order(cluster_vec)
result_ordered <- result_matrix_matched[, sample_order]
cluster_ordered <- cluster_vec[sample_order]

# Annotation
annotation_col <- data.frame(
  Cluster = factor(cluster_ordered, labels = c("Proliferative", "Immune-Inflammatory"))
)
rownames(annotation_col) <- colnames(result_ordered)

annotation_colors <- list(
  Cluster = c("Proliferative" = "#2E9FDF", "Immune-Inflammatory" = "#E7B800")
)

# Heatmap
pdf("04_Figures/15_Immune_Deconvolution/immune_cell_heatmap.pdf",
    width = 14, height = 8)
pheatmap(result_ordered,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Immune Cell Deconvolution (MCP-counter)",
         fontsize_row = 10)
dev.off()

cat("✓ Saved: immune_cell_heatmap.pdf\n")

# Boxplot for all cell types
result_long_plot <- result_long %>%
  mutate(cluster_label = factor(cluster, labels = c("Proliferative", "Immune-Inflammatory")))

p <- ggplot(result_long_plot, aes(x = reorder(cell_type, score, FUN = median),
                                    y = score, fill = cluster_label)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Proliferative" = "#2E9FDF",
                               "Immune-Inflammatory" = "#E7B800")) +
  labs(x = "Cell Type", y = "MCP-counter Score",
       title = "Immune Cell Type Abundance by Molecular Subtype",
       fill = "Molecular Subtype") +
  theme_bw() +
  theme(text = element_text(size = 11))

ggsave("04_Figures/15_Immune_Deconvolution/all_immune_cells_boxplot.pdf",
       p, width = 10, height = 6)

cat("✓ Saved: all_immune_cells_boxplot.pdf\n")

# Bar plot of mean differences
diff_plot_data <- diff_test_results %>%
  mutate(
    cell_type = reorder(cell_type, diff),
    sig = ifelse(fdr < 0.05, "FDR < 0.05", "NS")
  )

p2 <- ggplot(diff_plot_data, aes(x = cell_type, y = diff, fill = sig)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("FDR < 0.05" = "#E74C3C", "NS" = "#95A5A6")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Cell Type",
       y = "Mean Difference (Cluster 2 - Cluster 1)",
       title = "Differential Immune Cell Abundance",
       subtitle = "Positive = enriched in Immune-Inflammatory subtype",
       fill = "Significance") +
  theme_bw() +
  theme(text = element_text(size = 11))

ggsave("04_Figures/15_Immune_Deconvolution/differential_immune_abundance.pdf",
       p2, width = 10, height = 6)

cat("✓ Saved: differential_immune_abundance.pdf\n\n")

# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("IMMUNE CELL DECONVOLUTION ANALYSIS COMPLETE\n\n")

cat("Method used: MCP-counter\n")
cat("Cell types analyzed:", nrow(diff_test_results), "\n")
cat("Significant cell types (FDR < 0.05):", nrow(sig_cells), "\n\n")

if (nrow(sig_cells) > 0) {
  cat("Key findings:\n")
  cat("  Enriched in Cluster 2 (Immune-Inflammatory):\n")
  cluster2_enriched <- sig_cells %>% filter(diff > 0)
  if (nrow(cluster2_enriched) > 0) {
    for (i in 1:nrow(cluster2_enriched)) {
      cat("    -", cluster2_enriched$cell_type[i],
          ": mean diff =", round(cluster2_enriched$diff[i], 2),
          ", FDR =", format(cluster2_enriched$fdr[i], scientific = TRUE, digits = 2), "\n")
    }
  }

  cat("\n  Enriched in Cluster 1 (Proliferative):\n")
  cluster1_enriched <- sig_cells %>% filter(diff < 0)
  if (nrow(cluster1_enriched) > 0) {
    for (i in 1:nrow(cluster1_enriched)) {
      cat("    -", cluster1_enriched$cell_type[i],
          ": mean diff =", round(cluster1_enriched$diff[i], 2),
          ", FDR =", format(cluster1_enriched$fdr[i], scientific = TRUE, digits = 2), "\n")
    }
  }
} else {
  cat("No significant differences at FDR < 0.05\n")
  cat("Top trends (lowest p-values):\n")
  top_trends <- diff_test_results %>% head(5)
  for (i in 1:nrow(top_trends)) {
    cat("  ", top_trends$cell_type[i], ": p =",
        format(top_trends$p_value[i], scientific = TRUE, digits = 2),
        ", enriched in", top_trends$enriched_in[i], "\n")
  }
}

cat("\n### Immune Deconvolution Task COMPLETE ###\n")
