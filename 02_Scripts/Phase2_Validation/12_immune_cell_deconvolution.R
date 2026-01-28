# TASK: Immune Cell Deconvolution Analysis
# Validate "immune-inflammatory" label for Cluster 2 using deconvolution methods

library(tidyverse)
library(immunedeconv)
library(ggplot2)
library(pheatmap)

cat("=== IMMUNE CELL DECONVOLUTION ANALYSIS ===\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading expression data...\n")

# Load batch-corrected expression data
expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
cat("Expression data:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")

# Convert Ensembl IDs to gene symbols
# Load annotation from original expression file (just first 2 columns for speed)
cat("Loading gene annotations...\n")
original_expr <- read.delim("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
                            stringsAsFactors = FALSE, sep = "\t")
gene_annotations <- original_expr %>%
  select(Gene_ID = Gene.ID, Gene_Symbol = Gene.Symbol) %>%
  distinct()
cat("Gene annotations:", nrow(gene_annotations), "genes\n")

# Match Ensembl IDs with gene symbols
ensembl_ids <- rownames(expr_data)
matched_symbols <- gene_annotations$Gene_Symbol[match(ensembl_ids, gene_annotations$Gene_ID)]

# Remove genes without symbols
has_symbol <- !is.na(matched_symbols) & matched_symbols != ""
expr_data_symbol <- expr_data[has_symbol, ]
rownames(expr_data_symbol) <- matched_symbols[has_symbol]

# Remove duplicates by aggregating (take mean)
expr_aggregated <- aggregate(expr_data_symbol,
                             by = list(gene = rownames(expr_data_symbol)),
                             FUN = mean)
rownames(expr_aggregated) <- expr_aggregated$gene
expr_aggregated$gene <- NULL
expr_data <- as.matrix(expr_aggregated)

cat("After converting to gene symbols:", nrow(expr_data), "genes\n")

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
cat("Cluster assignments:", nrow(clusters), "samples\n\n")

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

# immunedeconv requires:
# - Gene symbols as rownames
# - Linear scale expression (not log-transformed)
# - Non-negative values

# Check if data is log-transformed (typically log2(TPM+1) or similar)
# If max value is < 30, likely log-transformed
if (max(expr_matched, na.rm = TRUE) < 30) {
  cat("Data appears log-transformed, converting to linear scale...\n")
  expr_linear <- 2^expr_matched - 1  # Assuming log2(x+1) transformation
  expr_linear[expr_linear < 0] <- 0  # Remove any negative values
} else {
  cat("Data appears to be in linear scale\n")
  expr_linear <- expr_matched
}

# Ensure non-negative
expr_linear[expr_linear < 0] <- 0

cat("Expression range after conversion:", round(min(expr_linear), 2),
    "to", round(max(expr_linear), 2), "\n\n")

# ============================================================================
# 3. RUN DECONVOLUTION METHODS
# ============================================================================

cat("=== RUNNING DECONVOLUTION METHODS ===\n\n")

# We'll use multiple methods available in immunedeconv
# Different methods work better for different tissues

deconv_results <- list()

# Method 1: EPIC (good for tumor samples)
cat("Running EPIC deconvolution...\n")
tryCatch({
  epic_res <- deconvolute(expr_linear, "epic")
  deconv_results$epic <- epic_res
  cat("✓ EPIC completed\n")
}, error = function(e) {
  cat("✗ EPIC failed:", e$message, "\n")
})

# Method 2: quanTIseq (designed for RNA-seq)
cat("Running quanTIseq deconvolution...\n")
tryCatch({
  quantiseq_res <- deconvolute(expr_linear, "quantiseq", tumor = TRUE)
  deconv_results$quantiseq <- quantiseq_res
  cat("✓ quanTIseq completed\n")
}, error = function(e) {
  cat("✗ quanTIseq failed:", e$message, "\n")
})

# Method 3: xCell (comprehensive cell types)
cat("Running xCell deconvolution...\n")
tryCatch({
  xcell_res <- deconvolute(expr_linear, "xcell")
  deconv_results$xcell <- xcell_res
  cat("✓ xCell completed\n")
}, error = function(e) {
  cat("✗ xCell failed:", e$message, "\n")
})

# Method 4: MCP-counter (simple, robust)
cat("Running MCP-counter deconvolution...\n")
tryCatch({
  mcp_res <- deconvolute(expr_linear, "mcp_counter")
  deconv_results$mcp_counter <- mcp_res
  cat("✓ MCP-counter completed\n")
}, error = function(e) {
  cat("✗ MCP-counter failed:", e$message, "\n")
})

cat("\nCompleted", length(deconv_results), "deconvolution methods\n\n")

# ============================================================================
# 4. ANALYZE RESULTS
# ============================================================================

cat("=== ANALYZING DECONVOLUTION RESULTS ===\n\n")

dir.create("03_Results/16_Immune_Deconvolution", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/15_Immune_Deconvolution", showWarnings = FALSE, recursive = TRUE)

# Process each method
for (method_name in names(deconv_results)) {
  cat("--- Method:", method_name, "---\n")

  result_df <- deconv_results[[method_name]]

  # Convert to long format
  # Check if cell_type column exists
  if ("cell_type" %in% colnames(result_df)) {
    result_long <- result_df %>%
      pivot_longer(-cell_type, names_to = "sample_id", values_to = "score") %>%
      left_join(cluster_matched %>% select(sample_id, cluster), by = "sample_id")
  } else {
    # cell_type is rownames
    result_df <- result_df %>%
      rownames_to_column("cell_type")
    result_long <- result_df %>%
      pivot_longer(-cell_type, names_to = "sample_id", values_to = "score") %>%
      left_join(cluster_matched %>% select(sample_id, cluster), by = "sample_id")
  }

  # Calculate mean scores by cluster
  cluster_summary <- result_long %>%
    group_by(cell_type, cluster) %>%
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      median_score = median(score, na.rm = TRUE),
      .groups = "drop"
    )

  # Test for differences
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
        method = method_name,
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
  diff_test_results <- diff_test_results %>% arrange(fdr)

  # Save results
  write.csv(diff_test_results,
            paste0("03_Results/16_Immune_Deconvolution/", method_name, "_cell_enrichment.csv"),
            row.names = FALSE)

  cat("Significant cell types (FDR < 0.05):\n")
  sig_cells <- diff_test_results %>% filter(fdr < 0.05)
  if (nrow(sig_cells) > 0) {
    print(sig_cells %>% select(cell_type, diff, p_value, fdr), row.names = FALSE)
  } else {
    cat("  None\n")
  }
  cat("\n")
}

# ============================================================================
# 5. COMBINE RESULTS ACROSS METHODS
# ============================================================================

cat("=== CONSENSUS ACROSS METHODS ===\n\n")

# Load all results
all_results <- list()
for (method_name in names(deconv_results)) {
  all_results[[method_name]] <- read.csv(
    paste0("03_Results/16_Immune_Deconvolution/", method_name, "_cell_enrichment.csv")
  )
}

# Combine
combined_results <- bind_rows(all_results)

# Identify cell types significant in multiple methods
consensus_cells <- combined_results %>%
  filter(fdr < 0.05) %>%
  group_by(cell_type) %>%
  summarise(
    n_methods = n(),
    mean_diff = mean(diff),
    min_fdr = min(fdr),
    methods = paste(method, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_methods), min_fdr)

cat("Cell types enriched in Cluster 2 (Immune-Inflammatory) across methods:\n")
print(consensus_cells %>% filter(mean_diff > 0), row.names = FALSE)

cat("\nCell types enriched in Cluster 1 (Proliferative) across methods:\n")
print(consensus_cells %>% filter(mean_diff < 0), row.names = FALSE)

# Save consensus
write.csv(consensus_cells,
          "03_Results/16_Immune_Deconvolution/consensus_immune_enrichment.csv",
          row.names = FALSE)

# ============================================================================
# 6. VISUALIZATIONS
# ============================================================================

cat("\n=== CREATING VISUALIZATIONS ===\n\n")

# Use the method with most results (usually xCell or quanTIseq)
best_method <- names(deconv_results)[1]
if ("quantiseq" %in% names(deconv_results)) {
  best_method <- "quantiseq"
} else if ("xcell" %in% names(deconv_results)) {
  best_method <- "xcell"
}

cat("Creating plots using", best_method, "results...\n")

result_df <- deconv_results[[best_method]]

# Convert to matrix for heatmap
result_matrix <- result_df %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# Match with clusters
result_matrix <- result_matrix[, common_samples]
cluster_vec <- cluster_matched$cluster[match(common_samples, cluster_matched$sample_id)]

# Order by cluster
sample_order <- order(cluster_vec)
result_ordered <- result_matrix[, sample_order]
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
         main = paste0("Immune Cell Deconvolution (", best_method, ")"),
         fontsize_row = 8)
dev.off()

cat("✓ Saved: immune_cell_heatmap.pdf\n")

# Box plots for top differential cell types
method_results <- read.csv(
  paste0("03_Results/16_Immune_Deconvolution/", best_method, "_cell_enrichment.csv")
)

top_cells <- method_results %>%
  arrange(fdr) %>%
  head(10) %>%
  pull(cell_type)

# Prepare data for plotting
result_long <- result_df %>%
  pivot_longer(-cell_type, names_to = "sample_id", values_to = "score") %>%
  left_join(cluster_matched %>% select(sample_id, cluster), by = "sample_id") %>%
  filter(cell_type %in% top_cells) %>%
  mutate(cluster_label = factor(cluster, labels = c("Proliferative", "Immune-Inflammatory")))

# Boxplot
p <- ggplot(result_long, aes(x = cell_type, y = score, fill = cluster_label)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Proliferative" = "#2E9FDF",
                               "Immune-Inflammatory" = "#E7B800")) +
  labs(x = "Cell Type", y = "Deconvolution Score",
       title = paste0("Top 10 Differential Immune Cell Types (", best_method, ")"),
       fill = "Molecular Subtype") +
  theme_bw() +
  theme(text = element_text(size = 11))

ggsave("04_Figures/15_Immune_Deconvolution/top10_immune_cells_boxplot.pdf",
       p, width = 10, height = 6)

cat("✓ Saved: top10_immune_cells_boxplot.pdf\n\n")

# ============================================================================
# 7. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("IMMUNE CELL DECONVOLUTION COMPLETE\n\n")

cat("Methods completed:", length(deconv_results), "\n")
cat("  Methods:", paste(names(deconv_results), collapse = ", "), "\n\n")

cat("Consensus findings (cell types enriched across multiple methods):\n")
cat("  Cell types enriched in Cluster 2 (n=", nrow(consensus_cells %>% filter(mean_diff > 0)), "):\n")
if (nrow(consensus_cells %>% filter(mean_diff > 0)) > 0) {
  top_cluster2 <- consensus_cells %>%
    filter(mean_diff > 0) %>%
    head(5)
  for (i in 1:nrow(top_cluster2)) {
    cat("    -", top_cluster2$cell_type[i], "(", top_cluster2$n_methods[i], "methods, FDR=",
        format(top_cluster2$min_fdr[i], scientific = TRUE, digits = 2), ")\n")
  }
}

cat("\n  Cell types enriched in Cluster 1 (n=", nrow(consensus_cells %>% filter(mean_diff < 0)), "):\n")
if (nrow(consensus_cells %>% filter(mean_diff < 0)) > 0) {
  top_cluster1 <- consensus_cells %>%
    filter(mean_diff < 0) %>%
    head(5)
  for (i in 1:nrow(top_cluster1)) {
    cat("    -", top_cluster1$cell_type[i], "(", top_cluster1$n_methods[i], "methods, FDR=",
        format(top_cluster1$min_fdr[i], scientific = TRUE, digits = 2), ")\n")
  }
}

cat("\n### Immune Deconvolution Task COMPLETE ###\n")
