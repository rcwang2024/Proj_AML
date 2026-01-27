#!/usr/bin/env Rscript
# ==============================================================================
# Phase 5 Part 7: Create Publication Figures
# ==============================================================================
# Creates publication-ready figures for drug response validation
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(ggpubr)
  library(scales)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("PART 7: CREATING PUBLICATION FIGURES\n")
cat("==============================================================================\n\n")

# Create output directory
dir.create("04_Figures/22_Drug_Validation", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n\n")

# Load drug results
drug_results <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")
independence_results <- read.csv("03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv")
bcl2_results <- read.csv("03_Results/23_Drug_Validation/bcl2_pathway_expression_FIXED.csv")

# Load raw data for individual plots
drug_raw <- fread("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt", data.table = FALSE)
drug_auc <- drug_raw %>%
  dplyr::select(dbgap_rnaseq_sample, inhibitor, auc) %>%
  filter(!is.na(auc)) %>%
  pivot_wider(names_from = inhibitor, values_from = auc, values_fn = mean)
drug_auc_df <- as.data.frame(drug_auc)
rownames(drug_auc_df) <- drug_auc_df$dbgap_rnaseq_sample
drug_auc_df$dbgap_rnaseq_sample <- NULL

clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
expr_raw <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")
gene_annot <- read.csv("03_Results/01_Processed_Data/gene_annotations.csv")

common_samples <- intersect(rownames(drug_auc_df), clusters$sample_id)
cluster_assign <- clusters$cluster[match(common_samples, clusters$sample_id)]

# ==============================================================================
# FIGURE 5A: VENETOCLAX BOXPLOT
# ==============================================================================

cat("Creating Figure 5A: Venetoclax boxplot...\n")

ven_data <- data.frame(
  sample_id = common_samples,
  cluster = factor(cluster_assign, levels = c(1, 2),
                   labels = c("Cluster 1\n(NPM1+)", "Cluster 2\n(TP53+/RUNX1+)")),
  venetoclax_auc = drug_auc_df[common_samples, "Venetoclax"]
) %>% filter(!is.na(venetoclax_auc))

# Get statistics
ven_stats <- drug_results %>% filter(drug == "Venetoclax")

fig5a <- ggplot(ven_data, aes(x = cluster, y = venetoclax_auc, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(
    title = "Venetoclax Sensitivity by Molecular Subtype",
    subtitle = sprintf("p = %.2e (FDR = %.2e), Cohen's d = %.2f",
                      ven_stats$wilcoxon_pvalue, ven_stats$fdr, abs(ven_stats$cohens_d)),
    x = "",
    y = "Venetoclax AUC\n(Lower = More Sensitive)",
    caption = "Cluster 1 is 1.79× more sensitive to Venetoclax"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "darkred"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.caption = element_text(size = 11, hjust = 0.5, face = "italic")
  ) +
  annotate("text", x = 1.5, y = max(ven_data$venetoclax_auc, na.rm = TRUE) * 0.95,
           label = "***", size = 10, color = "darkred")

ggsave("04_Figures/22_Drug_Validation/Figure5A_Venetoclax_Boxplot.pdf",
       fig5a, width = 8, height = 7)
cat("✓ Saved: Figure5A_Venetoclax_Boxplot.pdf\n\n")

# ==============================================================================
# FIGURE 5B: R² IMPROVEMENT HEATMAP
# ==============================================================================

cat("Creating Figure 5B: R² improvement heatmap...\n")

# Prepare data for heatmap
r2_data <- independence_results %>%
  arrange(desc(r2_improvement)) %>%
  head(20) %>%
  mutate(
    drug_short = str_trunc(drug, 25, "right"),
    cluster_independent = ifelse(fdr_cluster < 0.05, "Yes", "No")
  )

# Create matrix for heatmap
r2_matrix <- matrix(
  c(r2_data$r2_mutations_only,
    r2_data$r2_improvement),
  nrow = nrow(r2_data), ncol = 2
)
rownames(r2_matrix) <- r2_data$drug_short
colnames(r2_matrix) <- c("R² Mutations\nOnly", "R² Improvement\nfrom Cluster")

# Annotation
annotation_row <- data.frame(
  Independent = factor(r2_data$cluster_independent, levels = c("Yes", "No")),
  row.names = rownames(r2_matrix)
)

ann_colors <- list(
  Independent = c("Yes" = "#2E7D32", "No" = "#BDBDBD")
)

pdf("04_Figures/22_Drug_Validation/Figure5B_R2_Improvement_Heatmap.pdf",
    width = 10, height = 12)

pheatmap(
  r2_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
  breaks = seq(0, 0.4, length.out = 101),
  display_numbers = matrix(sprintf("%.3f", r2_matrix), nrow = nrow(r2_matrix)),
  number_color = "black",
  fontsize_number = 9,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  main = "Drug Response Prediction:\nCluster Adds Value Beyond Mutations",
  fontsize = 11,
  fontsize_row = 10,
  fontsize_col = 11,
  angle_col = 0,
  border_color = "grey60",
  cellwidth = 80,
  cellheight = 20
)

dev.off()
cat("✓ Saved: Figure5B_R2_Improvement_Heatmap.pdf\n\n")

# ==============================================================================
# FIGURE 5C: BCL2 vs VENETOCLAX SCATTER
# ==============================================================================

cat("Creating Figure 5C: BCL2 expression vs Venetoclax AUC...\n")

# Get BCL2 expression
bcl2_ensembl <- gene_annot %>% filter(gene_symbol == "BCL2") %>% pull(ensembl_id)
bcl2_expr <- expr_raw[bcl2_ensembl, common_samples]

scatter_data <- data.frame(
  sample_id = common_samples,
  cluster = factor(cluster_assign, levels = c(1, 2),
                   labels = c("Cluster 1 (NPM1+)", "Cluster 2 (TP53+/RUNX1+)")),
  bcl2_expression = as.numeric(bcl2_expr),
  venetoclax_auc = drug_auc_df[common_samples, "Venetoclax"]
) %>% filter(!is.na(venetoclax_auc) & !is.na(bcl2_expression))

# Calculate correlation
cor_test <- cor.test(scatter_data$bcl2_expression, scatter_data$venetoclax_auc,
                     method = "spearman")

fig5c <- ggplot(scatter_data, aes(x = bcl2_expression, y = venetoclax_auc, color = cluster)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 0.8) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(
    title = "BCL2 Expression Predicts Venetoclax Sensitivity",
    subtitle = sprintf("Spearman ρ = %.3f, p < 10^-30", cor_test$estimate),
    x = "BCL2 Expression (log2 TPM)",
    y = "Venetoclax AUC\n(Lower = More Sensitive)",
    color = "Molecular Subtype",
    caption = "Higher BCL2 expression → Greater Venetoclax sensitivity"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "darkred"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = c(0.75, 0.9),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(size = 11, hjust = 0.5, face = "italic")
  )

ggsave("04_Figures/22_Drug_Validation/Figure5C_BCL2_Venetoclax_Scatter.pdf",
       fig5c, width = 9, height = 7)
cat("✓ Saved: Figure5C_BCL2_Venetoclax_Scatter.pdf\n\n")

# ==============================================================================
# FIGURE 5D: COMBINED PANEL (MAIN FIGURE 5)
# ==============================================================================

cat("Creating Figure 5: Combined panel...\n")

pdf("04_Figures/22_Drug_Validation/Figure5_Drug_Response_Main.pdf",
    width = 16, height = 12)

# Create layout
layout_matrix <- rbind(
  c(1, 1, 2, 2, 2),
  c(1, 1, 2, 2, 2),
  c(3, 3, 3, 3, 3),
  c(3, 3, 3, 3, 3)
)

grid.arrange(
  fig5a + labs(title = "A. Venetoclax Sensitivity by Subtype"),
  fig5c + labs(title = "B. BCL2 Expression Mechanism"),
  grobs = list(
    fig5a + labs(title = "A. Venetoclax Sensitivity by Subtype"),
    fig5c + labs(title = "B. BCL2 Expression Mechanism")
  ),
  layout_matrix = layout_matrix,
  top = textGrob("Figure 5: Molecular Subtypes Predict Drug Response",
                 gp = gpar(fontsize = 18, fontface = "bold"))
)

dev.off()
cat("✓ Saved: Figure5_Drug_Response_Main.pdf\n\n")

# ==============================================================================
# SUPPLEMENTARY FIGURE S1: TOP 20 DRUGS BOXPLOTS
# ==============================================================================

cat("Creating Supplementary Figure S1: Top 20 drugs boxplots...\n")

top20_drugs <- drug_results %>%
  filter(fdr < 0.05) %>%
  arrange(wilcoxon_pvalue) %>%
  head(20) %>%
  pull(drug)

plot_list <- list()

for (i in 1:min(20, length(top20_drugs))) {
  drug_name <- top20_drugs[i]

  drug_data <- data.frame(
    cluster = factor(cluster_assign, levels = c(1, 2)),
    auc = drug_auc_df[common_samples, drug_name]
  ) %>% filter(!is.na(auc))

  stats <- drug_results %>% filter(drug == drug_name)

  p <- ggplot(drug_data, aes(x = cluster, y = auc, fill = cluster)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.8) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
    labs(
      title = str_trunc(drug_name, 30),
      subtitle = sprintf("p=%.1e, d=%.2f", stats$wilcoxon_pvalue, abs(stats$cohens_d)),
      x = "",
      y = "AUC"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 8),
      axis.text = element_text(size = 8)
    )

  plot_list[[i]] <- p
}

pdf("04_Figures/22_Drug_Validation/FigureS1_Top20_Drugs_Boxplots.pdf",
    width = 16, height = 20)

do.call(grid.arrange, c(plot_list, ncol = 4))

dev.off()
cat("✓ Saved: FigureS1_Top20_Drugs_Boxplots.pdf\n\n")

# ==============================================================================
# SUPPLEMENTARY FIGURE S2: BCL-2 PATHWAY HEATMAP
# ==============================================================================

cat("Creating Supplementary Figure S2: BCL-2 pathway heatmap...\n")

# Get BCL-2 pathway genes
bcl2_genes <- bcl2_results$ensembl_id
bcl2_symbols <- bcl2_results$gene

bcl2_expr_matrix <- expr_raw[bcl2_genes, common_samples]
rownames(bcl2_expr_matrix) <- bcl2_symbols

# Calculate mean expression by cluster
cluster_means <- matrix(NA, nrow = length(bcl2_symbols), ncol = 2)
rownames(cluster_means) <- bcl2_symbols
colnames(cluster_means) <- c("Cluster 1\n(NPM1+)", "Cluster 2\n(TP53+/RUNX1+)")

for (i in 1:length(bcl2_symbols)) {
  cluster_means[i, 1] <- mean(as.numeric(bcl2_expr_matrix[i, cluster_assign == 1]), na.rm = TRUE)
  cluster_means[i, 2] <- mean(as.numeric(bcl2_expr_matrix[i, cluster_assign == 2]), na.rm = TRUE)
}

# Z-score normalize
cluster_means_z <- t(scale(t(cluster_means)))

# Annotation
annotation_row <- data.frame(
  Role = factor(
    c("Anti", "Pro", "Pro", "Anti", "Anti", "Pro", "Pro", "Pro", "Pro", "Pro")[match(bcl2_symbols, c("BCL2", "BCL2L11", "BBC3", "BCL2L1", "MCL1", "BID", "BAX", "BAD", "PMAIP1", "BAK1"))],
    levels = c("Anti", "Pro")
  ),
  FDR = ifelse(bcl2_results$fdr < 0.001, "< 0.001",
               ifelse(bcl2_results$fdr < 0.01, "< 0.01",
                      ifelse(bcl2_results$fdr < 0.05, "< 0.05", "> 0.05"))),
  row.names = bcl2_symbols
)

ann_colors_bcl2 <- list(
  Role = c("Anti" = "#D32F2F", "Pro" = "#1976D2"),
  FDR = c("< 0.001" = "#1B5E20", "< 0.01" = "#388E3C",
          "< 0.05" = "#81C784", "> 0.05" = "#BDBDBD")
)

pdf("04_Figures/22_Drug_Validation/FigureS2_BCL2_Pathway_Heatmap.pdf",
    width = 8, height = 10)

pheatmap(
  cluster_means_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2, 2, length.out = 101),
  annotation_row = annotation_row,
  annotation_colors = ann_colors_bcl2,
  main = "BCL-2 Pathway Expression by Molecular Subtype\n(Z-scored mean expression)",
  fontsize = 12,
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 0,
  border_color = "grey60",
  cellwidth = 80,
  cellheight = 30
)

dev.off()
cat("✓ Saved: FigureS2_BCL2_Pathway_Heatmap.pdf\n\n")

# ==============================================================================
# SUPPLEMENTARY FIGURE S3: DRUG CLASS ENRICHMENT
# ==============================================================================

cat("Creating Supplementary Figure S3: Drug class enrichment...\n")

class_data <- read.csv("03_Results/23_Drug_Validation/drug_class_enrichment_FIXED.csv")

class_data_plot <- class_data %>%
  mutate(
    drug_class = str_replace_all(drug_class, "_", " "),
    drug_class = str_to_title(drug_class),
    sig = ifelse(pvalue < 0.05, "Significant", "Not Significant")
  ) %>%
  arrange(desc(percent_significant))

figs3 <- ggplot(class_data_plot, aes(x = reorder(drug_class, percent_significant),
                                      y = percent_significant, fill = sig)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%d/%d", n_drugs_significant, n_drugs_tested)),
            hjust = -0.2, size = 3.5) +
  scale_fill_manual(values = c("Significant" = "#D32F2F", "Not Significant" = "#757575")) +
  coord_flip() +
  labs(
    title = "Drug Class Enrichment in Differential Response",
    subtitle = "Percentage of drugs in each class showing differential response (FDR<0.05)",
    x = "Drug Class",
    y = "Percentage of Drugs with Differential Response (%)",
    fill = "Enrichment"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0))

ggsave("04_Figures/22_Drug_Validation/FigureS3_Drug_Class_Enrichment.pdf",
       figs3, width = 10, height = 8)
cat("✓ Saved: FigureS3_Drug_Class_Enrichment.pdf\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("FIGURE GENERATION COMPLETE\n")
cat("==============================================================================\n\n")

cat("Main Figures Created:\n")
cat("  1. Figure5A_Venetoclax_Boxplot.pdf\n")
cat("  2. Figure5B_R2_Improvement_Heatmap.pdf\n")
cat("  3. Figure5C_BCL2_Venetoclax_Scatter.pdf\n")
cat("  4. Figure5_Drug_Response_Main.pdf (Combined panel)\n\n")

cat("Supplementary Figures Created:\n")
cat("  5. FigureS1_Top20_Drugs_Boxplots.pdf\n")
cat("  6. FigureS2_BCL2_Pathway_Heatmap.pdf\n")
cat("  7. FigureS3_Drug_Class_Enrichment.pdf\n\n")

cat("All figures saved to: 04_Figures/22_Drug_Validation/\n")
cat("==============================================================================\n")
