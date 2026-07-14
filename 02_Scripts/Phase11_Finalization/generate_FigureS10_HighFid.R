# R script to generate Upgraded High-Fidelity Supplementary Figure S10 (Two-Panel Validation Suite)
# Panel A: Continuous VRS vs Venetoclax AUC validation scatter plot (BeatAML 2.0, Spearman p = 3.3e-37)
# Panel B: Discrete Predicted Cluster 1 vs Cluster 2 ex vivo Venetoclax AUC validation boxplot (Wilcoxon p = 2.7e-28)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)
library(randomForest)
library(gridExtra)
margin <- ggplot2::margin

cat("=== GENERATING DUAL-PANEL HIGH-FIDELITY SUPPLEMENTARY FIGURE S10 ===\n")

# ============================================================
# 1. LOAD CLASSIFIER & CLASSIFY VALIDATION SAMPLES
# ============================================================
cat("Predicting Cluster assignments for BeatAML 2.0 validation samples...\n")

# Load Classifier and Signature
rf_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")
signature_genes <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")

# Load Expression Data for BeatAML 2.0
expr_file <- "03_Results/29_ExternalValidation/beataml_waves1to4_norm_exp_dbgap.txt"
expr_raw <- read_tsv(expr_file, show_col_types = FALSE)
sample_cols <- grep("^BA", colnames(expr_raw), value = TRUE)

# Load VRS validation data
validation_df <- read_csv("03_Results/29_ExternalValidation/beataml2_VRS_validation.csv", show_col_types = FALSE)
common_samples <- intersect(sample_cols, validation_df$sample_id)

# Subset Expression to Signature Genes
expr_subset <- expr_raw %>%
    filter(stable_id %in% signature_genes$gene) %>%
    select(gene = stable_id, all_of(common_samples))

# Load BeatAML 1.0 (Discovery) to scale
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

# Scale BeatAML 2.0 to match BeatAML 1.0 distribution
common_genes <- intersect(signature_genes$gene, expr_subset$gene)
expr_subset_df <- expr_subset %>%
    filter(gene %in% common_genes) %>%
    column_to_rownames("gene")

# Scale each gene
for (gene in common_genes) {
    if (gene %in% rownames(beataml_expr)) {
        b1_mean <- mean(beataml_expr[gene, ], na.rm = TRUE)
        b1_sd <- sd(beataml_expr[gene, ], na.rm = TRUE)
        
        b2_mean <- mean(as.numeric(expr_subset_df[gene, ]), na.rm = TRUE)
        b2_sd <- sd(as.numeric(expr_subset_df[gene, ]), na.rm = TRUE)
        
        if (b2_sd > 0) {
            expr_subset_df[gene, ] <- ((as.numeric(expr_subset_df[gene, ]) - b2_mean) / b2_sd) * b1_sd + b1_mean
        }
    }
}

# Predict Clusters
pred_data <- as.data.frame(t(expr_subset_df))
missing_genes <- setdiff(signature_genes$gene, colnames(pred_data))
for (gene in missing_genes) {
    pred_data[[gene]] <- 0
}
pred_data <- pred_data[, signature_genes$gene]
colnames(pred_data) <- make.names(colnames(pred_data))

predictions <- predict(rf_classifier, newdata = pred_data, type = "response")
prediction_probs <- predict(rf_classifier, newdata = pred_data, type = "prob")

pred_df <- data.frame(
    sample_id = rownames(pred_data),
    predicted_cluster = predictions,
    stringsAsFactors = FALSE
)

# Merge predictions with validation data
merged_df <- inner_join(validation_df, pred_df, by = "sample_id")

# ============================================================
# 2. SET UP HIGH-FIDELITY THEME
# ============================================================
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 14, color = "darkblue", margin = margin(b=10)),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "plain", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "plain", color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    plot.margin = margin(15, 15, 15, 15)
  )

# ============================================================
# 3. PANEL A: CONTINUOUS VRS CORRELATION
# ============================================================
p_scatter <- ggplot(merged_df, aes(x = VRS, y = venetoclax_auc)) +
  geom_point(alpha = 0.5, color = "#2C3E50", size = 3) +
  geom_smooth(method = "lm", color = "#E74C3C", linewidth = 1.5, se = TRUE) +
  annotate("text", x = 70, y = 260, 
           label = "Spearman r = -0.60\np = 3.3e-37", 
           size = 4.0, fontface = "bold", color = "darkred", hjust = 0) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(
    title = "A. Continuous VRS Validation",
    subtitle = "Independent Cohort Validation (N = 367)",
    x = "Venetoclax Response Score (VRS)",
    y = "Ex Vivo Venetoclax AUC (Lower = Sensitive)"
  ) +
  theme_hf

# ============================================================
# 4. PANEL B: DISCRETE CLUSTER RESPONSE DIFFERENCE
# ============================================================
color_c1 <- "#3498DB"  # Soft Blue (Primitive/Cluster 1)
color_c2 <- "#E67E22"  # Soft Orange (Monocytic/Cluster 2)

p_boxplot <- ggplot(merged_df, aes(x = factor(predicted_cluster, labels = c("Cluster 1\n(Primitive)", "Cluster 2\n(Monocytic)")), 
                                  y = venetoclax_auc, 
                                  fill = factor(predicted_cluster))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.4, linewidth = 1.0, color = "black") +
  geom_jitter(aes(color = factor(predicted_cluster)), width = 0.15, alpha = 0.4, size = 1.8) +
  scale_fill_manual(values = c("1" = color_c1, "2" = color_c2)) +
  scale_color_manual(values = c("1" = color_c1, "2" = color_c2)) +
  scale_y_continuous(limits = c(0, 360)) +
  annotate("text", x = 1.5, y = 330, 
           label = "Wilcoxon p = 2.7e-28\nMean AUC: 102.3 vs 193.9", 
           size = 4.0, fontface = "bold", color = "black", hjust = 0.5) +
  labs(
    title = "B. Discrete Subtype Validation",
    subtitle = "Ex Vivo Venetoclax AUC by Predicted Cluster",
    x = "Predicted Molecular Subtype",
    y = "Ex Vivo Venetoclax AUC (Lower = Sensitive)"
  ) +
  theme_hf +
  theme(legend.position = "none")

# ============================================================
# 5. ASSEMBLE AND SAVE
# ============================================================
output_dir <- "05_Submission/Submission_Hub/03_Supplementary_Figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(output_dir, "FigureS10.pdf"), width = 14, height = 7)
grid.arrange(p_scatter, p_boxplot, ncol = 2)
dev.off()

png(file.path(output_dir, "FigureS10.png"), width = 14, height = 7, units = "in", res = 300)
grid.arrange(p_scatter, p_boxplot, ncol = 2)
dev.off()

cat("✓ Dual-panel High-Fidelity Supplementary Figure S10 generated successfully!\n")
