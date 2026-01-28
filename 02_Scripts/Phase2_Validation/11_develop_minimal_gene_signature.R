# TASK: Develop Minimal 50-Gene Signature Classifier
# Use LASSO + Random Forest to identify minimal gene set for clinical deployment

library(tidyverse)
library(glmnet)
library(randomForest)
library(caret)
library(pROC)

cat("=== DEVELOPING MINIMAL 50-GENE SIGNATURE ===\n\n")

set.seed(42)  # For reproducibility

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load batch-corrected expression data
expr_data <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
cat("Expression data:", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
cat("Cluster assignments:", nrow(clusters), "samples\n\n")

# Match samples
common_samples <- intersect(colnames(expr_data), clusters$sample_id)
expr_matched <- expr_data[, common_samples]
cluster_vec <- clusters$cluster[match(common_samples, clusters$sample_id)]

cat("Matched samples:", length(common_samples), "\n")
cat("Cluster 1:", sum(cluster_vec == 1), "\n")
cat("Cluster 2:", sum(cluster_vec == 2), "\n\n")

# ============================================================================
# 2. FEATURE SELECTION: DIFFERENTIAL EXPRESSION
# ============================================================================

cat("=== STEP 1: DIFFERENTIAL EXPRESSION ===\n\n")

# Wilcoxon test for each gene
cat("Running differential expression analysis...\n")

de_results <- data.frame()

for (i in 1:nrow(expr_matched)) {
  gene_name <- rownames(expr_matched)[i]

  cluster1_expr <- expr_matched[i, cluster_vec == 1]
  cluster2_expr <- expr_matched[i, cluster_vec == 2]

  # Wilcoxon test
  wilcox_result <- wilcox.test(cluster1_expr, cluster2_expr)

  # Fold change (log2)
  mean_c1 <- mean(cluster1_expr, na.rm = TRUE)
  mean_c2 <- mean(cluster2_expr, na.rm = TRUE)
  log2fc <- log2((mean_c1 + 0.01) / (mean_c2 + 0.01))

  de_results <- rbind(de_results, data.frame(
    gene = gene_name,
    mean_cluster1 = mean_c1,
    mean_cluster2 = mean_c2,
    log2fc = log2fc,
    p_value = wilcox_result$p.value,
    stringsAsFactors = FALSE
  ))

  if (i %% 2000 == 0) {
    cat("  Processed", i, "genes...\n")
  }
}

# FDR correction
de_results$fdr <- p.adjust(de_results$p_value, method = "BH")

# Sort by p-value
de_results <- de_results %>% arrange(p_value)

cat("\nTop 10 differentially expressed genes:\n")
print(de_results %>% select(gene, log2fc, p_value, fdr) %>% head(10), row.names = FALSE)

# Select top 200 DE genes for further analysis
top_de_genes <- de_results$gene[1:200]
cat("\nSelected top 200 DE genes for feature selection\n\n")

# ============================================================================
# 3. FEATURE SELECTION: LASSO REGULARIZATION
# ============================================================================

cat("=== STEP 2: LASSO FEATURE SELECTION ===\n\n")

# Prepare data for LASSO (use top 200 DE genes)
expr_subset <- t(expr_matched[top_de_genes, ])  # Transpose: samples x genes
y <- ifelse(cluster_vec == 1, 0, 1)  # Binary outcome

cat("Input dimensions:", nrow(expr_subset), "samples x", ncol(expr_subset), "genes\n")

# Cross-validated LASSO
cat("Running cross-validated LASSO...\n")
cv_lasso <- cv.glmnet(expr_subset, y,
                      family = "binomial",
                      alpha = 1,  # LASSO penalty
                      nfolds = 10,
                      type.measure = "class")

cat("Optimal lambda:", cv_lasso$lambda.min, "\n")
cat("1SE lambda:", cv_lasso$lambda.1se, "\n")

# Extract coefficients at lambda.1se (more parsimonious)
lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
lasso_genes <- rownames(lasso_coefs)[which(lasso_coefs != 0)][-1]  # Remove intercept

cat("LASSO selected", length(lasso_genes), "genes\n")
cat("Top 10 LASSO genes:\n")
print(head(lasso_genes, 10))
cat("\n")

# ============================================================================
# 4. FEATURE SELECTION: RANDOM FOREST IMPORTANCE
# ============================================================================

cat("=== STEP 3: RANDOM FOREST VARIABLE IMPORTANCE ===\n\n")

# Prepare data
rf_data <- data.frame(expr_subset)
rf_data$cluster <- as.factor(cluster_vec)

cat("Training Random Forest (this may take a few minutes)...\n")
rf_model <- randomForest(cluster ~ .,
                         data = rf_data,
                         ntree = 500,
                         importance = TRUE,
                         mtry = sqrt(ncol(expr_subset)))

cat("Random Forest OOB error rate:", round(rf_model$err.rate[500, "OOB"] * 100, 2), "%\n")

# Extract variable importance
importance_scores <- importance(rf_model, type = 1)  # Mean decrease in accuracy
importance_df <- data.frame(
  gene = rownames(importance_scores),
  importance = importance_scores[, 1],
  stringsAsFactors = FALSE
) %>% arrange(desc(importance))

cat("\nTop 10 most important genes by Random Forest:\n")
print(head(importance_df, 10), row.names = FALSE)

# Select top 50 by importance
rf_top_genes <- importance_df$gene[1:50]
cat("\nSelected top 50 genes by RF importance\n\n")

# ============================================================================
# 5. COMBINE METHODS: MINIMAL SIGNATURE
# ============================================================================

cat("=== STEP 4: CREATING MINIMAL 50-GENE SIGNATURE ===\n\n")

# Combine genes from LASSO and RF
combined_genes <- unique(c(lasso_genes, rf_top_genes))

cat("Combined gene pool:", length(combined_genes), "genes\n")

# If more than 50, prioritize by:
# 1. Present in both LASSO and RF
# 2. Highest RF importance
# 3. Lowest DE p-value

genes_in_both <- intersect(lasso_genes, rf_top_genes)
cat("Genes selected by both methods:", length(genes_in_both), "\n")

# Rank combined genes by multiple criteria
gene_scores <- data.frame(
  gene = combined_genes,
  in_lasso = combined_genes %in% lasso_genes,
  in_rf = combined_genes %in% rf_top_genes,
  stringsAsFactors = FALSE
) %>%
  left_join(importance_df, by = "gene") %>%
  left_join(de_results %>% select(gene, p_value, log2fc), by = "gene") %>%
  mutate(
    importance = ifelse(is.na(importance), 0, importance),
    score = (as.numeric(in_lasso) + as.numeric(in_rf)) * 1000 +
            importance - log10(p_value) * 10
  ) %>%
  arrange(desc(score))

# Select top 50
minimal_signature <- gene_scores$gene[1:50]

cat("\nFinal 50-gene signature created\n")
cat("  Genes in both LASSO and RF:", sum(minimal_signature %in% genes_in_both), "\n")
cat("  Genes only in LASSO:", sum(minimal_signature %in% lasso_genes & !(minimal_signature %in% rf_top_genes)), "\n")
cat("  Genes only in RF:", sum(minimal_signature %in% rf_top_genes & !(minimal_signature %in% lasso_genes)), "\n")

cat("\nTop 20 signature genes:\n")
print(head(gene_scores %>% select(gene, in_lasso, in_rf, importance, log2fc, p_value), 20),
      row.names = FALSE)

# ============================================================================
# 6. BUILD CLASSIFIER WITH 50-GENE SIGNATURE
# ============================================================================

cat("\n=== STEP 5: BUILDING CLASSIFIER WITH 50-GENE SIGNATURE ===\n\n")

# Prepare data with 50 genes only
expr_signature <- t(expr_matched[minimal_signature, ])
colnames(expr_signature) <- minimal_signature

# Train-test split (70-30)
train_indices <- createDataPartition(cluster_vec, p = 0.7, list = FALSE)
train_data <- data.frame(expr_signature[train_indices, ])
train_data$cluster <- as.factor(cluster_vec[train_indices])
test_data <- data.frame(expr_signature[-train_indices, ])
test_data$cluster <- as.factor(cluster_vec[-train_indices])

cat("Training set:", nrow(train_data), "samples\n")
cat("Test set:", nrow(test_data), "samples\n\n")

# Train Random Forest with 50 genes
cat("Training final Random Forest classifier...\n")
final_rf <- randomForest(cluster ~ .,
                         data = train_data,
                         ntree = 1000,
                         importance = TRUE)

cat("Training OOB error:", round(final_rf$err.rate[1000, "OOB"] * 100, 2), "%\n\n")

# Test set predictions
test_pred <- predict(final_rf, test_data, type = "prob")
test_pred_class <- predict(final_rf, test_data)

# Performance metrics
conf_matrix <- confusionMatrix(test_pred_class, test_data$cluster)
cat("=== TEST SET PERFORMANCE ===\n")
print(conf_matrix)

# ROC curve
roc_obj <- roc(test_data$cluster, test_pred[, 2])
auc_value <- auc(roc_obj)

cat("\nAUC:", round(auc_value, 3), "\n\n")

# ============================================================================
# 7. CROSS-VALIDATION PERFORMANCE
# ============================================================================

cat("=== STEP 6: 10-FOLD CROSS-VALIDATION ===\n\n")

# Prepare full dataset
full_data <- data.frame(expr_signature)
full_data$cluster <- as.factor(cluster_vec)

# 10-fold CV
cat("Running 10-fold cross-validation...\n")
cv_control <- trainControl(method = "cv",
                          number = 10,
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          savePredictions = TRUE)

# Rename levels to valid R names (required by caret)
levels(full_data$cluster) <- c("Cluster1", "Cluster2")

cv_model <- train(cluster ~ .,
                  data = full_data,
                  method = "rf",
                  trControl = cv_control,
                  ntree = 500,
                  metric = "ROC")

cat("\n10-Fold CV Results:\n")
print(cv_model$results)

cat("\nMean CV Accuracy:", round(mean(cv_model$resample$Accuracy), 3), "\n")
cat("Mean CV ROC:", round(mean(cv_model$resample$ROC), 3), "\n")
cat("Mean CV Sensitivity:", round(mean(cv_model$resample$Sens), 3), "\n")
cat("Mean CV Specificity:", round(mean(cv_model$resample$Spec), 3), "\n\n")

# ============================================================================
# 8. SAVE RESULTS
# ============================================================================

cat("=== SAVING RESULTS ===\n\n")

dir.create("03_Results/15_Gene_Signature", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/14_Gene_Signature", showWarnings = FALSE, recursive = TRUE)

# Save signature genes
signature_info <- gene_scores %>%
  filter(gene %in% minimal_signature) %>%
  select(gene, in_lasso, in_rf, importance, log2fc, p_value, score)

write.csv(signature_info,
          "03_Results/15_Gene_Signature/50_gene_signature.csv",
          row.names = FALSE)

# Save full DE results
write.csv(de_results,
          "03_Results/15_Gene_Signature/full_differential_expression.csv",
          row.names = FALSE)

# Save LASSO genes
write.csv(data.frame(gene = lasso_genes),
          "03_Results/15_Gene_Signature/lasso_selected_genes.csv",
          row.names = FALSE)

# Save RF importance
write.csv(importance_df,
          "03_Results/15_Gene_Signature/rf_variable_importance.csv",
          row.names = FALSE)

# Save classifier model
saveRDS(final_rf, "03_Results/15_Gene_Signature/final_rf_classifier.rds")

# Save performance metrics
performance_summary <- data.frame(
  metric = c("Test_Accuracy", "Test_Sensitivity", "Test_Specificity",
             "Test_AUC", "CV_Accuracy", "CV_ROC", "CV_Sensitivity", "CV_Specificity"),
  value = c(
    conf_matrix$overall["Accuracy"],
    conf_matrix$byClass["Sensitivity"],
    conf_matrix$byClass["Specificity"],
    auc_value,
    mean(cv_model$resample$Accuracy),
    mean(cv_model$resample$ROC),
    mean(cv_model$resample$Sens),
    mean(cv_model$resample$Spec)
  )
)

write.csv(performance_summary,
          "03_Results/15_Gene_Signature/classifier_performance.csv",
          row.names = FALSE)

cat("✓ Saved: 50_gene_signature.csv\n")
cat("✓ Saved: full_differential_expression.csv\n")
cat("✓ Saved: lasso_selected_genes.csv\n")
cat("✓ Saved: rf_variable_importance.csv\n")
cat("✓ Saved: final_rf_classifier.rds\n")
cat("✓ Saved: classifier_performance.csv\n\n")

# ============================================================================
# 9. CREATE VISUALIZATIONS
# ============================================================================

cat("=== CREATING VISUALIZATIONS ===\n\n")

# Plot 1: Variable importance for top 20 genes
top20_importance <- signature_info %>%
  arrange(desc(importance)) %>%
  head(20)

p1 <- ggplot(top20_importance, aes(x = reorder(gene, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "#2E9FDF") +
  coord_flip() +
  labs(x = "Gene", y = "Random Forest Importance",
       title = "Top 20 Genes in 50-Gene Signature",
       subtitle = "Ranked by Random Forest Variable Importance") +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("04_Figures/14_Gene_Signature/top20_gene_importance.pdf",
       p1, width = 8, height = 6)

# Plot 2: ROC curve
pdf("04_Figures/14_Gene_Signature/roc_curve_50gene_classifier.pdf",
    width = 7, height = 7)
plot(roc_obj,
     main = paste0("ROC Curve - 50-Gene Signature Classifier\nAUC = ", round(auc_value, 3)),
     col = "#2E9FDF", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")
dev.off()

# Plot 3: Heatmap of 50 genes
library(pheatmap)

# Order samples by cluster
sample_order <- order(cluster_vec)
expr_heatmap <- expr_matched[minimal_signature, sample_order]
cluster_ordered <- cluster_vec[sample_order]

# Annotation
annotation_col <- data.frame(
  Cluster = factor(cluster_ordered, labels = c("Proliferative", "Immune-Inflammatory"))
)
rownames(annotation_col) <- colnames(expr_heatmap)

annotation_colors <- list(
  Cluster = c("Proliferative" = "#2E9FDF", "Immune-Inflammatory" = "#E7B800")
)

pdf("04_Figures/14_Gene_Signature/heatmap_50gene_signature.pdf",
    width = 12, height = 10)
pheatmap(expr_heatmap,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "50-Gene Signature Expression Heatmap",
         fontsize_row = 6)
dev.off()

cat("✓ Saved: top20_gene_importance.pdf\n")
cat("✓ Saved: roc_curve_50gene_classifier.pdf\n")
cat("✓ Saved: heatmap_50gene_signature.pdf\n\n")

# ============================================================================
# 10. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("50-GENE SIGNATURE DEVELOPMENT COMPLETE\n\n")

cat("Signature Composition:\n")
cat("  Total genes:", length(minimal_signature), "\n")
cat("  Selected by both LASSO and RF:", sum(minimal_signature %in% genes_in_both), "\n")
cat("  LASSO only:", sum(minimal_signature %in% lasso_genes & !(minimal_signature %in% rf_top_genes)), "\n")
cat("  RF only:", sum(minimal_signature %in% rf_top_genes & !(minimal_signature %in% lasso_genes)), "\n\n")

cat("Classifier Performance:\n")
cat("  Test Set Accuracy:", round(conf_matrix$overall["Accuracy"] * 100, 1), "%\n")
cat("  Test Set Sensitivity:", round(conf_matrix$byClass["Sensitivity"] * 100, 1), "%\n")
cat("  Test Set Specificity:", round(conf_matrix$byClass["Specificity"] * 100, 1), "%\n")
cat("  Test Set AUC:", round(auc_value, 3), "\n\n")

cat("  10-Fold CV Accuracy:", round(mean(cv_model$resample$Accuracy) * 100, 1), "%\n")
cat("  10-Fold CV ROC:", round(mean(cv_model$resample$ROC), 3), "\n")
cat("  10-Fold CV Sensitivity:", round(mean(cv_model$resample$Sens) * 100, 1), "%\n")
cat("  10-Fold CV Specificity:", round(mean(cv_model$resample$Spec) * 100, 1), "%\n\n")

cat("Top 10 signature genes:\n")
print(head(signature_info %>% select(gene, log2fc, importance), 10), row.names = FALSE)

cat("\n### 50-Gene Signature Task COMPLETE ###\n")
