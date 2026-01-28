# TASK: Apply BeatAML Classifier to TCGA-LAML Data
# External validation of molecular subtypes on independent cohort

library(tidyverse)
library(randomForest)
library(survival)
library(survminer)
library(ggplot2)

cat("=== TCGA-LAML EXTERNAL VALIDATION ===\n\n")

# ============================================================================
# 1. LOAD BEATAML CLASSIFIER
# ============================================================================

cat("=== STEP 1: LOADING BEATAML CLASSIFIER ===\n\n")

# Load the trained 50-gene Random Forest classifier
rf_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")
cat("✓ Loaded Random Forest classifier\n")

# Load the 50-gene signature
signature_genes <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")
cat("Signature genes:", nrow(signature_genes), "\n")
cat("First 10 genes:\n")
print(head(signature_genes$gene, 10))
cat("\n")

# ============================================================================
# 2. LOAD TCGA DATA
# ============================================================================

cat("=== STEP 2: LOADING TCGA-LAML DATA ===\n\n")

# Load TCGA expression data
tcga_expr <- readRDS("01_Data/TCGA_LAML/tcga_laml_expression_normalized.rds")
cat("TCGA expression data:", nrow(tcga_expr), "genes x", ncol(tcga_expr), "samples\n")

# Load TCGA clinical data
tcga_clinical <- read.csv("01_Data/TCGA_LAML/tcga_laml_clinical.csv")
cat("TCGA clinical data:", nrow(tcga_clinical), "samples\n\n")

# ============================================================================
# 3. MATCH SIGNATURE GENES
# ============================================================================

cat("=== STEP 3: MATCHING SIGNATURE GENES ===\n\n")

# Check which signature genes are present in TCGA data
genes_in_tcga <- signature_genes$gene %in% rownames(tcga_expr)
cat("Signature genes found in TCGA:", sum(genes_in_tcga), "of", nrow(signature_genes), "\n")

if (sum(genes_in_tcga) < 40) {
  cat("WARNING: Less than 40 signature genes found in TCGA data\n")
  cat("Missing genes:\n")
  print(signature_genes$gene[!genes_in_tcga])
  cat("\n")
}

# Use only genes present in both datasets
common_genes <- intersect(signature_genes$gene, rownames(tcga_expr))
cat("Using", length(common_genes), "common genes for classification\n\n")

# Subset TCGA expression to signature genes
tcga_expr_subset <- tcga_expr[common_genes, ]

# Load BeatAML data to match normalization
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

# Scale TCGA data to match BeatAML distribution for each gene
cat("Scaling TCGA data to match BeatAML distribution...\n")
for (gene in common_genes) {
  if (gene %in% rownames(beataml_expr)) {
    # Get BeatAML gene statistics
    beataml_mean <- mean(beataml_expr[gene, ], na.rm = TRUE)
    beataml_sd <- sd(beataml_expr[gene, ], na.rm = TRUE)

    # Get TCGA gene statistics
    tcga_mean <- mean(tcga_expr_subset[gene, ], na.rm = TRUE)
    tcga_sd <- sd(tcga_expr_subset[gene, ], na.rm = TRUE)

    # Z-score normalize TCGA then scale to BeatAML distribution
    if (tcga_sd > 0) {
      tcga_expr_subset[gene, ] <- ((tcga_expr_subset[gene, ] - tcga_mean) / tcga_sd) * beataml_sd + beataml_mean
    }
  }
}
cat("✓ TCGA data scaled to BeatAML distribution\n\n")

# ============================================================================
# 4. APPLY CLASSIFIER
# ============================================================================

cat("=== STEP 4: APPLYING CLASSIFIER TO TCGA DATA ===\n\n")

# Prepare data for prediction
# Random Forest expects samples as rows, genes as columns
tcga_for_prediction <- as.data.frame(t(tcga_expr_subset))

# Add missing genes with zero values (not detected)
missing_genes <- setdiff(signature_genes$gene, colnames(tcga_for_prediction))
if (length(missing_genes) > 0) {
  cat("Adding", length(missing_genes), "missing genes with zero values (not detected):\n")
  print(missing_genes)
  cat("\n")

  # For each missing gene, add a column with 0 (gene not detected/expressed in TCGA)
  for (gene in missing_genes) {
    tcga_for_prediction[[gene]] <- 0
  }
}

# Ensure columns are in the same order as signature genes
tcga_for_prediction <- tcga_for_prediction[, signature_genes$gene]

# Make sure column names match (remove any special characters)
colnames(tcga_for_prediction) <- make.names(colnames(tcga_for_prediction))

cat("Prediction data dimensions:", nrow(tcga_for_prediction), "samples x",
    ncol(tcga_for_prediction), "genes\n")

# Apply classifier
cat("Applying Random Forest classifier...\n")
predictions <- predict(rf_classifier, newdata = tcga_for_prediction, type = "response")
prediction_probs <- predict(rf_classifier, newdata = tcga_for_prediction, type = "prob")

cat("✓ Classification complete\n\n")

# Create results dataframe
tcga_predictions <- data.frame(
  sample_id = rownames(tcga_for_prediction),
  predicted_cluster = predictions,
  prob_cluster1 = prediction_probs[, 1],
  prob_cluster2 = prediction_probs[, 2],
  confidence = apply(prediction_probs, 1, max),
  stringsAsFactors = FALSE
)

# Summary
cat("TCGA-LAML Classification Results:\n")
cat("  Cluster 1 (Proliferative):", sum(predictions == 1), "samples (",
    round(sum(predictions == 1) / length(predictions) * 100, 1), "%)\n")
cat("  Cluster 2 (Immune-Inflammatory):", sum(predictions == 2), "samples (",
    round(sum(predictions == 2) / length(predictions) * 100, 1), "%)\n")
cat("  Mean classification confidence:", round(mean(tcga_predictions$confidence), 3), "\n\n")

# ============================================================================
# 5. MERGE WITH CLINICAL DATA
# ============================================================================

cat("=== STEP 5: MERGING WITH CLINICAL DATA ===\n\n")

# Merge predictions with clinical data
tcga_merged <- tcga_clinical %>%
  left_join(tcga_predictions, by = "sample_id")

cat("Merged data:", nrow(tcga_merged), "samples\n")
cat("Samples with survival data:", sum(!is.na(tcga_merged$OS_days)), "\n\n")

# ============================================================================
# 6. SURVIVAL ANALYSIS
# ============================================================================

cat("=== STEP 6: SURVIVAL ANALYSIS ===\n\n")

dir.create("03_Results/17_TCGA_Validation", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/16_TCGA_Validation", showWarnings = FALSE, recursive = TRUE)

# Filter samples with survival data
tcga_survival <- tcga_merged %>%
  filter(!is.na(OS_days) & !is.na(predicted_cluster))

cat("Samples for survival analysis:", nrow(tcga_survival), "\n")
cat("  Cluster 1:", sum(tcga_survival$predicted_cluster == 1), "\n")
cat("  Cluster 2:", sum(tcga_survival$predicted_cluster == 2), "\n")
cat("  Events:", sum(tcga_survival$OS_event), "\n\n")

# Create survival object
surv_obj <- Surv(time = tcga_survival$OS_months, event = tcga_survival$OS_event)

# Kaplan-Meier analysis
cat("Running Kaplan-Meier analysis...\n")
km_fit <- survfit(surv_obj ~ predicted_cluster, data = tcga_survival)

# Log-rank test
km_test <- survdiff(surv_obj ~ predicted_cluster, data = tcga_survival)
km_pvalue <- 1 - pchisq(km_test$chisq, df = 1)

cat("Log-rank test p-value:", format(km_pvalue, scientific = TRUE, digits = 3), "\n")

# Median survival times
surv_summary <- summary(km_fit)$table
cat("\nMedian survival times:\n")
print(surv_summary[, c("median", "0.95LCL", "0.95UCL")])
cat("\n")

# Cox proportional hazards
cat("Running Cox proportional hazards model...\n")
cox_model <- coxph(surv_obj ~ predicted_cluster, data = tcga_survival)
cox_summary <- summary(cox_model)

hr <- cox_summary$conf.int[1, 1]
hr_lower <- cox_summary$conf.int[1, 3]
hr_upper <- cox_summary$conf.int[1, 4]
cox_pvalue <- cox_summary$coefficients[1, 5]

cat("Hazard Ratio (Cluster 2 vs Cluster 1):", round(hr, 3),
    "(95% CI:", round(hr_lower, 3), "-", round(hr_upper, 3), ")\n")
cat("Cox model p-value:", format(cox_pvalue, scientific = TRUE, digits = 3), "\n\n")

# Save results
survival_results <- data.frame(
  cohort = "TCGA-LAML",
  n_samples = nrow(tcga_survival),
  n_cluster1 = sum(tcga_survival$predicted_cluster == 1),
  n_cluster2 = sum(tcga_survival$predicted_cluster == 2),
  n_events = sum(tcga_survival$OS_event),
  logrank_pvalue = km_pvalue,
  hazard_ratio = hr,
  hr_95_lower = hr_lower,
  hr_95_upper = hr_upper,
  cox_pvalue = cox_pvalue
)

write.csv(survival_results,
          "03_Results/17_TCGA_Validation/tcga_survival_validation.csv",
          row.names = FALSE)
cat("✓ Saved: tcga_survival_validation.csv\n\n")

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

cat("=== STEP 7: CREATING VISUALIZATIONS ===\n\n")

# Kaplan-Meier plot
tcga_survival$cluster_label <- factor(tcga_survival$predicted_cluster,
                                      labels = c("Proliferative", "Immune-Inflammatory"))

p_km <- ggsurvplot(
  km_fit,
  data = tcga_survival,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.3,
  legend.title = "Molecular Subtype",
  legend.labs = c("Proliferative (Cluster 1)", "Immune-Inflammatory (Cluster 2)"),
  palette = c("#2E9FDF", "#E7B800"),
  title = "TCGA-LAML External Validation: Overall Survival",
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  ggtheme = theme_bw()
)

pdf("04_Figures/16_TCGA_Validation/tcga_kaplan_meier.pdf", width = 10, height = 8)
print(p_km)
dev.off()

cat("✓ Saved: tcga_kaplan_meier.pdf\n")

# Classification confidence distribution
p_conf <- ggplot(tcga_predictions, aes(x = factor(predicted_cluster), y = confidence,
                                        fill = factor(predicted_cluster))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("1" = "#2E9FDF", "2" = "#E7B800"),
                    labels = c("Proliferative", "Immune-Inflammatory")) +
  labs(x = "Predicted Cluster",
       y = "Classification Confidence",
       title = "TCGA-LAML Classification Confidence",
       fill = "Cluster") +
  theme_bw() +
  theme(text = element_text(size = 12))

ggsave("04_Figures/16_TCGA_Validation/classification_confidence.pdf",
       p_conf, width = 8, height = 6)

cat("✓ Saved: classification_confidence.pdf\n")

# Cluster proportion comparison
beataml_props <- data.frame(
  cohort = "BeatAML",
  cluster = c("Cluster 1", "Cluster 2"),
  proportion = c(320/707, 387/707)
)

tcga_props <- data.frame(
  cohort = "TCGA-LAML",
  cluster = c("Cluster 1", "Cluster 2"),
  proportion = c(
    sum(predictions == 1) / length(predictions),
    sum(predictions == 2) / length(predictions)
  )
)

props_combined <- rbind(beataml_props, tcga_props)

p_props <- ggplot(props_combined, aes(x = cohort, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Cluster 1" = "#2E9FDF", "Cluster 2" = "#E7B800")) +
  labs(x = "Cohort", y = "Proportion", title = "Cluster Proportions: BeatAML vs TCGA",
       fill = "Molecular Subtype") +
  theme_bw() +
  theme(text = element_text(size = 12)) +
  geom_text(aes(label = paste0(round(proportion * 100, 1), "%")),
            position = position_dodge(width = 0.9), vjust = -0.5)

ggsave("04_Figures/16_TCGA_Validation/cluster_proportions_comparison.pdf",
       p_props, width = 8, height = 6)

cat("✓ Saved: cluster_proportions_comparison.pdf\n\n")

# ============================================================================
# 8. SAVE ALL RESULTS
# ============================================================================

cat("=== STEP 8: SAVING RESULTS ===\n\n")

# Save predictions
write.csv(tcga_predictions,
          "03_Results/17_TCGA_Validation/tcga_sample_predictions.csv",
          row.names = FALSE)
cat("✓ Saved: tcga_sample_predictions.csv\n")

# Save merged clinical + predictions
write.csv(tcga_merged,
          "03_Results/17_TCGA_Validation/tcga_clinical_with_predictions.csv",
          row.names = FALSE)
cat("✓ Saved: tcga_clinical_with_predictions.csv\n\n")

# ============================================================================
# 9. SUMMARY
# ============================================================================

cat("=== SUMMARY ===\n\n")

cat("TCGA-LAML EXTERNAL VALIDATION COMPLETE\n\n")

cat("Cohort Comparison:\n")
cat("  BeatAML: 707 samples, Cluster 1: 320 (45.3%), Cluster 2: 387 (54.7%)\n")
cat("  TCGA:   ", length(predictions), "samples, Cluster 1:",
    sum(predictions == 1), "(",
    round(sum(predictions == 1) / length(predictions) * 100, 1), "%), Cluster 2:",
    sum(predictions == 2), "(",
    round(sum(predictions == 2) / length(predictions) * 100, 1), "%)\n\n")

cat("Classification Performance:\n")
cat("  Mean confidence:", round(mean(tcga_predictions$confidence), 3), "\n")
cat("  Genes used:", length(common_genes), "of", nrow(signature_genes), "\n\n")

cat("Survival Validation:\n")
cat("  Log-rank p-value:", format(km_pvalue, scientific = TRUE, digits = 3), "\n")
cat("  Hazard Ratio:", round(hr, 3), "(95% CI:", round(hr_lower, 3), "-", round(hr_upper, 3), ")\n")
if (km_pvalue < 0.05) {
  cat("  ✓ Subtypes significantly associated with survival in TCGA\n")
} else {
  cat("  ✗ No significant survival difference in TCGA\n")
}
cat("\n")

cat("Interpretation:\n")
if (km_pvalue < 0.05 && hr > 1) {
  cat("  The BeatAML molecular subtypes successfully validated in TCGA-LAML.\n")
  cat("  Cluster 2 (Immune-Inflammatory) shows worse survival, consistent with BeatAML.\n")
} else if (km_pvalue < 0.05 && hr < 1) {
  cat("  Subtypes are prognostic in TCGA but direction differs from BeatAML.\n")
  cat("  Further investigation needed.\n")
} else {
  cat("  No significant survival difference in TCGA-LAML.\n")
  cat("  Possible reasons: smaller cohort, different treatments, or population differences.\n")
}

cat("\n### TCGA Validation Task COMPLETE ###\n")
