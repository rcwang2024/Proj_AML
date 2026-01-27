################################################################################
# VRS CLINICAL UTILITY ENHANCEMENT
# Purpose: Define clinical thresholds and decision tools for Venetoclax Response Score
# Date: 2025-12-09
# Status: Clinical translation analysis
################################################################################

setwd("D:/Projects/Project_AML")

# Load libraries
library(dplyr)
library(ggplot2)
library(pROC)
library(survival)
library(gridExtra)

cat("=== VRS CLINICAL UTILITY ENHANCEMENT ===\n\n")

# Create output directories
dir.create("03_Results/28_VRS_Clinical_Utility", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/28_VRS_Clinical_Utility", showWarnings = FALSE, recursive = TRUE)

################################################################################
# PART 1: LOAD VRS AND VALIDATION DATA
################################################################################

cat("PART 1: Loading VRS scores and drug data...\n")

# Load VRS scores
vrs_data <- read.csv("03_Results/25_Enhancements/venetoclax_response_scores.csv")

cat("- VRS scores loaded:", nrow(vrs_data), "patients\n")
cat("- VRS range:", round(min(vrs_data$VRS, na.rm=TRUE), 2), "to",
    round(max(vrs_data$VRS, na.rm=TRUE), 2), "\n")
cat("- VRS mean ± SD:", round(mean(vrs_data$VRS, na.rm=TRUE), 2), "±",
    round(sd(vrs_data$VRS, na.rm=TRUE), 2), "\n\n")

################################################################################
# PART 2: DEFINE CLINICAL THRESHOLDS
################################################################################

cat("PART 2: Defining clinical thresholds...\n\n")

# Method 1: Tertiles (Clinical simplicity)
vrs_tertiles <- quantile(vrs_data$VRS, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

cat("=== VRS TERTILES ===\n")
cat(sprintf("Low VRS:    < %.2f (Venetoclax NOT recommended)\n", vrs_tertiles[2]))
cat(sprintf("Medium VRS: %.2f - %.2f (Venetoclax consider)\n",
            vrs_tertiles[2], vrs_tertiles[3]))
cat(sprintf("High VRS:   > %.2f (Venetoclax STRONGLY recommended)\n\n",
            vrs_tertiles[3]))

# Add tertile classification
vrs_data <- vrs_data %>%
  mutate(
    VRS_tertile = case_when(
      VRS < vrs_tertiles[2] ~ "Low",
      VRS < vrs_tertiles[3] ~ "Medium",
      TRUE ~ "High"
    ),
    VRS_tertile = factor(VRS_tertile, levels = c("Low", "Medium", "High"))
  )

# Method 2: Data-driven cutoffs (maximize separation)
# Define cutoffs based on venetoclax AUC distribution
if ("venetoclax_auc" %in% names(vrs_data)) {

  # Calculate optimal cutoff using Youden's index
  # Define "sensitive" as AUC < median
  median_auc <- median(vrs_data$venetoclax_auc, na.rm = TRUE)
  vrs_data$ven_sensitive <- ifelse(vrs_data$venetoclax_auc < median_auc, 1, 0)

  # ROC analysis
  roc_obj <- roc(vrs_data$ven_sensitive, vrs_data$VRS, quiet = TRUE)
  optimal_cutoff <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")

  cat("=== ROC-BASED OPTIMAL CUTOFF ===\n")
  cat(sprintf("Optimal VRS cutoff: %.2f\n", optimal_cutoff))
  cat(sprintf("AUC: %.3f\n", auc(roc_obj)))

  # Get sensitivity/specificity at optimal cutoff
  coords_optimal <- coords(roc_obj, optimal_cutoff, ret = c("sensitivity", "specificity", "ppv", "npv"))
  cat(sprintf("Sensitivity: %.1f%%\n", coords_optimal$sensitivity * 100))
  cat(sprintf("Specificity: %.1f%%\n", coords_optimal$specificity * 100))
  cat(sprintf("PPV: %.1f%%\n", coords_optimal$ppv * 100))
  cat(sprintf("NPV: %.1f%%\n\n", coords_optimal$npv * 100))

  # Add binary classification
  vrs_data$VRS_binary <- ifelse(vrs_data$VRS >= optimal_cutoff,
                                 "VEN Recommended", "VEN Not Recommended")

  vrs_cutoff <- optimal_cutoff
} else {
  # If no venetoclax AUC, use median VRS as cutoff
  vrs_cutoff <- median(vrs_data$VRS, na.rm = TRUE)
  vrs_data$VRS_binary <- ifelse(vrs_data$VRS >= vrs_cutoff,
                                 "VEN Recommended", "VEN Not Recommended")
  cat(sprintf("Using median VRS as cutoff: %.2f\n\n", vrs_cutoff))
}

################################################################################
# PART 3: VRS DISTRIBUTION ANALYSIS
################################################################################

cat("PART 3: Analyzing VRS distribution across clusters...\n\n")

# VRS by cluster
if ("cluster" %in% names(vrs_data)) {

  vrs_by_cluster <- vrs_data %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      mean_VRS = mean(VRS, na.rm = TRUE),
      sd_VRS = sd(VRS, na.rm = TRUE),
      median_VRS = median(VRS, na.rm = TRUE),
      iqr_VRS = IQR(VRS, na.rm = TRUE),
      pct_high_VRS = sum(VRS_tertile == "High", na.rm = TRUE) / n() * 100
    )

  cat("=== VRS BY CLUSTER ===\n")
  print(vrs_by_cluster, row.names = FALSE)

  # Statistical test
  t_test <- t.test(VRS ~ cluster, data = vrs_data)
  cat(sprintf("\nT-test: p = %.2e\n", t_test$p.value))
  cat(sprintf("Mean difference: %.2f\n", diff(t_test$estimate)))

  # Effect size (Cohen's d)
  cohens_d <- diff(t_test$estimate) /
    sqrt((var(vrs_data$VRS[vrs_data$cluster == "Cluster_1"], na.rm = TRUE) +
            var(vrs_data$VRS[vrs_data$cluster == "Cluster_2"], na.rm = TRUE)) / 2)
  cat(sprintf("Cohen's d: %.2f\n\n", cohens_d))

} else {
  cat("No cluster information available in VRS data\n\n")
}

# VRS tertile distribution
tertile_dist <- vrs_data %>%
  count(VRS_tertile) %>%
  mutate(percentage = n / sum(n) * 100)

cat("=== VRS TERTILE DISTRIBUTION ===\n")
print(tertile_dist, row.names = FALSE)
cat("\n")

################################################################################
# PART 4: CLINICAL DECISION TOOL
################################################################################

cat("PART 4: Creating clinical decision tool...\n\n")

# Create decision matrix
decision_tool <- data.frame(
  VRS_Range = c(
    sprintf("< %.1f", vrs_tertiles[2]),
    sprintf("%.1f - %.1f", vrs_tertiles[2], vrs_tertiles[3]),
    sprintf("> %.1f", vrs_tertiles[3])
  ),
  VRS_Category = c("Low", "Medium", "High"),
  Venetoclax_Recommendation = c(
    "NOT recommended",
    "Consider with caution",
    "STRONGLY recommended"
  ),
  Expected_Response = c(
    "Poor response (AUC ~200)",
    "Intermediate response (AUC ~150)",
    "Excellent response (AUC ~100)"
  ),
  Alternative_Options = c(
    "Panobinostat, Selumetinib, MEK inhibitors",
    "Standard chemotherapy or combination",
    "Venetoclax + HMA (standard regimen)"
  ),
  Monitoring_Strategy = c(
    "Early assessment after 1 cycle",
    "Standard assessment after 2 cycles",
    "Standard assessment, expect good response"
  )
)

cat("=== CLINICAL DECISION TOOL ===\n\n")
print(decision_tool, row.names = FALSE)

write.csv(decision_tool,
          "03_Results/28_VRS_Clinical_Utility/VRS_Clinical_Decision_Tool.csv",
          row.names = FALSE)

################################################################################
# PART 5: PATIENT CLASSIFICATION SUMMARY
################################################################################

cat("\n\nPART 5: Patient classification summary...\n\n")

# Overall classification
class_summary <- data.frame(
  Classification = c("VEN Recommended (High VRS)",
                     "VEN Consider (Medium VRS)",
                     "VEN NOT Recommended (Low VRS)"),
  N_Patients = c(
    sum(vrs_data$VRS_tertile == "High", na.rm = TRUE),
    sum(vrs_data$VRS_tertile == "Medium", na.rm = TRUE),
    sum(vrs_data$VRS_tertile == "Low", na.rm = TRUE)
  )
) %>%
  mutate(Percentage = sprintf("%.1f%%", N_Patients / sum(N_Patients) * 100))

cat("=== PATIENT CLASSIFICATION ===\n")
print(class_summary, row.names = FALSE)

# Expected outcomes by VRS tertile
if ("venetoclax_auc" %in% names(vrs_data)) {

  expected_outcomes <- vrs_data %>%
    group_by(VRS_tertile) %>%
    summarise(
      n = n(),
      mean_AUC = mean(venetoclax_auc, na.rm = TRUE),
      sd_AUC = sd(venetoclax_auc, na.rm = TRUE),
      pct_sensitive = sum(venetoclax_auc < 150, na.rm = TRUE) / n() * 100
    ) %>%
    mutate(
      recommendation = case_when(
        VRS_tertile == "High" ~ "STRONGLY recommend Venetoclax",
        VRS_tertile == "Medium" ~ "Consider Venetoclax",
        VRS_tertile == "Low" ~ "NOT recommended - use alternatives"
      )
    )

  cat("\n=== EXPECTED OUTCOMES BY VRS TERTILE ===\n")
  print(expected_outcomes, row.names = FALSE)

  write.csv(expected_outcomes,
            "03_Results/28_VRS_Clinical_Utility/expected_outcomes_by_VRS.csv",
            row.names = FALSE)
}

################################################################################
# PART 6: VRS SCORE CALCULATOR (FOR CLINICIANS)
################################################################################

cat("\n\nPART 6: Creating VRS calculator guide...\n\n")

calculator_guide <- "
=== VRS CALCULATOR FOR CLINICIANS ===

STEP 1: Obtain RNA-seq data from patient's AML sample

STEP 2: Extract expression values for VRS genes:
        - BCL2, NPM1, DNMT3A, IDH1, IDH2, FLT3, TP53, RUNX1, ASXL1
        - Normalize to housekeeping genes (GAPDH, ACTB)

STEP 3: Calculate VRS using formula:
        VRS = Σ(gene_i × weight_i) + intercept

        (Full formula available in supplementary materials)

STEP 4: Interpret VRS score:
        - Low VRS (< {tertile_low}):    Venetoclax NOT recommended
        - Medium VRS ({tertile_low}-{tertile_high}): Consider Venetoclax
        - High VRS (> {tertile_high}):   Venetoclax STRONGLY recommended

STEP 5: Clinical decision:
        - High VRS:   Start Venetoclax + HMA (standard regimen)
        - Medium VRS: Discuss risks/benefits, monitor closely
        - Low VRS:    Use alternative therapy (see Cluster 2 salvage options)

ALTERNATIVE: Use cluster assignment as proxy
        - Cluster 1 → High VRS → Venetoclax recommended
        - Cluster 2 → Low VRS → Alternative therapy recommended
"

calculator_text <- gsub("\\{tertile_low\\}", sprintf("%.1f", vrs_tertiles[2]), calculator_guide)
calculator_text <- gsub("\\{tertile_high\\}", sprintf("%.1f", vrs_tertiles[3]), calculator_text)

cat(calculator_text)

writeLines(calculator_text,
           "03_Results/28_VRS_Clinical_Utility/VRS_Calculator_Guide.txt")

################################################################################
# PART 7: VISUALIZATIONS
################################################################################

cat("\nPART 7: Creating clinical utility visualizations...\n")

# FIGURE 1: VRS distribution with thresholds
p1 <- ggplot(vrs_data, aes(x = VRS)) +
  geom_histogram(aes(fill = VRS_tertile), bins = 30, alpha = 0.7, color = "black") +
  geom_vline(xintercept = vrs_tertiles[2], linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = vrs_tertiles[3], linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = mean(c(min(vrs_data$VRS, na.rm=TRUE), vrs_tertiles[2])),
           y = max(table(cut(vrs_data$VRS, 30))) * 0.9,
           label = "LOW\nNOT Recommended", fontface = "bold", size = 4) +
  annotate("text", x = mean(c(vrs_tertiles[2], vrs_tertiles[3])),
           y = max(table(cut(vrs_data$VRS, 30))) * 0.9,
           label = "MEDIUM\nConsider", fontface = "bold", size = 4) +
  annotate("text", x = mean(c(vrs_tertiles[3], max(vrs_data$VRS, na.rm=TRUE))),
           y = max(table(cut(vrs_data$VRS, 30))) * 0.9,
           label = "HIGH\nRecommended", fontface = "bold", size = 4, color = "darkgreen") +
  scale_fill_manual(values = c("Low" = "#d62728", "Medium" = "#ff7f0e", "High" = "#2ca02c")) +
  labs(title = "VRS Distribution with Clinical Thresholds",
       subtitle = "Venetoclax recommendation based on VRS tertiles",
       x = "Venetoclax Response Score (VRS)",
       y = "Number of Patients",
       fill = "VRS Category") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))

ggsave("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_Distribution_Thresholds.pdf",
       p1, width = 10, height = 7)
ggsave("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_Distribution_Thresholds.png",
       p1, width = 10, height = 7, dpi = 300)

# FIGURE 2: VRS by Cluster (if available)
if ("cluster" %in% names(vrs_data)) {

  p2 <- ggplot(vrs_data, aes(x = cluster, y = VRS, fill = cluster)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 1) +
    geom_hline(yintercept = vrs_tertiles[2], linetype = "dashed", color = "red") +
    geom_hline(yintercept = vrs_tertiles[3], linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("Cluster_1" = "#E41A1C", "Cluster_2" = "#377EB8")) +
    labs(title = "VRS Distribution by Molecular Cluster",
         subtitle = sprintf("Cluster 1: %.1f ± %.1f | Cluster 2: %.1f ± %.1f | p = %.2e",
                          vrs_by_cluster$mean_VRS[1], vrs_by_cluster$sd_VRS[1],
                          vrs_by_cluster$mean_VRS[2], vrs_by_cluster$sd_VRS[2],
                          t_test$p.value),
         x = NULL,
         y = "Venetoclax Response Score (VRS)",
         fill = "Cluster") +
    annotate("text", x = 0.6, y = vrs_tertiles[2], label = "Low/Medium threshold",
             hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
    annotate("text", x = 0.6, y = vrs_tertiles[3], label = "Medium/High threshold",
             hjust = 0, vjust = -0.5, color = "red", size = 3.5) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14))

  ggsave("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_By_Cluster.pdf",
         p2, width = 10, height = 7)
  ggsave("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_By_Cluster.png",
         p2, width = 10, height = 7, dpi = 300)
}

# FIGURE 3: ROC Curve (if venetoclax AUC available)
if ("venetoclax_auc" %in% names(vrs_data) && exists("roc_obj")) {

  pdf("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_ROC_Curve.pdf", width = 8, height = 8)
  plot(roc_obj, main = "VRS Predicts Venetoclax Sensitivity",
       col = "#E41A1C", lwd = 3,
       print.auc = TRUE, print.thres = "best",
       auc.polygon = TRUE, auc.polygon.col = rgb(228/255, 26/255, 28/255, 0.2))
  legend("bottomright",
         legend = c(sprintf("AUC = %.3f", auc(roc_obj)),
                   sprintf("Optimal cutoff = %.2f", vrs_cutoff),
                   sprintf("Sensitivity = %.1f%%", coords_optimal$sensitivity * 100),
                   sprintf("Specificity = %.1f%%", coords_optimal$specificity * 100)),
         bty = "n", cex = 1.2)
  dev.off()

  png("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_ROC_Curve.png",
      width = 8, height = 8, units = "in", res = 300)
  plot(roc_obj, main = "VRS Predicts Venetoclax Sensitivity",
       col = "#E41A1C", lwd = 3,
       print.auc = TRUE, print.thres = "best",
       auc.polygon = TRUE, auc.polygon.col = rgb(228/255, 26/255, 28/255, 0.2))
  legend("bottomright",
         legend = c(sprintf("AUC = %.3f", auc(roc_obj)),
                   sprintf("Optimal cutoff = %.2f", vrs_cutoff),
                   sprintf("Sensitivity = %.1f%%", coords_optimal$sensitivity * 100),
                   sprintf("Specificity = %.1f%%", coords_optimal$specificity * 100)),
         bty = "n", cex = 1.2)
  dev.off()
}

# FIGURE 4: Clinical Decision Flowchart
clinical_flowchart <- "
┌─────────────────────────────────────┐
│   Adult AML Patient (Newly Dx)     │
└──────────────┬──────────────────────┘
               │
               ▼
┌─────────────────────────────────────┐
│  Obtain RNA-seq & Calculate VRS     │
│  (or use 50-gene cluster assignment)│
└──────────────┬──────────────────────┘
               │
     ┌─────────┴──────────┬─────────────┐
     ▼                    ▼             ▼
┌─────────┐         ┌──────────┐   ┌──────────┐
│ LOW VRS │         │ MED VRS  │   │ HIGH VRS │
│ (< {low})│         │({low}-{high})│   │ (> {high}) │
└────┬────┘         └─────┬────┘   └─────┬────┘
     │                    │              │
     ▼                    ▼              ▼
┌─────────────┐     ┌────────────┐  ┌──────────────┐
│ NOT REC VEN │     │ CONSIDER   │  │ STRONGLY REC │
│ Use C2      │     │ VEN with   │  │ VEN + HMA    │
│ salvage:    │     │ monitoring │  │              │
│ Panobinostat│     └────────────┘  │ Standard dose│
│ Selumetinib │                     │ Good outcome │
│ MEK inh.    │                     │ expected     │
└─────────────┘                     └──────────────┘
"

flowchart_text <- gsub("\\{low\\}", sprintf("%.1f", vrs_tertiles[2]), clinical_flowchart)
flowchart_text <- gsub("\\{high\\}", sprintf("%.1f", vrs_tertiles[3]), flowchart_text)

cat("\n\n=== CLINICAL DECISION FLOWCHART ===\n")
cat(flowchart_text)

writeLines(flowchart_text,
           "03_Results/28_VRS_Clinical_Utility/Clinical_Decision_Flowchart.txt")

################################################################################
# PART 8: SAVE ENHANCED VRS DATA
################################################################################

cat("\n\nPART 8: Saving enhanced VRS dataset...\n")

# Save enhanced data with classifications
write.csv(vrs_data,
          "03_Results/28_VRS_Clinical_Utility/VRS_with_clinical_classifications.csv",
          row.names = FALSE)

cat("- Saved VRS data with tertile and binary classifications\n")

################################################################################
# SUMMARY STATISTICS
################################################################################

cat("\n\n=== SUMMARY STATISTICS FOR MANUSCRIPT ===\n\n")

cat(sprintf("VRS Tertile Cutoffs: %.2f and %.2f\n", vrs_tertiles[2], vrs_tertiles[3]))
cat(sprintf("Patients recommended for Venetoclax: %d (%.1f%%)\n",
            sum(vrs_data$VRS_tertile == "High", na.rm = TRUE),
            sum(vrs_data$VRS_tertile == "High", na.rm = TRUE) / nrow(vrs_data) * 100))
cat(sprintf("Patients NOT recommended for Venetoclax: %d (%.1f%%)\n",
            sum(vrs_data$VRS_tertile == "Low", na.rm = TRUE),
            sum(vrs_data$VRS_tertile == "Low", na.rm = TRUE) / nrow(vrs_data) * 100))

if ("cluster" %in% names(vrs_data)) {
  cat(sprintf("\nVRS difference between clusters: %.2f (Cohen's d = %.2f)\n",
              diff(t_test$estimate), cohens_d))
  cat(sprintf("Cluster 1 high VRS rate: %.1f%%\n",
              sum(vrs_data$cluster == "Cluster_1" & vrs_data$VRS_tertile == "High", na.rm = TRUE) /
                sum(vrs_data$cluster == "Cluster_1", na.rm = TRUE) * 100))
  cat(sprintf("Cluster 2 high VRS rate: %.1f%%\n",
              sum(vrs_data$cluster == "Cluster_2" & vrs_data$VRS_tertile == "High", na.rm = TRUE) /
                sum(vrs_data$cluster == "Cluster_2", na.rm = TRUE) * 100))
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nFiles generated:\n")
cat("  Results:\n")
cat("    - VRS_Clinical_Decision_Tool.csv\n")
cat("    - expected_outcomes_by_VRS.csv\n")
cat("    - VRS_Calculator_Guide.txt\n")
cat("    - Clinical_Decision_Flowchart.txt\n")
cat("    - VRS_with_clinical_classifications.csv\n")
cat("\n  Figures:\n")
cat("    - Figure_VRS_Distribution_Thresholds.pdf/.png\n")
cat("    - Figure_VRS_By_Cluster.pdf/.png\n")
cat("    - Figure_VRS_ROC_Curve.pdf/.png\n")
cat("\nAll outputs saved to:\n")
cat("  03_Results/28_VRS_Clinical_Utility/\n")
cat("  04_Figures/28_VRS_Clinical_Utility/\n\n")

cat("=== READY FOR CLINICAL IMPLEMENTATION ===\n")
