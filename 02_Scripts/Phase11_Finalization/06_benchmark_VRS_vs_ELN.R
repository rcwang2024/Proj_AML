# R script to benchmark VRS against ELN 2017 clinical risk categories
setwd("d:/Proj_AML")
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)

cat("=== BENCHMARKING VRS VS ELN 2017 RISK CLASSIFICATION ===\n\n")

# Setup Directories
dir.create("03_Results/28_VRS_Clinical_Utility", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/Phase11_Finalization", recursive = TRUE, showWarnings = FALSE)

# 1. Load data
cat("Loading clinical, cluster, and VRS data...\n")
clinical <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
vrs_data <- read_csv("03_Results/25_Enhancements/vrs_9gene_scores.csv", show_col_types = FALSE)
monocytic_results <- read_csv("03_Results/Phase10_Analysis/10_1_Monocytic_Mapping_Results.csv", show_col_types = FALSE)

# 2. Merge and filter
merged_df <- clinical %>%
  filter(!is.na(dbgap_rnaseq_sample)) %>%
  inner_join(vrs_data, by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  inner_join(monocytic_results %>% select(sample_id, auc, cluster), by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  select(sample_id = dbgap_rnaseq_sample,
         ELN2017,
         VRS9,
         ven_auc = auc,
         cluster) %>%
  filter(ELN2017 %in% c("Favorable", "Intermediate", "Adverse")) %>%  # Keep only standard ELN categories
  filter(!is.na(VRS9)) %>%      # Remove samples without VRS
  filter(!is.na(ven_auc))       # Remove samples without drug response

cat("Matched samples with clinical, VRS, and drug response:", nrow(merged_df), "\n\n")

# Check distribution of ELN risk
cat("Samples by ELN 2017 category:\n")
print(table(merged_df$ELN2017))
cat("\n")

# === 3. COMPARE VRS SCORES ACROSS ELN GROUPS ===
cat("=== 3. COMPARISON OF VRS SCORES ACROSS ELN GROUPS ===\n")
kruskal_vrs_eln <- kruskal.test(VRS9 ~ ELN2017, data = merged_df)
cat(sprintf("Kruskal-Wallis p-value for VRS9 by ELN risk: %.2e\n\n", kruskal_vrs_eln$p.value))

# === 4. HIERARCHICAL LINEAR REGRESSION ===
cat("=== 4. HIERARCHICAL REGRESSION FOR INDEPENDENT PREDICTIVE VALUE ===\n")
# Model 1: ELN alone predicting Venetoclax AUC
fit_eln_only <- lm(ven_auc ~ factor(ELN2017), data = merged_df)
r2_eln_only <- summary(fit_eln_only)$r.squared
cat(sprintf("Model 1 (ELN only): R² = %.4f\n", r2_eln_only))

# Model 2: ELN + VRS9 predicting Venetoclax AUC
fit_combined <- lm(ven_auc ~ factor(ELN2017) + VRS9, data = merged_df)
r2_combined <- summary(fit_combined)$r.squared
cat(sprintf("Model 2 (ELN + VRS9): R² = %.4f\n", r2_combined))

# F-test for model comparison
anova_res <- anova(fit_eln_only, fit_combined)
lrt_p <- anova_res$`Pr(>F)`[2]
r2_diff <- r2_combined - r2_eln_only
pct_improve <- (r2_diff / r2_eln_only) * 100

cat(sprintf("VRS9 Incremental Explained Variance (ΔR²): +%.4f\n", r2_diff))
cat(sprintf("VRS9 Percentage R² Improvement beyond ELN: +%.1f%%\n", pct_improve))
cat(sprintf("F-test p-value: %.2e\n\n", lrt_p))

# Save regression stats
write_csv(
  data.frame(
    metric = c("ELN_R2", "Combined_R2", "Delta_R2", "Percent_Improvement", "F_test_p"),
    value = c(r2_eln_only, r2_combined, r2_diff, pct_improve, lrt_p)
  ),
  "03_Results/28_VRS_Clinical_Utility/eln_benchmarking_regression.csv"
)

# === 5. STRATIFIED SUBGROUP CORRELATION ANALYSIS ===
cat("=== 5. SUBGROUP CORRELATION ANALYSIS (VRS9 VS VEN_AUC WITHIN ELN CATEGORIES) ===\n")
subgroup_correlations <- merged_df %>%
  group_by(ELN2017) %>%
  summarise(
    N = n(),
    spearman_rho = cor(VRS9, ven_auc, method = "spearman"),
    spearman_p = cor.test(VRS9, ven_auc, method = "spearman")$p.value,
    pearson_r = cor(VRS9, ven_auc, method = "pearson"),
    pearson_p = cor.test(VRS9, ven_auc, method = "pearson")$p.value,
    .groups = "drop"
  )

print(subgroup_correlations)
cat("\n")

write_csv(subgroup_correlations, "03_Results/28_VRS_Clinical_Utility/eln_subgroup_correlations.csv")

# === 6. VISUALIZATIONS ===
cat("=== 6. GENERATING PLOTS ===\n")

# Colors
color_c1 <- "#3498DB" # Soft Blue
color_c2 <- "#E67E22" # Warm Orange

# Custom Theme
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 13, color = "darkblue", margin = margin(b=10)),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    axis.text = element_text(size = 11, face = "plain", color = "black"),
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 11, face = "plain", color = "black"),
    strip.text = element_text(size = 13, face = "bold", color = "black"),
    strip.background = element_rect(fill = "gray95", color = "gray80", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    plot.margin = margin(15, 15, 15, 15)
  )

# Plot 1: Boxplot of VRS9 by ELN Category
p_box <- ggplot(merged_df, aes(x = factor(ELN2017, levels=c("Favorable", "Intermediate", "Adverse")), y = VRS9)) +
  geom_boxplot(aes(fill = ELN2017), outlier.shape = NA, alpha = 0.7, width = 0.4, color = "black") +
  geom_jitter(color = "gray30", width = 0.1, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c("Favorable" = "#2ECC71", "Intermediate" = "#F1C40F", "Adverse" = "#E74C3C")) +
  labs(
    title = "VRS Distribution across ELN 2017 Risk Groups",
    subtitle = sprintf("Kruskal-Wallis p = %.2e", kruskal_vrs_eln$p.value),
    x = "ELN 2017 Risk Classification",
    y = "9-Gene Venetoclax Response Score (VRS)"
  ) +
  theme_hf +
  theme(legend.position = "none")

ggsave("04_Figures/Phase11_Finalization/VRS_distribution_by_ELN.pdf", p_box, width = 6, height = 5.5, device = cairo_pdf)
ggsave("04_Figures/Phase11_Finalization/VRS_distribution_by_ELN.png", p_box, width = 6, height = 5.5, dpi = 300)

# Plot 2: Scatter plot stratified by ELN Category
# Create formatted labels for facets with correlation stats
facet_labels <- subgroup_correlations %>%
  mutate(label = sprintf("%s\n(n = %d, ρ = %.2f\np = %.1e)", ELN2017, N, spearman_rho, spearman_p)) %>%
  select(ELN2017, label)

plot_data <- merged_df %>%
  inner_join(facet_labels, by = "ELN2017") %>%
  mutate(label = factor(label, levels = facet_labels$label[order(match(facet_labels$ELN2017, c("Favorable", "Intermediate", "Adverse")))]))

p_scatter <- ggplot(plot_data, aes(x = VRS9, y = ven_auc)) +
  geom_point(aes(color = factor(cluster, levels=c(1, 2), labels=c("Cluster 1", "Cluster 2"))), size = 2.0, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", size = 1.0, se = TRUE, fill = "gray90") +
  scale_color_manual(values = c("Cluster 1" = color_c1, "Cluster 2" = color_c2), name = "Molecular Subtype") +
  facet_wrap(~ label, scales = "free_x") +
  labs(
    title = "Venetoclax Sensitivity Prediction within ELN 2017 Risk Categories",
    subtitle = sprintf("Incremental Explained Variance ΔR² = +%.4f (F-test p = %.2e)", r2_diff, lrt_p),
    x = "Venetoclax Response Score (VRS, 0-100)",
    y = "Venetoclax ex vivo AUC\n(Lower = More Sensitive)"
  ) +
  theme_hf +
  theme(legend.position = "bottom")

ggsave("04_Figures/Phase11_Finalization/VRS_vs_AUC_stratified_by_ELN.pdf", p_scatter, width = 10, height = 5.5, device = cairo_pdf)
ggsave("04_Figures/Phase11_Finalization/VRS_vs_AUC_stratified_by_ELN.png", p_scatter, width = 10, height = 5.5, dpi = 300)

cat("✓ Saved: 04_Figures/Phase11_Finalization/VRS_distribution_by_ELN.pdf\n")
cat("✓ Saved: 04_Figures/Phase11_Finalization/VRS_vs_AUC_stratified_by_ELN.pdf\n")
cat("\n### ELN Benchmarking Complete ###\n")
