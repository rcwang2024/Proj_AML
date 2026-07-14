# R script to generate Figure 5 Panel E (VIALE-A Clinical Trial Validation)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)

cat("=== GENERATING FIGURE 5 PANEL E (VIALE-A VALIDATION) ===\n")

# Load data
trial_results <- read_csv("03_Results/Phase10_Analysis/10_4_Trial_Validation_Results.csv")
vrs_scores <- read_csv("03_Results/25_Enhancements/vrs_9gene_scores.csv")

# Merge data
merged_data <- inner_join(trial_results, vrs_scores, by="sample_id") %>%
  mutate(cluster = factor(cluster, levels=c(1,2), labels=c("Cluster 1 (Sensitive)", "Cluster 2 (Resistant)")))

# Calculate correlations
spearman_res <- cor.test(merged_data$VRS9, merged_data$TRIAL_SCORE, method="spearman")
pearson_res <- cor.test(merged_data$VRS9, merged_data$TRIAL_SCORE, method="pearson")

cat(sprintf("Spearman correlation: rho = %.3f (p-value = %.2e)\n", spearman_res$estimate, spearman_res$p.value))
cat(sprintf("Pearson correlation: r = %.3f (p-value = %.2e)\n", pearson_res$estimate, pearson_res$p.value))

# Save stats to CSV for record
write_csv(
  data.frame(
    metric = c("Spearman_rho", "Spearman_p", "Pearson_r", "Pearson_p"),
    value = c(spearman_res$estimate, spearman_res$p.value, pearson_res$estimate, pearson_res$p.value)
  ),
  "03_Results/Phase10_Analysis/10_4_TSR_VRS_Correlation.csv"
)

# Colors matching Figure 1 and Panel A
color_c1 <- "#3498DB" # Soft Blue
color_c2 <- "#E67E22" # Warm Orange

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

# Build Scatter Plot
p5e <- ggplot(merged_data, aes(x = VRS9, y = TRIAL_SCORE)) +
  geom_point(aes(fill = cluster), shape = 21, color = "black", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 1.2, se = TRUE, fill = "gray80") +
  scale_fill_manual(values = c("Cluster 1 (Sensitive)" = color_c1, "Cluster 2 (Resistant)" = color_c2)) +
  scale_color_manual(values = c("Cluster 1 (Sensitive)" = color_c1, "Cluster 2 (Resistant)" = color_c2)) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Venetoclax Response Score (VRS, 0-100)",
    y = "VIALE-A Trial Gene Signature Score\n(Resistant = Higher)"
  ) +
  theme_hf +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size=4)))

# Save as vector PDF
output_pdf <- "04_Figures/29_ExternalValidation/Figure5E_VIALE_Validation.pdf"
ggsave(output_pdf, p5e, width = 5.8, height = 6.9, device = cairo_pdf)

cat("✓ Successfully saved Panel E to:", output_pdf, "\n")
