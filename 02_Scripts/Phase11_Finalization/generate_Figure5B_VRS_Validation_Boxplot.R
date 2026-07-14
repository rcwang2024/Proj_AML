# R script to generate Figure 5 Panel B (Stratified Boxplot for Ex Vivo Validation)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)
library(ggpubr)

cat("=== GENERATING FIGURE 5 PANEL B (STRATIFIED BOXPLOT) ===\n")

# Load validation data
merged_data <- read_csv("03_Results/29_ExternalValidation/beataml2_VRS_validation.csv", show_col_types = FALSE) %>%
  mutate(VRS_tertile = case_when(
    VRS < 41.8 ~ "Low",
    VRS > 71.0 ~ "High",
    TRUE ~ "Medium"
  )) %>%
  mutate(VRS_tertile = factor(VRS_tertile, levels=c("Low", "Medium", "High"))) %>%
  rename(ven_auc = venetoclax_auc)

# Calculate counts for labels
counts <- merged_data %>%
  group_by(VRS_tertile) %>%
  summarise(n = n(), .groups = "drop")

# Create custom labels with sample counts
label_low <- sprintf("Low VRS\n(n=%d)", counts$n[counts$VRS_tertile == "Low"])
label_med <- sprintf("Medium VRS\n(n=%d)", counts$n[counts$VRS_tertile == "Medium"])
label_high <- sprintf("High VRS\n(n=%d)", counts$n[counts$VRS_tertile == "High"])

merged_data <- merged_data %>%
  mutate(x_label = case_when(
    VRS_tertile == "Low" ~ label_low,
    VRS_tertile == "Medium" ~ label_med,
    VRS_tertile == "High" ~ label_high
  )) %>%
  mutate(x_label = factor(x_label, levels = c(label_low, label_med, label_high)))

# Run Kruskal-Wallis test
kruskal_res <- kruskal.test(ven_auc ~ VRS_tertile, data = merged_data)
cat(sprintf("Kruskal-Wallis test p-value: %.2e\n", kruskal_res$p.value))

# Perform pairwise Wilcoxon tests
w_low_high <- wilcox.test(ven_auc ~ VRS_tertile, data = merged_data %>% filter(VRS_tertile %in% c("Low", "High")))
w_low_med <- wilcox.test(ven_auc ~ VRS_tertile, data = merged_data %>% filter(VRS_tertile %in% c("Low", "Medium")))
w_med_high <- wilcox.test(ven_auc ~ VRS_tertile, data = merged_data %>% filter(VRS_tertile %in% c("Medium", "High")))

cat(sprintf("Wilcoxon Low vs High p-value: %.2e\n", w_low_high$p.value))
cat(sprintf("Wilcoxon Low vs Medium p-value: %.2e\n", w_low_med$p.value))
cat(sprintf("Wilcoxon Medium vs High p-value: %.2e\n", w_med_high$p.value))

# Boxplot Colors (matching Figure 5A)
color_low <- "#E67E22"  # Warm Orange (Resistant)
color_med <- "#95A5A6"  # Soft Gray (Moderate)
color_high <- "#3498DB" # Soft Blue (Sensitive)

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

# Plot Boxplot
p5b <- ggplot(merged_data, aes(x = x_label, y = ven_auc)) +
  geom_boxplot(
    aes(fill = VRS_tertile),
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.5,
    linewidth = 1.0,
    color = "black"
  ) +
  geom_jitter(
    color = "gray20",
    width = 0.15,
    alpha = 0.4,
    size = 1.5
  ) +
  scale_fill_manual(values = c("Low" = color_low, "Medium" = color_med, "High" = color_high)) +
  scale_color_manual(values = c("Low" = color_low, "Medium" = color_med, "High" = color_high)) +
  labs(
    title = "B. Ex Vivo Functional Validation",
    subtitle = "VRS stratifies venetoclax AUC (BeatAML cohort)",
    x = "VRS Clinical Tier",
    y = "Venetoclax AUC\n(Lower = Sensitive)"
  ) +
  theme_hf +
  theme(legend.position = "none") +
  # Add significance bars using ggpubr
  stat_compare_means(
    comparisons = list(c(label_low, label_high), c(label_low, label_med), c(label_med, label_high)),
    method = "wilcox.test",
    label = "p.signif",
    step_increase = 0.08,
    size = 3.5,
    bracket.size = 0.6
  )

# Save as vector PDF overwriting VRS_validation_beataml2_scatter.pdf
output_pdf <- "04_Figures/29_ExternalValidation/VRS_validation_beataml2_scatter.pdf"
output_png <- "04_Figures/29_ExternalValidation/VRS_validation_beataml2_scatter.png"

ggsave(output_pdf, p5b, width = 7.4, height = 6.9, device = cairo_pdf)
ggsave(output_png, p5b, width = 6, height = 5, dpi = 300)

cat("✓ Successfully saved Panel B Boxplot to:", output_pdf, "\n")
