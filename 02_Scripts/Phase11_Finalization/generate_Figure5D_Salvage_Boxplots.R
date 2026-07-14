# R script to generate Figure 5 Panel D (Stratified Boxplots for Salvage Targets)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

cat("=== GENERATING FIGURE 5 PANEL D (SALVAGE TARGETS BOXPLOTS) ===\n")

# Load data
vrs_class <- read_csv("03_Results/28_VRS_Clinical_Utility/VRS_with_clinical_classifications.csv", show_col_types = FALSE)
drug_auc <- read_tsv("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt", show_col_types = FALSE)

# Filter Panobinostat and Selumetinib
salvage_data <- drug_auc %>%
  filter(inhibitor %in% c("Panobinostat", "Selumetinib (AZD6244)")) %>%
  select(sample_id = dbgap_rnaseq_sample, inhibitor, auc) %>%
  group_by(sample_id, inhibitor) %>%
  summarise(auc = mean(auc, na.rm = TRUE), .groups = "drop")

# Merge
merged_data <- inner_join(vrs_class, salvage_data, by="sample_id") %>%
  mutate(VRS_tertile = factor(VRS_tertile, levels=c("Low", "Medium", "High")))

# Colors matching Figure 5A/B
color_low <- "#E67E22"  # Warm Orange (Ven-Resistant / Salvage-Sensitive)
color_med <- "#95A5A6"  # Soft Gray (Moderate)
color_high <- "#3498DB" # Soft Blue (Ven-Sensitive / Salvage-Resistant)

# Custom theme
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

# --- Plot 1: Panobinostat Boxplot ---
pan_df <- merged_data %>% filter(inhibitor == "Panobinostat")
counts_pan <- pan_df %>% group_by(VRS_tertile) %>% summarise(n = n(), .groups = "drop")
lbl_low_pan <- sprintf("Low VRS\n(n=%d)", counts_pan$n[counts_pan$VRS_tertile == "Low"])
lbl_med_pan <- sprintf("Medium VRS\n(n=%d)", counts_pan$n[counts_pan$VRS_tertile == "Medium"])
lbl_high_pan <- sprintf("High VRS\n(n=%d)", counts_pan$n[counts_pan$VRS_tertile == "High"])

pan_df <- pan_df %>%
  mutate(x_label = case_when(
    VRS_tertile == "Low" ~ lbl_low_pan,
    VRS_tertile == "Medium" ~ lbl_med_pan,
    VRS_tertile == "High" ~ lbl_high_pan
  )) %>%
  mutate(x_label = factor(x_label, levels = c(lbl_low_pan, lbl_med_pan, lbl_high_pan)))

p_pan <- ggplot(pan_df, aes(x = x_label, y = auc)) +
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
    x = "VRS Clinical Tier",
    y = "Panobinostat AUC\n(Lower = Sensitive)",
    title = "Panobinostat (HDACi)",
    subtitle = NULL
  ) +
  theme_hf + theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(c(lbl_low_pan, lbl_high_pan), c(lbl_low_pan, lbl_med_pan), c(lbl_med_pan, lbl_high_pan)),
    method = "wilcox.test",
    label = "p.signif",
    step_increase = 0.08,
    size = 3.5,
    bracket.size = 0.6
  )

# --- Plot 2: Selumetinib Boxplot ---
sel_df <- merged_data %>% filter(inhibitor == "Selumetinib (AZD6244)")
counts_sel <- sel_df %>% group_by(VRS_tertile) %>% summarise(n = n(), .groups = "drop")
lbl_low_sel <- sprintf("Low VRS\n(n=%d)", counts_sel$n[counts_sel$VRS_tertile == "Low"])
lbl_med_sel <- sprintf("Medium VRS\n(n=%d)", counts_sel$n[counts_sel$VRS_tertile == "Medium"])
lbl_high_sel <- sprintf("High VRS\n(n=%d)", counts_sel$n[counts_sel$VRS_tertile == "High"])

sel_df <- sel_df %>%
  mutate(x_label = case_when(
    VRS_tertile == "Low" ~ lbl_low_sel,
    VRS_tertile == "Medium" ~ lbl_med_sel,
    VRS_tertile == "High" ~ lbl_high_sel
  )) %>%
  mutate(x_label = factor(x_label, levels = c(lbl_low_sel, lbl_med_sel, lbl_high_sel)))

p_sel <- ggplot(sel_df, aes(x = x_label, y = auc)) +
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
    title = "Selumetinib (AZD6244)",
    subtitle = NULL,
    x = "",
    y = ""
  ) +
  theme_hf + theme(legend.position = "none") +
  stat_compare_means(
    comparisons = list(c(lbl_low_sel, lbl_high_sel), c(lbl_low_sel, lbl_med_sel), c(lbl_med_sel, lbl_high_sel)),
    method = "wilcox.test",
    label = "p.signif",
    step_increase = 0.08,
    size = 3.5,
    bracket.size = 0.6
  )

# Combine using patchwork
p_combined <- p_pan + p_sel + plot_layout(ncol = 2) +
  plot_annotation(
    title = "D. Therapeutic Vulnerabilities",
    subtitle = "VRS-High blasts remain sensitive to salvage targets",
    theme = theme(
      plot.title = element_text(face = "bold", size = 12.7, color = "darkblue", margin = margin(b=8)),
      plot.subtitle = element_text(face = "plain", size = 11.1, color = "darkblue", margin = margin(b=10))
    )
  )

# Save as vector PDF and PNG
output_pdf <- "04_Figures/29_ExternalValidation/Figure5D_Salvage_Boxplots.pdf"
output_png <- "04_Figures/29_ExternalValidation/Figure5D_Salvage_Boxplots.png"

dir.create("04_Figures/29_ExternalValidation", showWarnings = FALSE, recursive = TRUE)
ggsave(output_pdf, p_combined, width = 12, height = 5.5, device = cairo_pdf)
ggsave(output_png, p_combined, width = 12, height = 5.5, dpi = 300)

cat("✓ Successfully saved Panel D Boxplots to:", output_pdf, "\n")
