# R script to generate Figure 6 Panels E and F (Patient-level Boxplots for Top Salvage Candidates)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)
library(ggpubr)

cat("=== GENERATING FIGURE 6 PANELS E & F (PATIENT-LEVEL BOXPLOTS) ===\n")

# 1. Load Data
clusters <- read_csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv", show_col_types = FALSE)
drug_auc <- read_tsv("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt", show_col_types = FALSE)

# 2. Filter for Panobinostat and Selumetinib
salvage_data <- drug_auc %>%
  filter(inhibitor %in% c("Panobinostat", "Selumetinib (AZD6244)")) %>%
  select(sample_id = dbgap_rnaseq_sample, inhibitor, auc) %>%
  group_by(sample_id, inhibitor) %>%
  summarise(auc = mean(auc, na.rm = TRUE), .groups = "drop")

# 3. Merge
merged_data <- inner_join(clusters, salvage_data, by="sample_id") %>%
  mutate(cluster = factor(sprintf("Cluster %d", cluster), levels = c("Cluster 1", "Cluster 2")))

# Define standard cluster colors (Blue for C1, Orange for C2)
color_c1 <- "#3498DB"
color_c2 <- "#E67E22"

# Theme matching salvage comprehensive script
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 18.2, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 15.9, color = "darkblue", margin = margin(b=10)),
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
counts_pan <- pan_df %>% group_by(cluster) %>% summarise(n = n(), .groups = "drop")
lbl_c1_pan <- sprintf("Cluster 1\n(n=%d)", counts_pan$n[counts_pan$cluster == "Cluster 1"])
lbl_c2_pan <- sprintf("Cluster 2\n(n=%d)", counts_pan$n[counts_pan$cluster == "Cluster 2"])

pan_df <- pan_df %>%
  mutate(x_label = case_when(
    cluster == "Cluster 1" ~ lbl_c1_pan,
    cluster == "Cluster 2" ~ lbl_c2_pan
  )) %>%
  mutate(x_label = factor(x_label, levels = c(lbl_c1_pan, lbl_c2_pan)))

p_pan <- ggplot(pan_df, aes(x = x_label, y = auc)) +
  geom_boxplot(
    aes(fill = cluster),
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.5,
    linewidth = 1.0,
    color = "black"
  ) +
  geom_jitter(
    aes(color = cluster),
    width = 0.15,
    alpha = 0.4,
    size = 1.5
  ) +
  scale_fill_manual(values = c("Cluster 1" = color_c1, "Cluster 2" = color_c2)) +
  scale_color_manual(values = c("Cluster 1" = color_c1, "Cluster 2" = color_c2)) +
  labs(
    x = NULL,
    y = "Panobinostat AUC\n(Lower = Sensitive)",
    title = "E. Panobinostat Patient Response",
    subtitle = "HDAC inhibitor, FDA approved"
  ) +
  theme_hf +
  stat_compare_means(
    comparisons = list(c(lbl_c1_pan, lbl_c2_pan)),
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.5,
    size = 4,
    bracket.size = 0.6
  )

# Save Panobinostat
ggsave("04_Figures/27_Cluster2_Salvage/Figure_Panobinostat_Boxplot.pdf", p_pan, width = 10, height = 8)
ggsave("04_Figures/27_Cluster2_Salvage/Figure_Panobinostat_Boxplot.png", p_pan, width = 10, height = 8, dpi = 300)

# --- Plot 2: Selumetinib Boxplot ---
sel_df <- merged_data %>% filter(inhibitor == "Selumetinib (AZD6244)")
counts_sel <- sel_df %>% group_by(cluster) %>% summarise(n = n(), .groups = "drop")
lbl_c1_sel <- sprintf("Cluster 1\n(n=%d)", counts_sel$n[counts_sel$cluster == "Cluster 1"])
lbl_c2_sel <- sprintf("Cluster 2\n(n=%d)", counts_sel$n[counts_sel$cluster == "Cluster 2"])

sel_df <- sel_df %>%
  mutate(x_label = case_when(
    cluster == "Cluster 1" ~ lbl_c1_sel,
    cluster == "Cluster 2" ~ lbl_c2_sel
  )) %>%
  mutate(x_label = factor(x_label, levels = c(lbl_c1_sel, lbl_c2_sel)))

p_sel <- ggplot(sel_df, aes(x = x_label, y = auc)) +
  geom_boxplot(
    aes(fill = cluster),
    outlier.shape = NA,
    alpha = 0.7,
    width = 0.5,
    linewidth = 1.0,
    color = "black"
  ) +
  geom_jitter(
    aes(color = cluster),
    width = 0.15,
    alpha = 0.4,
    size = 1.5
  ) +
  scale_fill_manual(values = c("Cluster 1" = color_c1, "Cluster 2" = color_c2)) +
  scale_color_manual(values = c("Cluster 1" = color_c1, "Cluster 2" = color_c2)) +
  labs(
    x = NULL,
    y = "Selumetinib AUC\n(Lower = Sensitive)",
    title = "F. Selumetinib Patient Response",
    subtitle = "MEK inhibitor, FDA approved"
  ) +
  theme_hf +
  stat_compare_means(
    comparisons = list(c(lbl_c1_sel, lbl_c2_sel)),
    method = "wilcox.test",
    label = "p.format",
    label.x = 1.5,
    size = 4,
    bracket.size = 0.6
  )

# Save Selumetinib
ggsave("04_Figures/27_Cluster2_Salvage/Figure_Selumetinib_Boxplot.pdf", p_sel, width = 10, height = 8)
ggsave("04_Figures/27_Cluster2_Salvage/Figure_Selumetinib_Boxplot.png", p_sel, width = 10, height = 8, dpi = 300)

cat("✓ Successfully saved patient-level validation boxplots!\n")
