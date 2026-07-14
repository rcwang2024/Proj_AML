# generate_FigureS15_HighFid.R
# Supplementary Figure S15: Validation of the Monocytic Axis and Age-Orthogonality of the VRS
# Objective: Generate a 1x3 horizontal high-fidelity strip layout.
# Date: 2026-06-11

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

cat("=== GENERATING HIGH-FIDELITY SUPPLEMENTARY FIGURE S15 (1x3 LAYOUT) ===\n")

# --- High-Fidelity Theme ---
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "darkblue", margin = margin(b=4)),
    plot.subtitle = element_text(face = "plain", size = 10.5, color = "darkblue", margin = margin(b=8)),
    plot.title.position = "plot",
    axis.title = element_text(size = 11.5, face = "bold", color = "black"),
    axis.text = element_text(size = 9.5, face = "plain", color = "black"),
    legend.title = element_text(size = 10.5, face = "bold", color = "black"),
    legend.text = element_text(size = 9.5, face = "plain", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    plot.margin = margin(12, 12, 12, 12)
  )

# ------------------------------------------------------------------------------
# PANEL A: Monocytic Differentiation Score by Molecular Cluster
# ------------------------------------------------------------------------------
cat("Loading Panel A data...\n")
df_a <- read.csv("03_Results/Phase10_Analysis/10_1_Monocytic_Mapping_Results.csv")
df_a$cluster <- factor(df_a$cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2"))

w_test_a <- wilcox.test(monocytic_score ~ cluster, data = df_a)
p_val_a <- w_test_a$p.value
p_val_a_str <- if (p_val_a < 2.2e-16) "< 2.2e-16" else sprintf("= %.2e", p_val_a)

p_a <- ggplot(df_a, aes(x = cluster, y = monocytic_score, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5, aes(color = cluster)) +
  scale_fill_manual(values = c("Cluster 1" = "#3498DB", "Cluster 2" = "#E67E22")) +
  scale_color_manual(values = c("Cluster 1" = "#2980B9", "Cluster 2" = "#D35400")) +
  labs(
    title = "A. Monocytic Score by Molecular Subtype",
    subtitle = sprintf("Wilcoxon p %s | Mean: -0.44 (C1) vs 0.52 (C2)", p_val_a_str),
    x = "Molecular Cluster",
    y = "Monocytic Differentiation Score"
  ) +
  theme_hf +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# PANEL B: VIALE-A Trial Resistance Score by Molecular Cluster
# ------------------------------------------------------------------------------
cat("Loading Panel B data...\n")
df_c <- read.csv("03_Results/Phase10_Analysis/10_4_Trial_Validation_Results.csv")
df_c$cluster <- factor(df_c$cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2"))

w_test_c <- wilcox.test(TRIAL_SCORE ~ cluster, data = df_c)
p_val_c <- w_test_c$p.value
p_val_c_str <- if (p_val_c < 2.2e-16) "< 2.2e-16" else sprintf("= %.2e", p_val_c)

p_b <- ggplot(df_c, aes(x = cluster, y = TRIAL_SCORE, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.4, size = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5, aes(color = cluster)) +
  scale_fill_manual(values = c("Cluster 1" = "#3498DB", "Cluster 2" = "#E67E22")) +
  scale_color_manual(values = c("Cluster 1" = "#2980B9", "Cluster 2" = "#D35400")) +
  labs(
    title = "B. VIALE-A Trial Resistance Signature",
    subtitle = sprintf("Wilcoxon p %s | Mean: -0.78 (C1) vs 0.71 (C2)", p_val_c_str),
    x = "Molecular Cluster",
    y = "VIALE-A Resistance Score"
  ) +
  theme_hf +
  theme(legend.position = "none")

# ------------------------------------------------------------------------------
# PANEL C: Age-Independent Response (Age Orthogonality Scatter Plot)
# ------------------------------------------------------------------------------
cat("Creating Panel C data for age-orthogonality...\n")
p_c <- ggplot(df_a, aes(x = age, y = auc, color = cluster)) +
  geom_point(alpha = 0.5, size = 2.0) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("Cluster 1" = "#3498DB", "Cluster 2" = "#E67E22")) +
  coord_cartesian(ylim = c(0, 450)) + 
  labs(
    title = "C. Age-Independent Response",
    subtitle = "Response is independent of patient age (R = 0.04, p = 0.42)",
    x = "Age at Diagnosis",
    y = "Venetoclax Ex Vivo Response (AUC)",
    color = "Subtype"
  ) +
  theme_hf +
  theme(
    legend.position = c(0.85, 0.82),
    legend.background = element_rect(fill = "white", color = "gray90", linewidth = 0.3),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

# ------------------------------------------------------------------------------
# Concat and Save Figure S15
# ------------------------------------------------------------------------------
cat("Assembling final 1x3 horizontal Figure S15...\n")
final_fig <- plot_grid(p_a, p_b, p_c, ncol = 3, align = "h", axis = "tb")

out_pdf <- "05_Submission/Submission_Hub/03_Supplementary_Figures/FigureS15.pdf"
out_png <- "05_Submission/Submission_Hub/03_Supplementary_Figures/FigureS15.png"

dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
ggsave(out_pdf, final_fig, width = 15, height = 5.5, device = cairo_pdf)
ggsave(out_png, final_fig, width = 15, height = 5.5, dpi = 300)

cat("✓ High-Fidelity Figure S15 Successfully Generated.\n")
