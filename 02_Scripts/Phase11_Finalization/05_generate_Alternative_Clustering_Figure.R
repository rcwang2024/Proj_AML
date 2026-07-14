# ==============================================================================
# HIGH-FIDELITY SUPPLEMENTARY FIGURE S1: ALTERNATIVE CLUSTERING COMPARISON
# ==============================================================================
# Purpose: Generates a premium, publication-quality multi-panel visualization of
#          cluster stability and survival metrics for different values of k.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

setwd("d:/Proj_AML")

# 1. LOAD DATA
# ==============================================================================
quality_file <- "03_Results/11_Survival_Reanalysis/07_alternative_clustering/cluster_quality_comparison.csv"
survival_file <- "03_Results/11_Survival_Reanalysis/07_alternative_clustering/survival_comparison.csv"

# Load the generated datasets
if (!file.exists(quality_file) || !file.exists(survival_file)) {
    stop("Input data not found. Please run 09_alternative_clustering.R first.")
}

quality_data <- read.csv(quality_file)
survival_data <- read.csv(survival_file)

# Merge datasets
alt_cluster_data <- merge(quality_data, survival_data, by = "k")

# 2. THEME DEFINITION
# ==============================================================================
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

# Color palettes
line_col <- "#2B8CBE"
point_col <- "#2B8CBE"
highlight_col <- "#E34A33"
text_highlight <- "#B30000"

# 3. GENERATE PANELS
# ==============================================================================

# Panel A: Consensus Score by k
p_A <- ggplot(alt_cluster_data, aes(x = k, y = mean_consensus)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_line(color = line_col, linewidth = 1.2) +
  geom_point(aes(color = factor(k == 2), size = factor(k == 2))) +
  scale_color_manual(values = c("FALSE" = point_col, "TRUE" = highlight_col)) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5)) +
  scale_x_continuous(breaks = 2:5) +
  scale_y_continuous(limits = c(0.7, 1.0), breaks = seq(0.7, 1.0, 0.1)) +
  annotate("text", x = 2, y = alt_cluster_data$mean_consensus[alt_cluster_data$k == 2] - 0.05, 
           label = "k=2\n(Selected)", fontface = "bold", color = text_highlight, size = 4) +
  labs(title = "A. Mean Consensus Score",
       x = "Number of Clusters (k)",
       y = "Mean Consensus Score") +
  theme_hf +
  theme(legend.position = "none")

# Panel B: Silhouette Score by k
p_B <- ggplot(alt_cluster_data, aes(x = k, y = mean_silhouette)) +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_line(color = line_col, linewidth = 1.2) +
  geom_point(aes(color = factor(k == 2), size = factor(k == 2))) +
  scale_color_manual(values = c("FALSE" = point_col, "TRUE" = highlight_col)) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5)) +
  scale_x_continuous(breaks = 2:5) +
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
  annotate("text", x = 2, y = alt_cluster_data$mean_silhouette[alt_cluster_data$k == 2] + 0.02, 
           label = "k=2\n(Selected)", fontface = "bold", color = text_highlight, size = 4) +
  labs(title = "B. Mean Silhouette Width",
       x = "Number of Clusters (k)",
       y = "Silhouette Score") +
  theme_hf +
  theme(legend.position = "none")

# Panel C: Survival Significance
p_C <- ggplot(alt_cluster_data, aes(x = k, y = -log10(logrank_p))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 4.5, y = -log10(0.05) + 0.15, label = "p = 0.05", color = "red", fontface = "italic") +
  geom_line(color = line_col, linewidth = 1.2) +
  geom_point(aes(color = factor(k == 2), size = factor(k == 2))) +
  scale_color_manual(values = c("FALSE" = point_col, "TRUE" = highlight_col)) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5)) +
  scale_x_continuous(breaks = 2:5) +
  scale_y_continuous(limits = c(0, max(-log10(alt_cluster_data$logrank_p)) * 1.2)) +
  labs(title = "C. Survival Separation Significance",
       x = "Number of Clusters (k)",
       y = expression(-log[10](italic(p)-value))) +
  theme_hf +
  theme(legend.position = "none")

# Panel D: Cluster Size Variability (CV)
p_D <- ggplot(alt_cluster_data, aes(x = k, y = size_cv)) +
  geom_line(color = line_col, linewidth = 1.2) +
  geom_point(aes(color = factor(k == 2), size = factor(k == 2))) +
  scale_color_manual(values = c("FALSE" = point_col, "TRUE" = highlight_col)) +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5)) +
  scale_x_continuous(breaks = 2:5) +
  scale_y_continuous(limits = c(0, max(alt_cluster_data$size_cv) * 1.2)) +
  labs(title = "D. Cluster Size Variation",
       x = "Number of Clusters (k)",
       y = "Size CV (lower = more balanced)") +
  theme_hf +
  theme(legend.position = "none")

# 4. ASSEMBLE AND SAVE
# ==============================================================================
composite <- (p_A | p_B) / (p_C | p_D) + 
  plot_annotation(
    title = "Optimization of Transcriptomic Lineage State Number (k)",
    subtitle = "Evaluating stability, cohesion, and prognostic separation across clustering resolutions.",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40", margin = margin(b = 20))
    )
  )

out_dir <- "05_Submission/Submission_Hub/03_Supplementary_Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "FigureS1.pdf"), composite, width = 12, height = 10, dpi = 300)
ggsave(file.path(out_dir, "FigureS1.png"), composite, width = 12, height = 10, dpi = 600, bg="white")

cat("SUCCESS: High-fidelity Figure S1 has been created and saved to Submission_Hub.\n")
