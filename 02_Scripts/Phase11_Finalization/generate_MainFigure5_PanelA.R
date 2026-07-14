# R script to generate clean high-fidelity Figure 5 Panel A (VRS Histogram)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)

cat("=== GENERATING FIGURE 5 PANEL A ===\n")

# Load data
survival_data <- read_csv("D:/Proj_AML/03_Results/08_Survival_Analysis/survival_data_with_clusters.csv") %>%
  mutate(cluster = factor(cluster, levels=c(1,2), labels=c("Cluster 1", "Cluster 2"))) %>%
  left_join(read_csv("D:/Proj_AML/03_Results/25_Enhancements/vrs_9gene_scores.csv"), by="sample_id")

color_c1 <- "#3498DB" # Blue (Sensitive / Cluster 1)
color_c2 <- "#E67E22" # Orange (Resistant / Cluster 2)

# Professional Theme
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

# Build Histogram
p5a <- ggplot(survival_data, aes(x = VRS9)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 30, color = "black", linewidth = 0.3) +
    geom_vline(xintercept = c(41.8, 71.0), linetype = "dashed", linewidth = 1.5, color = "red") +
    scale_fill_gradient2(low = color_c2, mid = "gray90", high = color_c1, midpoint = 50) +
    annotate("text", x = 20, y = 40, label = "Resistant\n(Low VRS)", color=color_c2, fontface="bold", size=3.5) +
    annotate("text", x = 85, y = 40, label = "Sensitive\n(High VRS)", color=color_c1, fontface="bold", size=3.5) +
    labs(
      title = "A. Venetoclax Response Score (VRS)",
      subtitle = "Continuous gradient of clinical resistance",
      x = "Venetoclax Response Score (0-100)", 
      y = "Patient Count"
    ) +
    theme_hf + theme(legend.position = "none")

# Save as vector PDF
dir.create("04_Figures/29_ExternalValidation", showWarnings = FALSE, recursive = TRUE)
output_pdf <- "04_Figures/29_ExternalValidation/Figure5A_VRS_Histogram.pdf"
ggsave(output_pdf, p5a, width = 7.4, height = 6.9, device = cairo_pdf)

cat("✓ Successfully saved Panel A to:", output_pdf, "\n")
