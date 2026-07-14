# R script to generate High-Fidelity Multivariate Forest Plot
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)

cat("=== GENERATING MULTIVARIATE FOREST PLOT ===\n")

# Data from Table 2 (Multivariate Cox Regression)
df <- data.frame(
  Variable = c("Cluster 2 (vs Cluster 1)", "Age (per 10 years)", "Sex (Male vs Female)",
               "TP53 mutation", "TET2 mutation", "RUNX1 mutation", "ASXL1 mutation"),
  HR = c(1.06, 1.03^10, 1.12, 2.96, 1.42, 1.13, 1.21),
  Lower = c(0.83, 1.02^10, 0.91, 2.10, 1.03, 0.78, 0.82),
  Upper = c(1.36, 1.04^10, 1.38, 4.17, 1.94, 1.64, 1.79),
  P_Value = c("0.649", "7.3e-12", "0.278", "5.6e-10", "0.031", "0.518", "0.331"),
  Significant = c("No", "Yes", "No", "Yes", "Yes", "No", "No")
)

# Convert Age (per year) to Age (per 10 years) for better visualization scale
# HR for 10 years = 1.03^10 = 1.34, CI: 1.02^10 - 1.04^10

# Order variables (Cluster 2 at the top, followed by clinical, then mutations)
df$Variable <- factor(df$Variable, levels = rev(c(
  "Cluster 2 (vs Cluster 1)",
  "Age (per 10 years)",
  "Sex (Male vs Female)",
  "TP53 mutation",
  "TET2 mutation",
  "RUNX1 mutation",
  "ASXL1 mutation"
)))

# Shared theme to match Panel A-E
theme_hf <- theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
    axis.title    = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 12, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

p <- ggplot(df, aes(x = HR, y = Variable, color = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c("No" = "gray50", "Yes" = "black")) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 4, 5), limits = c(0.7, 5.5)) +
  geom_text(aes(label = paste0("HR=", sprintf("%.2f", HR), " (p=", P_Value, ")"), 
                x = Upper * 1.1), hjust = 0, size = 3.5, fontface = "plain", show.legend = FALSE) +
  labs(title = "F. Multivariate Analysis for Overall Survival",
       subtitle = "Subtype prognostic effect is driven entirely by genomic co-linearity",
       x = "Hazard Ratio (95% CI, log scale)", y = "") +
  theme_hf +
  expand_limits(x = 6.5)

dir.create("05_Submission/Submission_Hub/03_Supplementary_Figures", showWarnings = FALSE, recursive = TRUE)

pdf("05_Submission/Submission_Hub/03_Supplementary_Figures/FigureS12.pdf", width = 10, height = 4)
print(p)
dev.off()

png("05_Submission/Submission_Hub/03_Supplementary_Figures/FigureS12.png", width = 1500, height = 600, res = 150)
print(p)
dev.off()

cat("✓ Multivariate Forest Plot generated successfully.\n")
