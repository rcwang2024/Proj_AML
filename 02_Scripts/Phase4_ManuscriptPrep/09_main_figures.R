#!/usr/bin/env Rscript
# Phase 4 Part 9: Main Manuscript Figures
# Purpose: Generate publication-ready main figures

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(gridExtra)
  library(grid)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("MAIN MANUSCRIPT FIGURES\n")
cat("==============================================================================\n\n")

# Create output directory
dir.create("04_Figures/21_Main_Figures", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FIGURE 1: MUTATION LANDSCAPE
# ==============================================================================

cat("Creating Figure 1: Mutation Landscape...\n")

# Load data
mutation_enrich <- read.csv("03_Results/11_Survival_Reanalysis/05_mutation_enrichment.csv")
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutations <- read.csv("03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv", row.names = 1)

# Merge
merged <- survival_data %>%
  left_join(mutations %>% tibble::rownames_to_column("sample_id"), by = "sample_id") %>%
  filter(!is.na(cluster))

# Top enriched mutations
top_mutations <- mutation_enrich %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue) %>%
  head(8) %>%
  pull(gene)

# Calculate frequencies
freq_data <- merged %>%
  select(cluster, all_of(top_mutations)) %>%
  tidyr::pivot_longer(cols = all_of(top_mutations), names_to = "Gene", values_to = "Mutated") %>%
  group_by(cluster, Gene) %>%
  summarise(Frequency = mean(Mutated, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(cluster = factor(cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2")))

# Create figure
pdf("04_Figures/21_Main_Figures/Figure1_mutation_landscape.pdf", width = 10, height = 7)

p1 <- ggplot(freq_data, aes(x = reorder(Gene, -Frequency), y = Frequency, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8"),
                    name = "Molecular Subtype") +
  labs(x = "", y = "Mutation Frequency (%)",
       title = "Mutation Enrichment by Molecular Subtype") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "top",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16),
        panel.grid.major.x = element_blank())

print(p1)
dev.off()

cat("  ✓ Figure1_mutation_landscape.pdf\n")

# ==============================================================================
# FIGURE 2: SURVIVAL META-ANALYSIS
# ==============================================================================

cat("Creating Figure 2: Survival Meta-Analysis...\n")

# Load survival data
surv_obj <- Surv(time = survival_data$OS_months, event = survival_data$OS_event)
fit_beataml <- survfit(surv_obj ~ cluster, data = survival_data)

# Create figure with 2 panels
pdf("04_Figures/21_Main_Figures/Figure2_survival_meta_analysis.pdf", width = 12, height = 10)
layout(matrix(c(1, 2), nrow = 2), heights = c(3, 2))

# Panel A: BeatAML KM curve
par(mar = c(5, 5, 4, 2))
plot(fit_beataml, col = c("#E41A1C", "#377EB8"), lwd = 3,
     xlab = "Overall Survival (months)",
     ylab = "Survival Probability",
     main = "A. BeatAML Discovery Cohort (n=671)",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)

legend("topright",
       legend = c("Cluster 1 (better prognosis)", "Cluster 2 (worse prognosis)"),
       col = c("#E41A1C", "#377EB8"),
       lwd = 3, cex = 1.2, bty = "n")

# Add HR text
text(100, 0.2, "HR = 1.39 (1.13-1.68), p = 0.001", cex = 1.2, font = 2)

# Panel B: Forest plot meta-analysis
par(mar = c(5, 8, 4, 2))

meta_data <- data.frame(
  Cohort = c("BeatAML", "TCGA-LAML", "Adult Pooled"),
  N = c(671, 151, 822),
  Events = c(398, 97, 495),
  HR = c(1.38, 1.24, 1.35),
  Lower = c(1.13, 0.80, 1.13),
  Upper = c(1.68, 1.94, 1.62),
  Type = c("Individual", "Individual", "Pooled")
)

y_pos <- seq(nrow(meta_data), 1, -1)
colors <- ifelse(meta_data$Type == "Pooled", "darkred", "navy")

plot(meta_data$HR, y_pos,
     xlim = c(0.5, 3),
     ylim = c(0.5, nrow(meta_data) + 0.5),
     pch = ifelse(meta_data$Type == "Pooled", 18, 19),
     cex = ifelse(meta_data$Type == "Pooled", 3, 2),
     col = colors,
     xlab = "Hazard Ratio (95% CI)",
     ylab = "",
     yaxt = "n",
     main = "B. Meta-Analysis: Adult Cohorts",
     cex.lab = 1.3, cex.main = 1.4)

# Add CI lines
segments(meta_data$Lower, y_pos, meta_data$Upper, y_pos,
         lwd = ifelse(meta_data$Type == "Pooled", 3, 2), col = colors)

# Add reference line
abline(v = 1, lty = 2, col = "gray40", lwd = 2)

# Add cohort labels
axis(2, at = y_pos, labels = paste0(meta_data$Cohort, " (n=", meta_data$N, ")"),
     las = 1, cex.axis = 1.1)

# Add HR values
text(3, y_pos, sprintf("%.2f (%.2f-%.2f)", meta_data$HR, meta_data$Lower, meta_data$Upper),
     cex = 1.0)

# Add heterogeneity text
text(0.6, 0.7, "I² = 0%, p = 0.67\n(No heterogeneity)", cex = 1.1, pos = 4)

dev.off()

cat("  ✓ Figure2_survival_meta_analysis.pdf\n")

# ==============================================================================
# FIGURE 3: AGE-SPECIFIC HETEROGENEITY
# ==============================================================================

cat("Creating Figure 3: Age Heterogeneity...\n")

pdf("04_Figures/21_Main_Figures/Figure3_age_heterogeneity.pdf", width = 12, height = 10)
layout(matrix(c(1, 2), nrow = 2), heights = c(3, 2))

# Panel A: BeatAML KM (same as Figure 2)
par(mar = c(5, 5, 4, 2))
plot(fit_beataml, col = c("#E41A1C", "#377EB8"), lwd = 3,
     xlab = "Overall Survival (months)",
     ylab = "Survival Probability",
     main = "A. BeatAML (Adult, median age 62y) - HR=1.38, p<0.001",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4)

legend("topright",
       legend = c("Cluster 1", "Cluster 2"),
       col = c("#E41A1C", "#377EB8"),
       lwd = 3, cex = 1.2, bty = "n")

# Panel B: Forest plot showing heterogeneity
par(mar = c(5, 10, 4, 2))

all_cohorts <- data.frame(
  Cohort = c("BeatAML (Adult)", "TCGA (Adult)", "TARGET (Pediatric)"),
  Age_Group = c("Adult", "Adult", "Pediatric"),
  N = c(671, 151, 1713),
  HR = c(1.38, 1.24, 0.81),
  Lower = c(1.13, 0.80, 0.66),
  Upper = c(1.68, 1.94, 1.00)
)

y_pos <- seq(nrow(all_cohorts), 1, -1)
colors <- ifelse(all_cohorts$Age_Group == "Pediatric", "darkgreen", "navy")

plot(all_cohorts$HR, y_pos,
     xlim = c(0.5, 2.5),
     ylim = c(0.5, nrow(all_cohorts) + 0.5),
     pch = 19, cex = 2.5,
     col = colors,
     xlab = "Hazard Ratio (95% CI)",
     ylab = "",
     yaxt = "n",
     main = "B. Cross-Cohort Comparison: Age-Specific Heterogeneity",
     cex.lab = 1.3, cex.main = 1.4)

# Add CI lines
segments(all_cohorts$Lower, y_pos, all_cohorts$Upper, y_pos,
         lwd = 3, col = colors)

# Add reference line
abline(v = 1, lty = 2, col = "gray40", lwd = 2)

# Add labels
axis(2, at = y_pos, labels = paste0(all_cohorts$Cohort, "\n(n=", all_cohorts$N, ")"),
     las = 1, cex.axis = 1.0)

# Add HR text
text(2.3, y_pos, sprintf("%.2f\n(%.2f-%.2f)", all_cohorts$HR, all_cohorts$Lower, all_cohorts$Upper),
     cex = 0.9)

# Add legend
legend("topright",
       legend = c("Adult cohort (C2 worse)", "Pediatric cohort (C2 better!)", "No effect"),
       col = c("navy", "darkgreen", "gray40"),
       pch = c(19, 19, NA),
       lty = c(NA, NA, 2),
       lwd = c(NA, NA, 2),
       cex = 1.1, bty = "n")

# Add heterogeneity text
text(0.6, 1, "Heterogeneity:\nI² = 84.8%, p = 0.001", cex = 1.1, pos = 4, font = 2, col = "red")

dev.off()

cat("  ✓ Figure3_age_heterogeneity.pdf\n")

# ==============================================================================
# FIGURE 4: MULTIVARIATE ANALYSIS
# ==============================================================================

cat("Creating Figure 4: Multivariate Analysis...\n")

# Load multivariate results
multivar_coef <- read.csv("03_Results/11_Survival_Reanalysis/05_full_model_coefficients.csv")

# Clean variable names
multivar_plot <- multivar_coef %>%
  mutate(
    Variable_clean = case_when(
      variable == "cluster_assignmentCluster2" ~ "Cluster 2 (vs 1)",
      variable == "AGE" ~ "Age (per year)",
      variable == "SEXM" ~ "Male (vs Female)",
      variable == "TP53" ~ "TP53 mutation",
      variable == "TET2" ~ "TET2 mutation",
      variable == "RUNX1" ~ "RUNX1 mutation",
      variable == "ASXL1" ~ "ASXL1 mutation",
      TRUE ~ variable
    ),
    Significant = pvalue < 0.05
  )

pdf("04_Figures/21_Main_Figures/Figure4_multivariate_analysis.pdf", width = 10, height = 7)
par(mar = c(5, 8, 4, 2))

y_pos <- seq(nrow(multivar_plot), 1, -1)
colors <- ifelse(multivar_plot$Significant, "darkred", "gray50")

plot(multivar_plot$HR, y_pos,
     xlim = c(0.5, 5),
     ylim = c(0.5, nrow(multivar_plot) + 0.5),
     pch = 19, cex = 2,
     col = colors,
     xlab = "Hazard Ratio (95% CI)",
     ylab = "",
     yaxt = "n",
     log = "x",
     main = "Multivariate Cox Regression: Full Model",
     cex.lab = 1.3, cex.main = 1.4)

# Add CI lines
segments(multivar_plot$HR_lower, y_pos, multivar_plot$HR_upper, y_pos,
         lwd = 2, col = colors)

# Add reference line
abline(v = 1, lty = 2, col = "gray40", lwd = 2)

# Add variable labels
axis(2, at = y_pos, labels = multivar_plot$Variable_clean,
     las = 1, cex.axis = 1.1)

# Add HR and p-value text
text(4.5, y_pos,
     sprintf("%.2f\np=%.3f", multivar_plot$HR, multivar_plot$pvalue),
     cex = 0.85)

# Add legend
legend("topright",
       legend = c("Significant (p<0.05)", "Not significant"),
       col = c("darkred", "gray50"),
       pch = 19, pt.cex = 2,
       cex = 1.1, bty = "n")

# Add note about cluster
text(0.55, 2, "Cluster NOT independent\nof mutations (p=0.649)", cex = 1.1, pos = 4, col = "red", font = 2)

# Add model info
mtext("n=459, 282 events | Concordance=0.685", side = 1, line = 4, cex = 1.0)

dev.off()

cat("  ✓ Figure4_multivariate_analysis.pdf\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("MAIN FIGURES COMPLETE\n")
cat("==============================================================================\n\n")

cat("Generated 4 publication-ready figures:\n")
cat("  1. Figure1_mutation_landscape.pdf - Mutation enrichment by cluster\n")
cat("  2. Figure2_survival_meta_analysis.pdf - BeatAML + meta-analysis\n")
cat("  3. Figure3_age_heterogeneity.pdf - Age-specific validation\n")
cat("  4. Figure4_multivariate_analysis.pdf - Independence testing\n\n")

cat("All figures saved to: 04_Figures/21_Main_Figures/\n\n")

cat("KEY MESSAGES:\n")
cat("  Figure 1: Distinct mutation profiles (NPM1+ vs RUNX1/TP53+)\n")
cat("  Figure 2: Robust prognostic effect in adults (meta-analysis HR=1.35, p=0.001)\n")
cat("  Figure 3: Age-specific heterogeneity (opposite effect in children)\n")
cat("  Figure 4: NOT independent of mutations (p=0.649)\n\n")
