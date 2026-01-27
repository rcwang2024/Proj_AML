################################################################################
# CREATE SUPPLEMENTARY FIGURES S1-S3
# Date: 2025-12-09
################################################################################

setwd("D:/Projects/Project_AML")
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(png)

cat("=== CREATING SUPPLEMENTARY FIGURES ===\n\n")

# Create output directory
dir.create("05_Manuscript/Supplementary_Figures", showWarnings = FALSE, recursive = TRUE)

################################################################################
# FIGURE S1: ALTERNATIVE CLUSTERING COMPARISON
################################################################################

cat("Creating Figure S1: Alternative Clustering Comparison...\n")

# Load data
quality_data <- read.csv("03_Results/11_Survival_Reanalysis/07_alternative_clustering/cluster_quality_comparison.csv")
survival_data <- read.csv("03_Results/11_Survival_Reanalysis/07_alternative_clustering/survival_comparison.csv")

# Merge datasets
alt_cluster_data <- merge(quality_data, survival_data, by = "k")

# Create multi-panel figure
pdf("05_Manuscript/Supplementary_Figures/Figure_S1_Alternative_Clustering.pdf",
    width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

# Panel A: Consensus Score by k
plot(alt_cluster_data$k, alt_cluster_data$mean_consensus,
     type = "b", pch = 19, cex = 2, lwd = 2, col = "#2166AC",
     xlab = "Number of Clusters (k)", ylab = "Mean Consensus Score",
     main = "A. Consensus Score by k",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4,
     ylim = c(0.7, 1.0), xaxt = "n")
axis(1, at = 2:5, labels = 2:5, cex.axis = 1.2)
abline(h = 0.8, lty = 2, col = "gray50")
points(2, alt_cluster_data$mean_consensus[1], pch = 19, cex = 3, col = "#D6604D")
text(2, alt_cluster_data$mean_consensus[1] - 0.05, "k=2\n(selected)", cex = 1.1, font = 2)

# Panel B: Silhouette Score by k
plot(alt_cluster_data$k, alt_cluster_data$mean_silhouette,
     type = "b", pch = 19, cex = 2, lwd = 2, col = "#2166AC",
     xlab = "Number of Clusters (k)", ylab = "Mean Silhouette Score",
     main = "B. Silhouette Score by k",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4,
     ylim = c(0, 0.15), xaxt = "n")
axis(1, at = 2:5, labels = 2:5, cex.axis = 1.2)
abline(h = 0.0, lty = 2, col = "gray50")
points(2, alt_cluster_data$mean_silhouette[1], pch = 19, cex = 3, col = "#D6604D")
text(2, alt_cluster_data$mean_silhouette[1] + 0.015, "k=2\n(selected)", cex = 1.1, font = 2)

# Panel C: Survival Significance
plot(alt_cluster_data$k, -log10(alt_cluster_data$logrank_p),
     type = "b", pch = 19, cex = 2, lwd = 2, col = "#2166AC",
     xlab = "Number of Clusters (k)", ylab = "-log10(p-value)",
     main = "C. Survival Significance by k",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4,
     xaxt = "n")
axis(1, at = 2:5, labels = 2:5, cex.axis = 1.2)
abline(h = -log10(0.05), lty = 2, col = "red")
text(4.5, -log10(0.05) + 0.1, "p=0.05", cex = 1.0, col = "red")
points(2, -log10(alt_cluster_data$logrank_p[1]), pch = 19, cex = 3, col = "#D6604D")

# Panel D: Cluster Size Variability
plot(alt_cluster_data$k, alt_cluster_data$size_cv,
     type = "b", pch = 19, cex = 2, lwd = 2, col = "#2166AC",
     xlab = "Number of Clusters (k)", ylab = "Size CV (lower = better balance)",
     main = "D. Cluster Size Balance by k",
     cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.4,
     xaxt = "n")
axis(1, at = 2:5, labels = 2:5, cex.axis = 1.2)
points(2, alt_cluster_data$size_cv[1], pch = 19, cex = 3, col = "#D6604D")

dev.off()

cat("✓ Figure S1 created\n\n")

################################################################################
# FIGURE S2: PROPORTIONAL HAZARDS DIAGNOSTICS
################################################################################

cat("Creating Figure S2: Proportional Hazards Diagnostics...\n")

# Copy and combine existing PH diagnostic figures
pdf("05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics.pdf",
    width = 12, height = 10)

# Create a layout with 2 rows, 2 columns
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

# Note: Since we can't directly combine PDFs in R without additional packages,
# we'll document what should be included and create a placeholder

plot.new()
text(0.5, 0.7, "Figure S2: Proportional Hazards Diagnostics",
     cex = 2, font = 2)
text(0.5, 0.5, "Panel A: Schoenfeld Residuals", cex = 1.5)
text(0.5, 0.45, "Source: 04_Figures/11_Survival_Reanalysis/01_schoenfeld_residuals.pdf",
     cex = 1.0, col = "blue")

plot.new()
text(0.5, 0.7, "", cex = 2, font = 2)
text(0.5, 0.5, "Panel B: Hazard Ratio Over Time", cex = 1.5)
text(0.5, 0.45, "Source: 04_Figures/11_Survival_Reanalysis/02_hr_over_time.pdf",
     cex = 1.0, col = "blue")

plot.new()
text(0.5, 0.7, "", cex = 2, font = 2)
text(0.5, 0.5, "Panel C: Landmark Analysis HRs", cex = 1.5)
text(0.5, 0.45, "Source: 04_Figures/11_Survival_Reanalysis/03_hr_by_landmark.pdf",
     cex = 1.0, col = "blue")

plot.new()
text(0.5, 0.7, "", cex = 2, font = 2)
text(0.5, 0.5, "Panel D: RMST Difference", cex = 1.5)
text(0.5, 0.45, "Source: 04_Figures/11_Survival_Reanalysis/04_rmst_difference.pdf",
     cex = 1.0, col = "blue")

dev.off()

cat("✓ Figure S2 placeholder created\n")
cat("  NOTE: Combine individual PDF panels manually or use Adobe Acrobat/pdftk\n\n")

################################################################################
# FIGURE S3: META-ANALYSIS WITH PEDIATRIC COHORT
################################################################################

cat("Creating Figure S3: Meta-Analysis Including Pediatric...\n")

# Copy the existing forest plot
if (file.exists("04_Figures/18_TARGET_Validation/forest_plot_all_cohorts.pdf")) {
  file.copy("04_Figures/18_TARGET_Validation/forest_plot_all_cohorts.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S3_Meta_Analysis_All_Cohorts.pdf",
           overwrite = TRUE)
  cat("✓ Figure S3 created (forest plot with all cohorts including TARGET)\n\n")
} else {
  cat("⚠ Forest plot not found, skipping Figure S3\n\n")
}

################################################################################
# ORGANIZE EXISTING SUPPLEMENTARY FIGURES
################################################################################

cat("Organizing existing supplementary figures...\n\n")

# Figure S4: Drug Class Enrichment
if (file.exists("04_Figures/22_Drug_Validation/FigureS3_Drug_Class_Enrichment.pdf")) {
  file.copy("04_Figures/22_Drug_Validation/FigureS3_Drug_Class_Enrichment.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S4_Drug_Class_Enrichment.pdf",
           overwrite = TRUE)
  cat("✓ Copied Figure S4 (Drug Class Enrichment)\n")
}

# Figure S5: Top 20 Drugs Boxplots
if (file.exists("04_Figures/22_Drug_Validation/FigureS1_Top20_Drugs_Boxplots.pdf")) {
  file.copy("04_Figures/22_Drug_Validation/FigureS1_Top20_Drugs_Boxplots.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S5_Top20_Drugs_Boxplots.pdf",
           overwrite = TRUE)
  cat("✓ Copied Figure S5 (Top 20 Drugs)\n")
}

# Figure S6: BCL-2 Pathway Heatmap
if (file.exists("04_Figures/22_Drug_Validation/FigureS2_BCL2_Pathway_Heatmap.pdf")) {
  file.copy("04_Figures/22_Drug_Validation/FigureS2_BCL2_Pathway_Heatmap.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S6_BCL2_Pathway_Heatmap.pdf",
           overwrite = TRUE)
  cat("✓ Copied Figure S6 (BCL-2 Pathway)\n")
}

# Figure S7: Cluster 2 Drug Profile
if (file.exists("04_Figures/27_Cluster2_Salvage/Figure_Cluster_Comparison.pdf")) {
  file.copy("04_Figures/27_Cluster2_Salvage/Figure_Cluster_Comparison.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S7_Cluster2_Drug_Profile.pdf",
           overwrite = TRUE)
  cat("✓ Copied Figure S7 (Cluster 2 Drug Profile)\n")
}

# Figure S8: VRS Distribution
if (file.exists("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_Distribution_Thresholds.pdf")) {
  file.copy("04_Figures/28_VRS_Clinical_Utility/Figure_VRS_Distribution_Thresholds.pdf",
           "05_Manuscript/Supplementary_Figures/Figure_S8_VRS_Distribution.pdf",
           overwrite = TRUE)
  cat("✓ Copied Figure S8 (VRS Distribution)\n")
}

################################################################################
# VERIFY ALL FIGURES
################################################################################

cat("\n=== VERIFYING ALL SUPPLEMENTARY FIGURES ===\n\n")

figures_list <- c(
  "Figure_S1_Alternative_Clustering.pdf",
  "Figure_S2_PH_Diagnostics.pdf",
  "Figure_S3_Meta_Analysis_All_Cohorts.pdf",
  "Figure_S4_Drug_Class_Enrichment.pdf",
  "Figure_S5_Top20_Drugs_Boxplots.pdf",
  "Figure_S6_BCL2_Pathway_Heatmap.pdf",
  "Figure_S7_Cluster2_Drug_Profile.pdf",
  "Figure_S8_VRS_Distribution.pdf"
)

status <- data.frame(
  Figure = paste0("S", 1:8),
  Filename = figures_list,
  Exists = "Checking...",
  Size_KB = 0
)

for (i in 1:nrow(status)) {
  filepath <- file.path("05_Manuscript/Supplementary_Figures", figures_list[i])
  if (file.exists(filepath)) {
    size_kb <- round(file.info(filepath)$size / 1024, 1)
    status$Exists[i] <- "✓ Yes"
    status$Size_KB[i] <- size_kb
  } else {
    status$Exists[i] <- "✗ No"
    status$Size_KB[i] <- NA
  }
}

print(status, row.names = FALSE)

n_complete <- sum(status$Exists == "✓ Yes")
cat(sprintf("\n\nCOMPLETE: %d / %d figures (%.0f%%)\n",
           n_complete, nrow(status),
           n_complete/nrow(status)*100))

if (n_complete == nrow(status)) {
  cat("\n🎉 ALL SUPPLEMENTARY FIGURES COMPLETE!\n\n")
} else {
  cat("\n⚠ Some figures are missing. Check file locations.\n\n")
}

################################################################################
# MANUAL STEPS REQUIRED
################################################################################

cat("=== MANUAL STEPS REQUIRED ===\n\n")

cat("Figure S2 (PH Diagnostics) requires manual assembly:\n")
cat("1. Open Adobe Acrobat or use pdftk\n")
cat("2. Combine these 4 PDFs into a 2x2 panel layout:\n")
cat("   - Panel A: 04_Figures/11_Survival_Reanalysis/01_schoenfeld_residuals.pdf\n")
cat("   - Panel B: 04_Figures/11_Survival_Reanalysis/02_hr_over_time.pdf\n")
cat("   - Panel C: 04_Figures/11_Survival_Reanalysis/03_hr_by_landmark.pdf\n")
cat("   - Panel D: 04_Figures/11_Survival_Reanalysis/04_rmst_difference.pdf\n")
cat("3. Save as: 05_Manuscript/Supplementary_Figures/Figure_S2_PH_Diagnostics.pdf\n\n")

cat("Alternative: Use ImageMagick convert command:\n")
cat("convert -density 300 panel_A.pdf panel_B.pdf panel_C.pdf panel_D.pdf -append Figure_S2.pdf\n\n")

cat("All figures saved to: 05_Manuscript/Supplementary_Figures/\n")
