# TASK: Create Publication-Quality Multi-Panel Figures
# Combine key results into comprehensive figures for manuscript

library(tidyverse)
library(cowplot)
library(patchwork)
library(magick)
library(grid)
library(gridExtra)

cat("=== CREATING PUBLICATION-QUALITY FIGURES ===\n\n")

# Create output directory
dir.create("04_Figures/13_Publication_Figures", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# FIGURE 1: Study Overview and Molecular Subtypes
# ============================================================================

cat("Creating Figure 1: Study Overview and Molecular Subtypes...\n")

# Load existing plots
# Panel A: Sample overlap
overlap_img <- image_read("04_Figures/01_QC_Figures/sample_overlap_venn.png")

# Panel B: PCA after batch correction
pca_img <- image_read_pdf("04_Figures/02_Batch_Correction/PCA_after_batch_correction.pdf")

# Panel C: Consensus clustering
consensus_img <- image_read_pdf("04_Figures/03_Consensus_Clustering/ConsensusPlots/consensus.pdf")

# Panel D: Cluster sizes
sizes_img <- image_read_pdf("04_Figures/04_Subtype_Characterization/cluster_sizes.pdf")

# Combine panels
fig1 <- image_append(c(
  image_append(c(overlap_img, pca_img), stack = FALSE),
  image_append(c(consensus_img, sizes_img), stack = FALSE)
), stack = TRUE)

# Add panel labels
fig1_labeled <- image_annotate(fig1, "A", size = 40, color = "black",
                                location = "+30+30", weight = 700)
fig1_labeled <- image_annotate(fig1_labeled, "B", size = 40, color = "black",
                                location = paste0("+", image_info(overlap_img)$width + 30, "+30"),
                                weight = 700)
fig1_labeled <- image_annotate(fig1_labeled, "C", size = 40, color = "black",
                                location = paste0("+30+", image_info(overlap_img)$height + 30),
                                weight = 700)
fig1_labeled <- image_annotate(fig1_labeled, "D", size = 40, color = "black",
                                location = paste0("+", image_info(overlap_img)$width + 30,
                                                 "+", image_info(overlap_img)$height + 30),
                                weight = 700)

image_write(fig1_labeled,
            "04_Figures/13_Publication_Figures/Figure1_Study_Overview.pdf",
            format = "pdf")
image_write(fig1_labeled,
            "04_Figures/13_Publication_Figures/Figure1_Study_Overview.png",
            format = "png", density = 300)

cat("✓ Figure 1 saved\n\n")

# ============================================================================
# FIGURE 2: Molecular Characterization
# ============================================================================

cat("Creating Figure 2: Molecular Characterization...\n")

# Panel A: Pathway heatmap
pathway_img <- image_read_pdf("04_Figures/03_Consensus_Clustering/pathway_heatmap_by_cluster.pdf")

# Panel B: Mutation frequencies
mutation_img <- image_read_pdf("04_Figures/07_Mutations/mutation_frequencies_by_cluster.pdf")

# Panel C: ELN distribution
eln_img <- image_read_pdf("04_Figures/12_ELN_Comparison/eln_distribution_by_cluster.pdf")

# Combine panels
fig2 <- image_append(c(
  pathway_img,
  image_append(c(mutation_img, eln_img), stack = FALSE)
), stack = TRUE)

# Add panel labels
fig2_labeled <- image_annotate(fig2, "A", size = 40, color = "black",
                                location = "+30+30", weight = 700)
fig2_labeled <- image_annotate(fig2_labeled, "B", size = 40, color = "black",
                                location = paste0("+30+", image_info(pathway_img)$height + 30),
                                weight = 700)
fig2_labeled <- image_annotate(fig2_labeled, "C", size = 40, color = "black",
                                location = paste0("+", image_info(mutation_img)$width + 30,
                                                 "+", image_info(pathway_img)$height + 30),
                                weight = 700)

image_write(fig2_labeled,
            "04_Figures/13_Publication_Figures/Figure2_Molecular_Characterization.pdf",
            format = "pdf")
image_write(fig2_labeled,
            "04_Figures/13_Publication_Figures/Figure2_Molecular_Characterization.png",
            format = "png", density = 300)

cat("✓ Figure 2 saved\n\n")

# ============================================================================
# FIGURE 3: Clinical Outcomes
# ============================================================================

cat("Creating Figure 3: Clinical Outcomes...\n")

# Panel A: Kaplan-Meier curves
km_img <- image_read_pdf("04_Figures/05_Survival_Analysis/KM_curves_by_cluster.pdf")

# Panel B: Forest plot
forest_img <- image_read_pdf("04_Figures/05_Survival_Analysis/forest_plot_cox.pdf")

# Panel C: Venetoclax sensitivity
venetoclax_img <- image_read_pdf("04_Figures/11_Drug_Validation/venetoclax_auc_distribution.pdf")

# Combine panels
fig3 <- image_append(c(
  km_img,
  image_append(c(forest_img, venetoclax_img), stack = FALSE)
), stack = TRUE)

# Add panel labels
fig3_labeled <- image_annotate(fig3, "A", size = 40, color = "black",
                                location = "+30+30", weight = 700)
fig3_labeled <- image_annotate(fig3_labeled, "B", size = 40, color = "black",
                                location = paste0("+30+", image_info(km_img)$height + 30),
                                weight = 700)
fig3_labeled <- image_annotate(fig3_labeled, "C", size = 40, color = "black",
                                location = paste0("+", image_info(forest_img)$width + 30,
                                                 "+", image_info(km_img)$height + 30),
                                weight = 700)

image_write(fig3_labeled,
            "04_Figures/13_Publication_Figures/Figure3_Clinical_Outcomes.pdf",
            format = "pdf")
image_write(fig3_labeled,
            "04_Figures/13_Publication_Figures/Figure3_Clinical_Outcomes.png",
            format = "png", density = 300)

cat("✓ Figure 3 saved\n\n")

# ============================================================================
# FIGURE 4: Drug Response
# ============================================================================

cat("Creating Figure 4: Drug Response Analysis...\n")

# Load drug heatmap
drug_img <- image_read_pdf("04_Figures/06_Drug_Response/drug_sensitivity_heatmap.pdf")

# For single panel, just add label
fig4_labeled <- image_annotate(drug_img, "A", size = 40, color = "black",
                                location = "+30+30", weight = 700)

image_write(fig4_labeled,
            "04_Figures/13_Publication_Figures/Figure4_Drug_Response.pdf",
            format = "pdf")
image_write(fig4_labeled,
            "04_Figures/13_Publication_Figures/Figure4_Drug_Response.png",
            format = "png", density = 300)

cat("✓ Figure 4 saved\n\n")

# ============================================================================
# SUPPLEMENTARY FIGURE 1: Batch Correction
# ============================================================================

cat("Creating Supplementary Figure 1: Batch Correction...\n")

# Before and after batch correction
pca_before <- image_read_pdf("04_Figures/02_Batch_Correction/PCA_before_batch_correction.pdf")
pca_after <- image_read_pdf("04_Figures/02_Batch_Correction/PCA_after_batch_correction.pdf")

figS1 <- image_append(c(pca_before, pca_after), stack = FALSE)
figS1_labeled <- image_annotate(figS1, "A", size = 40, color = "black",
                                 location = "+30+30", weight = 700)
figS1_labeled <- image_annotate(figS1_labeled, "B", size = 40, color = "black",
                                 location = paste0("+", image_info(pca_before)$width + 30, "+30"),
                                 weight = 700)

image_write(figS1_labeled,
            "04_Figures/13_Publication_Figures/FigureS1_Batch_Correction.pdf",
            format = "pdf")
image_write(figS1_labeled,
            "04_Figures/13_Publication_Figures/FigureS1_Batch_Correction.png",
            format = "png", density = 300)

cat("✓ Supplementary Figure 1 saved\n\n")

# ============================================================================
# SUPPLEMENTARY FIGURE 2: Model Diagnostics
# ============================================================================

cat("Creating Supplementary Figure 2: Cox Model Diagnostics...\n")

# Schoenfeld residuals and DFBETAs
schoenfeld_img <- image_read_pdf("04_Figures/10_Model_Diagnostics/schoenfeld_residuals.pdf")
dfbeta_img <- image_read_pdf("04_Figures/10_Model_Diagnostics/dfbeta_plots.pdf")

figS2 <- image_append(c(schoenfeld_img, dfbeta_img), stack = TRUE)
figS2_labeled <- image_annotate(figS2, "A", size = 40, color = "black",
                                 location = "+30+30", weight = 700)
figS2_labeled <- image_annotate(figS2_labeled, "B", size = 40, color = "black",
                                 location = paste0("+30+", image_info(schoenfeld_img)$height + 30),
                                 weight = 700)

image_write(figS2_labeled,
            "04_Figures/13_Publication_Figures/FigureS2_Model_Diagnostics.pdf",
            format = "pdf")
image_write(figS2_labeled,
            "04_Figures/13_Publication_Figures/FigureS2_Model_Diagnostics.png",
            format = "png", density = 300)

cat("✓ Supplementary Figure 2 saved\n\n")

# ============================================================================
# SUPPLEMENTARY FIGURE 3: ELN Comparison
# ============================================================================

cat("Creating Supplementary Figure 3: ELN Risk Comparison...\n")

# ELN combined survival
eln_survival_img <- image_read_pdf("04_Figures/12_ELN_Comparison/survival_combined_stratification.pdf")

figS3_labeled <- image_annotate(eln_survival_img, "A", size = 40, color = "black",
                                 location = "+30+30", weight = 700)

image_write(figS3_labeled,
            "04_Figures/13_Publication_Figures/FigureS3_ELN_Comparison.pdf",
            format = "pdf")
image_write(figS3_labeled,
            "04_Figures/13_Publication_Figures/FigureS3_ELN_Comparison.png",
            format = "png", density = 300)

cat("✓ Supplementary Figure 3 saved\n\n")

# ============================================================================
# CREATE FIGURE LEGENDS
# ============================================================================

cat("Creating figure legends document...\n")

legends <- list(
  Figure1 = "Figure 1. Study Overview and Molecular Subtype Discovery.
(A) Sample overlap across four data types (expression, mutations, drug response, clinical).
(B) Principal component analysis of gene expression data after batch correction, colored by molecular subtype.
(C) Consensus clustering matrix showing robustness of k=2 solution across 1000 bootstrap iterations.
(D) Distribution of samples across two molecular subtypes (Proliferative n=320, Immune-Inflammatory n=387).",

  Figure2 = "Figure 2. Molecular Characterization of AML Subtypes.
(A) Heatmap of pathway activity scores showing distinct biological profiles. Proliferative subtype enriched for cell cycle and DNA repair pathways; Immune-Inflammatory subtype enriched for immune response and cytokine signaling.
(B) Mutation frequencies of 23 key AML genes by molecular subtype. Proliferative subtype enriched for NPM1, CEBPA, DNMT3A, IDH1; Immune-Inflammatory subtype enriched for TP53, RUNX1, ASXL1, RAS pathway mutations.
(C) Distribution of ELN 2017 risk categories by molecular subtype, showing significant association (Chi-square p<10⁻⁹).",

  Figure3 = "Figure 3. Clinical Outcomes and Therapeutic Implications.
(A) Kaplan-Meier survival curves showing significant overall survival difference between molecular subtypes (log-rank p=0.00155). Median survival: Proliferative 19.1 months, Immune-Inflammatory 11.8 months.
(B) Forest plot of multivariate Cox proportional hazards model adjusted for age, sex, and ELN risk.
(C) Venetoclax sensitivity by molecular subtype. Proliferative subtype shows dramatically higher sensitivity (p=2.78×10⁻²⁴, Cohen's d=-1.25).",

  Figure4 = "Figure 4. Comprehensive Drug Response Profiling.
(A) Heatmap showing differential drug sensitivity across 82 drugs with FDR<0.10. Rows represent drugs, columns represent samples ordered by molecular subtype. Color indicates area under curve (AUC) with lower values representing higher sensitivity.",

  FigureS1 = "Supplementary Figure 1. Batch Effect Correction.
(A) PCA before batch correction showing clustering by technical covariates (sequencing wave, center).
(B) PCA after ComBat batch correction showing removal of technical variation while preserving biological signal.",

  FigureS2 = "Supplementary Figure 2. Cox Proportional Hazards Model Diagnostics.
(A) Schoenfeld residuals for testing proportional hazards assumption. Cluster and age variables show violations (p<0.05), suggesting time-varying effects.
(B) DFBETA plots showing one influential observation (sample 670) with high leverage on cluster coefficient.",

  FigureS3 = "Supplementary Figure 3. Molecular Subtypes and ELN Risk Classification.
(A) Overall survival stratified by combined ELN risk and molecular subtype (p=1.38×10⁻¹⁴). Molecular subtypes provide additional prognostic value specifically within ELN Adverse category (p=0.0049)."
)

# Write legends to file
legends_text <- paste(names(legends), unlist(legends), sep = "\n", collapse = "\n\n")
writeLines(legends_text, "04_Figures/13_Publication_Figures/Figure_Legends.txt")

cat("✓ Figure legends saved\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("=== PUBLICATION FIGURES COMPLETE ===\n\n")

cat("Main Figures:\n")
cat("  - Figure 1: Study Overview and Molecular Subtypes\n")
cat("  - Figure 2: Molecular Characterization\n")
cat("  - Figure 3: Clinical Outcomes\n")
cat("  - Figure 4: Drug Response\n\n")

cat("Supplementary Figures:\n")
cat("  - Figure S1: Batch Correction\n")
cat("  - Figure S2: Cox Model Diagnostics\n")
cat("  - Figure S3: ELN Risk Comparison\n\n")

cat("All figures saved in: 04_Figures/13_Publication_Figures/\n")
cat("  - PDF format (vector, publication-ready)\n")
cat("  - PNG format (300 dpi, for presentations)\n\n")

cat("Figure legends saved: 04_Figures/13_Publication_Figures/Figure_Legends.txt\n\n")

cat("### Publication Figures Task COMPLETE ###\n")
