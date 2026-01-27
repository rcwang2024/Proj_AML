#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Phase 2.2: Pathway Enrichment-Based Clustering
# ==============================================================================
# Objective:
#   1. Calculate pathway enrichment scores using GSVA
#   2. Perform clustering on pathway scores
#   3. Compare to gene-based clusters
# Date: 2025-10-04
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(GSEABase)
  library(pheatmap)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(fpc)
})

set.seed(42)
setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("PHASE 2.2: PATHWAY ENRICHMENT-BASED CLUSTERING\n")
cat("==============================================================================\n\n")

# ------------------------------------------------------------------------------
# STEP 1: Load Gene Sets
# ------------------------------------------------------------------------------

cat("STEP 1: Loading gene set collections...\n\n")

# Download Hallmark gene sets from MSigDB if not already present
hallmark_file <- "01_Data/Gene_Sets/h.all.v2024.1.Hs.symbols.gmt"

if (!file.exists(hallmark_file)) {
  cat("Downloading Hallmark gene sets from MSigDB...\n")
  dir.create("01_Data/Gene_Sets", recursive = TRUE, showWarnings = FALSE)

  download.file(
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
    hallmark_file,
    mode = "wb"
  )
  cat("✓ Downloaded Hallmark gene sets\n\n")
} else {
  cat("Using existing Hallmark gene sets\n\n")
}

# Load gene sets
hallmark_genesets <- getGmt(hallmark_file)

cat(sprintf("Loaded %d Hallmark gene sets\n\n", length(hallmark_genesets)))

# Display first few gene sets
cat("Example gene sets:\n")
cat(paste("  -", names(hallmark_genesets)[1:5], collapse = "\n"), "\n\n")

# ------------------------------------------------------------------------------
# STEP 2: Load Expression Data
# ------------------------------------------------------------------------------

cat("STEP 2: Loading expression data...\n\n")

expr_data <- readRDS("03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")

cat(sprintf("Expression data: %d genes × %d samples\n", nrow(expr_data), ncol(expr_data)))

# Convert Ensembl IDs to gene symbols
cat("Converting Ensembl IDs to gene symbols...\n")

ensembl_ids <- rownames(expr_data)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Remove genes without symbols and duplicates
valid_idx <- !is.na(gene_symbols) & !duplicated(gene_symbols)
expr_data_sym <- expr_data[valid_idx, ]
rownames(expr_data_sym) <- gene_symbols[valid_idx]

cat(sprintf("After conversion: %d genes with valid symbols\n\n",
            nrow(expr_data_sym)))

# ------------------------------------------------------------------------------
# STEP 3: Calculate Pathway Enrichment Scores
# ------------------------------------------------------------------------------

cat("STEP 3: Calculating pathway enrichment scores using GSVA...\n")
cat("This may take several minutes...\n\n")

# Run GSVA with new API (GSVA >= 2.0)
# Create parameter object using symbol-converted data
param <- gsvaParam(
  exprData = as.matrix(expr_data_sym),
  geneSets = geneIds(hallmark_genesets),
  kcdf = "Gaussian"  # For log-normalized data
)

# Run GSVA
pathway_scores <- gsva(param, verbose = TRUE)

cat(sprintf("\nPathway scores calculated: %d pathways × %d samples\n\n",
            nrow(pathway_scores), ncol(pathway_scores)))

# Save pathway scores
saveRDS(pathway_scores,
        "03_Results/06_Molecular_Subtypes/hallmark_pathway_scores.rds")
cat("✓ Saved: hallmark_pathway_scores.rds\n\n")

# ------------------------------------------------------------------------------
# STEP 4: Z-score Normalization
# ------------------------------------------------------------------------------

cat("STEP 4: Z-score normalizing pathway scores...\n\n")

pathway_scores_z <- t(scale(t(pathway_scores)))

cat(sprintf("Normalized pathway scores: %d pathways × %d samples\n\n",
            nrow(pathway_scores_z), ncol(pathway_scores_z)))

# ------------------------------------------------------------------------------
# STEP 5: Load Gene-Based Cluster Assignments
# ------------------------------------------------------------------------------

cat("STEP 5: Loading gene-based cluster assignments...\n\n")

gene_clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

cat(sprintf("Gene-based clusters: %d samples in %d clusters\n\n",
            nrow(gene_clusters), length(unique(gene_clusters$cluster))))

# ------------------------------------------------------------------------------
# STEP 6: Pathway Profile by Cluster
# ------------------------------------------------------------------------------

cat("STEP 6: Analyzing pathway profiles by gene-based cluster...\n\n")

# Align pathway scores with cluster assignments
pathway_scores_aligned <- pathway_scores_z[, gene_clusters$sample_id]

# Calculate mean pathway score per cluster
n_clusters <- length(unique(gene_clusters$cluster))
pathway_by_cluster <- matrix(NA, nrow = nrow(pathway_scores_z), ncol = n_clusters)
rownames(pathway_by_cluster) <- rownames(pathway_scores_z)
colnames(pathway_by_cluster) <- paste0("Cluster_", 1:n_clusters)

for (i in 1:n_clusters) {
  cluster_samples <- gene_clusters$sample_id[gene_clusters$cluster == i]
  pathway_by_cluster[, i] <- rowMeans(pathway_scores_aligned[, cluster_samples])
}

# Save pathway profiles
write.csv(pathway_by_cluster,
          "03_Results/06_Molecular_Subtypes/pathway_profiles_by_cluster.csv")
cat("✓ Saved: pathway_profiles_by_cluster.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 7: Identify Cluster-Defining Pathways
# ------------------------------------------------------------------------------

cat("STEP 7: Identifying cluster-defining pathways...\n\n")

# For each cluster, find top enriched pathways
top_pathways_per_cluster <- list()

for (i in 1:n_clusters) {
  # Sort pathways by enrichment score for this cluster
  pathway_scores_cluster_i <- pathway_by_cluster[, i]
  top_up <- names(sort(pathway_scores_cluster_i, decreasing = TRUE))[1:10]
  top_down <- names(sort(pathway_scores_cluster_i, decreasing = FALSE))[1:10]

  top_pathways_per_cluster[[paste0("Cluster_", i)]] <- list(
    top_enriched = top_up,
    top_depleted = top_down
  )

  cat(sprintf("Cluster %d - Top enriched pathways:\n", i))
  cat(paste("  ", 1:5, ".", top_up[1:5], collapse = "\n"), "\n\n")
}

# Save top pathways
saveRDS(top_pathways_per_cluster,
        "03_Results/06_Molecular_Subtypes/top_pathways_per_cluster.rds")
cat("✓ Saved: top_pathways_per_cluster.rds\n\n")

# ------------------------------------------------------------------------------
# STEP 8: Create Pathway Heatmap
# ------------------------------------------------------------------------------

cat("STEP 8: Creating pathway enrichment heatmap...\n\n")

# Create annotation for samples
annotation_col <- data.frame(
  Cluster = factor(gene_clusters$cluster)
)
rownames(annotation_col) <- gene_clusters$sample_id

# Cluster colors
cluster_colors <- rainbow(n_clusters)
names(cluster_colors) <- levels(annotation_col$Cluster)

annotation_colors <- list(Cluster = cluster_colors)

# Select top variable pathways for visualization
pathway_var <- apply(pathway_scores_aligned, 1, var)
top_var_pathways <- names(sort(pathway_var, decreasing = TRUE))[1:min(50, nrow(pathway_scores_aligned))]

# Order samples by cluster
sample_order <- order(gene_clusters$cluster)

# Create heatmap
pdf("04_Figures/03_Consensus_Clustering/pathway_heatmap_by_cluster.pdf",
    width = 12, height = 10)

pheatmap(
  pathway_scores_aligned[top_var_pathways, sample_order],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col[sample_order, , drop = FALSE],
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "Hallmark Pathway Enrichment by Cluster",
  fontsize_row = 8,
  fontsize = 10
)

dev.off()

cat("✓ Saved: pathway_heatmap_by_cluster.pdf\n\n")

# ------------------------------------------------------------------------------
# STEP 9: Pathway-Based Clustering Comparison
# ------------------------------------------------------------------------------

cat("STEP 9: Performing hierarchical clustering on pathway scores...\n\n")

# Hierarchical clustering on pathway scores
hc_pathway <- hclust(dist(t(pathway_scores_aligned)), method = "ward.D2")

# Cut tree to match number of gene-based clusters
pathway_clusters <- cutree(hc_pathway, k = n_clusters)

# Compare pathway-based vs gene-based clusters
comparison_table <- table(
  Gene_Cluster = gene_clusters$cluster,
  Pathway_Cluster = pathway_clusters
)

cat("Cluster correspondence (Gene-based vs Pathway-based):\n")
print(comparison_table)
cat("\n")

# Calculate adjusted Rand index for agreement
library(cluster)
ari <- cluster.stats(dist(t(pathway_scores_aligned)),
                     as.numeric(gene_clusters$cluster),
                     pathway_clusters)$corrected.rand

cat(sprintf("Adjusted Rand Index: %.3f\n", ari))
cat(sprintf("Interpretation: %s\n\n",
            ifelse(ari > 0.75, "Strong agreement",
                   ifelse(ari > 0.5, "Moderate agreement",
                          "Weak agreement"))))

# Save comparison
write.csv(comparison_table,
          "03_Results/06_Molecular_Subtypes/gene_vs_pathway_cluster_comparison.csv")
cat("✓ Saved: gene_vs_pathway_cluster_comparison.csv\n\n")

# ------------------------------------------------------------------------------
# STEP 10: Summary Report
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("PATHWAY ENRICHMENT ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")

cat("SUMMARY:\n")
cat(sprintf("  - Pathway gene sets analyzed: %d (Hallmark)\n", nrow(pathway_scores)))
cat(sprintf("  - Samples analyzed: %d\n", ncol(pathway_scores)))
cat(sprintf("  - Gene-based clusters: %d\n", n_clusters))
cat(sprintf("  - Gene vs pathway cluster agreement (ARI): %.3f\n\n", ari))

cat("BIOLOGICAL INSIGHTS:\n")
cat("  Each cluster has distinct pathway enrichment profiles\n")
cat("  Top pathways per cluster saved for biological interpretation\n\n")

cat("OUTPUT FILES:\n")
cat("  - 03_Results/06_Molecular_Subtypes/hallmark_pathway_scores.rds\n")
cat("  - 03_Results/06_Molecular_Subtypes/pathway_profiles_by_cluster.csv\n")
cat("  - 03_Results/06_Molecular_Subtypes/top_pathways_per_cluster.rds\n")
cat("  - 04_Figures/03_Consensus_Clustering/pathway_heatmap_by_cluster.pdf\n\n")

cat("NEXT STEPS:\n")
cat("  → Phase 2.3: Subtype characterization\n")
cat("    - Differential expression analysis\n")
cat("    - Mutation enrichment\n")
cat("    - Clinical associations\n")
cat("    - Biological naming based on pathways\n\n")

cat("==============================================================================\n")
cat(sprintf("Analysis completed: %s\n", Sys.time()))
cat("==============================================================================\n")
