# Install packages for Phase 2 clustering analysis

cat("Installing packages for Phase 2 analysis...\n\n")

# Install ConsensusClusterPlus
if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  cat("Installing ConsensusClusterPlus...\n")
  BiocManager::install("ConsensusClusterPlus", update = FALSE, ask = FALSE)
} else {
  cat("✓ ConsensusClusterPlus already installed\n")
}

# Install GSVA for pathway analysis
if (!requireNamespace("GSVA", quietly = TRUE)) {
  cat("Installing GSVA...\n")
  BiocManager::install("GSVA", update = FALSE, ask = FALSE)
} else {
  cat("✓ GSVA already installed\n")
}

# Install ComplexHeatmap
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  cat("Installing ComplexHeatmap...\n")
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
} else {
  cat("✓ ComplexHeatmap already installed\n")
}

# Install cluster for silhouette
if (!requireNamespace("cluster", quietly = TRUE)) {
  cat("Installing cluster...\n")
  install.packages("cluster", repos = "https://cloud.r-project.org")
} else {
  cat("✓ cluster already installed\n")
}

cat("\n✓ All Phase 2 packages ready!\n")
