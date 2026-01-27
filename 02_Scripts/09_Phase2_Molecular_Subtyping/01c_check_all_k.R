#!/usr/bin/env Rscript
# Quick script to check cluster sizes for all k values

results <- readRDS("04_Figures/03_Consensus_Clustering/ConsensusPlots/consensus_results.rds")

cat("Cluster sizes for all k values:\n\n")

for (k in 2:10) {
  sizes <- table(results[[k]]$consensusClass)
  pcts <- sizes / sum(sizes) * 100

  cat(sprintf("k=%d: ", k))
  cat(paste(sprintf("%d (%.1f%%)", sizes, pcts), collapse = ", "))
  cat("\n")
}

cat("\n")

# Find k with best balance (smallest cluster >5%)
cat("Finding best balanced k:\n\n")

for (k in 2:10) {
  sizes <- table(results[[k]]$consensusClass)
  min_pct <- min(sizes) / sum(sizes) * 100

  if (min_pct >= 5) {
    cat(sprintf("k=%d has minimum cluster size %.1f%% - GOOD\n", k, min_pct))
  } else {
    cat(sprintf("k=%d has minimum cluster size %.1f%% - too small\n", k, min_pct))
  }
}
