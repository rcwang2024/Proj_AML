#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - MASTER SEQUENTIAL ANALYSIS SCRIPT
# ==============================================================================
# This script runs all analyses in sequence:
#   Phase 2.1: Consensus clustering (IF NOT ALREADY RUN)
#   Phase 2.2: Pathway enrichment
#   Phase 2.3: Subtype characterization
#   Phase 3.1: Survival analysis
#   Phase 4.1: Drug response integration
# ==============================================================================

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("BEATML PROJECT - SEQUENTIAL ANALYSIS PIPELINE\n")
cat("==============================================================================\n\n")

start_time <- Sys.time()

# Check if consensus clustering results exist
if (file.exists("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")) {
  cat("✓ Phase 2.1 (Consensus Clustering) already complete\n")
  cat("  Skipping to Phase 2.2...\n\n")
} else {
  cat("Phase 2.1: Consensus Clustering not yet complete\n")
  cat("Please wait for 01_consensus_clustering.R to finish\n\n")
  stop("Consensus clustering must complete first")
}

# ------------------------------------------------------------------------------
# Phase 2.2: Pathway Enrichment
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("PHASE 2.2: PATHWAY ENRICHMENT\n")
cat("==============================================================================\n\n")

source("02_Scripts/09_Phase2_Molecular_Subtyping/02_pathway_enrichment_clustering.R")

cat("\n✓ Phase 2.2 complete\n\n")

# ------------------------------------------------------------------------------
# Phase 2.3: Subtype Characterization
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("PHASE 2.3: SUBTYPE CHARACTERIZATION\n")
cat("==============================================================================\n\n")

source("02_Scripts/09_Phase2_Molecular_Subtyping/03_subtype_characterization.R")

cat("\n✓ Phase 2.3 complete\n\n")

# ------------------------------------------------------------------------------
# Phase 3.1: Survival Analysis
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("PHASE 3.1: SURVIVAL ANALYSIS\n")
cat("==============================================================================\n\n")

source("02_Scripts/10_Phase3_Survival_Analysis/01_survival_by_subtype.R")

cat("\n✓ Phase 3.1 complete\n\n")

# ------------------------------------------------------------------------------
# Phase 4.1: Drug Response
# ------------------------------------------------------------------------------

cat("==============================================================================\n")
cat("PHASE 4.1: DRUG RESPONSE INTEGRATION\n")
cat("==============================================================================\n\n")

source("02_Scripts/11_Phase4_Drug_Response/01_drug_response_by_subtype.R")

cat("\n✓ Phase 4.1 complete\n\n")

# ------------------------------------------------------------------------------
# Final Summary
# ------------------------------------------------------------------------------

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("==============================================================================\n")
cat("ALL ANALYSES COMPLETE!\n")
cat("==============================================================================\n\n")

cat(sprintf("Total analysis time: %.1f minutes\n\n", as.numeric(total_time)))

cat("COMPLETED PHASES:\n")
cat("  ✓ Phase 2.2: Pathway Enrichment\n")
cat("  ✓ Phase 2.3: Subtype Characterization\n")
cat("  ✓ Phase 3.1: Survival Analysis\n")
cat("  ✓ Phase 4.1: Drug Response Integration\n\n")

cat("OUTPUT LOCATIONS:\n")
cat("  - Results: 03_Results/\n")
cat("  - Figures: 04_Figures/\n")
cat("  - Documentation: 06_Documentation/\n\n")

cat("NEXT STEPS:\n")
cat("  → Review all figures and results\n")
cat("  → Generate final integrated summary report\n")
cat("  → Prepare manuscript figures\n\n")

cat("==============================================================================\n")
cat(sprintf("Pipeline completed: %s\n", Sys.time()))
cat("==============================================================================\n")
