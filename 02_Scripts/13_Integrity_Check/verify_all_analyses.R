#!/usr/bin/env Rscript
# ==============================================================================
# BeatAML Project - Comprehensive Integrity Check
# ==============================================================================
# Verifies the integrity and correctness of all analyses
# Date: 2025-10-09
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("COMPREHENSIVE ANALYSIS INTEGRITY CHECK\n")
cat("==============================================================================\n\n")

cat("Starting integrity verification...\n")
cat(sprintf("Date: %s\n\n", Sys.time()))

# Initialize results tracker
integrity_results <- list()
issues_found <- character()
warnings_found <- character()

# ==============================================================================
# PHASE 1: BATCH CORRECTION
# ==============================================================================
cat("==============================================================================\n")
cat("PHASE 1: BATCH CORRECTION VERIFICATION\n")
cat("==============================================================================\n\n")

phase1_pass <- TRUE

# Check batch correction file
batch_file <- "03_Results/04_Batch_Correction/beataml_expression_batchcorrected.rds"
if (file.exists(batch_file)) {
  expr_corrected <- readRDS(batch_file)
  cat(sprintf("✓ Batch-corrected expression matrix: %d genes × %d samples\n",
              nrow(expr_corrected), ncol(expr_corrected)))

  # Check for missing values
  n_missing <- sum(is.na(expr_corrected))
  if (n_missing > 0) {
    issues_found <- c(issues_found, sprintf("Batch-corrected data has %d missing values", n_missing))
    phase1_pass <- FALSE
  } else {
    cat("✓ No missing values in batch-corrected data\n")
  }

  # Check for infinite values
  n_inf <- sum(is.infinite(expr_corrected))
  if (n_inf > 0) {
    issues_found <- c(issues_found, sprintf("Batch-corrected data has %d infinite values", n_inf))
    phase1_pass <- FALSE
  } else {
    cat("✓ No infinite values in batch-corrected data\n")
  }

  # Check expected dimensions
  if (nrow(expr_corrected) < 10000) {
    warnings_found <- c(warnings_found, sprintf("Low number of genes (%d), expected >10,000", nrow(expr_corrected)))
  }

  if (ncol(expr_corrected) < 500) {
    warnings_found <- c(warnings_found, sprintf("Low number of samples (%d), expected >500", ncol(expr_corrected)))
  }

} else {
  issues_found <- c(issues_found, "Batch-corrected expression file not found")
  phase1_pass <- FALSE
}

# Check validation results
validation_file <- "03_Results/04_Batch_Correction/batch_correction_summary.csv"
if (file.exists(validation_file)) {
  validation <- read.csv(validation_file)
  cat(sprintf("✓ Batch correction validation file found (%d metrics)\n", nrow(validation)))

  # Check p-values
  after_pvals <- validation %>% filter(timepoint == "After") %>% pull(p_value)
  if (any(after_pvals < 0.05)) {
    warnings_found <- c(warnings_found, "Some batch effects remain after correction (p<0.05)")
  } else {
    cat("✓ Batch effects successfully removed (all p>0.05)\n")
  }
} else {
  warnings_found <- c(warnings_found, "Batch correction validation file not found")
}

integrity_results$phase1 <- list(pass = phase1_pass, status = ifelse(phase1_pass, "PASS", "FAIL"))
cat(sprintf("\nPhase 1 Status: %s\n\n", integrity_results$phase1$status))

# ==============================================================================
# PHASE 2.1: CONSENSUS CLUSTERING
# ==============================================================================
cat("==============================================================================\n")
cat("PHASE 2.1: CONSENSUS CLUSTERING VERIFICATION\n")
cat("==============================================================================\n\n")

phase2_1_pass <- TRUE

# Check cluster assignments
cluster_file <- "03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv"
if (file.exists(cluster_file)) {
  clusters <- read.csv(cluster_file)
  cat(sprintf("✓ Cluster assignments file found: %d samples\n", nrow(clusters)))

  # Check for required columns
  required_cols <- c("sample_id", "cluster")
  if (all(required_cols %in% colnames(clusters))) {
    cat("✓ Required columns present\n")
  } else {
    issues_found <- c(issues_found, "Missing required columns in cluster assignments")
    phase2_1_pass <- FALSE
  }

  # Check cluster distribution
  cluster_counts <- table(clusters$cluster)
  cat(sprintf("✓ Cluster distribution: %s\n", paste(paste0("Cluster ", names(cluster_counts), ": ", cluster_counts), collapse = ", ")))

  # Check for balanced clusters
  if (length(cluster_counts) == 2) {
    prop1 <- as.numeric(cluster_counts[1]) / sum(cluster_counts)
    if (prop1 < 0.20 || prop1 > 0.80) {
      warnings_found <- c(warnings_found, sprintf("Imbalanced clusters: %.1f%% vs %.1f%%", prop1*100, (1-prop1)*100))
    } else {
      cat(sprintf("✓ Balanced clusters: %.1f%% vs %.1f%%\n", prop1*100, (1-prop1)*100))
    }
  }

  # Check for missing cluster assignments
  n_missing_cluster <- sum(is.na(clusters$cluster))
  if (n_missing_cluster > 0) {
    issues_found <- c(issues_found, sprintf("%d samples with missing cluster assignments", n_missing_cluster))
    phase2_1_pass <- FALSE
  } else {
    cat("✓ No missing cluster assignments\n")
  }

} else {
  issues_found <- c(issues_found, "Cluster assignments file not found")
  phase2_1_pass <- FALSE
}

# Check consensus results
consensus_results_file <- "03_Results/06_Molecular_Subtypes/consensus_clustering_results.rds"
if (file.exists(consensus_results_file)) {
  cat("✓ Consensus clustering results file found\n")
} else {
  warnings_found <- c(warnings_found, "Consensus clustering results RDS not found")
}

integrity_results$phase2_1 <- list(pass = phase2_1_pass, status = ifelse(phase2_1_pass, "PASS", "FAIL"))
cat(sprintf("\nPhase 2.1 Status: %s\n\n", integrity_results$phase2_1$status))

# ==============================================================================
# PHASE 2.2: PATHWAY ENRICHMENT
# ==============================================================================
cat("==============================================================================\n")
cat("PHASE 2.2: PATHWAY ENRICHMENT VERIFICATION\n")
cat("==============================================================================\n\n")

phase2_2_pass <- TRUE

# Check pathway scores
pathway_file <- "03_Results/07_Pathway_Enrichment/pathway_scores_by_cluster.csv"
if (file.exists(pathway_file)) {
  pathways <- read.csv(pathway_file)
  cat(sprintf("✓ Pathway scores file found: %d pathways\n", nrow(pathways)))

  # Check for required columns
  if ("pathway" %in% colnames(pathways)) {
    cat("✓ Pathway column present\n")
  } else {
    issues_found <- c(issues_found, "Missing pathway column")
    phase2_2_pass <- FALSE
  }

  # Check for cluster columns
  cluster_cols <- grep("cluster|Cluster", colnames(pathways), value = TRUE)
  if (length(cluster_cols) >= 2) {
    cat(sprintf("✓ Found %d cluster columns\n", length(cluster_cols)))
  } else {
    issues_found <- c(issues_found, "Missing cluster score columns")
    phase2_2_pass <- FALSE
  }

  # Check for missing values
  n_missing_pathway <- sum(is.na(pathways))
  if (n_missing_pathway > 0) {
    warnings_found <- c(warnings_found, sprintf("%d missing values in pathway scores", n_missing_pathway))
  } else {
    cat("✓ No missing values in pathway scores\n")
  }

  # Check expected number of pathways
  if (nrow(pathways) < 30) {
    warnings_found <- c(warnings_found, sprintf("Low number of pathways (%d), expected ~50", nrow(pathways)))
  } else {
    cat(sprintf("✓ Adequate number of pathways (%d)\n", nrow(pathways)))
  }

} else {
  issues_found <- c(issues_found, "Pathway scores file not found")
  phase2_2_pass <- FALSE
}

integrity_results$phase2_2 <- list(pass = phase2_2_pass, status = ifelse(phase2_2_pass, "PASS", "FAIL"))
cat(sprintf("\nPhase 2.2 Status: %s\n\n", integrity_results$phase2_2$status))

# ==============================================================================
# PHASE 2.3: SUBTYPE CHARACTERIZATION
# ==============================================================================
cat("==============================================================================\n")
cat("PHASE 2.3: SUBTYPE CHARACTERIZATION VERIFICATION\n")
cat("==============================================================================\n\n")

phase2_3_pass <- TRUE

# Check DEG results
deg_file <- "03_Results/06_Molecular_Subtypes/cluster_deg_results.csv"
if (file.exists(deg_file)) {
  degs <- read.csv(deg_file)
  cat(sprintf("✓ DEG results file found: %d genes tested\n", nrow(degs)))

  # Check for required columns
  required_deg_cols <- c("gene", "logFC", "adj.P.Val")
  if (all(required_deg_cols %in% colnames(degs))) {
    cat("✓ Required DEG columns present\n")
  } else {
    issues_found <- c(issues_found, "Missing required DEG columns")
    phase2_3_pass <- FALSE
  }

  # Check number of significant DEGs
  sig_degs <- degs %>% filter(adj.P.Val < 0.05)
  cat(sprintf("✓ Significant DEGs (FDR<0.05): %d genes\n", nrow(sig_degs)))

  if (nrow(sig_degs) < 100) {
    warnings_found <- c(warnings_found, sprintf("Low number of significant DEGs (%d)", nrow(sig_degs)))
  }

  # Check for NA p-values
  n_na_pvals <- sum(is.na(degs$adj.P.Val))
  if (n_na_pvals > 0) {
    warnings_found <- c(warnings_found, sprintf("%d genes with NA p-values", n_na_pvals))
  } else {
    cat("✓ No missing p-values\n")
  }

} else {
  issues_found <- c(issues_found, "DEG results file not found")
  phase2_3_pass <- FALSE
}

# Check clinical associations
clinical_file <- "03_Results/06_Molecular_Subtypes/cluster_clinical_associations.csv"
if (file.exists(clinical_file)) {
  cat("✓ Clinical associations file found\n")
  clinical <- read.csv(clinical_file)
  cat(sprintf("  %d clinical variables tested\n", nrow(clinical)))
} else {
  warnings_found <- c(warnings_found, "Clinical associations file not found")
}

integrity_results$phase2_3 <- list(pass = phase2_3_pass, status = ifelse(phase2_3_pass, "PASS", "FAIL"))
cat(sprintf("\nPhase 2.3 Status: %s\n\n", integrity_results$phase2_3$status))

# ==============================================================================
# PHASE 3: SURVIVAL ANALYSIS
# ==============================================================================
cat("==============================================================================\n")
cat("PHASE 3: SURVIVAL ANALYSIS VERIFICATION\n")
cat("==============================================================================\n\n")

phase3_pass <- TRUE

# Check survival results
survival_file <- "03_Results/08_Survival_Analysis/survival_results_summary.csv"
if (file.exists(survival_file)) {
  survival <- read.csv(survival_file)
  cat(sprintf("✓ Survival results file found: %d clusters analyzed\n", nrow(survival)))

  # Check for required columns
  required_surv_cols <- c("cluster", "n", "median_survival_months")
  if (all(required_surv_cols %in% colnames(survival))) {
    cat("✓ Required survival columns present\n")
  } else {
    issues_found <- c(issues_found, "Missing required survival columns")
    phase3_pass <- FALSE
  }

  # Print median survival
  for (i in 1:nrow(survival)) {
    cat(sprintf("  Cluster %s: n=%d, median survival=%.1f months\n",
                survival$cluster[i], survival$n[i], survival$median_survival_months[i]))
  }

  # Check for reasonable survival times
  if (any(survival$median_survival_months < 0, na.rm = TRUE)) {
    issues_found <- c(issues_found, "Negative survival times detected")
    phase3_pass <- FALSE
  }

  if (any(survival$median_survival_months > 120, na.rm = TRUE)) {
    warnings_found <- c(warnings_found, "Very long survival times detected (>10 years)")
  }

} else {
  issues_found <- c(issues_found, "Survival results file not found")
  phase3_pass <- FALSE
}

# Check Cox regression results
cox_file <- "03_Results/08_Survival_Analysis/cox_regression_results.csv"
if (file.exists(cox_file)) {
  cox <- read.csv(cox_file)
  cat(sprintf("✓ Cox regression file found: %d variables\n", nrow(cox)))

  # Check for cluster variable
  if (any(grepl("cluster", cox$variable, ignore.case = TRUE))) {
    cluster_row <- cox[grepl("cluster", cox$variable, ignore.case = TRUE), ][1, ]
    cat(sprintf("  Cluster HR=%.2f, p=%.4f\n", cluster_row$hazard_ratio, cluster_row$p_value))

    if (cluster_row$p_value >= 0.05) {
      warnings_found <- c(warnings_found, "Cluster not significantly associated with survival (p>=0.05)")
    } else {
      cat("✓ Cluster significantly associated with survival\n")
    }
  }
} else {
  warnings_found <- c(warnings_found, "Cox regression results file not found")
}

integrity_results$phase3 <- list(pass = phase3_pass, status = ifelse(phase3_pass, "PASS", "FAIL"))
cat(sprintf("\nPhase 3 Status: %s\n\n", integrity_results$phase3$status))

# ==============================================================================
# PHASE 4: DRUG RESPONSE
# ==============================================================================
cat("==============================================================================\n")
cat("PHASE 4: DRUG RESPONSE VERIFICATION\n")
cat("==============================================================================\n\n")

phase4_pass <- TRUE

# Check drug association results
drug_file <- "03_Results/09_Drug_Response/drug_cluster_associations.csv"
if (file.exists(drug_file)) {
  drugs <- read.csv(drug_file)
  cat(sprintf("✓ Drug association file found: %d drugs tested\n", nrow(drugs)))

  # Check for required columns
  required_drug_cols <- c("drug", "n_samples", "kruskal_pvalue", "adj_pvalue")
  if (all(required_drug_cols %in% colnames(drugs))) {
    cat("✓ Required drug columns present\n")
  } else {
    issues_found <- c(issues_found, "Missing required drug columns")
    phase4_pass <- FALSE
  }

  # Check significant drugs
  sig_drugs <- drugs %>% filter(adj_pvalue < 0.10)
  cat(sprintf("✓ Significant drugs (FDR<0.10): %d drugs\n", nrow(sig_drugs)))

  if (nrow(sig_drugs) == 0) {
    warnings_found <- c(warnings_found, "No significant drug associations found")
  } else {
    # Show top drugs
    top_drugs <- sig_drugs %>% arrange(adj_pvalue) %>% head(5)
    cat("\nTop 5 significant drugs:\n")
    for (i in 1:min(5, nrow(top_drugs))) {
      cat(sprintf("  %d. %s (p=%.2e, FDR=%.2e)\n",
                  i, top_drugs$drug[i], top_drugs$kruskal_pvalue[i], top_drugs$adj_pvalue[i]))
    }
  }

  # Check for NA p-values
  n_na_drug_pvals <- sum(is.na(drugs$adj_pvalue))
  if (n_na_drug_pvals > 0) {
    warnings_found <- c(warnings_found, sprintf("%d drugs with NA p-values", n_na_drug_pvals))
  } else {
    cat("\n✓ No missing p-values\n")
  }

} else {
  issues_found <- c(issues_found, "Drug association file not found")
  phase4_pass <- FALSE
}

integrity_results$phase4 <- list(pass = phase4_pass, status = ifelse(phase4_pass, "PASS", "FAIL"))
cat(sprintf("\nPhase 4 Status: %s\n\n", integrity_results$phase4$status))

# ==============================================================================
# DATA CONSISTENCY CHECKS
# ==============================================================================
cat("==============================================================================\n")
cat("DATA CONSISTENCY CHECKS\n")
cat("==============================================================================\n\n")

consistency_pass <- TRUE

# Check sample IDs match across phases
if (exists("clusters") && exists("expr_corrected")) {
  n_overlap <- length(intersect(clusters$sample_id, colnames(expr_corrected)))
  cat(sprintf("Sample overlap (clustering vs expression): %d samples\n", n_overlap))

  if (n_overlap < nrow(clusters) * 0.9) {
    issues_found <- c(issues_found, "Low sample overlap between clustering and expression data")
    consistency_pass <- FALSE
  } else {
    cat("✓ Good sample overlap (>90%)\n")
  }
}

# Check figure files exist
cat("\nFigure file verification:\n")
expected_figures <- c(
  "04_Figures/02_Batch_Correction/PCA_before_batch_correction.pdf",
  "04_Figures/02_Batch_Correction/PCA_after_batch_correction.pdf",
  "04_Figures/03_Consensus_Clustering/ConsensusPlots/consensus.pdf",
  "04_Figures/03_Consensus_Clustering/pathway_heatmap_by_cluster.pdf",
  "04_Figures/04_Subtype_Characterization/cluster_sizes.pdf",
  "04_Figures/05_Survival_Analysis/KM_curves_by_cluster.pdf",
  "04_Figures/06_Drug_Response/drug_sensitivity_heatmap.pdf"
)

n_figures_found <- 0
for (fig in expected_figures) {
  if (file.exists(fig)) {
    n_figures_found <- n_figures_found + 1
  } else {
    warnings_found <- c(warnings_found, sprintf("Figure missing: %s", basename(fig)))
  }
}
cat(sprintf("✓ Found %d of %d expected figures\n", n_figures_found, length(expected_figures)))

integrity_results$consistency <- list(pass = consistency_pass, status = ifelse(consistency_pass, "PASS", "FAIL"))
cat(sprintf("\nConsistency Status: %s\n\n", integrity_results$consistency$status))

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("==============================================================================\n")
cat("INTEGRITY CHECK SUMMARY\n")
cat("==============================================================================\n\n")

all_pass <- all(sapply(integrity_results, function(x) x$pass))

cat("Phase Results:\n")
cat(sprintf("  Phase 1 (Batch Correction):      %s\n", integrity_results$phase1$status))
cat(sprintf("  Phase 2.1 (Consensus Clustering): %s\n", integrity_results$phase2_1$status))
cat(sprintf("  Phase 2.2 (Pathway Enrichment):   %s\n", integrity_results$phase2_2$status))
cat(sprintf("  Phase 2.3 (Characterization):     %s\n", integrity_results$phase2_3$status))
cat(sprintf("  Phase 3 (Survival Analysis):      %s\n", integrity_results$phase3$status))
cat(sprintf("  Phase 4 (Drug Response):          %s\n", integrity_results$phase4$status))
cat(sprintf("  Data Consistency:                 %s\n\n", integrity_results$consistency$status))

if (length(issues_found) > 0) {
  cat("CRITICAL ISSUES FOUND:\n")
  for (i in 1:length(issues_found)) {
    cat(sprintf("  ✗ %d. %s\n", i, issues_found[i]))
  }
  cat("\n")
}

if (length(warnings_found) > 0) {
  cat("WARNINGS:\n")
  for (i in 1:length(warnings_found)) {
    cat(sprintf("  ⚠ %d. %s\n", i, warnings_found[i]))
  }
  cat("\n")
}

cat("OVERALL STATUS: ")
if (all_pass && length(issues_found) == 0) {
  cat("✓ ALL CHECKS PASSED\n\n")
  cat("The analyses appear to be complete and correct.\n")
  cat("All expected files are present with reasonable values.\n")
} else if (length(issues_found) > 0) {
  cat("✗ CRITICAL ISSUES DETECTED\n\n")
  cat("Please review and address the critical issues listed above.\n")
} else {
  cat("⚠ PASSED WITH WARNINGS\n\n")
  cat("The analyses are complete but some warnings were noted.\n")
  cat("Review the warnings to ensure they are acceptable.\n")
}

# Save summary
summary_df <- data.frame(
  phase = c("Batch Correction", "Consensus Clustering", "Pathway Enrichment",
            "Subtype Characterization", "Survival Analysis", "Drug Response", "Data Consistency"),
  status = c(integrity_results$phase1$status, integrity_results$phase2_1$status,
             integrity_results$phase2_2$status, integrity_results$phase2_3$status,
             integrity_results$phase3$status, integrity_results$phase4$status,
             integrity_results$consistency$status),
  pass = c(integrity_results$phase1$pass, integrity_results$phase2_1$pass,
           integrity_results$phase2_2$pass, integrity_results$phase2_3$pass,
           integrity_results$phase3$pass, integrity_results$phase4$pass,
           integrity_results$consistency$pass)
)

write.csv(summary_df, "03_Results/INTEGRITY_CHECK_SUMMARY.csv", row.names = FALSE)
cat("\n✓ Integrity check summary saved: 03_Results/INTEGRITY_CHECK_SUMMARY.csv\n")

cat("\n==============================================================================\n")
cat(sprintf("Integrity check completed: %s\n", Sys.time()))
cat("==============================================================================\n")
