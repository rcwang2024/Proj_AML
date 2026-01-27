#!/usr/bin/env Rscript
# Phase 4 Part 4: Multiple Testing Correction Catalog
# Purpose: Document ALL statistical tests and apply appropriate corrections

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

setwd("D:/Projects/Project_AML")

cat("==============================================================================\n")
cat("MULTIPLE TESTING CATALOG\n")
cat("==============================================================================\n\n")

# Create comprehensive catalog of all tests
all_tests <- data.frame(
  Analysis_Part = character(),
  Test_Type = character(),
  Comparison = character(),
  Raw_P_value = numeric(),
  Test_Category = character(),
  Priority = character(),
  stringsAsFactors = FALSE
)

# ==============================================================================
# PART 1: MUTATION ENRICHMENT TESTS
# ==============================================================================

cat("Loading mutation enrichment results...\n")
if (file.exists("03_Results/11_Survival_Reanalysis/05_mutation_enrichment.csv")) {
  mutation_enrich <- read.csv("03_Results/11_Survival_Reanalysis/05_mutation_enrichment.csv")

  mutation_tests <- mutation_enrich %>%
    mutate(
      Analysis_Part = "Mutation Enrichment",
      Test_Type = "Fisher's Exact Test",
      Comparison = paste(gene, "enrichment by cluster"),
      Raw_P_value = pvalue,
      Test_Category = "Exploratory",
      Priority = "Secondary"
    ) %>%
    select(Analysis_Part, Test_Type, Comparison, Raw_P_value, Test_Category, Priority)

  all_tests <- rbind(all_tests, mutation_tests)
  cat(sprintf("  Added %d mutation enrichment tests\n", nrow(mutation_tests)))
}

# ==============================================================================
# PART 2: SURVIVAL ANALYSES (Phase 3)
# ==============================================================================

cat("Adding survival analysis tests...\n")

survival_tests <- data.frame(
  Analysis_Part = c(
    "Survival - Stratified Cox",
    "Survival - Time-varying",
    "Survival - Landmark (6m)",
    "Survival - Landmark (12m)",
    "Survival - Landmark (18m)",
    "Survival - Landmark (24m)",
    "Survival - Landmark (36m)",
    "Survival - RMST (24m)",
    "Survival - RMST (60m)"
  ),
  Test_Type = c(rep("Log-rank / Cox", 7), rep("RMST", 2)),
  Comparison = "Cluster survival difference",
  Raw_P_value = c(0.00155, 0.058, 0.017, 0.011, 0.015, 0.018, 0.020, 0.029, 0.007),
  Test_Category = "Confirmatory",
  Priority = "Primary",
  stringsAsFactors = FALSE
)

all_tests <- rbind(all_tests, survival_tests)
cat(sprintf("  Added %d survival analysis tests\n", nrow(survival_tests)))

# ==============================================================================
# PART 3: MULTIVARIATE MODELS
# ==============================================================================

cat("Adding multivariate analysis tests...\n")

multivar_tests <- data.frame(
  Analysis_Part = "Multivariate Analysis",
  Test_Type = "Likelihood Ratio Test",
  Comparison = c(
    "Cluster vs Clinical only",
    "Cluster in full model (vs Clinical+Mutations)"
  ),
  Raw_P_value = c(0.052, 0.649),
  Test_Category = "Confirmatory",
  Priority = "Primary",
  stringsAsFactors = FALSE
)

all_tests <- rbind(all_tests, multivar_tests)
cat(sprintf("  Added %d multivariate tests\n", nrow(multivar_tests)))

# ==============================================================================
# PART 4: EXTERNAL VALIDATION
# ==============================================================================

cat("Adding external validation tests...\n")

validation_tests <- data.frame(
  Analysis_Part = c("TCGA Validation", "TARGET Validation", "Meta-analysis (Adult)", "Meta-analysis (All cohorts)"),
  Test_Type = c("Cox regression", "Cox regression", "Random effects", "Random effects"),
  Comparison = "Cluster survival difference",
  Raw_P_value = c(0.353, 0.052, 0.001, 0.0014),  # Last is heterogeneity test
  Test_Category = c("Confirmatory", "Exploratory", "Confirmatory", "Exploratory"),
  Priority = c("Primary", "Secondary", "Primary", "Secondary"),
  stringsAsFactors = FALSE
)

all_tests <- rbind(all_tests, validation_tests)
cat(sprintf("  Added %d validation tests\n", nrow(validation_tests)))

# ==============================================================================
# PART 5: UNIVARIATE ANALYSES (if available)
# ==============================================================================

cat("Loading univariate analyses...\n")
if (file.exists("03_Results/11_Survival_Reanalysis/05_univariate_results.csv")) {
  univar <- read.csv("03_Results/11_Survival_Reanalysis/05_univariate_results.csv")

  univar_tests <- univar %>%
    mutate(
      Analysis_Part = "Univariate Survival",
      Test_Type = "Univariate Cox",
      Comparison = paste(variable, "prognostic value"),
      Raw_P_value = pvalue,
      Test_Category = "Exploratory",
      Priority = "Secondary"
    ) %>%
    select(Analysis_Part, Test_Type, Comparison, Raw_P_value, Test_Category, Priority)

  all_tests <- rbind(all_tests, univar_tests)
  cat(sprintf("  Added %d univariate tests\n", nrow(univar_tests)))
}

# ==============================================================================
# APPLY MULTIPLE TESTING CORRECTIONS
# ==============================================================================

cat("\n")
cat("Applying multiple testing corrections...\n\n")

# Within-analysis FDR correction
all_tests <- all_tests %>%
  group_by(Analysis_Part) %>%
  mutate(
    FDR_within_analysis = p.adjust(Raw_P_value, method = "BH"),
    Bonferroni_within_analysis = p.adjust(Raw_P_value, method = "bonferroni")
  ) %>%
  ungroup()

# Study-wide correction for PRIMARY tests only
primary_tests <- all_tests %>% filter(Priority == "Primary")
primary_tests$FDR_studywide <- p.adjust(primary_tests$Raw_P_value, method = "BH")
primary_tests$Bonferroni_studywide <- p.adjust(primary_tests$Raw_P_value, method = "bonferroni")

# Merge back
all_tests <- all_tests %>%
  left_join(
    primary_tests %>% select(Comparison, Analysis_Part, FDR_studywide, Bonferroni_studywide),
    by = c("Comparison", "Analysis_Part")
  )

# Flag significant results
all_tests <- all_tests %>%
  mutate(
    Sig_raw = Raw_P_value < 0.05,
    Sig_FDR_within = FDR_within_analysis < 0.05,
    Sig_FDR_studywide = !is.na(FDR_studywide) & FDR_studywide < 0.05,
    Interpretation = case_when(
      is.na(Raw_P_value) ~ "Not tested",
      Sig_FDR_studywide ~ "Significant (study-wide FDR<0.05)",
      Sig_FDR_within ~ "Significant (within-analysis FDR<0.05)",
      Sig_raw ~ "Nominally significant (p<0.05)",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(Priority, Analysis_Part, Raw_P_value)

# Save full catalog
write.csv(all_tests,
          "03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv",
          row.names = FALSE)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

cat("=== MULTIPLE TESTING SUMMARY ===\n\n")
cat(sprintf("Total tests performed: %d\n", nrow(all_tests)))
cat(sprintf("  Primary confirmatory tests: %d\n", sum(all_tests$Priority == "Primary")))
cat(sprintf("  Secondary/exploratory tests: %d\n", sum(all_tests$Priority == "Secondary")))
cat("\n")

cat("Significant tests:\n")
cat(sprintf("  Raw p<0.05: %d (%.1f%%)\n",
            sum(all_tests$Sig_raw, na.rm = TRUE),
            sum(all_tests$Sig_raw, na.rm = TRUE) / nrow(all_tests) * 100))
cat(sprintf("  FDR<0.05 (within analysis): %d (%.1f%%)\n",
            sum(all_tests$Sig_FDR_within, na.rm = TRUE),
            sum(all_tests$Sig_FDR_within, na.rm = TRUE) / nrow(all_tests) * 100))
cat(sprintf("  Study-wide FDR<0.05 (primary only): %d (%.1f%%)\n\n",
            sum(all_tests$Sig_FDR_studywide, na.rm = TRUE),
            sum(all_tests$Sig_FDR_studywide, na.rm = TRUE) / sum(all_tests$Priority == "Primary") * 100))

# Key findings that survive correction
cat("=== KEY FINDINGS (Study-wide FDR < 0.05) ===\n\n")
key_findings <- all_tests %>%
  filter(Sig_FDR_studywide) %>%
  arrange(FDR_studywide) %>%
  select(Analysis_Part, Comparison, Raw_P_value, FDR_studywide)

if (nrow(key_findings) > 0) {
  print(key_findings, row.names = FALSE)
} else {
  cat("  No primary tests survive study-wide FDR correction\n")
}

cat("\n")

# Findings by category
cat("=== SIGNIFICANT FINDINGS BY CATEGORY ===\n\n")

for (cat_name in unique(all_tests$Analysis_Part)) {
  cat_tests <- all_tests %>% filter(Analysis_Part == cat_name)
  sig_within <- sum(cat_tests$Sig_FDR_within, na.rm = TRUE)

  cat(sprintf("%s:\n", cat_name))
  cat(sprintf("  Total tests: %d\n", nrow(cat_tests)))
  cat(sprintf("  Significant (FDR<0.05 within): %d\n", sig_within))

  if (sig_within > 0) {
    sig_tests <- cat_tests %>%
      filter(Sig_FDR_within) %>%
      head(3) %>%
      select(Comparison, Raw_P_value, FDR_within_analysis)

    print(sig_tests, row.names = FALSE)
  }
  cat("\n")
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

cat("Creating visualization...\n")

# P-value distribution
pdf("04_Figures/20_Manuscript_Prep/pvalue_distribution.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# Panel A: Histogram of raw p-values
hist(all_tests$Raw_P_value, breaks = 20,
     main = "A. Distribution of Raw P-values",
     xlab = "P-value",
     ylab = "Frequency",
     col = "lightblue",
     border = "navy")
abline(v = 0.05, col = "red", lwd = 2, lty = 2)
text(0.05, par("usr")[4] * 0.9, "p=0.05", col = "red", pos = 4)

# Panel B: FDR vs raw p-values (primary tests)
primary <- all_tests %>% filter(Priority == "Primary" & !is.na(FDR_studywide))
plot(primary$Raw_P_value, primary$FDR_studywide,
     main = "B. Raw P-values vs Study-wide FDR (Primary Tests)",
     xlab = "Raw P-value",
     ylab = "FDR-adjusted P-value",
     pch = 19,
     col = ifelse(primary$Sig_FDR_studywide, "red", "gray50"),
     cex = 1.5)
abline(0, 1, lty = 2, col = "gray30")
abline(h = 0.05, lty = 2, col = "red")
abline(v = 0.05, lty = 2, col = "red")
legend("topleft",
       legend = c("Significant (FDR<0.05)", "Not significant"),
       col = c("red", "gray50"),
       pch = 19,
       bty = "n")

dev.off()

cat("\n✅ Multiple testing catalog complete\n\n")
cat("Files generated:\n")
cat("  - 03_Results/21_Manuscript_Prep/all_statistical_tests_catalog.csv\n")
cat("  - 04_Figures/20_Manuscript_Prep/pvalue_distribution.pdf\n\n")

cat("==============================================================================\n")
cat("RECOMMENDATION FOR MANUSCRIPT\n")
cat("==============================================================================\n\n")

n_primary <- sum(all_tests$Priority == "Primary")
n_primary_sig <- sum(all_tests$Sig_FDR_studywide, na.rm = TRUE)

cat(sprintf("Primary confirmatory tests: %d\n", n_primary))
cat(sprintf("Significant after study-wide FDR correction: %d\n", n_primary_sig))
cat("\n")

if (n_primary_sig > 0) {
  cat("✓ Multiple key findings survive study-wide correction\n")
  cat("✓ Report study-wide FDR for primary tests\n")
  cat("✓ Report within-analysis FDR for exploratory tests\n")
} else {
  cat("⚠ No primary tests survive study-wide FDR correction\n")
  cat("  Recommendation: Report raw p-values with clear labeling as exploratory\n")
  cat("  Consider pre-specifying single primary endpoint for future studies\n")
}

cat("\n")
