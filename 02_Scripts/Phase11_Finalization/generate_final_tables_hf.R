# R script to generate Definitive High-Fidelity Tables (N=320)
setwd("d:/Proj_AML")
library(tidyverse)

cat("=== GENERATING DEFINITIVE HIGH-FIDELITY TABLES (N=320) ===\n")

# 1. Load Gold Standard Cohort (N=320)
survival_data <- read_csv("D:/Proj_AML/03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutation_data <- read_csv("D:/Proj_AML/03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv") %>% rename(sample_id = 1)
eln_data <- read_csv("D:/Proj_AML/03_Results/12_ELN_Comparison/samples_with_eln_and_cluster.csv")

# Ensure N=320 (this is the gold standard cohort)
gs_samples <- survival_data$sample_id
cat("Gold Standard N =", length(gs_samples), "\n")

# --- Table S1: Baseline Characteristics (N=320) ---
char_df <- survival_data %>%
  filter(sample_id %in% mutation_data$sample_id) %>%
  left_join(mutation_data, by = "sample_id") %>%
  left_join(eln_data %>% select(sample_id, ELN2017), by = "sample_id") %>%
  mutate(cluster = factor(cluster, levels=c(1,2), labels=c("Cluster 1", "Cluster 2")))

# Summary function
get_summary <- function(df) {
  res <- list()
  res$N <- nrow(df)
  res$Age_Median <- paste0(round(median(df$age, na.rm=T), 1), " (", min(df$age, na.rm=T), "-", max(df$age, na.rm=T), ")")
  res$Male_N <- paste0(sum(df$sex == "M", na.rm=T), " (", round(mean(df$sex == "M", na.rm=T)*100, 1), "%)")
  res$NPM1 <- paste0(sum(df$NPM1 == 1, na.rm=T), " (", round(mean(df$NPM1 == 1, na.rm=T)*100, 1), "%)")
  res$TP53 <- paste0(sum(df$TP53 == 1, na.rm=T), " (", round(mean(df$TP53 == 1, na.rm=T)*100, 1), "%)")
  res$ELN_Adverse <- paste0(sum(df$ELN2017 == "Adverse", na.rm=T), " (", round(mean(df$ELN2017 == "Adverse", na.rm=T)*100, 1), "%)")
  return(as.data.frame(res))
}

table_s1 <- char_df %>%
  group_by(cluster) %>%
  do(get_summary(.)) %>%
  ungroup()

dir.create("05_Submission/Submission_Hub/04_Supplementary_Tables", showWarnings = FALSE, recursive = TRUE)
write_csv(table_s1, "05_Submission/Submission_Hub/04_Supplementary_Tables/TableS1_Sample_Characteristics.csv")
write_csv(table_s1, "05_Submission/Submission_Hub/01_Manuscript_Source/tables/TableS1_Sample_Characteristics.csv")

# --- Table S6: Full Multivariate Cox (N=320) ---
coefs <- read_csv("03_Results/11_Survival_Reanalysis/05_full_model_coefficients.csv", show_col_types = FALSE)
name_map <- c(
  "cluster_assignmentCluster2" = "Cluster 2",
  "AGE" = "Age",
  "SEXM" = "Sex (Male)",
  "TP53" = "TP53",
  "TET2" = "TET2",
  "RUNX1" = "RUNX1",
  "ASXL1" = "ASXL1"
)
table_s6 <- coefs %>%
  mutate(
    Variable = name_map[variable],
    CI_95 = sprintf("%.2f-%.2f", HR_lower, HR_upper),
    HR = round(HR, 2),
    P_value = signif(pvalue, 2)
  ) %>%
  select(Variable, HR, CI_95, P_value)

write_csv(table_s6, "05_Submission/Submission_Hub/04_Supplementary_Tables/TableS6_Multivariate_Cox.csv")
write_csv(table_s6, "05_Submission/Submission_Hub/01_Manuscript_Source/tables/TableS6_Multivariate_Cox.csv")

cat("✓ Definitive High-Fidelity Tables generated and synced to N=320.\n")
