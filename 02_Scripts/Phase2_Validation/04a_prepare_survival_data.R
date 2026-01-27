# Prepare survival data with clusters for Cox model analysis

library(tidyverse)
library(readxl)

cat("=== PREPARING SURVIVAL DATA WITH CLUSTERS ===\n\n")

# Load clinical data
clinical <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
cat("Loaded clinical data:", nrow(clinical), "samples\n")

# Load cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
cat("Loaded cluster assignments:", nrow(clusters), "samples\n\n")

# Merge (use dbgap_rnaseq_sample for RNA IDs)
survival_data <- clinical %>%
  left_join(clusters, by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  filter(!is.na(cluster)) %>%
  select(sample_id = dbgap_rnaseq_sample, cluster, age = ageAtDiagnosis, sex = consensus_sex,
         OS_months = overallSurvival, OS_event = vitalStatus) %>%
  mutate(
    OS_event = ifelse(OS_event == "Dead", 1, 0),
    sex = ifelse(sex == "Male", "M", "F")
  ) %>%
  filter(!is.na(OS_months) & !is.na(OS_event))

cat("Merged survival data:", nrow(survival_data), "samples\n")
cat("Cluster distribution:\n")
print(table(survival_data$cluster))
cat("\nEvents:", sum(survival_data$OS_event), "\n\n")

# Save
write.csv(survival_data, "03_Results/08_Survival_Analysis/survival_data_with_clusters.csv",
          row.names = FALSE)

cat("✓ Saved: 03_Results/08_Survival_Analysis/survival_data_with_clusters.csv\n")
cat("Columns:", paste(colnames(survival_data), collapse = ", "), "\n")
