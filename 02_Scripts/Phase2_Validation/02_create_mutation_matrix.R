# TASK 1.2: Create Mutation Matrix with Matched Samples
# Uses the ID mapping created in Task 1.1

library(tidyverse)

cat("=== TASK 1.2: CREATE MUTATION MATRIX ===\n\n")

# Load ID mapping
if (!file.exists("03_Results/sample_id_mapping.csv")) {
  cat("ERROR: No ID mapping found. Run Task 1.1 first.\n")
  quit(save = "no", status = 1)
}

id_mapping <- read.csv("03_Results/sample_id_mapping.csv")
cat("Loaded ID mapping:", nrow(id_mapping), "matched samples\n\n")

# Load mutation data
cat("Loading mutation data...\n")
mutation_data <- read.delim("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt",
                            stringsAsFactors = FALSE, sep = "\t")

cat("Total mutations:", nrow(mutation_data), "\n")

# Filter to matched samples only
mutation_data_matched <- mutation_data %>%
  filter(dbgap_sample_id %in% id_mapping$mutation_id)

cat("Mutations for matched samples:", nrow(mutation_data_matched), "\n\n")

# Define key AML genes
key_aml_genes <- c("NPM1", "FLT3", "DNMT3A", "TET2", "IDH1", "IDH2",
                   "TP53", "RUNX1", "CEBPA", "NRAS", "KRAS", "PTPN11",
                   "ASXL1", "SRSF2", "SF3B1", "U2AF1", "ZRSR2", "EZH2",
                   "BCOR", "STAG2", "RAD21", "SMC1A", "SMC3")

cat("Focusing on", length(key_aml_genes), "key AML genes\n\n")

# Create binary mutation matrix
# Keep only non-synonymous mutations
non_syn_classes <- c("missense_variant", "nonsense_mutation",
                     "frameshift_variant", "frame_shift_del", "frame_shift_ins",
                     "splice_site", "in_frame_del", "in_frame_ins",
                     "stop_gained", "start_lost")

mutation_data_filtered <- mutation_data_matched %>%
  filter(symbol %in% key_aml_genes) %>%
  filter(variant_classification %in% non_syn_classes) %>%
  distinct(symbol, dbgap_sample_id)

cat("Non-synonymous mutations in key genes:", nrow(mutation_data_filtered), "\n\n")

# Create binary matrix (0 = wild-type, 1 = mutated)
mutation_matrix_long <- mutation_data_filtered %>%
  mutate(mutated = 1) %>%
  select(symbol, dbgap_sample_id, mutated)

# Pivot to wide format
mutation_matrix_wide <- mutation_matrix_long %>%
  pivot_wider(names_from = dbgap_sample_id,
              values_from = mutated,
              values_fill = 0)

# Convert to matrix
genes <- mutation_matrix_wide$symbol
mut_matrix <- as.matrix(mutation_matrix_wide[, -1])
rownames(mut_matrix) <- genes

cat("Mutation matrix dimensions BEFORE mapping:\n")
cat("  Genes:", nrow(mut_matrix), "\n")
cat("  Samples:", ncol(mut_matrix), "\n\n")

# Map column names back to expression IDs (originalid with R suffix)
cat("Mapping sample IDs back to expression format...\n")
colnames_mapped <- id_mapping$original_id[match(colnames(mut_matrix),
                                                 id_mapping$mutation_id)]

# Check for any NAs in mapping
if (any(is.na(colnames_mapped))) {
  cat("WARNING:", sum(is.na(colnames_mapped)), "samples could not be mapped\n")
}

colnames(mut_matrix) <- colnames_mapped

cat("\nMutation matrix dimensions AFTER mapping:\n")
cat("  Genes:", nrow(mut_matrix), "\n")
cat("  Samples:", ncol(mut_matrix), "\n\n")

# Calculate mutation frequencies
mut_freq <- rowSums(mut_matrix) / ncol(mut_matrix) * 100
mut_freq_df <- data.frame(
  gene = names(mut_freq),
  n_mutated = rowSums(mut_matrix),
  n_total = ncol(mut_matrix),
  frequency_pct = mut_freq,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(frequency_pct))

cat("=== MUTATION FREQUENCIES ===\n")
print(mut_freq_df, row.names = FALSE)

cat("\nTop 10 most frequently mutated genes:\n")
print(head(mut_freq_df, 10), row.names = FALSE)

# Save results
saveRDS(mut_matrix, "03_Results/10_Mutations/mutation_matrix.rds")
write.csv(mut_freq_df, "03_Results/10_Mutations/mutation_frequencies.csv",
          row.names = FALSE)

cat("\n✓ Mutation matrix saved: 03_Results/10_Mutations/mutation_matrix.rds\n")
cat("✓ Mutation frequencies saved: 03_Results/10_Mutations/mutation_frequencies.csv\n")

cat("\n### Task 1.2 COMPLETE ###\n")
