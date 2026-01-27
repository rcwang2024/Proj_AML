# Synergy Prediction Analysis
# Objective: Identify drugs effective in Venetoclax-resistant AML (Cluster 2)
# Focus: Synthetic Lethality & Inflammation Targets (MEK/JAK)
# Date: Dec 8, 2025

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
    library(data.table)
})

BASE_DIR <- "d:/Projects/Project_AML"
setwd(BASE_DIR)

# Output Paths
dir.create("03_Results/25_Enhancements", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/25_Enhancements", showWarnings = FALSE, recursive = TRUE)

# 1. Load Data
drug_data <- readRDS("03_Results/01_Processed_Data/drug_response_auc.rds")
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# 2. Define Venetoclax Resistance
# Get Venetoclax AUCs
venetoclax <- drug_data %>%
    filter(inhibitor == "Venetoclax") %>%
    select(sample_id = dbgap_rnaseq_sample, ven_auc = auc)

# Merge with Cluster info
cohort <- inner_join(venetoclax, clusters, by = "sample_id")

# Define "Resistant" Group (Cluster 2 OR High AUC)
# We strictly want to solve Cluster 2, so let's filter for Cluster 2 samples.
cluster2_samples <- cohort %>%
    filter(cluster == 2) %>%
    pull(sample_id)
message(paste("Analyzing", length(cluster2_samples), "Cluster 2 (Resistant) patients..."))

# 3. Screen for "Salvage" Drugs
# Find drugs that are effective (Low AUC) specifically in Cluster 2
# Strategy: Calculate mean AUC for every drug in Cluster 2.
# Compare to Venetoclax mean AUC in Cluster 2.

# Filter drug data to Cluster 2
c2_drug_data <- drug_data %>%
    filter(dbgap_rnaseq_sample %in% cluster2_samples) %>%
    rename(sample_id = dbgap_rnaseq_sample)

# Calculate Stats per Drug
drug_stats <- c2_drug_data %>%
    group_by(inhibitor) %>%
    summarise(
        mean_auc = mean(auc, na.rm = TRUE),
        median_auc = median(auc, na.rm = TRUE),
        n = n(),
        # Effect Size vs Venetoclax (Is it better?)
        # We compare this drug's AUC distribution to Venetoclax's AUC distribution IN CLUSTER 2.
        p_vs_ven = tryCatch(wilcox.test(auc, cohort$ven_auc[cohort$cluster == 2], alternative = "less")$p.value, error = function(e) NA)
    ) %>%
    filter(n >= 20) %>% # Minimum sample size
    arrange(mean_auc)

# 4. Mechanism Matching
# Highlight targets of interest: MEK, JAK, MCL1, BCL2, FLT3
targets_of_interest <- c(
    "Trametinib", "Cobimetinib", "Selumetinib", # MEK
    "Ruxolitinib", "Tofacitinib", # JAK
    "MCL1", "S63845", "A-1210477", # MCL1 (Check names)
    "Sorafenib", "Midostaurin", # FLT3 (often co-mutated)
    "Idelalisib", "Duvelisib", # PI3K/mTOR (Myeloid signaling)
    "Panobinostat"
) # HDAC (Broad)

# Annotate
drug_stats <- drug_stats %>%
    mutate(
        mechanism = case_when(
            grepl("tinib", inhibitor) ~ "Kinase Inhibitor",
            grepl("MCL|S63845", inhibitor, ignore.case = TRUE) ~ "MCL1 Inhibitor",
            inhibitor %in% targets_of_interest ~ "Priority Candidate",
            TRUE ~ "Other"
        ),
        is_priority = inhibitor %in% targets_of_interest
    ) %>%
    arrange(desc(is_priority), mean_auc)

# Print Top Hits
message("Top Candidates for Cluster 2 (Ranked by Efficacy/Low AUC):")
print(head(drug_stats %>% select(inhibitor, mean_auc, p_vs_ven), 15))

# Save
write_csv(drug_stats, "03_Results/25_Enhancements/cluster2_salvage_drugs.csv")

# 5. Synergy Plot (Venetoclax AUC vs Candidate AUC)
# Select top candidate (lowest AUC among priorities)
top_cand <- drug_stats %>%
    filter(is_priority) %>%
    head(1) %>%
    pull(inhibitor)

if (length(top_cand) > 0) {
    cand_data <- drug_data %>%
        filter(inhibitor == top_cand) %>%
        select(sample_id = dbgap_rnaseq_sample, cand_auc = auc)

    plot_df <- inner_join(cohort, cand_data, by = "sample_id")

    p <- ggplot(plot_df, aes(x = ven_auc, y = cand_auc, color = as.factor(cluster))) +
        geom_point(alpha = 0.6) +
        geom_vline(xintercept = median(cohort$ven_auc), linetype = "dashed", color = "gray") +
        geom_hline(yintercept = median(plot_df$cand_auc), linetype = "dashed", color = "gray") +
        scale_color_manual(values = c("#2E9FDF", "#E7B800"), labels = c("Sensitive", "Resistant")) +
        labs(
            title = paste("Synergy Potential: Venetoclax vs", top_cand),
            subtitle = "Goal: Find drugs active (Low Y) where Venetoclax fails (High X)",
            x = "Venetoclax AUC", y = paste(top_cand, "AUC"),
            color = "Cluster"
        ) +
        theme_minimal()

    ggsave(paste0("04_Figures/25_Enhancements/synergy_ven_", top_cand, ".pdf"), p, width = 6, height = 5)
}
