# Immune-Drug Correlation Analysis
# Objective: Correlate immune microenvironment with Venetoclax response
# Date: Dec 8, 2025

# 1. Setup
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
    library(immunedeconv)
    library(tibble)
    library(data.table)
})

BASE_DIR <- "d:/Projects/Project_AML"
setwd(BASE_DIR)

# Output Paths
dir.create("03_Results/25_Enhancements", showWarnings = FALSE, recursive = TRUE)
dir.create("04_Figures/25_Enhancements", showWarnings = FALSE, recursive = TRUE)

# 2. Load Data
# Expression (Batch Corrected)
expr_file <- "03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds"
expr_data <- readRDS(expr_file)

# Drug Data
drug_file <- "03_Results/01_Processed_Data/drug_response_auc.rds"
drug_data <- readRDS(drug_file)

# 3. Process Expression for Deconvolution
# Convert Ensembl to Symbols (MCP-counter needs symbols)
annot <- fread("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt",
    select = c("stable_id", "display_label")
)

gene_map <- annot %>%
    distinct(stable_id, display_label) %>%
    filter(display_label != "") %>%
    rename(ensembl = stable_id, symbol = display_label)

# Map matrix rownames
current_ids <- rownames(expr_data)
matched_syms <- gene_map$symbol[match(current_ids, gene_map$ensembl)]
keep <- !is.na(matched_syms)
expr_sym <- expr_data[keep, ]
rownames(expr_sym) <- matched_syms[keep]

# Aggregate duplicates (mean)
expr_sym_df <- as.data.frame(expr_sym) %>%
    rownames_to_column("gene") %>%
    group_by(gene) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("gene") %>%
    as.matrix()

# Linearize (if log)
if (max(expr_sym_df, na.rm = TRUE) < 50) {
    expr_lin <- 2^expr_sym_df - 1
} else {
    expr_lin <- expr_sym_df
}
# Ensure non-negative
expr_lin[expr_lin < 0] <- 0

# 4. Run MCP-counter
message("Running MCP-counter...")
mcp_res <- deconvolute(expr_lin, "mcp_counter")

# 5. Correlate with Venetoclax
# Prepare Drug Data
venetoclax <- drug_data %>%
    filter(inhibitor == "Venetoclax") %>%
    select(sample_id = dbgap_rnaseq_sample, auc)

# Reshape Immune Scores
immune_long <- mcp_res %>%
    pivot_longer(-cell_type, names_to = "sample_id", values_to = "score")

# Join
analysis_df <- inner_join(immune_long, venetoclax, by = "sample_id")

# Calculate Correlations
cor_results <- analysis_df %>%
    group_by(cell_type) %>%
    summarise(
        spearman_rho = cor(score, auc, method = "spearman", use = "complete.obs"),
        p_value = cor.test(score, auc, method = "spearman")$p.value,
        n = n()
    ) %>%
    arrange(spearman_rho)

# FDR
cor_results$fdr <- p.adjust(cor_results$p_value, method = "BH")

# Print
print(cor_results)

# Save
write_csv(cor_results, "03_Results/25_Enhancements/immune_drug_correlation.csv")

# 6. Plot Top Hits
if (nrow(cor_results) > 0) {
    # Plot top 2 significant hits
    top_hits <- head(cor_results %>% arrange(p_value), 2)

    for (i in 1:nrow(top_hits)) {
        cell <- top_hits$cell_type[i]
        rho <- top_hits$spearman_rho[i]

        p <- ggplot(analysis_df %>% filter(cell_type == cell), aes(x = score, y = auc)) +
            geom_point(alpha = 0.6) +
            geom_smooth(method = "lm", color = ifelse(rho < 0, "blue", "red")) +
            labs(
                title = paste("Immune vs Venetoclax:", cell),
                subtitle = paste("Rho =", round(rho, 3), "| p =", signif(top_hits$p_value[i], 3)),
                x = "MCP-counter Score", y = "Venetoclax AUC"
            ) +
            theme_minimal()

        ggsave(paste0("04_Figures/25_Enhancements/corr_", cell, "_venetoclax.pdf"), p, width = 5, height = 5)
    }
}
