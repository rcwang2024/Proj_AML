# VRS Development Script
# Objective: Develop continuous Venetoclax Response Score (0-100)
# Date: Dec 8, 2025

# VRS Development Script
# Objective: Develop continuous Venetoclax Response Score (0-100)
# Date: Dec 8, 2025

# 1. Setup
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
    library(data.table) # Faster reading
})

# Set paths
BASE_DIR <- "d:/Projects/Project_AML"
setwd(BASE_DIR)

# Create output dirs
dir.create("03_Results/25_Enhancements", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/25_Enhancements", recursive = TRUE, showWarnings = FALSE)

# 2. Load Data
# Signature Weights
sig_file <- "03_Results/15_Gene_Signature/50_gene_signature.csv"
if (!file.exists(sig_file)) stop("Signature file not found")
signature <- read_csv(sig_file, show_col_types = FALSE)

# Expression Data (Optimized with fread)
expr_file <- "03_Results/05_Analysis_Ready_Data/expression_gold_standard.csv"
if (!file.exists(expr_file)) stop("Expression file not found")
# Read with data.table
expression <- fread(expr_file)

# Drug Data
drug_file <- "03_Results/01_Processed_Data/drug_response_auc.rds"
if (!file.exists(drug_file)) stop("Drug file not found")
drug_data <- readRDS(drug_file)

# 3. Calculate Scores

# A. Prepare Weights
# Filter expression to signature genes
sig_genes <- signature$gene
# Check overlap
valid_genes <- intersect(sig_genes, colnames(expression))
message(paste("Found", length(valid_genes), "signature genes in expression data"))

expr_subset <- expression %>% select(matches(paste0("^(", paste(valid_genes, collapse = "|"), ")$")))
patient_ids <- expression$samples # Check for sample column name
if ("sample_id" %in% colnames(expression)) {
    patient_ids <- expression$sample_id
} else if ("...1" %in% colnames(expression)) {
    patient_ids <- expression$...1
} else {
    # If no ID column found, assume first column is ID if string, or rownames
    # data.table doesn't leverage rownames.
    # Check first column name
    first_col <- colnames(expression)[1]
    message(paste("Using first column as ID:", first_col))
    patient_ids <- expression[[1]]
}

# Standardize Expression (Z-score per gene across cohort)
expr_z <- scale(expr_subset)

# Get Venetoclax Reference for Directionality Check
venetoclax_ref <- drug_data %>%
    filter(inhibitor == "Venetoclax") %>%
    select(sample_id = dbgap_rnaseq_sample, auc)

# Create a temporary df for correlation
temp_df <- data.frame(sample_id = patient_ids, expr_z) %>%
    inner_join(venetoclax_ref, by = "sample_id")

# Calculate gene-specific correlation with AUC
common_genes <- colnames(expr_subset)
gene_cors <- sapply(common_genes, function(g) {
    cor(temp_df[[g]], temp_df$auc, use = "complete.obs", method = "spearman")
})

# Define Weights:
# We want High Score = High Sensitivity (Low AUC).
# Weight = Importance * -1 * Sign(Correlation)
sig_subset <- signature %>% filter(gene %in% common_genes)
weights_map <- setNames(sig_subset$importance, sig_subset$gene)

final_weights <- sapply(common_genes, function(g) {
    imp <- weights_map[[g]]
    dir <- sign(gene_cors[[g]])
    return(imp * -1 * dir)
})

# Calculate Raw Score
raw_scores <- as.matrix(expr_z) %*% final_weights
vrs_df <- data.frame(
    sample_id = patient_ids,
    raw_score = as.numeric(raw_scores)
)

# Rescale 0-100
min_s <- min(vrs_df$raw_score, na.rm = TRUE)
max_s <- max(vrs_df$raw_score, na.rm = TRUE)
vrs_df$VRS <- ((vrs_df$raw_score - min_s) / (max_s - min_s)) * 100

# 4. Validate
validation <- inner_join(vrs_df, venetoclax_ref, by = "sample_id")

# Statistics
cor_res <- cor.test(validation$VRS, validation$auc, method = "spearman")
message(paste("Spearman Correlation:", round(cor_res$estimate, 3)))
message(paste("P-value:", signif(cor_res$p.value, 3)))

# 5. Plots
p <- ggplot(validation, aes(x = VRS, y = auc)) +
    geom_point(alpha = 0.6, color = "#2C3E50") +
    geom_smooth(method = "lm", color = "#E74C3C") +
    theme_minimal() +
    labs(
        title = "Venetoclax Response Score (VRS) Validation",
        subtitle = paste("Spearman rho =", round(cor_res$estimate, 3), "| p =", signif(cor_res$p.value, 3)),
        x = "VRS (0 = Resistant, 100 = Sensitive)",
        y = "Actual Venetoclax AUC (Lower = Better)"
    )

ggsave("04_Figures/25_Enhancements/VRS_validation_plot.pdf", p, width = 6, height = 5)
ggsave("04_Figures/25_Enhancements/VRS_validation_plot.png", p, width = 6, height = 5)

# Save Scores
write_csv(vrs_df, "03_Results/25_Enhancements/venetoclax_response_scores.csv")

# 6. Output Summary (Corrected)
summary_file <- "03_Results/25_Enhancements/VRS_summary.txt"
sink(summary_file)
cat("Venetoclax Response Score (VRS) Stats\n")
cat("-------------------------------------\n")
cat("N samples scored:", nrow(vrs_df), "\n")
cat("N validated:", nrow(validation), "\n")
cat("Correlation (Rho):", round(cor_res$estimate, 3), "\n")
cat("P-value:", signif(cor_res$p.value, 3), "\n")
cat("\nInterpretation:\n")
if (cor_res$estimate < -0.3) {
    cat("SUCCESS: Strong negative correlation. High VRS predicts Low AUC (Sensitivity).\n")
} else {
    cat("WARNING: Correlation weak. Check weight directions.\n")
}
sink()
