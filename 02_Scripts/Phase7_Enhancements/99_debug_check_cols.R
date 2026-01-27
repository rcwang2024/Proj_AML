# Diagnostic
library(readr)
library(dplyr)

drug_file <- "03_Results/01_Processed_Data/drug_response_auc.rds"
if (file.exists(drug_file)) {
    d <- readRDS(drug_file)
    print("Drug Data Cols:")
    print(colnames(d))
    print(head(d))
}

expr_file <- "03_Results/05_Analysis_Ready_Data/expression_gold_standard.csv"
if (file.exists(expr_file)) {
    e <- read_csv(expr_file, n_max = 5, show_col_types = FALSE)
    print("Expression Data Cols:")
    print(colnames(e)[1:5])
}
