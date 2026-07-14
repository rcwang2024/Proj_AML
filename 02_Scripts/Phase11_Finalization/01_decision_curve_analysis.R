# Phase 11: Task 11.1 - Decision Curve Analysis (DCA)
# Purpose: Quantify the Net Clinical Benefit of using VRS/Cluster vs standard mutations

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# 0. Setup Directories
dir.create("03_Results/Phase11_Finalization", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/Phase11_Finalization", recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
message(">>> Loading data for DCA...")
# We'll use the monocytic mapping results as a base (contains AUC and Cluster)
data <- read_csv("03_Results/Phase10_Analysis/10_1_Monocytic_Mapping_Results.csv")

# Load mutations (we need NPM1 and TP53 for the 'Standard of Care' model)
# I'll use the mutation matrix from Phase 10
mut_matrix <- readRDS("03_Results/01_Processed_Data/mutation_matrix_binary.rds")
mut_df <- as.data.frame(mut_matrix)
mut_df$sample_id <- gsub("D$", "R", rownames(mut_matrix))

# Merge
final_df <- data %>%
  inner_join(mut_df %>% select(sample_id, NPM1, TP53, FLT3, DNMT3A), by = "sample_id") %>%
  inner_join(read_csv("03_Results/25_Enhancements/vrs_9gene_scores.csv", show_col_types = FALSE) %>% select(sample_id, VRS9), by = "sample_id") %>%
  drop_na()

# 2. Define Outcome: "Venetoclax Sensitivity"
# Threshold based on the mid-point between cluster means (~150)
final_df$Response <- ifelse(final_df$auc < 150, 1, 0)
message("Response rate: ", mean(final_df$Response))

# 3. Decision Curve Functions
calculate_net_benefit <- function(probs, outcomes, thresholds) {
  sapply(thresholds, function(pt) {
    preds <- ifelse(probs >= pt, 1, 0)
    tp <- sum(preds == 1 & outcomes == 1)
    fp <- sum(preds == 1 & outcomes == 0)
    n <- length(outcomes)
    (tp / n) - (fp / n) * (pt / (1 - pt))
  })
}

# 4. Models
# Model 1: Genomic Only (NPM1 + TP53 + FLT3 + DNMT3A)
m1 <- glm(Response ~ NPM1 + TP53 + FLT3 + DNMT3A, data = final_df, family = "binomial")
prob_genomic <- predict(m1, type = "response")

# Model 2: Genomic + VRS9
m2 <- glm(Response ~ NPM1 + TP53 + FLT3 + DNMT3A + VRS9, data = final_df, family = "binomial")
prob_full <- predict(m2, type = "response")

# Threshold range (10% to 80%)
thresholds <- seq(0.01, 0.8, by = 0.01)

# 5. Calculate NB
nb_genomic <- calculate_net_benefit(prob_genomic, final_df$Response, thresholds)
nb_full <- calculate_net_benefit(prob_full, final_df$Response, thresholds)
nb_all <- (sum(final_df$Response) / nrow(final_df)) - ((nrow(final_df) - sum(final_df$Response)) / nrow(final_df)) * (thresholds / (1 - thresholds))
nb_none <- rep(0, length(thresholds))

# 6. Plotting
plot_df <- data.frame(
  Threshold = thresholds,
  Genomic = nb_genomic,
  Full = nb_full,
  All = nb_all,
  None = nb_none
) %>%
  pivot_longer(-Threshold, names_to = "Strategy", values_to = "NetBenefit")

p1 <- ggplot(plot_df, aes(x = Threshold, y = NetBenefit, color = Strategy)) +
  geom_line(size = 1.2) +
  coord_cartesian(ylim = c(-0.05, 0.6)) +
  scale_color_manual(values = c("Full" = "red", "Genomic" = "blue", "All" = "gray", "None" = "black"),
                     labels = c("Genomic + Cluster (VRS)", "Genomic Only (NPM1/TP53)", "Treat All", "Treat None")) +
  theme_minimal(base_size = 14) +
  labs(title = "Decision Curve Analysis: Clinical Utility of VRS",
       subtitle = "Net Benefit for predicting Venetoclax Sensitivity",
       x = "Threshold Probability (Pt)",
       y = "Net Benefit",
       color = "Prediction Strategy") +
  theme(legend.position = "bottom")

# 7. Save
ggsave("04_Figures/Phase11_Finalization/Figure11_1_DCA_NetBenefit.pdf", p1, width = 8, height = 6)
write_csv(plot_df, "03_Results/Phase11_Finalization/11_1_DCA_Results.csv")

message(">>> DCA Analysis Complete. High-impact evidence generated.")
