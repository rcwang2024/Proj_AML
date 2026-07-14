# R script to run Noise Sensitivity Analysis on VRS
# Testing model stability under technical variation for Nature Medicine
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)

cat("=== STARTING VRS NOISE SENSITIVITY ANALYSIS ===\n\n")

# Setup Directories
dir.create("03_Results/28_VRS_Clinical_Utility", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/Phase11_Finalization", recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(42)

# 1. Load VRS data
vrs_data <- read_csv("03_Results/25_Enhancements/vrs_9gene_scores.csv", show_col_types = FALSE)
true_scores <- vrs_data$VRS9
true_sd <- sd(true_scores)

# Function to classify VRS into clinical tiers
classify_vrs <- function(scores) {
  case_when(
    scores < 41.8 ~ "Low",
    scores > 71.0 ~ "High",
    TRUE ~ "Medium"
  )
}

true_classes <- classify_vrs(true_scores)

# 2. Run simulation over noise levels (0% to 50% of true SD, in 5% increments)
noise_levels <- seq(0, 0.50, by = 0.05)
n_iterations <- 50

results_list <- list()

for (noise_pct in noise_levels) {
  concordances <- c()
  correlations <- c()
  
  noise_sd <- noise_pct * true_sd
  
  for (iter in 1:n_iterations) {
    if (noise_pct == 0) {
      concordances <- c(concordances, 1.0)
      correlations <- c(correlations, 1.0)
      next
    }
    
    # Add Gaussian noise
    noisy_scores <- true_scores + rnorm(length(true_scores), mean = 0, sd = noise_sd)
    # Clamp between 0 and 100
    noisy_scores <- pmax(pmin(noisy_scores, 100), 0)
    
    # Reclassify
    noisy_classes <- classify_vrs(noisy_scores)
    
    # Calculate concordance
    concordance <- mean(noisy_classes == true_classes)
    concordances <- c(concordances, concordance)
    
    # Calculate correlation
    corr <- cor(true_scores, noisy_scores, method = "pearson")
    correlations <- c(correlations, corr)
  }
  
  results_list[[paste0("N_", noise_pct)]] <- tibble(
    Noise_Level = noise_pct * 100, # as percentage
    Mean_Concordance = mean(concordances) * 100, # as percentage
    SD_Concordance = sd(concordances) * 100,
    Mean_Correlation = mean(correlations),
    SD_Correlation = sd(correlations)
  )
}

results_df <- bind_rows(results_list)
print(results_df)
cat("\n")

# Save results
write_csv(results_df, "03_Results/28_VRS_Clinical_Utility/VRS_noise_sensitivity_results.csv")
cat("✓ Saved data to: 03_Results/28_VRS_Clinical_Utility/VRS_noise_sensitivity_results.csv\n")

# === 3. PLOT RESULTS ===
cat("Generating noise sensitivity plots...\n")

# Melt data for ggplot
plot_df <- results_df %>%
  select(Noise_Level, Concordance = Mean_Concordance, Correlation = Mean_Correlation) %>%
  mutate(Correlation = Correlation * 100) %>%  # Scale to percentage
  pivot_longer(-Noise_Level, names_to = "Metric", values_to = "Value")

# Custom Theme
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 13, color = "darkblue", margin = margin(b=10)),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    axis.text = element_text(size = 11, face = "plain", color = "black"),
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 11, face = "plain", color = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    plot.margin = margin(15, 15, 15, 15),
    legend.position = "bottom"
  )

p_noise <- ggplot(plot_df, aes(x = Noise_Level, y = Value, color = Metric, group = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(
    values = c("Concordance" = "#3498DB", "Correlation" = "#E67E22"),
    labels = c("Classification Concordance (%)", "Pearson Correlation (r × 100)")
  ) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  scale_y_continuous(breaks = seq(50, 100, by = 10), limits = c(50, 105)) +
  labs(
    title = "9-Gene VRS Noise Sensitivity Analysis",
    subtitle = "Robustness of clinical stratification under simulated technical noise",
    x = "Simulated Noise Level (% of True SD)",
    y = "Metric Score (%)",
    color = "Evaluation Metric"
  ) +
  theme_hf +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 20.5, y = 55, label = "Typical Inter-Batch Variance (~20%)", 
           color = "red", hjust = 0, size = 3.5, fontface = "italic")

ggsave("04_Figures/Phase11_Finalization/VRS_noise_sensitivity.pdf", p_noise, width = 7.5, height = 5.5, device = cairo_pdf)
ggsave("04_Figures/Phase11_Finalization/VRS_noise_sensitivity.png", p_noise, width = 7.5, height = 5.5, dpi = 300)

cat("✓ Saved plots to: 04_Figures/Phase11_Finalization/VRS_noise_sensitivity.pdf\n")
cat("\n=== VRS NOISE SENSITIVITY ANALYSIS COMPLETE ===\n")
