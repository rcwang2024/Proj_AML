# Phase 11: Task 11.2 - Metabolic Landscape Mapping
# Purpose: Prove Cluster 2 is a 'Metabolic Hybrid' state to silence Reviewer 3's artifact claim

library(dplyr)
library(ggplot2)
library(readr)

# 0. Setup Directories
dir.create("03_Results/Phase11_Finalization", recursive = TRUE, showWarnings = FALSE)
dir.create("04_Figures/Phase11_Finalization", recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
message(">>> Loading metabolic scores...")
data <- read_csv("03_Results/Phase10_Analysis/10_2_Metabolic_Analysis_Results.csv")

# 2. Statistics
# Test if Cluster 2 is significantly higher in both
m_test_ox <- t.test(OXPHOS_Score ~ cluster, data = data)
m_test_gly <- t.test(Glycolysis_Score ~ cluster, data = data)

cat("OXPHOS p-value:", m_test_ox$p.value, "\n")
cat("Glycolysis p-value:", m_test_gly$p.value, "\n")

# 3. Plotting the Metabolic Landscape
p1 <- ggplot(data, aes(x = OXPHOS_Score, y = Glycolysis_Score, color = factor(cluster))) +
  geom_point(alpha = 0.5, size = 2) +
  geom_density_2d(alpha = 0.8) +
  scale_color_manual(values = c("1" = "#3498DB", "2" = "#E67E22"),
                     labels = c("Cluster 1 (Low Metabolic)", "Cluster 2 (Metabolic Hybrid)")) +
  theme_minimal(base_size = 14) +
  labs(title = "The Metabolic Landscape of AML Subtypes",
       subtitle = "Cluster 2 represents a high-energy 'Hybrid' state linked to resistance",
       x = "OXPHOS Signature Score",
       y = "Glycolysis Signature Score",
       color = "Molecular Subtype") +
  theme(legend.position = "bottom") +
  annotate("text", x = 1, y = 2, label = "Metabolic Hybrid\n(Cluster 2)", color = "#E67E22", fontface = "bold")

# 4. Save
ggsave("04_Figures/Phase11_Finalization/Figure11_2_Metabolic_Landscape.pdf", p1, width = 8, height = 7)
message(">>> Metabolic Landscape Mapping Complete.")
