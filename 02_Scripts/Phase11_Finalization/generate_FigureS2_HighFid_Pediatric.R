# R script to generate High-Fidelity Supplementary Figure S2 (Pediatric TARGET-AML OS Validation)
setwd("d:/Proj_AML")
library(tidyverse)
library(survival)
library(survminer)
library(gridExtra)

cat("=== GENERATING HIGH-FIDELITY SUPPLEMENTARY FIGURE S2 (PEDIATRIC KM) ===\n")

# Canonical color palette
color_c1 <- "#3498DB" # Blue
color_c2 <- "#E67E22" # Orange

# Shared theme for KM plots
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 14, color = "darkblue", margin = margin(b=10)),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    axis.text = element_text(size = 12, face = "plain", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "plain", color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    plot.margin = margin(15, 15, 15, 15)
  )

# Load data
target_data <- read_csv("03_Results/18_TARGET_Validation/target_cluster_assignments.csv", show_col_types = FALSE) %>%
  filter(!is.na(OS_MONTHS) & !is.na(OS_STATUS) & !is.na(cluster)) %>%
  mutate(
    predicted_cluster = factor(cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2"))
  )

# Compute HR annotation string
surv_obj <- Surv(target_data$OS_MONTHS, target_data$OS_STATUS)
cox_fit <- coxph(surv_obj ~ target_data$predicted_cluster)
s <- summary(cox_fit)
hr  = round(s$conf.int[1, 1], 2)
lo  = round(s$conf.int[1, 3], 2)
hi  = round(s$conf.int[1, 4], 2)
p   = s$coefficients[1, 5]
p_str = ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                formatC(p, format = "f", digits = 3))
hr_label <- sprintf("HR = %s (95%% CI: %s\u2013%s)\nLog-rank p = %s", hr, lo, hi, p_str)

# Compute median OS string
km_fit <- survfit(surv_obj ~ predicted_cluster, data = target_data)
med <- summary(km_fit)$table
m1 <- ifelse(is.na(med[1, "median"]), "NR", round(med[1, "median"], 1))
m2 <- ifelse(is.na(med[2, "median"]), "NR", round(med[2, "median"], 1))
median_label <- sprintf("Median OS: %s vs %s mo", m1, m2)

# Build Kaplan-Meier Plot
p <- ggsurvplot(
  km_fit, data = target_data,
  palette           = c(color_c1, color_c2),
  linetype          = "solid",
  conf.int          = TRUE,
  conf.int.alpha    = 0.15,
  pval              = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.y.text = FALSE,
  risk.table.fontsize = 4.0,
  tables.theme      = theme_cleantable(),
  legend.title      = "Subtype",
  legend.labs       = c(sprintf("Cluster 1 (n=%d)", sum(target_data$cluster == 1)),
                        sprintf("Cluster 2 (n=%d)", sum(target_data$cluster == 2))),
  title             = "Pediatric TARGET-AML Validation Cohort (N=1,713)",
  xlab              = "Time (months)",
  ylab              = "Overall Survival",
  ggtheme           = theme_hf
)

p$plot <- p$plot +
  annotate("text", x = Inf, y = 0.97, hjust = 1.03, vjust = 1,
           label = hr_label,
           size = 4.5, fontface = "bold", color = "grey20") +
  annotate("text", x = Inf, y = 0.80, hjust = 1.03, vjust = 1,
           label = median_label,
           size = 4.0, fontface = "italic", color = "grey40")

# Save outputs
dir.create("04_Figures/18_TARGET_Validation/HighFid", showWarnings = FALSE, recursive = TRUE)

# Save PDF and PNG
pdf("04_Figures/18_TARGET_Validation/HighFid/FigureS2_Pediatric_KM.pdf", width = 8, height = 7)
print(p, newpage = FALSE)
dev.off()

png("04_Figures/18_TARGET_Validation/HighFid/FigureS2_Pediatric_KM.png", width = 800, height = 700, res = 120)
print(p, newpage = FALSE)
dev.off()

cat("✓ High-fidelity TARGET-AML Kaplan-Meier plot generated.\n")
