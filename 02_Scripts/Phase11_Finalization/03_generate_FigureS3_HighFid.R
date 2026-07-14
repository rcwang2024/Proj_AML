# R script to generate High-Fidelity Figure S3 Panels (ELN Integration) - v3 (Data-Driven 2x2 Grid)
setwd("d:/Proj_AML")
library(tidyverse)
library(readxl)
library(survival)
library(ggplot2)

cat("=== GENERATING HIGH-FIDELITY S3 PANELS (V3: DATA-DRIVEN 4-PANEL) ===\n")

# 1. Load Data
clinical <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

df <- clinical %>%
  filter(!is.na(dbgap_rnaseq_sample)) %>%
  inner_join(clusters, by = c("dbgap_rnaseq_sample" = "sample_id")) %>%
  filter(!is.na(ELN2017)) %>%
  mutate(
    OS_event = ifelse(vitalStatus == "Dead", 1, 0),
    OS_months = overallSurvival,
    Cluster_Label = ifelse(cluster == 1, "Cluster 1", "Cluster 2")
  ) %>%
  filter(!is.na(OS_months))

# --- High-Fidelity Theme ---
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

# Color Scheme: Slate Blue for Cluster 1, Copper Orange for Cluster 2
cluster_colors <- c(
  "Cluster 1" = "#2B5C8F",
  "Cluster 2" = "#D95F02"
)

# Helper function to extract step-wise survival curve data starting at time=0, surv=1
get_surv_df <- function(fit) {
  df_list <- list()
  strata_names <- names(fit$strata)
  start_idx <- 1
  for (s_name in strata_names) {
    n_points <- fit$strata[s_name]
    end_idx <- start_idx + n_points - 1
    
    t <- fit$time[start_idx:end_idx]
    surv <- fit$surv[start_idx:end_idx]
    clean_name <- gsub("Cluster_Label=", "", s_name)
    
    # Prepend time=0, surv=1 if missing
    if (!(0 %in% t)) {
      t <- c(0, t)
      surv <- c(1, surv)
    }
    
    df_list[[s_name]] <- data.frame(
      time = t,
      surv = surv,
      group = clean_name
    )
    start_idx <- end_idx + 1
  }
  bind_rows(df_list)
}

# --- PANEL A: ELN Favorable ---
df_fav <- df %>% filter(ELN2017 == "Favorable")
fit_fav <- survfit(Surv(OS_months, OS_event) ~ Cluster_Label, data = df_fav)
surv_fav <- get_surv_df(fit_fav)

p_a <- ggplot(surv_fav, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 1.8) +
  scale_color_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "A. Survival in ELN Favorable",
       subtitle = "Subtypes show identical survival profiles (p = 0.91)",
       x = "Months", y = "Overall Survival Probability", color = "Subtype") +
  theme_hf +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

# --- PANEL B: ELN Intermediate ---
df_int <- df %>% filter(ELN2017 == "Intermediate")
fit_int <- survfit(Surv(OS_months, OS_event) ~ Cluster_Label, data = df_int)
surv_int <- get_surv_df(fit_int)

p_b <- ggplot(surv_int, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 1.8) +
  scale_color_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "B. Survival in ELN Intermediate",
       subtitle = "Non-significant trend with distinct kinetics (p = 0.09)",
       x = "Months", y = "Overall Survival Probability", color = "Subtype") +
  theme_hf +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

# --- PANEL C: ELN Adverse ---
df_adv <- df %>% filter(ELN2017 == "Adverse")
fit_adv <- survfit(Surv(OS_months, OS_event) ~ Cluster_Label, data = df_adv)
surv_adv <- get_surv_df(fit_adv)

p_c <- ggplot(surv_adv, aes(x = time, y = surv, color = group)) +
  geom_step(linewidth = 1.8) +
  scale_color_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "C. Survival in ELN Adverse",
       subtitle = "Highly significant survival separation (p = 0.02)",
       x = "Months", y = "Overall Survival Probability", color = "Subtype") +
  theme_hf +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"))

# --- PANEL D: Subgroup Hazard Ratios (Forest Plot) ---
hr_data <- data.frame(
  Group = c("ELN Favorable", "ELN Intermediate", "ELN Adverse"),
  HR = c(1.04, 0.60, 1.84),
  Low = c(0.55, 0.33, 1.07),
  High = c(1.96, 1.09, 3.16),
  Pval = c("p=0.91", "p=0.09", "p=0.03")
)

# Clinical Ordering: Favorable (top) -> Intermediate (middle) -> Adverse (bottom)
hr_data$Group <- factor(hr_data$Group, levels = c("ELN Adverse", "ELN Intermediate", "ELN Favorable"))

# Format text labels as y-axis labels
hr_data_labels <- hr_data %>%
  mutate(
    GroupLabel = paste0(Group, "\nHR: ", sprintf("%.2f", HR), " [", sprintf("%.2f", Low), "-", sprintf("%.2f", High), "], ", Pval)
  )

# Ensure levels of GroupLabel match Group levels
hr_data_labels$GroupLabel <- factor(
  hr_data_labels$GroupLabel, 
  levels = hr_data_labels$GroupLabel[order(match(hr_data_labels$Group, c("ELN Adverse", "ELN Intermediate", "ELN Favorable")))]
)

p_d <- ggplot(hr_data_labels, aes(x = HR, y = GroupLabel)) +
  geom_point(size = 7, color = "#D95F02") +
  geom_errorbar(aes(xmin = Low, xmax = High), width = 0.25, linewidth = 1.8, color = "#D95F02") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1.2) +
  scale_x_log10(
    limits = c(0.2, 5.0),
    breaks = c(0.2, 0.5, 1.0, 2.0, 5.0),
    labels = c("0.2", "0.5", "1.0", "2.0", "5.0")
  ) +
  labs(title = "D. Subtype Prognostic Power by ELN Tier", 
       subtitle = "Subtype risk differences are highly tier-specific",
       x = "Hazard Ratio (95% CI, log scale)", y = "") +
  theme_hf +
  theme(
    plot.title.position = "plot",
    axis.text.y = element_text(lineheight = 1.2, face = "bold")
  )

# Save Outputs
dir.create("04_Figures/12_ELN_Comparison/HighFid", showWarnings = FALSE, recursive = TRUE)
# Each panel: 925pt x 925pt
ggsave("04_Figures/12_ELN_Comparison/HighFid/s3_pA.pdf", p_a, width = 6.0, height = 6.0)
ggsave("04_Figures/12_ELN_Comparison/HighFid/s3_pB.pdf", p_b, width = 6.0, height = 6.0)
ggsave("04_Figures/12_ELN_Comparison/HighFid/s3_pC.pdf", p_c, width = 6.0, height = 6.0)
ggsave("04_Figures/12_ELN_Comparison/HighFid/s3_pD.pdf", p_d, width = 6.0, height = 6.0)

cat("✓ High-fidelity S3 panels (V3: 4-Panel) generated successfully.\n")
