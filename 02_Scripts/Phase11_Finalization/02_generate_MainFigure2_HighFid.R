# R script to generate High-Fidelity Main Figure 2 (Survival) with full significance
setwd("d:/Proj_AML")
library(tidyverse)
library(survival)
library(survminer)
library(gridExtra)
library(grid)

cat("=== GENERATING HIGH-FIDELITY MAIN FIGURE 2 (WITH SIGNIFICANCE) ===\n")

# Canonical color palette
color_c1 <- "#3498DB" # Blue
color_c2 <- "#E67E22" # Orange

# Helper: compute HR annotation string
compute_hr_label <- function(surv_data, time_var, event_var) {
  surv_obj <- Surv(surv_data[[time_var]], surv_data[[event_var]])
  cox_fit <- coxph(surv_obj ~ surv_data$predicted_cluster)
  s <- summary(cox_fit)
  hr  <- round(s$conf.int[1, 1], 2)
  lo  <- round(s$conf.int[1, 3], 2)
  hi  <- round(s$conf.int[1, 4], 2)
  p   <- s$coefficients[1, 5]
  p_str <- ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                  formatC(p, format = "f", digits = 3))
  sprintf("HR = %s (95%% CI: %s\u2013%s)\nLog-rank p = %s", hr, lo, hi, p_str)
}

# Helper: compute median OS string
compute_median_os <- function(km_fit) {
  med <- summary(km_fit)$table
  if (is.matrix(med)) {
    m1 <- ifelse(is.na(med[1, "median"]), "NR", round(med[1, "median"], 1))
    m2 <- ifelse(is.na(med[2, "median"]), "NR", round(med[2, "median"], 1))
  } else {
    m1 <- "NR"; m2 <- "NR"
  }
  sprintf("Median OS: %s vs %s mo", m1, m2)
}

# Shared theme for KM plots
theme_km <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 14, color = "darkblue", margin = margin(b=10)),
    axis.title.x = element_text(size = 13, face = "bold", color = "black", margin = margin(t=4)),
    axis.title.y = element_text(size = 13, face = "bold", color = "black", margin = margin(r=2, l=0)),
    axis.text = element_text(size = 12, face = "plain", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "plain", color = "black"),
    panel.grid.minor = element_blank()
  )

# Shared ggsurvplot generator
make_km_panel <- function(data, time_var, event_var, title_str) {
  # Rename columns temporarily for ggsurvplot compatibility
  data$.time <- data[[time_var]]
  data$.event <- data[[event_var]]

  km_fit <- survfit(Surv(.time, .event) ~ predicted_cluster, data = data)
  hr_label <- compute_hr_label(data, time_var, event_var)
  median_label <- compute_median_os(km_fit)

  p <- ggsurvplot(
    km_fit, data = data,
    palette         = c(color_c1, color_c2),
    linetype        = "solid",
    conf.int        = TRUE,
    conf.int.alpha  = 0.15,
    pval            = TRUE,
    pval.method     = TRUE,
    pval.size       = 3.8,
    risk.table      = TRUE,
    risk.table.height = 0.22,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 3.2,
    tables.theme    = theme_cleantable(),
    ncensor.plot    = FALSE,
    legend.title    = "Subtype",
    legend.labs     = c("Cluster 1", "Cluster 2"),
    title           = title_str,
    xlab            = "Time (months)",
    ylab            = "Overall Survival",
    ggtheme         = theme_km
  )

  # Add HR + median OS annotation
  p$plot <- p$plot +
    annotate("text", x = Inf, y = 0.15, hjust = 1.05, vjust = 0,
             label = hr_label,
             size = 3.3, fontface = "bold", color = "grey30") +
    annotate("text", x = Inf, y = 0.03, hjust = 1.05, vjust = 0,
             label = median_label,
             size = 3.0, fontface = "italic", color = "grey40")

  return(p)
}

# ===========================================================================
# Panel A: BeatAML Discovery
# ===========================================================================
beat_data <- read_csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv",
                      show_col_types = FALSE) %>%
  filter(!is.na(OS_months) & !is.na(OS_event) & !is.na(cluster)) %>%
  mutate(
    # NOTE: The 'OS_months' column in this file is actually in DAYS
    OS_months = OS_months / 30.4375,
    predicted_cluster = factor(cluster, levels = c(1, 2),
                               labels = c("Cluster 1", "Cluster 2"))
  )
p2a <- make_km_panel(beat_data, "OS_months", "OS_event",
                     "A. BeatAML Discovery (N=320)")

# ===========================================================================
# Panel B: TCGA-LAML Validation
# ===========================================================================
tcga_data <- read_csv("03_Results/17_TCGA_Validation/tcga_clinical_with_predictions.csv",
                      show_col_types = FALSE) %>%
  filter(!is.na(OS_months) & !is.na(OS_event) & !is.na(predicted_cluster)) %>%
  mutate(predicted_cluster = factor(predicted_cluster, levels = c(1, 2),
                                    labels = c("Cluster 1", "Cluster 2")))
p2b <- make_km_panel(tcga_data, "OS_months", "OS_event",
                     "B. TCGA-LAML Validation (N=151)")

# ===========================================================================
# Panel C: Prognostic Meta-Analysis Forest Plot
# ===========================================================================
meta_df <- data.frame(
  Cohort = c("A. BeatAML Discovery (N=320)", "B. TCGA-LAML Validation (N=151)", "C. Prognostic Meta-Analysis (Pooled Adult)"),
  HR = c(1.03, 1.24, 1.26),
  Lower = c(0.80, 0.80, 1.02),
  Upper = c(1.33, 1.94, 1.55),
  P_Value = c("0.830", "0.342", "0.030"),
  Type = c("Individual", "Individual", "Pooled"),
  y = c(3, 2, 1)
)

# Soft blue (#2E5B88) and soft red (#C0392B)
p_forest <- ggplot(meta_df, aes(x = HR, y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper, color = Type), width = 0.15, linewidth = 0.8) +
  geom_point(data = filter(meta_df, Type == "Individual"), aes(color = Type), shape = 16, size = 4) +
  geom_point(data = filter(meta_df, Type == "Pooled"), aes(color = Type), shape = 18, size = 5.5) +
  scale_color_manual(values = c("Individual" = "#2E5B88", "Pooled" = "#C0392B")) +
  scale_y_continuous(
    breaks = meta_df$y,
    labels = meta_df$Cohort,
    limits = c(0.5, 3.5)
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.8, 1.0, 1.26, 1.5, 2.0, 2.5),
    labels = c("0.5", "0.8", "1.0", "1.26", "1.5", "2.0", "2.5"),
    limits = c(0.5, 3.0)
  ) +
  geom_text(
    aes(x = 2.1, label = sprintf("HR = %.2f (95%% CI: %.2f\u2013%.2f), p = %s", HR, Lower, Upper, P_Value)),
    hjust = 0, size = 3.2, fontface = "bold", color = "grey20"
  ) +
  xlab("Hazard Ratio (95% CI, Log Scale)") +
  ylab("") +
  ggtitle("C. Prognostic Meta-Analysis: Adult General Cohorts") +
  theme_classic(base_size = 11) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 8)),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "black"),
    legend.position = "none",
    plot.margin = margin(t = 15, r = 25, b = 15, l = 25)
  )

# ===========================================================================
# Compile 3-panel grid (Panels A, B, C)
# ===========================================================================
km_grid <- arrange_ggsurvplots(
  list(p2a, p2b),
  print = FALSE,
  ncol = 2, nrow = 1,
  risk.table.height = 0.22
)
km_grob <- do.call(arrangeGrob, km_grid)

# Save outputs
dir.create("05_Submission/Submission_Hub/02_Main_Figures", showWarnings = FALSE, recursive = TRUE)

pdf("05_Submission/Submission_Hub/02_Main_Figures/Deprecated_Survival_Figure2.pdf",
    width = 12, height = 9)
grid.arrange(km_grob, p_forest, ncol = 1, heights = c(2, 1.2))
dev.off()

png("05_Submission/Submission_Hub/02_Main_Figures/Deprecated_Survival_Figure2.png",
    width = 1200, height = 900, res = 150)
grid.arrange(km_grob, p_forest, ncol = 1, heights = c(2, 1.2))
dev.off()

cat("✓ High-fidelity Main Figure 2 with 3 panels (Prognostic Meta-Analysis Forest Plot) generated successfully.\n")
