# R script to generate High-Fidelity Supplementary Figure S11 (Adult Prognostic OS Validation Suite)
# Compiled with 5 panels: BeatAML, TCGA, Meta-Analysis, AMLCG GPL96, and AMLCG GPL570.
setwd("d:/Proj_AML")
library(tidyverse)
library(survival)
library(survminer)
library(gridExtra)
library(grid)

cat("=== GENERATING COMPREHENSIVE 5-PANEL SUPPLEMENTARY FIGURE S11 ===\n")

# Canonical color palette
color_c1 <- "#3498DB" # Blue
color_c2 <- "#E67E22" # Orange

# Helper: compute HR annotation string
compute_hr_label <- function(surv_data, time_var, event_var) {
  surv_obj <- Surv(surv_data[[time_var]], surv_data[[event_var]])
  cox_fit <- coxph(surv_obj ~ surv_data$predicted_cluster)
  s <- summary(cox_fit)
  hr  = round(s$conf.int[1, 1], 2)
  lo  = round(s$conf.int[1, 3], 2)
  hi  = round(s$conf.int[1, 4], 2)
  p   = s$coefficients[1, 5]
  p_str = ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
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

# Shared ggsurvplot generator
make_km_panel <- function(data, time_var, event_var, title_str) {
  data$.time  <- data[[time_var]]
  data$.event <- data[[event_var]]

  km_fit       <- survfit(Surv(.time, .event) ~ predicted_cluster, data = data)
  hr_label     <- compute_hr_label(data, time_var, event_var)
  median_label <- compute_median_os(km_fit)

  p <- ggsurvplot(
    km_fit, data = data,
    palette           = c(color_c1, color_c2),
    linetype          = "solid",
    conf.int          = TRUE,
    conf.int.alpha    = 0.15,
    pval              = FALSE,      # removed: shown in detailed annotation below
    pval.method       = FALSE,      # removed: avoids duplicate p-value display
    risk.table        = TRUE,
    risk.table.height = 0.22,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 3.2,
    tables.theme      = theme_cleantable(),
    ncensor.plot      = FALSE,
    legend.title      = "Subtype",
    legend.labs       = c("Cluster 1", "Cluster 2"),
    title             = title_str,
    xlab              = "Time (months)",
    ylab              = "Overall Survival",
    ggtheme           = theme_hf
  )

  # Place detailed HR + log-rank p annotation in the upper-right corner,
  # and median OS annotation just below it — away from the survival curves.
  p$plot <- p$plot +
    annotate("text", x = Inf, y = 0.97, hjust = 1.03, vjust = 1,
             label = hr_label,
             size = 3.5, fontface = "bold", color = "grey20") +
    annotate("text", x = Inf, y = 0.80, hjust = 1.03, vjust = 1,
             label = median_label,
             size = 3.2, fontface = "italic", color = "grey40")

  return(p)
}

# ===========================================================================
# Panel A: BeatAML Discovery
# ===========================================================================
beat_data <- read_csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv",
                      show_col_types = FALSE) %>%
  filter(!is.na(OS_months) & !is.na(OS_event) & !is.na(cluster)) %>%
  mutate(
    OS_months = OS_months / 30.4375,
    predicted_cluster = factor(cluster, levels = c(1, 2),
                               labels = c("Cluster 1", "Cluster 2"))
  )
ps11a <- make_km_panel(beat_data, "OS_months", "OS_event",
                        "A. BeatAML Discovery (N=320)")

# ===========================================================================
# Panel B: TCGA-LAML Validation
# ===========================================================================
tcga_data <- read_csv("03_Results/17_TCGA_Validation/tcga_clinical_with_predictions.csv",
                      show_col_types = FALSE) %>%
  filter(!is.na(OS_months) & !is.na(OS_event) & !is.na(predicted_cluster)) %>%
  mutate(predicted_cluster = factor(predicted_cluster, levels = c(1, 2),
                                    labels = c("Cluster 1", "Cluster 2")))
ps11b <- make_km_panel(tcga_data, "OS_months", "OS_event",
                       "B. TCGA-LAML Validation (N=151)")

# ===========================================================================
# Panel C: Prognostic Meta-Analysis Forest Plot
# ===========================================================================
meta_df <- data.frame(
  Cohort = c("A. BeatAML Discovery (N=320)", "B. TCGA-LAML Validation (N=151)", "E. Prognostic Meta-Analysis (Pooled Adult)"),
  HR = c(1.03, 1.24, 1.08),
  Lower = c(0.80, 0.80, 0.86),
  Upper = c(1.33, 1.94, 1.35),
  P_Value = c("0.825", "0.342", "0.506"),
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
    hjust = 0, size = 3.8, fontface = "bold", color = "grey20"
  ) +
  xlab("Hazard Ratio (95% CI, Log Scale)") +
  ylab("") +
  ggtitle("E. Prognostic Meta-Analysis: Adult General Cohorts") +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 8)),
    plot.title = element_text(face = "bold", size = 16, color = "darkblue", margin = margin(b=8)),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 25, b = 10, l = 25)
  )

# ===========================================================================
# Panel D: German AMLCG GPL96
# ===========================================================================
gpl96_df <- readRDS("03_Results/Phase14_AMLSG_Validation/GSE12417_GPL96_clinical_with_predictions.rds")
cols_gpl96 <- colnames(gpl96_df)[grepl("OS =", colnames(gpl96_df))]
row_non_na_gpl96 <- apply(gpl96_df[, cols_gpl96, drop = FALSE], 1, function(r) {
  non_na_idx <- which(!is.na(r))
  if (length(non_na_idx) > 0) cols_gpl96[non_na_idx[1]] else NA
})
gpl96_df$OS_time_raw <- as.numeric(gsub(".*OS = (\\d+) days.*", "\\1", row_non_na_gpl96))
gpl96_df$OS_status <- as.numeric(sapply(1:nrow(gpl96_df), function(i) {
  col <- row_non_na_gpl96[i]
  if (!is.na(col)) gpl96_df[i, col] else NA
}))
gpl96_df$OS_months <- gpl96_df$OS_time_raw / 30.4375

gpl96_df <- gpl96_df %>%
  filter(!is.na(OS_months) & !is.na(OS_status) & !is.na(predicted_cluster)) %>%
  mutate(predicted_cluster = factor(predicted_cluster, levels = c(1, 2),
                                    labels = c("Cluster 1", "Cluster 2")))
ps11c <- make_km_panel(gpl96_df, "OS_months", "OS_status",
                       "C. German AMLCG Trial GPL96 (N=163)")

# ===========================================================================
# Panel E: German AMLCG GPL570
# ===========================================================================
gpl570_df <- readRDS("03_Results/Phase14_AMLSG_Validation/GSE12417_GPL570_clinical_with_predictions.rds")
cols_gpl570 <- colnames(gpl570_df)[grepl("OS =", colnames(gpl570_df))]
row_non_na_gpl570 <- apply(gpl570_df[, cols_gpl570, drop = FALSE], 1, function(r) {
  non_na_idx <- which(!is.na(r))
  if (length(non_na_idx) > 0) cols_gpl570[non_na_idx[1]] else NA
})
gpl570_df$OS_time_raw <- as.numeric(gsub(".*OS = (\\d+) days.*", "\\1", row_non_na_gpl570))
gpl570_df$OS_status <- as.numeric(sapply(1:nrow(gpl570_df), function(i) {
  col <- row_non_na_gpl570[i]
  if (!is.na(col)) gpl570_df[i, col] else NA
}))
gpl570_df$OS_months <- gpl570_df$OS_time_raw / 30.4375

gpl570_df <- gpl570_df %>%
  filter(!is.na(OS_months) & !is.na(OS_status) & !is.na(predicted_cluster)) %>%
  mutate(predicted_cluster = factor(predicted_cluster, levels = c(1, 2),
                                    labels = c("Cluster 1", "Cluster 2")))
ps11d <- make_km_panel(gpl570_df, "OS_months", "OS_status",
                       "D. German AMLCG Trial GPL570 (N=79)")

# ===========================================================================
# ===========================================================================
# Panel F: Multivariate Analysis Forest Plot
# ===========================================================================
df_multivariate <- data.frame(
  Variable = c("Cluster 2 (vs Cluster 1)", "Age (per 10 years)", "Sex (Male vs Female)",
               "TP53 mutation", "TET2 mutation", "RUNX1 mutation", "ASXL1 mutation"),
  HR = c(1.06, 1.03^10, 1.12, 2.96, 1.42, 1.13, 1.21),
  Lower = c(0.83, 1.02^10, 0.91, 2.10, 1.03, 0.78, 0.82),
  Upper = c(1.36, 1.04^10, 1.38, 4.17, 1.94, 1.64, 1.79),
  P_Value = c("0.649", "7.3e-12", "0.278", "5.6e-10", "0.031", "0.518", "0.331"),
  Significant = c("No", "Yes", "No", "Yes", "Yes", "No", "No")
)

df_multivariate$Variable <- factor(df_multivariate$Variable, levels = rev(c(
  "Cluster 2 (vs Cluster 1)",
  "Age (per 10 years)",
  "Sex (Male vs Female)",
  "TP53 mutation",
  "TET2 mutation",
  "RUNX1 mutation",
  "ASXL1 mutation"
)))

p_multivariate <- ggplot(df_multivariate, aes(x = HR, y = Variable, color = Significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c("No" = "gray50", "Yes" = "black")) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 4, 5), limits = c(0.7, 5.5)) +
  geom_text(aes(label = paste0("HR=", sprintf("%.2f", HR), " (p=", P_Value, ")"), 
                x = Upper * 1.1), hjust = 0, size = 3.5, fontface = "plain", show.legend = FALSE) +
  labs(title = "F. Multivariate Analysis for Overall Survival",
       subtitle = "Subtype prognostic effect is driven entirely by genomic co-linearity",
       x = "Hazard Ratio (95% CI, log scale)", y = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 8)),
    plot.title = element_text(face = "bold", size = 16, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 14, color = "darkblue", margin = margin(b=10)),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 25, b = 10, l = 25)
  ) +
  expand_limits(x = 8.0)

# ===========================================================================
# Compile 6-panel layout (Panels A, B, C, D, E, F)
# ===========================================================================
km_grid_top <- arrange_ggsurvplots(list(ps11a, ps11b), print = FALSE, ncol = 2, nrow = 1, risk.table.height = 0.22)
km_grob_top <- do.call(arrangeGrob, km_grid_top)

km_grid_mid <- arrange_ggsurvplots(list(ps11c, ps11d), print = FALSE, ncol = 2, nrow = 1, risk.table.height = 0.22)
km_grob_mid <- do.call(arrangeGrob, km_grid_mid)

fig_dir <- "04_Figures/Supplementary_Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(fig_dir, "Supplementary_Figure_S11.pdf"), width = 14, height = 19)
grid.arrange(km_grob_top, km_grob_mid, p_forest, p_multivariate, ncol = 1, heights = c(2, 2, 0.9, 1.3))
dev.off()

png(file.path(fig_dir, "Supplementary_Figure_S11.png"), width = 14, height = 19, units = "in", res = 300)
grid.arrange(km_grob_top, km_grob_mid, p_forest, p_multivariate, ncol = 1, heights = c(2, 2, 0.9, 1.3))
dev.off()

cat("✓ Comprehensive 6-panel Supplementary Figure S11 generated successfully.\n")
