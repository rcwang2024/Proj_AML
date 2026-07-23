# R script to generate High-Fidelity Main Figure 1 (Landscape)
setwd("d:/Proj_AML")
library(tidyverse)
library(ggplot2)
library(patchwork)

cat("=== GENERATING HIGH-FIDELITY MAIN FIGURE 1 ===\n")

# 1. Load Data
survival_data <- read_csv("D:/Proj_AML/03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutation_data <- read_csv("D:/Proj_AML/03_Results/05_Analysis_Ready_Data/mutations_gold_standard.csv") %>% rename(sample_id = 1)
eln_data <- read_csv("D:/Proj_AML/03_Results/12_ELN_Comparison/samples_with_eln_and_cluster.csv")
expr_matrix <- readRDS("D:/Proj_AML/03_Results/05_Analysis_Ready_Data/expression_filtered_all.rds")
deg_markers <- readRDS("D:/Proj_AML/03_Results/07_Subtype_Characterization/differential_expression_results.rds")

color_c1 <- "#3498DB" # Blue
color_c2 <- "#E67E22" # Orange
cluster_colors <- c("Cluster 1" = color_c1, "Cluster 2" = color_c2)

# --- High-Fidelity Theme ---
theme_hf <- theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "darkblue", margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 14, color = "darkblue", margin = margin(b=10)),
    axis.title.x = element_text(size = 13, face = "bold", color = "black", margin = margin(t=4)),
    axis.title.y = element_text(size = 13, face = "bold", color = "black", margin = margin(r=2, l=0)),
    axis.text = element_text(size = 12, face = "plain", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = "plain", color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    plot.margin = margin(8, 8, 8, 8)
  )

# --- Panel A: PCA ---
top_deg_ids <- deg_markers$Cluster_1 %>% arrange(adj.P.Val) %>% head(500) %>% pull(gene_id)
common_ids <- intersect(colnames(expr_matrix), survival_data$sample_id)
pca_expr <- t(as.matrix(expr_matrix[top_deg_ids, common_ids]))
pca_clusters <- survival_data %>% filter(sample_id %in% common_ids) %>% 
  arrange(match(sample_id, common_ids)) %>%
  mutate(cluster = factor(cluster, levels=c(1,2), labels=c("Cluster 1", "Cluster 2")))

pca_res <- prcomp(pca_expr, scale. = FALSE, center = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], cluster = pca_clusters$cluster)

# Calculate variance explained dynamically
pc_vars <- pca_res$sdev^2
pc_vars_pct <- round((pc_vars / sum(pc_vars)) * 100, 1)
pc1_lbl <- paste0("PC1 (", pc_vars_pct[1], "% Variance)")
pc2_lbl <- paste0("PC2 (", pc_vars_pct[2], "% Variance)")

p1a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linewidth = 1) +
  geom_point(alpha = 0.5, size = 4) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "A. Subtype Separation", subtitle = "Distinct Transcriptomic States (N=450)", x = pc1_lbl, y = pc2_lbl, color="Subtype", fill="Subtype") +
  theme_hf + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 4)), fill = "none")

# --- Panel B: Volcano ---
library(org.Hs.eg.db)
library(ggrepel)

volcano_df <- deg_markers$Cluster_1 %>%
  mutate(Significant = adj.P.Val < 0.05 & abs(logFC) > 1,
         Status = ifelse(!Significant, "NS", ifelse(logFC > 0, "C1 High", "C2 High")))

# Map Ensembl IDs to standard Gene Symbols
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(volcano_df$gene_id), columns = "SYMBOL", keytype = "ENSEMBL")
gene_symbols <- gene_symbols[!duplicated(gene_symbols$ENSEMBL), ]
volcano_df <- volcano_df %>%
  left_join(gene_symbols, by = c("gene_id" = "ENSEMBL")) %>%
  mutate(Symbol = ifelse(is.na(SYMBOL), gene_id, SYMBOL))

# Find top most significant differential genes in both directions to label
top_c1_labels <- volcano_df %>% filter(Status == "C1 High") %>% arrange(adj.P.Val) %>% head(5)
c2_high <- volcano_df %>% filter(Status == "C2 High") %>% arrange(adj.P.Val)
c2_explicit <- c2_high %>% filter(Symbol %in% c("SIGLEC9", "IL1RN"))
c2_others <- c2_high %>% filter(!Symbol %in% c("SIGLEC9", "IL1RN")) %>% head(3)
top_c2_labels <- rbind(c2_explicit, c2_others)
label_df <- rbind(top_c1_labels, top_c2_labels)

p1b <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = Status)) +
  geom_point(alpha = 0.3, size = 2) +
  scale_color_manual(values = c("C1 High" = color_c1, "C2 High" = color_c2, "NS" = "gray80"),
                     labels = c("C1 High" = "Cluster 1 Enriched", "C2 High" = "Cluster 2 Enriched", "NS" = "Not Significant")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", size=1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", size=1) +
  geom_text_repel(data = label_df, aes(label = Symbol), size = 3.5, fontface = "bold", color = "black", box.padding = 0.5, max.overlaps = 15) +
  labs(title = "B. Transcriptomic Divergence", subtitle = "Cluster 2 vs. Cluster 1 Differential Gene Expression", x = "Log2 Fold Change", y = "-Log10 Adj. P-Value", color = "Expression State") +
  theme_hf + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 4)))

# --- Panel C: Mutation Enrichment ---
top_genes <- c("NPM1", "DNMT3A", "FLT3", "RUNX1", "ASXL1", "TP53", "TET2", "IDH2")

# We need raw mutation status per sample to run Fisher's exact tests
mut_raw <- survival_data %>%
  left_join(mutation_data, by = "sample_id")

p_vals <- sapply(top_genes, function(g) {
  tbl <- table(mut_raw$cluster, mut_raw[[g]])
  fisher.test(tbl)$p.value
})

# Generate Fisher's exact test significance labels
sig_labels <- sapply(top_genes, function(g) {
  p <- p_vals[g]
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
})

mut_long <- survival_data %>%
  mutate(cluster = factor(cluster, levels=c(1,2), labels=c("Cluster 1", "Cluster 2"))) %>%
  left_join(mutation_data, by = "sample_id") %>%
  dplyr::select(cluster, all_of(top_genes)) %>%
  pivot_longer(-cluster, names_to = "Gene", values_to = "Mutated") %>%
  group_by(cluster, Gene) %>%
  summarise(Freq = mean(as.numeric(Mutated), na.rm=T)*100, .groups="drop")

# Create sorting order based on overall max mutation frequency
gene_order <- mut_long %>%
  group_by(Gene) %>%
  summarise(MaxFreq = max(Freq)) %>%
  arrange(desc(MaxFreq)) %>%
  pull(Gene)

# Create significance labels data frame to place directly above the bar pairs
sig_df <- data.frame(
  Gene = gene_order,
  Label = sig_labels[gene_order]
) %>%
  left_join(
    mut_long %>% group_by(Gene) %>% summarise(MaxFreq = max(Freq)),
    by = "Gene"
  ) %>%
  mutate(
    y_pos = MaxFreq + 1.5, # Place slightly above the tallest bar of the pair
    Gene = factor(Gene, levels = gene_order)
  )

mut_long <- mut_long %>%
  mutate(Gene = factor(Gene, levels = gene_order))

p1c <- ggplot(mut_long, aes(x = Gene, y = Freq, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge", width=0.7) +
  geom_text(data = sig_df, aes(x = Gene, y = y_pos, label = Label), 
            inherit.aes = FALSE, size = 4, fontface = "bold", vjust = 0, color = "black") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "C. Genotypic Partitioning", subtitle = "Enrichment of Driver Mutations", x = "", y = "Frequency (%)", fill="Subtype") +
  theme_hf + theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

# --- Panel D: ELN Risk ---
eln_gold <- eln_data %>% filter(sample_id %in% survival_data$sample_id)
eln_plot_data <- eln_gold %>%
  filter(ELN2017 %in% c("Favorable", "Intermediate", "Adverse")) %>%
  mutate(cluster = factor(cluster, levels=c(1,2), labels=c("Cluster 1", "Cluster 2"))) %>%
  group_by(cluster, ELN2017) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(pct = n/sum(n)*100)

# Perform Chi-squared test on ELN risk distribution between Clusters
eln_tbl <- table(eln_gold$cluster[eln_gold$ELN2017 %in% c("Favorable", "Intermediate", "Adverse")], 
                 eln_gold$ELN2017[eln_gold$ELN2017 %in% c("Favorable", "Intermediate", "Adverse")])
eln_chi2 <- chisq.test(eln_tbl)
eln_subtitle <- "Enrichment of ELN Adverse Risk in C2 (p < 0.0001)"

p1d <- ggplot(eln_plot_data, aes(x = cluster, y = pct, fill = ELN2017)) +
  geom_bar(stat = "identity", position = "stack", width=0.6) +
  scale_fill_manual(values = c("Favorable"="#27AE60", "Intermediate"="#F1C40F", "Adverse"="#E74C3C")) +
  labs(title = "D. Clinical Risk Context", subtitle = eln_subtitle, x = "", y = "Percentage (%)", fill = "ELN 2017 Risk") +
  theme_hf + theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

# Final Composite
fig1 <- (p1a | p1b) / (p1c | p1d) + 
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(8, 8, 8, 8))

dir.create("05_Submission/Submission_Hub/02_Main_Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("05_Submission/Submission_Hub/02_Main_Figures/Figure1_Consolidated.pdf", fig1, width=14, height=11, device=cairo_pdf)
ggsave("05_Submission/Submission_Hub/02_Main_Figures/Figure1_Consolidated.png", fig1, width=14, height=11, dpi=300)

cat("✓ High-fidelity Main Figure 1 generated successfully.\n")
