# Supplementary Figure S1: Wide-Format Batch Correction Audit (HIGH-FIDELITY VERSION)
# Objective: Demonstrate technical normalization using Wide Side-by-Side Gene Plots
library(tidyverse)
library(sva)
library(ggplot2)
library(cowplot)
library(readxl)
library(org.Hs.eg.db)

cat("=== GENERATING HIGH-FIDELITY SUPPLEMENTARY FIGURE S1 ===\n")

# 1. Load Datasets
cat("Loading datasets...\n")
beataml_expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
tcga_expr <- readRDS("01_Data/TCGA_LAML/tcga_laml_expression_normalized.rds")
beataml_clinical <- read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")

# 2. Extract Batch Info
batch_map <- data.frame(
  sample_id = beataml_clinical$dbgap_rnaseq_sample,
  batch = beataml_clinical$centerID,
  stringsAsFactors = FALSE
) %>% filter(!is.na(sample_id) & !is.na(batch))

metadata_beataml <- data.frame(sample_id = colnames(beataml_expr)) %>%
  left_join(batch_map, by = "sample_id") %>%
  mutate(batch = ifelse(is.na(batch), "BeatAML_Other", batch), cohort = "BeatAML")

metadata_tcga <- data.frame(sample_id = colnames(tcga_expr), batch = "TCGA_LAML", cohort = "TCGA")
metadata_final <- rbind(metadata_beataml, metadata_tcga)

# Align
common_genes <- intersect(rownames(beataml_expr), rownames(tcga_expr))
common_samples <- intersect(metadata_final$sample_id, c(colnames(beataml_expr), colnames(tcga_expr)))
expr_joint <- cbind(beataml_expr[common_genes, intersect(colnames(beataml_expr), common_samples)], 
                    tcga_expr[common_genes, intersect(colnames(tcga_expr), common_samples)])
metadata_final <- metadata_final %>% filter(sample_id %in% colnames(expr_joint)) %>% arrange(match(sample_id, colnames(expr_joint)))

# 3. Identify top 10 batch genes
means_beataml <- rowMeans(beataml_expr[common_genes,])
means_tcga <- rowMeans(tcga_expr[common_genes,])
diffs <- abs(means_beataml - means_tcga)
top10_ensg <- names(sort(diffs, decreasing=TRUE))[1:10]

symbols <- select(org.Hs.eg.db, keys=top10_ensg, columns=c("SYMBOL"), keytype="ENSEMBL")
top10_names <- symbols$SYMBOL[match(top10_ensg, symbols$ENSEMBL)]
top10_names[is.na(top10_names)] <- top10_ensg[is.na(top10_names)]

# 4. Global PCA Audit
var_genes <- apply(expr_joint, 1, var)
top_pca_genes <- names(sort(var_genes, decreasing = TRUE))[1:2000]

pca_before <- prcomp(t(expr_joint[top_pca_genes, ]), scale. = TRUE)
expr_corrected <- ComBat(dat=as.matrix(expr_joint), batch=metadata_final$batch, mod=NULL)
pca_after <- prcomp(t(expr_corrected[top_pca_genes, ]), scale. = TRUE)

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

p_pca_before <- ggplot(data.frame(PC1=pca_before$x[,1], PC2=pca_before$x[,2], cohort=metadata_final$cohort), 
                       aes(x=PC1, y=PC2, color=cohort, shape=cohort)) +
  geom_point(size=4, alpha=0.5) + scale_shape_manual(values=c("BeatAML"=16, "TCGA"=17)) +
  scale_color_manual(values=c("BeatAML"="#3498db", "TCGA"="#e74c3c")) +
  labs(title="A. Global PCA (Raw)", x="PC1", y="PC2") + theme_hf + theme(legend.position="none")

p_pca_after <- ggplot(data.frame(PC1=pca_after$x[,1], PC2=pca_after$x[,2], cohort=metadata_final$cohort), 
                      aes(x=PC1, y=PC2, color=cohort, shape=cohort)) +
  geom_point(size=4, alpha=0.5) + scale_shape_manual(values=c("BeatAML"=16, "TCGA"=17)) +
  scale_color_manual(values=c("BeatAML"="#3498db", "TCGA"="#e74c3c")) +
  labs(title="B. Global PCA (Corrected)", x="PC1", y="PC2") + theme_hf + theme(legend.position="bottom")

# Boxplots
create_wide_gene_plot <- function(expr_mat, genes, names, title) {
  plot_data <- data.frame(t(expr_mat[genes, ]))
  colnames(plot_data) <- names
  plot_data$cohort <- metadata_final$cohort
  long_df <- plot_data %>% pivot_longer(cols = -cohort, names_to = "Gene", values_to = "Expression") %>%
    mutate(Gene = factor(Gene, levels = names))
  
  ggplot(long_df, aes(x=Gene, y=Expression, fill=cohort)) +
    geom_boxplot(alpha=0.7, outlier.shape = NA) +
    scale_fill_manual(values=c("BeatAML"="#3498db", "TCGA"="#e74c3c")) +
    labs(title=title, x="", y="Log-Expression") + theme_hf + 
    theme(axis.text.x = element_text(face="bold", size=12, angle=45, hjust=1))
}

p_wide_before <- create_wide_gene_plot(expr_joint, top10_ensg, top10_names, "C. Top 10 Batch Genes (Raw)")
p_wide_after <- create_wide_gene_plot(expr_corrected, top10_ensg, top10_names, "D. Top 10 Batch Genes (Corrected)")

# Assembly
top_row <- plot_grid(p_pca_before, p_pca_after, rel_widths = c(1, 1))
bottom_row <- plot_grid(p_wide_before, p_wide_after, ncol = 1)
final_fig <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 2))

ggsave("05_Submission/Submission_Hub/03_Supplementary_Figures/FigureS1.pdf", final_fig, width=10, height=12.5)
ggsave("05_Submission/Submission_Hub/03_Supplementary_Figures/FigureS1.png", final_fig, width=10, height=12.5, dpi=300)

cat("✓ High-Fidelity Figure S1 Successfully Generated.\n")
