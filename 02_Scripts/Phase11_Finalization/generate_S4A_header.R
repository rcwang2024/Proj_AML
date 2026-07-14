# R script to generate IDENTICAL header for Figure S4-A
setwd("d:/Proj_AML")
library(ggplot2)

# Same Theme as other plots
theme_hf <- theme_void() + 
  theme(
    plot.title = element_text(face = "bold", size = 33, color = "darkblue", hjust = 0, margin = margin(b=8)),
    plot.subtitle = element_text(face = "plain", size = 29, color = "darkblue", hjust = 0, margin = margin(b=0)),
    plot.margin = margin(15, 15, 15, 15)
  )

p_header <- ggplot() + 
  labs(title = "A. Immune Landscape Deconvolution (xCell/MCP-counter)", 
       subtitle = "Cluster 2 is significantly enriched for myelomonocytic and cytotoxic signatures") +
  theme_hf

# Generate at exact target size (1900pt x 120pt = 26.39in x 1.67in)
ggsave("04_Figures/15_Immune_Deconvolution/FigureS4_HeaderA.pdf", p_header, width=26.39, height=1.67, device=cairo_pdf)
