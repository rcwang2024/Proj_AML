# Install required packages for Phase 1 analysis

cat("Installing required R packages...\n\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Required CRAN packages
cran_packages <- c(
  "tidyverse",
  "data.table",
  "readxl",
  "pheatmap",
  "survival",
  "survminer",
  "ggplot2",
  "ggrepel",
  "RColorBrewer"
)

# Required Bioconductor packages
bioc_packages <- c(
  "sva",
  "DESeq2",
  "limma"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("  Installing ", pkg, "...\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(paste0("  ✓ ", pkg, " already installed\n"))
  }
}

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("  Installing ", pkg, "...\n"))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(paste0("  ✓ ", pkg, " already installed\n"))
  }
}

cat("\n✓ All packages installed successfully!\n")
