if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
}
cat("org.Hs.eg.db installed successfully\n")
