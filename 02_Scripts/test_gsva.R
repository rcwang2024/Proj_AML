#!/usr/bin/env Rscript
library(GSVA)
library(GSEABase)

cat("GSVA package version:", as.character(packageVersion("GSVA")), "\n\n")

# Try to see what methods are available
cat("Available gsva methods:\n")
print(methods("gsva"))
cat("\n")

# Check if gsvaParam exists
if (exists("gsvaParam")) {
  cat("gsvaParam function exists - use new API\n")
} else {
  cat("gsvaParam does not exist - use old API\n")
}

# Load small test data
expr <- matrix(rnorm(100), nrow=10)
rownames(expr) <- paste0("gene", 1:10)
colnames(expr) <- paste0("sample", 1:10)

geneset <- list(set1 = paste0("gene", 1:5), set2 = paste0("gene", 6:10))

cat("\nTesting GSVA call...\n")

# Try new API
tryCatch({
  param <- gsvaParam(expr, geneset, kcdf="Gaussian")
  result <- gsva(param)
  cat("SUCCESS with new API (gsvaParam)\n")
  print(dim(result))
}, error = function(e) {
  cat("New API failed:", conditionMessage(e), "\n")

  # Try old API
  tryCatch({
    result <- gsva(expr, geneset, kcdf="Gaussian", verbose=FALSE)
    cat("SUCCESS with old API (direct call)\n")
    print(dim(result))
  }, error = function(e2) {
    cat("Old API also failed:", conditionMessage(e2), "\n")
  })
})
