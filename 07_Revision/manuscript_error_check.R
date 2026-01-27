# Automated Manuscript Error Checking Script
# For AML_Molecular_Subtypes_Manuscript.tex
# Date: 2025-12-10

# Read manuscript
manuscript_file <- "D:/Projects/Project_AML/05_Manuscript/AML_Molecular_Subtypes_Manuscript.tex"
manuscript <- readLines(manuscript_file)

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("           MANUSCRIPT ERROR CHECK REPORT\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

# Initialize error counter
errors_found <- 0
warnings_found <- 0

# ============================================================================
# CHECK 1: Common Grammar Errors
# ============================================================================

cat("CHECK 1: Common Grammar Errors\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Check for "it's" (should be "its" in scientific writing)
its_lines <- grep("it's", manuscript, ignore.case = FALSE)
if (length(its_lines) > 0) {
  cat("⚠️  WARNING: Found \"it's\" (contraction) - should be \"its\" (possessive)\n")
  cat("   Lines:", paste(its_lines, collapse=", "), "\n")
  warnings_found <- warnings_found + length(its_lines)
} else {
  cat("✓ No \"it's\" contractions found\n")
}

# Check for "effect" vs "affect" misuse
affect_lines <- grep("affect", manuscript, ignore.case = TRUE)
if (length(affect_lines) > 0) {
  cat("⚠️  REVIEW: Found \"affect\" - verify correct usage (affect=verb, effect=noun)\n")
  cat("   Lines:", paste(affect_lines, collapse=", "), "\n")
  cat("   (Manual review needed)\n")
}

# Check for double spaces
double_space_lines <- grep("  ", manuscript, perl = TRUE)
if (length(double_space_lines) > 10) {
  cat("⚠️  WARNING: Found", length(double_space_lines), "lines with double spaces\n")
  warnings_found <- warnings_found + 1
} else {
  cat("✓ No major double space issues\n")
}

cat("\n")

# ============================================================================
# CHECK 2: Number Consistency
# ============================================================================

cat("CHECK 2: Critical Number Consistency\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Define expected numbers
expected_numbers <- list(
  "2,535" = "Total patients (671+151+1713)",
  "671" = "BeatAML patients",
  "520" = "Patients with drug data",
  "459" = "Multivariate analysis patients",
  "155" = "Drugs tested",
  "72" = "Drugs differential (FDR<0.05)",
  "1.25" = "Cohen's d for Venetoclax"
)

number_check <- data.frame(
  number = names(expected_numbers),
  description = unlist(expected_numbers),
  count = sapply(names(expected_numbers), function(x) {
    length(grep(gsub(",", "", x), manuscript, fixed = TRUE))
  }),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(number_check)) {
  if (number_check$count[i] == 0) {
    cat("❌ ERROR:", number_check$number[i], "-", number_check$description[i], "NOT FOUND\n")
    errors_found <- errors_found + 1
  } else if (number_check$count[i] < 2) {
    cat("⚠️  WARNING:", number_check$number[i], "- mentioned only", number_check$count[i], "time\n")
    warnings_found <- warnings_found + 1
  } else {
    cat("✓", number_check$number[i], "- mentioned", number_check$count[i], "times\n")
  }
}

cat("\n")

# ============================================================================
# CHECK 3: P-value Formatting
# ============================================================================

cat("CHECK 3: P-value Formatting Consistency\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Find all p-value patterns
pval_patterns <- c(
  "p=" = length(grep("p=", manuscript)),
  "p =" = length(grep("p =", manuscript)),
  "p<" = length(grep("p<", manuscript)),
  "p <" = length(grep("p <", manuscript)),
  "pval{" = length(grep("\\\\pval\\{", manuscript)),
  "P=" = length(grep("P=", manuscript))
)

cat("P-value format counts:\n")
for (pattern in names(pval_patterns)) {
  cat(" ", pattern, ":", pval_patterns[pattern], "\n")
}

if (length(unique(pval_patterns[pval_patterns > 0])) > 2) {
  cat("⚠️  WARNING: Multiple p-value formats used - consider standardizing\n")
  warnings_found <- warnings_found + 1
} else {
  cat("✓ P-value formatting relatively consistent\n")
}

cat("\n")

# ============================================================================
# CHECK 4: Citation Formatting
# ============================================================================

cat("CHECK 4: Citation Formatting\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Count citation styles
cite_latex <- length(grep("\\\\cite\\{", manuscript))
cite_bracket <- length(grep("\\[\\d+\\]", manuscript))

cat("Citation formats:\n")
cat("  LaTeX \\cite{}: ", cite_latex, "\n")
cat("  Bracket [1]: ", cite_bracket, "\n")

if (cite_latex > 0 && cite_bracket > 0) {
  cat("⚠️  WARNING: Mixed citation formats detected\n")
  warnings_found <- warnings_found + 1
} else {
  cat("✓ Consistent citation format\n")
}

cat("\n")

# ============================================================================
# CHECK 5: Figure and Table References
# ============================================================================

cat("CHECK 5: Figure and Table References\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Find figure references
fig_refs <- grep("Figure\\s+\\d|figure\\s+\\d|Fig\\.\\s+\\d|\\\\ref\\{fig:", manuscript, ignore.case = FALSE)
table_refs <- grep("Table\\s+\\d|table\\s+\\d|\\\\ref\\{tab:", manuscript, ignore.case = FALSE)

cat("Figure references found:", length(fig_refs), "\n")
cat("Table references found:", length(table_refs), "\n")

# Check for missing \\ref{}
fig_noref <- grep("Figure\\s+[1-9]", manuscript)
fig_withref <- grep("\\\\ref\\{fig:", manuscript)

if (length(fig_noref) > 0) {
  cat("⚠️  WARNING: Some figures referenced without \\ref{} (", length(fig_noref), "instances)\n")
  warnings_found <- warnings_found + 1
} else {
  cat("✓ All figures use \\ref{} format\n")
}

cat("\n")

# ============================================================================
# CHECK 6: Common LaTeX Errors
# ============================================================================

cat("CHECK 6: Common LaTeX Errors\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Check for unmatched braces (simple check)
open_braces <- sum(sapply(manuscript, function(x) length(gregexpr("\\{", x)[[1]])))
close_braces <- sum(sapply(manuscript, function(x) length(gregexpr("\\}", x)[[1]])))

if (open_braces != close_braces) {
  cat("❌ ERROR: Unmatched braces - { count:", open_braces, ", } count:", close_braces, "\n")
  errors_found <- errors_found + 1
} else {
  cat("✓ Brace matching OK (", open_braces, "pairs)\n")
}

# Check for common LaTeX typos
latex_issues <- list(
  "\\textit{" = "Italic text (check all closed)",
  "\\textbf{" = "Bold text (check all closed)",
  "$" = "Math mode (must be even count)"
)

for (cmd in names(latex_issues)) {
  count <- sum(sapply(manuscript, function(x) length(gregexpr(cmd, x, fixed = TRUE)[[1]])))
  if (cmd == "$" && count %% 2 != 0) {
    cat("❌ ERROR: Unmatched $ signs (", count, "found, must be even)\n")
    errors_found <- errors_found + 1
  } else {
    cat("✓", latex_issues[cmd], "-", count, "instances\n")
  }
}

cat("\n")

# ============================================================================
# CHECK 7: Supplementary Material References
# ============================================================================

cat("CHECK 7: Supplementary Material References\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Check for Table S1-S9
for (i in 1:9) {
  ref <- paste0("Table S", i)
  count <- length(grep(ref, manuscript, fixed = TRUE))
  if (count == 0) {
    cat("⚠️  WARNING:", ref, "not referenced in manuscript\n")
    warnings_found <- warnings_found + 1
  } else {
    cat("✓", ref, "- referenced", count, "time(s)\n")
  }
}

# Check for Figure S1-S8
for (i in 1:8) {
  ref <- paste0("Figure S", i)
  count <- length(grep(ref, manuscript, fixed = TRUE))
  if (count == 0) {
    cat("⚠️  WARNING:", ref, "not referenced in manuscript\n")
    warnings_found <- warnings_found + 1
  } else {
    cat("✓", ref, "- referenced", count, "time(s)\n")
  }
}

cat("\n")

# ============================================================================
# CHECK 8: Word Usage
# ============================================================================

cat("CHECK 8: Word Usage and Style\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Check for passive voice indicators (just flagging for review)
passive_indicators <- c("was measured", "were measured", "was performed", "were performed")
passive_count <- sum(sapply(passive_indicators, function(x) length(grep(x, manuscript, ignore.case = TRUE))))

cat("Passive voice instances:", passive_count, "\n")
if (passive_count > 20) {
  cat("⚠️  WARNING: High passive voice usage - consider active voice for clarity\n")
  warnings_found <- warnings_found + 1
} else {
  cat("✓ Acceptable passive voice usage\n")
}

# Check for abbreviations defined
abbrevs <- c("AML", "OS", "HR", "CI", "FDR", "AUC", "BCL-2", "PH", "RMST")
for (abbrev in abbrevs) {
  first_mention <- grep(abbrev, manuscript, fixed = TRUE)[1]
  if (!is.na(first_mention) && first_mention > 100) {
    # Check if it's defined in abstract/intro
    cat("⚠️  Review:", abbrev, "first mentioned at line", first_mention, "- check if defined\n")
  }
}

cat("\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("═══════════════════════════════════════════════════════════════════════\n")
cat("                        SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════\n\n")

cat("Total lines checked:", length(manuscript), "\n")
cat("❌ ERRORS found:", errors_found, "\n")
cat("⚠️  WARNINGS found:", warnings_found, "\n\n")

if (errors_found == 0 && warnings_found < 10) {
  cat("✓✓✓ MANUSCRIPT QUALITY: EXCELLENT\n")
  cat("Ready for submission with minor review of warnings\n\n")
} else if (errors_found == 0 && warnings_found < 20) {
  cat("✓✓ MANUSCRIPT QUALITY: GOOD\n")
  cat("Address warnings before submission\n\n")
} else if (errors_found > 0) {
  cat("⚠️  MANUSCRIPT QUALITY: NEEDS REVISION\n")
  cat("Fix errors before submission\n\n")
}

cat("═══════════════════════════════════════════════════════════════════════\n")

# Save report
sink("07_Revision/MANUSCRIPT_ERROR_CHECK_REPORT.txt")
cat("Manuscript Error Check Report\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("File:", manuscript_file, "\n\n")
cat("Errors found:", errors_found, "\n")
cat("Warnings found:", warnings_found, "\n\n")
sink()

cat("\nReport saved to: 07_Revision/MANUSCRIPT_ERROR_CHECK_REPORT.txt\n")
