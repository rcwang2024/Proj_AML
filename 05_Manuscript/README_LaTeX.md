# LaTeX Manuscript Compilation Guide

**Date**: 2025-12-09
**File**: `AML_Molecular_Subtypes_Manuscript.tex`

---

## Files Included

### Main Files
- **AML_Molecular_Subtypes_Manuscript.tex** - Main LaTeX manuscript (complete)
- **references.bib** - BibTeX references (18 citations)
- **compile_manuscript.bat** - Windows compilation script
- **README_LaTeX.md** - This file

### Required Figures (referenced in manuscript)
The manuscript references figures from:
- `04_Figures/21_Main_Figures/` - Main figures 1-5
- `04_Figures/22_Drug_Validation/` - Drug validation figures
- `04_Figures/11_Survival_Reanalysis/` - Survival analysis figures

---

## Prerequisites

### LaTeX Distribution

**Windows:**
- **MiKTeX** (recommended): https://miktex.org/download
- **TeX Live**: https://www.tug.org/texlive/

**Mac:**
- **MacTeX**: https://www.tug.org/mactex/

**Linux:**
```bash
sudo apt-get install texlive-full  # Debian/Ubuntu
sudo yum install texlive-scheme-full  # RHEL/CentOS
```

### Required LaTeX Packages

The manuscript uses these packages (usually included in full distributions):
- `inputenc`, `graphicx`, `amsmath`, `amssymb`
- `booktabs`, `multirow`, `longtable`, `array`
- `float`, `caption`, `subcaption`, `geometry`
- `setspace`, `lineno`, `natbib`, `hyperref`, `xcolor`

If you get "missing package" errors, install them via:
```bash
# MiKTeX (Windows)
mpm --install=<package-name>

# TeX Live (Linux/Mac)
tlmgr install <package-name>
```

---

## Compilation Instructions

### Option 1: Windows Batch Script (Easiest)

1. Double-click `compile_manuscript.bat`
2. The script will:
   - Run pdflatex three times (for references)
   - Clean up auxiliary files
   - Open the PDF automatically
3. Output: `AML_Molecular_Subtypes_Manuscript.pdf`

### Option 2: Command Line (All Platforms)

```bash
# Navigate to manuscript folder
cd 05_Manuscript/

# Compile (run 3 times for references)
pdflatex AML_Molecular_Subtypes_Manuscript.tex
pdflatex AML_Molecular_Subtypes_Manuscript.tex
pdflatex AML_Molecular_Subtypes_Manuscript.tex

# Clean up (optional)
rm *.aux *.log *.out *.toc *.bbl *.blg
```

### Option 3: LaTeX Editor

**Recommended Editors:**
- **TeXstudio** (cross-platform, beginner-friendly): https://www.texstudio.org/
- **Overleaf** (online, no installation): https://www.overleaf.com/
- **TeXmaker** (cross-platform): https://www.xm1math.net/texmaker/
- **VS Code** with LaTeX Workshop extension

**Steps:**
1. Open `AML_Molecular_Subtypes_Manuscript.tex` in your editor
2. Click "Build" or "Compile" (usually F5 or F6)
3. View PDF output

### Option 4: Overleaf (Online, No Installation)

1. Go to https://www.overleaf.com/
2. Create a free account
3. Click "New Project" → "Upload Project"
4. Upload `AML_Molecular_Subtypes_Manuscript.tex` and `references.bib`
5. Click "Recompile" to generate PDF

---

## Manuscript Structure

### Document Contents

```
Title Page
  - Title, Authors, Affiliations
  - Correspondence

Abstract (250 words)
  - Background, Methods, Results, Conclusions
  - Keywords

Main Text (~8,000 words)
  1. Introduction
  2. Methods (abbreviated, refers to Supplementary)
  3. Results (5 subsections)
  4. Discussion

References (18 citations in BibTeX format)

Tables (4 main tables)
  - Table 1: Baseline Characteristics
  - Table 2: Multivariate Cox Regression
  - Table 3: Top 10 Differential Drugs
  - Table 4: Cluster Independence (R² improvement)

Figure Legends (5 main figures)
  - Figure 1: Molecular Subtyping Overview
  - Figure 2: Survival Meta-Analysis
  - Figure 3: Drug Response and BCL-2 Validation
  - Figure 4: Multivariate Analysis
  - Figure 5: Clinical Decision Tool

Supplementary Information (reference to 9 tables + 8 figures)
```

### Line Numbers and Spacing

- **Line numbers**: Enabled (for review)
- **Spacing**: Double-spaced (for review)

To disable for final version, comment out these lines in the .tex file:
```latex
% \linenumbers      % Remove line numbers
% \doublespacing    % Change to \singlespacing
```

---

## Including Figures

The manuscript currently has **placeholders** for figure references. To include actual figures:

### Method 1: Uncomment Figure Blocks

At the end of the .tex file, uncomment and modify:

```latex
\begin{figure}[H]
\centering
\includegraphics[width=0.9\textwidth]{04_Figures/21_Main_Figures/Figure1_Overview.pdf}
\caption{See Figure 1 legend above}
\label{fig:overview}
\end{figure}
```

### Method 2: Create Symbolic Links (if figures are in different directory)

**Windows:**
```cmd
mklink /D figures ..\..\04_Figures
```

**Linux/Mac:**
```bash
ln -s ../../04_Figures figures
```

Then change paths in .tex to:
```latex
\includegraphics{figures/21_Main_Figures/Figure1_Overview.pdf}
```

### Method 3: Copy Figures to Manuscript Folder

```bash
# Create local figures folder
mkdir figures

# Copy main figures
cp ../04_Figures/21_Main_Figures/*.pdf figures/
cp ../04_Figures/22_Drug_Validation/Figure5*.pdf figures/
```

---

## Customization

### Journal-Specific Formatting

The manuscript uses a generic article format. To convert to journal templates:

**Nature Medicine:**
- Download template: https://www.nature.com/documents/nature-template.zip
- Copy content into their template
- Follow their reference style (numbered)

**Blood:**
- Download template: https://ashpublications.org/blood/pages/authors
- Use their LaTeX class file
- Convert references to Blood style

**JCO:**
- Use standard article class (current format is close)
- May need to adjust figure/table formatting

### Changing Citation Style

Current style: `natbib` with numbered citations

To change to author-year:
```latex
\usepackage[authoryear]{natbib}
```

To use numbered superscripts (Nature style):
```latex
\usepackage[super,sort&compress]{natbib}
```

### Adjusting Margins

Current: 1 inch all sides

To change:
```latex
\usepackage[margin=0.75in]{geometry}  % Smaller margins
\usepackage[top=1in,bottom=1in,left=1.25in,right=1.25in]{geometry}  % Custom
```

---

## Troubleshooting

### Common Errors

**Error: "File not found"**
- Solution: Check file paths, ensure all files are in correct directories
- Check that figure paths match actual file locations

**Error: "Missing package"**
- Solution: Install missing package via package manager (MiKTeX/TeX Live)
- Or install full TeX distribution

**Error: "Undefined references"**
- Solution: Run pdflatex 2-3 times (references need multiple passes)

**Warning: "Overfull \hbox"**
- Solution: LaTeX is having trouble fitting text (usually OK for draft)
- To fix: Adjust wording or use `\sloppy` command

**Tables/Figures not appearing**
- Solution: Check `\label{}` and `\ref{}` commands match
- Ensure figure files exist at specified paths

### Getting Help

**LaTeX Documentation:**
- Comprehensive Guide: https://en.wikibooks.org/wiki/LaTeX
- Package documentation: https://www.ctan.org/

**Forums:**
- TeX Stack Exchange: https://tex.stackexchange.com/
- LaTeX Community: https://latex.org/forum/

**Quick Reference:**
- LaTeX Symbols: https://www.rpi.edu/dept/arc/training/latex/LaTeX_symbols.pdf
- Math symbols: https://oeis.org/wiki/List_of_LaTeX_mathematical_symbols

---

## Manuscript Statistics

- **Word count**: ~8,000 words (main text)
- **Page count**: ~25-30 pages (with figures, double-spaced)
- **Tables**: 4 main + 9 supplementary
- **Figures**: 5 main + 8 supplementary
- **References**: 18 citations
- **Compilation time**: ~10-20 seconds (3 passes)

---

## Converting to Other Formats

### LaTeX → Word (for journals requiring .docx)

**Method 1: Pandoc** (best quality)
```bash
pandoc AML_Molecular_Subtypes_Manuscript.tex -o manuscript.docx
```

**Method 2: latex2rtf** (good for complex math)
```bash
latex2rtf AML_Molecular_Subtypes_Manuscript.tex
```

**Method 3: PDF → Word** (last resort, poor quality)
- Adobe Acrobat: File → Export To → Microsoft Word
- Online converters: https://pdf2doc.com/

### LaTeX → HTML

```bash
htlatex AML_Molecular_Subtypes_Manuscript.tex
```

### LaTeX → Plain Text

```bash
pdftotext AML_Molecular_Subtypes_Manuscript.pdf manuscript.txt
```

---

## File Checklist for Submission

### Minimum Required Files
- [ ] `AML_Molecular_Subtypes_Manuscript.pdf` (compiled PDF)
- [ ] `AML_Molecular_Subtypes_Manuscript.tex` (source)
- [ ] `references.bib` (if using BibTeX)
- [ ] All figure files (PDFs in `04_Figures/`)

### Journal Submission Package
- [ ] Main manuscript PDF
- [ ] LaTeX source files (if accepted)
- [ ] Individual figure files (separate PDFs, 300-600 DPI)
- [ ] Individual table files (Excel or CSV)
- [ ] Supplementary materials PDF
- [ ] Cover letter
- [ ] Author statements

---

## Next Steps

### For Review/Editing
1. Compile the manuscript to PDF
2. Review output for formatting issues
3. Check that all tables/figures are referenced correctly
4. Proofread for typos and grammar
5. Share PDF with co-authors for feedback

### For Submission
1. Add author names and affiliations (currently placeholders)
2. Include actual figure files (currently referenced but not embedded)
3. Write acknowledgments and author contributions
4. Add funding information
5. Declare conflicts of interest
6. Convert to journal-specific format (if required)
7. Generate final high-resolution PDF (600 DPI for figures)

---

## Support

For issues with this manuscript:
- Check compilation errors in `.log` file
- Verify all packages are installed
- Ensure figure paths are correct
- Run pdflatex 3 times for proper reference resolution

**Created**: 2025-12-09
**Status**: Ready for compilation
**Format**: LaTeX with standard article class
**Output**: PDF (ready for submission or further editing)
