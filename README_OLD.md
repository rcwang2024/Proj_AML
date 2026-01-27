# AML Multi-Omics Integration Project

**Project Title:** Comprehensive Multi-Omics Analysis of Acute Myeloid Leukemia
**Data Source:** Beat AML Database (Public Genomics Data)
**Project Start Date:** September 2025
**Last Updated:** October 2, 2025
**Status:** Phase 2 - Data Structure Analysis (In Progress)

---

## 📋 Table of Contents

1. [Project Overview](#project-overview)
2. [Project Objectives](#project-objectives)
3. [Data Description](#data-description)
4. [Project Structure](#project-structure)
5. [Progress Summary](#progress-summary)
6. [Key Findings](#key-findings)
7. [Analysis Pipeline](#analysis-pipeline)
8. [Requirements](#requirements)
9. [Usage Guide](#usage-guide)
10. [Important Notes](#important-notes)

---

## 🎯 Project Overview

This project focuses on building a comprehensive multi-omics analysis pipeline for **Acute Myeloid Leukemia (AML)** using publicly available data from the **Beat AML database**. The goal is to integrate multiple data types (genomics, transcriptomics, pharmacogenomics, and clinical data) to understand disease mechanisms, predict drug responses, and identify prognostic biomarkers.

### Project Context

**Previous Attempt:** Initially attempted analysis including MPN (Myeloproliferative Neoplasm) data but decided to **restart and focus exclusively on AML samples** for a more focused and clinically relevant analysis.

**Current Approach:** True multi-omics integration with mutation data as a critical component for understanding disease mechanisms and drug response.

---

## 🎯 Project Objectives

### Primary Objectives
1. **Data Integration:** Integrate genomics, transcriptomics, pharmacogenomics, and clinical data
2. **Molecular Subtyping:** Identify AML molecular subtypes using integrated multi-omics
3. **Drug Response Prediction:** Predict drug sensitivity from molecular features
4. **Biomarker Discovery:** Identify prognostic and predictive biomarkers
5. **Survival Analysis:** Build prognostic models using multi-omics features

### Secondary Objectives
- Understand mutation-expression relationships
- Identify drug-mutation associations
- Refine risk stratification beyond ELN2017
- Generate clinically actionable insights

---

## 📊 Data Description

### Data Source
**Beat AML Database:** A comprehensive dataset from the Beat AML clinical trial, including molecular profiling and drug sensitivity testing of AML patient samples.

### Available Data Types

#### 1. **Genomics (Mutations)**
- **File:** `beataml_mutations.txt` (3.5 MB)
- **Content:** Somatic mutation calls from WES/targeted sequencing
- **Samples:** 871 unique DNA samples
- **Mutations:** 11,720 variant calls
- **Annotations:** SIFT, PolyPhen, ExAC frequencies, functional predictions
- **Quality:** 6.04% missing data (Excellent)

#### 2. **Transcriptomics (Gene Expression)**
- **File:** `beataml_expression.txt` (268.3 MB)
- **Content:** Gene expression data (FPKM/TPM normalized)
- **Samples:** 707 RNA samples
- **Genes:** 22,843 genes profiled
- **Quality:** 0.53% missing data (Excellent)

#### 3. **Pharmacogenomics (Drug Response)**
- **Main File:** `beataml_drug_auc.txt` (18.2 MB)
- **Raw File:** `beataml_raw_inhibitor.txt` (47.1 MB)
- **Content:** Drug sensitivity testing results
- **Samples:** 542 samples tested
- **Compounds:** 166 unique drugs/inhibitors
- **Measurements:** 63,395 sample-drug combinations
- **Metrics:** AUC, IC50, IC10, IC25, IC75, IC90
- **Quality:** 0.63% missing data (Excellent)

#### 4. **Clinical Data**
- **File:** `beataml_clinical.xlsx` (0.47 MB)
- **Content:** Patient demographics, disease characteristics, outcomes
- **Samples:** 698 unique RNA samples, 903 DNA samples
- **Variables:** 95 clinical features
- **Key Data:**
  - Demographics (age, sex, race/ethnicity)
  - Disease classification (FAB, ELN2017 risk)
  - Survival outcomes (vital status, overall survival)
  - Molecular markers (FLT3-ITD, NPM1, TP53, etc.)
  - Laboratory values (blood counts, chemistry)
  - Treatment information
- **Quality:** 18.6% missing (Expected for clinical data)

### Integration Statistics

#### **Gold Standard Multi-Omics Dataset**
**478 samples** with complete 4-way integration:
- ✅ Gene Expression (22,843 genes)
- ✅ Drug Response (166 compounds)
- ✅ Clinical Data (95 variables)
- ✅ Somatic Mutations (annotated variants)

**Integration Success Rate:** 91.9% RNA-DNA sample mapping

#### **Extended Dataset**
**520 samples** with Expression + Drug + Clinical (3-way integration)

---

## 📁 Project Structure

```
Project_AML/
│
├── 01_Data/                                    # Raw and processed data
│   ├── BeatAML_Downloaded_Data/               # Original downloaded files
│   │   ├── beataml_expression.txt             # Gene expression (268 MB)
│   │   ├── beataml_drug_auc.txt               # Drug response AUC (18 MB)
│   │   ├── beataml_clinical.xlsx              # Clinical data (0.5 MB)
│   │   ├── beataml_mutations.txt              # Mutation calls (3.5 MB)
│   │   ├── beataml_raw_inhibitor.txt          # Raw drug data (47 MB)
│   │   └── beataml_drug_families.xlsx         # Drug classifications
│   │
│   └── [Previous analysis folders - archived]
│
├── 02_Scripts/                                 # Analysis scripts
│   ├── 01_Data_Processing/                    # Data processing scripts
│   │   ├── 01_verify_files.py                 # File verification (Python)
│   │   ├── 01_verify_files.R                  # File verification (R)
│   │   ├── 01_verify_files.ipynb              # File verification (Jupyter)
│   │   ├── 02_inspect_data.py                 # Data inspection (Python)
│   │   ├── 02_inspect_data.R                  # Data inspection (R)
│   │   ├── 02_inspect_data.ipynb              # Data inspection (Jupyter)
│   │   ├── 03_sample_inventory.py             # Sample inventory (Python)
│   │   ├── 03_sample_inventory.R              # Sample inventory (R)
│   │   ├── 04_drug_response_analysis.py       # Drug analysis (Phase 2)
│   │   ├── 05_clinical_data_analysis.py       # Clinical analysis (Phase 2)
│   │   └── README.md                          # Scripts documentation
│   │
│   ├── 02_QC_Analysis/                        # Quality control scripts
│   ├── 03_Molecular_Subtyping/                # Subtyping analyses
│   ├── 04_Clinical_Analysis/                  # Clinical analyses
│   ├── 05_Drug_Response/                      # Drug response modeling
│   ├── 06_Integration/                        # Multi-omics integration
│   └── 07_Visualization/                      # Visualization scripts
│
├── 03_Results/                                 # Analysis results
│   ├── 01_Processed_Data/                     # Processed data files
│   │   ├── sample_inventory.csv               # Sample inventory
│   │   ├── sample_integration_table.csv       # Integration matrix (520 samples)
│   │   ├── core_multi_omics_samples.txt       # Core sample list (478 samples)
│   │   ├── drug_response_summary.csv          # Drug statistics
│   │   ├── samples_drug_counts.csv            # Drug testing per sample
│   │   ├── clinical_data_summary.csv          # Clinical variable summary
│   │   ├── demographics_table.csv             # Demographics summary
│   │   ├── SAMPLE_INVENTORY_REPORT.md         # Detailed inventory report
│   │   ├── UPDATED_INVENTORY_WITH_CLINICAL.md # Updated with clinical data
│   │   ├── CLINICAL_VARIABLES_SUMMARY.md      # Clinical variables catalog
│   │   ├── PHASE2_TASKS_1_2_SUMMARY.md        # Phase 2 progress report
│   │   └── figures/                           # Generated plots
│   │       ├── auc_distribution.png           # Drug AUC distribution
│   │       ├── auc_boxplot.png                # AUC boxplot
│   │       ├── age_distribution.png           # Age histogram
│   │       └── clinical_missing_data.png      # Missing data visualization
│   │
│   └── 02_QC_Reports/                         # Quality control reports
│       ├── data_inspection_summary.txt        # Detailed data inspection
│       ├── INSPECTION_SUMMARY.md              # QC summary
│       └── clinical_completeness.csv          # Clinical data completeness
│
├── 04_Figures/                                 # Publication-ready figures
│
├── 05_Tables/                                  # Publication-ready tables
│
├── 06_Documentation/                           # Project documentation
│   ├── Data_Analysis_Log.txt                  # File verification log
│   ├── PHASE1_COMPLETION_REPORT.md            # Phase 1 report
│   ├── PHASE1_FINAL_SUMMARY.md                # Phase 1 final summary
│   └── README_files.md                        # Documentation of all files
│
└── README.md                                   # This file

```

---

## 📈 Progress Summary

### ✅ Phase 1: Initial Data Assessment (COMPLETE)

#### Task 1.1: File Verification ✅
- **Status:** Complete
- **Scripts:** Python, R, Jupyter
- **Outcome:** All 6 files verified and intact with MD5 checksums
- **Files Generated:**
  - `01_verify_files.py`, `01_verify_files.R`, `01_verify_files.ipynb`
  - `Data_Analysis_Log.txt`

#### Task 1.2: Data Inspection ✅
- **Status:** Complete
- **Scripts:** Python, R, Jupyter
- **Outcome:** Comprehensive analysis of all datasets
- **Key Results:**
  - Expression: 22,843 genes × 707 samples (0.53% missing)
  - Drug response: 63,395 measurements (0.63% missing)
  - Clinical: 942 patients × 95 variables (18.6% missing)
  - Mutations: 11,720 calls (6.04% missing)
- **Files Generated:**
  - `02_inspect_data.py`, `02_inspect_data.R`, `02_inspect_data.ipynb`
  - `data_inspection_summary.txt`
  - `INSPECTION_SUMMARY.md`

#### Task 1.3: Sample Inventory ✅
- **Status:** Complete
- **Scripts:** Python, R
- **Outcome:** Sample overlap analysis and integration matrix
- **Key Results:**
  - 520 samples with Expression + Drug + Clinical
  - 478 samples with complete 4-way multi-omics (GOLD STANDARD)
  - 91.9% RNA-DNA mapping success
- **Files Generated:**
  - `03_sample_inventory.py`, `03_sample_inventory.R`
  - `sample_inventory.csv`
  - `sample_integration_table.csv`
  - `core_multi_omics_samples.txt`
  - `SAMPLE_INVENTORY_REPORT.md`

---

### 🔄 Phase 2: Data Structure Analysis (IN PROGRESS)

#### Task 2.1: Drug Response Data Analysis ✅
- **Status:** Complete
- **Script:** `04_drug_response_analysis.py`
- **Key Findings:**
  - 603 samples, 166 drugs, 63,395 measurements
  - Mean 100.7 drugs tested per sample
  - AUC: mean 196.1, median 207.5
  - High quality: 0.63% missing overall
- **Files Generated:**
  - `drug_response_summary.csv`
  - `samples_drug_counts.csv`
  - `auc_distribution.png`, `auc_boxplot.png`

#### Task 2.2: Clinical Data Analysis ✅
- **Status:** Complete
- **Script:** `05_clinical_data_analysis.py`
- **Key Findings:**
  - 942 patients, 95 clinical variables
  - Age: mean 57.2 years, median 61 years
  - 60% deceased, median OS: 361 days
  - FLT3-ITD: 23.3% positive, NPM1: 25.6% positive
  - ELN2017 risk: 24.8% adverse, 20.5% favorable
- **Files Generated:**
  - `clinical_data_summary.csv`
  - `demographics_table.csv`
  - `clinical_completeness.csv`
  - `age_distribution.png`, `clinical_missing_data.png`

#### Task 2.3: Expression Data Analysis ⏳
- **Status:** Pending
- **Planned:** Comprehensive gene expression QC and structure analysis

---

### 🔮 Phase 3: Detailed QC and Preprocessing (PLANNED)
- Expression normalization and batch effect assessment
- Drug response QC and filtering
- Mutation filtering and annotation
- Sample quality filtering

### 🔮 Phase 4: Multi-Omics Integration (PLANNED)
- Molecular subtyping
- Drug response prediction models
- Survival analysis
- Pathway analysis

---

## 🔍 Key Findings

### Data Quality Summary

| Data Type | Samples | Features | Missing % | Quality |
|-----------|---------|----------|-----------|---------|
| **Expression** | 707 | 22,843 genes | 0.53% | ⭐⭐⭐⭐⭐ Excellent |
| **Drug Response** | 603 | 166 compounds | 0.63% | ⭐⭐⭐⭐⭐ Excellent |
| **Clinical** | 698/903 | 95 variables | 18.6% | ⭐⭐⭐⭐☆ Good |
| **Mutations** | 871 | 11,720 calls | 6.04% | ⭐⭐⭐⭐⭐ Excellent |

### Sample Integration

```
Total RNA Samples: 707
├─ With Clinical Data: 671 (95%)
├─ With Drug Response: 542 (77%)
│   └─ With Mutations: 478 (68% of all, 88% of drug-tested)
└─ Expression Only: 165 samples

🌟 Complete 4-way Multi-Omics: 478 samples
   ✓ Gene Expression
   ✓ Drug Response
   ✓ Clinical Data
   ✓ Somatic Mutations
```

### Clinical Characteristics

**Demographics:**
- Mean age: 57.2 years (range: 0-88)
- Sex: 56.1% male, 43.9% female

**Disease:**
- De novo AML: 42.7%
- Secondary AML: 54.9%
- Relapse: 6.6%

**Survival:**
- Deaths: 60%
- Median OS: 361 days (~12 months)
- Cause: 75% disease-related

**Molecular Markers:**
- FLT3-ITD positive: 23.3%
- NPM1 mutated: 25.6%
- ELN2017 Adverse risk: 24.8%

### Drug Response Characteristics

**Testing Coverage:**
- Mean drugs per sample: 100.7
- Median: 114 drugs
- Most comprehensive drugs: Imatinib, Sorafenib, Dasatinib

**AUC Distribution:**
- Mean: 196.1
- Median: 207.5
- Range: [0, 286.3]
- Outliers: 0.98%

---

## 🔬 Analysis Pipeline

### Current Workflow

```
1. DATA ACQUISITION
   └─ Download from Beat AML database ✅

2. FILE VERIFICATION
   └─ Integrity checks, MD5 checksums ✅

3. DATA INSPECTION
   └─ Structure, dimensions, quality ✅

4. SAMPLE INVENTORY
   └─ Cross-dataset mapping ✅

5. DATA STRUCTURE ANALYSIS
   ├─ Drug response ✅
   ├─ Clinical data ✅
   └─ Expression data ⏳

6. QUALITY CONTROL (Next)
   ├─ Expression QC
   ├─ Drug response QC
   ├─ Mutation filtering
   └─ Sample filtering

7. MULTI-OMICS INTEGRATION
   ├─ Molecular subtyping
   ├─ Drug response prediction
   ├─ Survival analysis
   └─ Pathway analysis

8. RESULTS & VISUALIZATION
   └─ Publication-ready outputs
```

### Analysis Scripts

All scripts created in **multiple formats** for reproducibility:
- **Python** (.py) - Primary analysis language
- **R** (.R) - Statistical analysis
- **Jupyter Notebooks** (.ipynb) - Interactive exploration

---

## 💻 Requirements

### Software Requirements

#### Python Environment
```bash
# Core packages
pandas >= 1.3.0
numpy >= 1.20.0
matplotlib >= 3.4.0
seaborn >= 0.11.0
scipy >= 1.7.0
openpyxl >= 3.0.0  # For Excel file reading

# Bioinformatics packages (for future analyses)
scikit-learn >= 0.24.0
statsmodels >= 0.12.0
lifelines  # For survival analysis
```

#### R Environment
```R
# Core packages
readr
readxl
dplyr
tibble
ggplot2

# Bioinformatics packages (for future analyses)
DESeq2
limma
survival
```

### Hardware Requirements
- **RAM:** Minimum 16 GB (32 GB recommended for expression analysis)
- **Storage:** ~10 GB for data and results
- **Processor:** Multi-core CPU recommended

---

## 📖 Usage Guide

### Running Analysis Scripts

#### 1. File Verification
```bash
# Python
cd D:\Projects\Project_AML
python 02_Scripts/01_Data_Processing/01_verify_files.py

# R
Rscript 02_Scripts/01_Data_Processing/01_verify_files.R
```

#### 2. Data Inspection
```bash
# Python
python 02_Scripts/01_Data_Processing/02_inspect_data.py

# R
Rscript 02_Scripts/01_Data_Processing/02_inspect_data.R
```

#### 3. Sample Inventory
```bash
# Python
python 02_Scripts/01_Data_Processing/03_sample_inventory.py

# R (if needed after installing openpyxl)
Rscript 02_Scripts/01_Data_Processing/03_sample_inventory.R
```

#### 4. Drug Response Analysis
```bash
python 02_Scripts/01_Data_Processing/04_drug_response_analysis.py
```

#### 5. Clinical Data Analysis
```bash
python 02_Scripts/01_Data_Processing/05_clinical_data_analysis.py
```

### Accessing Results

All outputs are organized in `03_Results/`:
- **Processed Data:** `03_Results/01_Processed_Data/`
- **QC Reports:** `03_Results/02_QC_Reports/`
- **Figures:** `03_Results/01_Processed_Data/figures/`

### Key Output Files

**Sample Lists:**
- `sample_integration_table.csv` - Integration status for 520 samples
- `core_multi_omics_samples.txt` - List of 478 complete multi-omics samples

**Data Summaries:**
- `drug_response_summary.csv` - Drug-level statistics
- `clinical_data_summary.csv` - All 95 clinical variables
- `demographics_table.csv` - Patient demographics

**Reports:**
- `SAMPLE_INVENTORY_REPORT.md` - Comprehensive sample analysis
- `CLINICAL_VARIABLES_SUMMARY.md` - Clinical data catalog
- `PHASE2_TASKS_1_2_SUMMARY.md` - Current progress summary

---

## 📌 Important Notes

### Critical Project Decisions

1. **AML-Only Focus:** Excluded MPN samples for cleaner, more focused analysis
2. **Multi-Omics Integration:** Prioritized samples with all 4 data types
3. **Mutation Data:** Included as critical component (not excluded despite initial consideration)
4. **Gold Standard Cohort:** 478 samples with complete multi-omics

### Data Considerations

#### Sample ID Formats
- **RNA samples:** BA####R (e.g., BA2000R)
- **DNA samples:** BA####D (e.g., BA2000D)
- **Mapping:** 91.9% success rate (478/520 RNA samples have DNA data)

#### Missing Data Patterns
- **Expression:** Minimal (0.53%) - highest quality
- **Drug response:** Minimal (0.63%) - structured missingness (not all drugs on all samples)
- **Clinical:** Moderate (18.6%) - expected for clinical data
- **Mutations:** Low (6.04%) - well-annotated

#### Quality Flags
- All datasets passed integrity checks
- MD5 checksums calculated and stored
- No corrupted files detected

### Clinical Data Note
**Important:** openpyxl package required to read clinical Excel file. Installed successfully after initial setup.

### Future Considerations

1. **Batch Effects:** Will assess in Phase 3 (expression data shows multiple "waves")
2. **Survival Analysis:** Complete data available (vital status: 100%, OS: 100%)
3. **Drug Filtering:** May filter drugs with <50 samples tested
4. **Gene Filtering:** Will filter low-expression genes in Phase 3

---

## 🎓 Project Significance

### Scientific Impact
- **Comprehensive Dataset:** One of largest multi-omics AML datasets
- **Complete Integration:** 478 samples with 4 data types
- **Clinical Relevance:** Real patient data with outcomes
- **Drug Testing:** 166 compounds including targeted therapies

### Potential Applications
1. **Precision Medicine:** Predict drug response from molecular features
2. **Biomarker Discovery:** Identify prognostic/predictive markers
3. **Subtype Refinement:** Improve molecular classification
4. **Treatment Optimization:** Match patients to optimal therapies
5. **Resistance Mechanisms:** Understand drug resistance

---

## 📞 Contact & Collaboration

**Project Lead:** [Your Name]
**Institution:** [Your Institution]
**Email:** [Your Email]

For questions, collaboration opportunities, or access to analysis results, please contact via email.

---

## 📚 References

### Data Source
- **Beat AML Database:** [https://biodev.github.io/BeatAML2/](https://biodev.github.io/BeatAML2/)
- Tyner JW, et al. "Functional genomic landscape of acute myeloid leukaemia." *Nature* 2018.

### Risk Stratification
- Döhner H, et al. "Diagnosis and management of AML in adults: 2017 ELN recommendations." *Blood* 2017.

---

## 📄 License

This project uses publicly available data from the Beat AML database. All analysis code is available for research purposes. Please cite appropriately if using this work.

---

## 🔄 Version History

- **v1.0 (Oct 2, 2025):** Phase 1 complete, Phase 2 in progress
  - File verification complete
  - Data inspection complete
  - Sample inventory complete
  - Drug response analysis complete
  - Clinical data analysis complete
  - 478 complete multi-omics samples identified

---

## 📊 Current Status

```
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
         PROJECT STATUS: Phase 2 (In Progress)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Phase 1: Initial Data Assessment       ✅ 100% Complete
  ├─ Task 1.1: File Verification       ✅ Complete
  ├─ Task 1.2: Data Inspection         ✅ Complete
  └─ Task 1.3: Sample Inventory        ✅ Complete

Phase 2: Data Structure Analysis       🔄 67% Complete
  ├─ Task 2.1: Drug Response           ✅ Complete
  ├─ Task 2.2: Clinical Data           ✅ Complete
  └─ Task 2.3: Expression Data         ⏳ Pending

Phase 3: QC and Preprocessing          ⏸️ Not Started
Phase 4: Multi-Omics Integration       ⏸️ Not Started

Total Scripts Created: 11
Total Documentation Files: 15+
Total Output Files: 20+

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

---

**Last Updated:** October 2, 2025
**Next Milestone:** Complete Task 2.3 (Expression Data Analysis)
**Overall Progress:** Phase 2 of 4 (50% of project phases complete)

---

*This project represents a comprehensive effort to understand AML through integrated multi-omics analysis, with the ultimate goal of improving patient outcomes through precision medicine approaches.*
