# AML Multi-Omics Project - Complete Data Structure Documentation

**Last Updated**: October 25, 2025
**Project**: Molecular Subtyping in Adult AML with Drug Response Validation

---

## 📁 PROJECT DIRECTORY STRUCTURE

```
D:/Projects/Project_AML/
│
├── 01_Data/                          # Raw and processed input data
│   ├── BeatAML_Downloaded_Data/      # Original BeatAML files (6 files, ~338 MB)
│   ├── BeatAML_MultiOmics/           # Organized multi-omics data
│   ├── TCGA_LAML/                    # TCGA validation cohort
│   └── Gene_Sets/                    # Reference gene sets
│
├── 02_Scripts/                       # Analysis scripts (R and Python)
│   ├── 01_Data_Processing/
│   ├── 09_Phase2_Molecular_Subtyping/
│   ├── 10_Phase3_Survival_Analysis/
│   ├── Phase3_CriticalValidation/
│   ├── Phase4_ManuscriptPrep/
│   └── [Other phase-specific folders]
│
├── 03_Results/                       # All analysis outputs
│   ├── 01_Processed_Data/            # ⭐ Clean, analysis-ready data
│   ├── 06_Molecular_Subtypes/        # ⭐ Cluster assignments
│   ├── 11_Survival_Reanalysis/       # Phase 3 survival analyses
│   ├── 21_Manuscript_Prep/           # Phase 4 materials
│   ├── 23_Drug_Validation/           # ⭐⭐⭐ Phase 5 drug response
│   └── [25 total result directories]
│
├── 04_Figures/                       # Publication-quality figures
│   ├── 21_Main_Figures/              # Main manuscript figures (1-4)
│   ├── 22_Drug_Validation/           # ⭐ Phase 5 figures (Figure 5, S1-S3)
│   └── [20+ figure directories]
│
├── 05_Reports/                       # QC and analysis reports
├── 06_Documentation/                 # Project documentation
│
└── [Root-level summary documents]    # ⭐ Start here for project overview
    ├── COMPLETE_PROJECT_SUMMARY_ALL_PHASES_V3.md
    ├── EXECUTIVE_SUMMARY_WITH_PHASE5.md
    ├── MANUSCRIPT_ABSTRACT_WITH_PHASE5.md
    ├── KEY_STATISTICS_PHASE5_QUICK_REFERENCE.md
    └── CLAUDE.md                     # Project guide for Claude Code
```

---

## 📊 DATA FLOW ARCHITECTURE

```
RAW DATA (01_Data/)
    ↓
PHASE 1: DATA PROCESSING & INTEGRATION
    ├─ Expression: 22,843 genes × 671 samples → Log2, batch-corrected
    ├─ Mutations: 11,720 variants → Binary matrix
    ├─ Drug Response: 166 compounds × 520 samples → AUC values
    └─ Clinical: 95 variables × 671 samples → Survival, demographics
    ↓
PROCESSED DATA (03_Results/01_Processed_Data/)
    ├─ expression_matrix.rds (99 MB)
    ├─ drug_response_auc.rds (4.8 MB)
    ├─ clinical_data.rds (150 KB)
    └─ master_sample_id_mapping.csv (73 KB) ⭐
    ↓
PHASE 2: MOLECULAR SUBTYPING
    ├─ Consensus Clustering → k=2 optimal
    ├─ 50-gene Classifier → 92.9% accuracy
    └─ Cluster Assignments → 671 samples
    ↓
CLUSTER ASSIGNMENTS (03_Results/06_Molecular_Subtypes/)
    └─ sample_cluster_assignments.csv ⭐⭐⭐
    ↓
PHASE 3-4: VALIDATION & MANUSCRIPT PREP
    ├─ Survival Analysis
    ├─ Multivariate Testing
    ├─ TCGA/TARGET Validation
    └─ Figure Generation
    ↓
PHASE 5: DRUG RESPONSE VALIDATION ⭐⭐⭐
    ├─ 155 drugs tested → 72 differential (FDR<0.05)
    ├─ Cluster independence → 19/20 drugs (mean +42% R²)
    ├─ BCL-2 pathway → 9/10 genes validated
    └─ Clinical actionability → Venetoclax biomarker
    ↓
PUBLICATION MATERIALS
    ├─ 7 Main Figures (Figure 1-5)
    ├─ 9 Supplementary Tables (S1-S9)
    └─ Manuscript with Phase 5 integration
```

---

## 🗂️ KEY DATA FILES

### **⭐ MOST IMPORTANT FILES (Start Here)**

#### 1. Sample Mapping & Integration
**File**: `03_Results/01_Processed_Data/master_sample_id_mapping.csv`
- **Purpose**: Master index connecting all data types
- **Dimensions**: 1,067 rows × 11 columns
- **Key Columns**:
  - `unified_sample_id`: Primary key (e.g., BA2000R)
  - `expression_id`: RNA-seq sample ID
  - `drug_response_id`: Drug screening sample ID (DNA-based, e.g., BA2000D)
  - `clinical_id`: Clinical data sample ID
  - `mutation_id`: Mutation data sample ID
  - `has_expression`, `has_drug_response`, `has_clinical`, `has_mutations`: Boolean flags
  - `n_data_types`: Count of available data types (1-4)
  - `cohort_category`: Classification (e.g., "complete_quad_omics", "triple_expr_clin_mut")

**Sample Cohorts**:
- **Complete quad-omics**: 478 samples (all 4 data types)
- **Drug response cohort**: 520 samples (expression + drug + clinical)
- **Survival cohort**: 671 samples (expression + clinical)

**Sample ID Convention**:
- RNA samples: `BA####R` (e.g., BA2000R)
- DNA samples: `BA####D` (e.g., BA2000D)
- Mapping rule: Replace 'R' with 'D' for RNA→DNA conversion

---

#### 2. Cluster Assignments
**File**: `03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv`
- **Purpose**: Molecular subtype labels for all samples
- **Dimensions**: 671 rows × 3 columns
- **Columns**:
  - `sample_id`: RNA sample ID (BA####R format)
  - `cluster`: Cluster assignment (1 or 2)
  - `consensus_score`: Clustering stability score (0-1, higher = more confident)

**Cluster Characteristics**:
- **Cluster 1** (n~220, 42%): NPM1+/DNMT3A+, favorable-risk, Venetoclax-sensitive
- **Cluster 2** (n~300, 58%): TP53+/RUNX1+/ASXL1+, adverse-risk, Venetoclax-resistant

---

#### 3. Drug Response Data (Phase 5 - Primary Output)
**File**: `03_Results/23_Drug_Validation/all_drugs_differential_response.csv`
- **Purpose**: Differential drug sensitivity between clusters (155 drugs tested)
- **Dimensions**: 155 rows × 15 columns
- **Key Columns**:
  - `drug`: Drug name
  - `n_samples`: Total samples with data for this drug
  - `n_cluster1`, `n_cluster2`: Samples per cluster
  - `mean_auc_cluster1`, `mean_auc_cluster2`: Mean AUC values (lower = more sensitive)
  - `sd_auc_cluster1`, `sd_auc_cluster2`: Standard deviations
  - `median_auc_cluster1`, `median_auc_cluster2`: Median AUC values
  - `wilcoxon_pvalue`: Differential sensitivity test p-value
  - `cohens_d`: Effect size
  - `fold_difference`: Ratio of mean AUCs
  - `cluster_more_sensitive`: Which cluster responds better
  - `fdr`: FDR-corrected p-value

**Key Statistics**:
- 72/155 drugs: FDR<0.05 (significant differential response)
- Top drug: Venetoclax (p=2.78×10⁻²⁴, Cohen's d=1.25)

---

#### 4. Cluster Independence Analysis (Phase 5 - Critical Finding)
**File**: `03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv`
- **Purpose**: Tests if clusters predict drug response BEYOND mutations
- **Dimensions**: 20 rows × 11 columns
- **Key Columns**:
  - `drug`: Drug name
  - `n_samples`: Sample size for this analysis
  - `n_mutations`: Number of mutations tested (11: NPM1, FLT3, DNMT3A, IDH1, IDH2, TET2, TP53, RUNX1, ASXL1, NRAS, KRAS)
  - `mutations_tested`: List of mutations included
  - `r2_mutations_only`: R² from model with only mutations
  - `r2_cluster_only`: R² from model with only cluster
  - `r2_mutations_plus_cluster`: R² from model with mutations + cluster
  - `r2_improvement`: Absolute R² increase (mutations+cluster vs mutations only)
  - `pct_improvement`: Percent R² increase
  - `p_cluster_adds_value`: ANOVA p-value testing if cluster significantly improves model
  - `fdr_cluster`: FDR-corrected p-value

**BREAKTHROUGH FINDING**:
- **19/20 drugs**: FDR<0.05 (clusters add independent predictive value)
- **Mean R² improvement**: +0.049 (+42% relative increase)
- **Venetoclax**: +0.225 R² improvement (+161%)

**Interpretation**: Clusters capture therapeutic vulnerabilities beyond genomic alterations

---

### **RAW DATA FILES (01_Data/BeatAML_Downloaded_Data/)**

#### Original BeatAML Files (6 files, 338 MB total)

1. **beataml_expression.txt** (269 MB)
   - **Format**: Tab-delimited
   - **Dimensions**: 22,843 genes × 671 samples
   - **Gene IDs**: Ensembl IDs
   - **Values**: FPKM (Fragments Per Kilobase Million)
   - **Processing**: Log2-transformed, batch-corrected (ComBat)

2. **beataml_drug_auc.txt** (19 MB)
   - **Format**: Tab-delimited, long format
   - **Dimensions**: ~86,000 rows (520 samples × 166 drugs, sparse)
   - **Columns**:
     - `dbgap_dnaseq_sample`: DNA sample ID (BA####D)
     - `inhibitor`: Drug name
     - `auc`: Area Under Curve (drug sensitivity metric)
   - **AUC Interpretation**: **Lower AUC = MORE sensitive**
   - **Coverage**: 166 compounds, ex vivo screening assay
   - **Missingness**: ~6% (expected for drug screening)

3. **beataml_mutations.txt** (3.5 MB)
   - **Format**: Tab-delimited, variant-level
   - **Dimensions**: ~11,720 variant calls
   - **Columns**:
     - `dbgap_dnaseq_sample`: DNA sample ID
     - `Hugo_Symbol`: Gene symbol (e.g., NPM1, TP53)
     - `Variant_Classification`: Mutation type (e.g., Missense, Nonsense, Frame_Shift)
     - `HGVSp_Short`: Protein change notation
     - `Chromosome`, `Start_Position`, `End_Position`: Genomic coordinates
   - **Processing**: Converted to binary mutation matrix (gene × sample)

4. **beataml_clinical.xlsx** (477 KB)
   - **Format**: Excel spreadsheet
   - **Dimensions**: 671 patients × 95 variables
   - **Key Variables**:
     - `dbgap_rnaseq_sample`: RNA sample ID (primary key)
     - `overallSurvival`: Survival time in days
     - `vitalStatus`: 0=alive, 1=deceased
     - `ageAtDiagnosis`: Age in years
     - `sex`: M/F
     - `specificDxAtAcquisition`: AML subtype
     - `priorMalignancy`: Prior cancer history
     - `cytogeneticComplexity`: Normal, Complex, Other
     - `ELN2017`: ELN risk group (Favorable, Intermediate, Adverse)

5. **beataml_drug_families.xlsx** (58 KB)
   - **Format**: Excel spreadsheet
   - **Purpose**: Drug classification by mechanism of action
   - **Content**:
     - Drug names
     - Target pathways (BCL-2, MEK, HDAC, FLT3, etc.)
     - FDA approval status

6. **beataml_raw_inhibitor.txt** (48 MB)
   - **Format**: Raw drug screening data
   - **Purpose**: Full dose-response curves (not used in primary analysis)

---

### **PROCESSED DATA FILES (03_Results/01_Processed_Data/)**

#### Expression Data

1. **expression_matrix.rds** (99 MB)
   - **Format**: R data structure (matrix)
   - **Dimensions**: 22,843 genes × 671 samples
   - **Values**: Log2-transformed FPKM
   - **Row names**: Ensembl gene IDs
   - **Column names**: RNA sample IDs (BA####R)
   - **Processing**: Quality-filtered, no batch correction at this stage

2. **03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds**
   - **Same dimensions as above**
   - **Processing**: ComBat batch correction applied
   - **Use**: Preferred for clustering and differential expression

#### Drug Response Data

3. **drug_response_auc.rds** (4.8 MB)
   - **Format**: R data structure (data.frame or matrix)
   - **Dimensions**: 520 samples × 166 drugs (wide format)
   - **Values**: AUC (lower = more sensitive)
   - **Missing data**: ~6% (varies by drug)

4. **drug_response_summary.csv** (3 KB)
   - **Summary statistics**: Drugs tested, coverage, missingness

#### Clinical Data

5. **clinical_data.rds** (150 KB)
   - **Format**: R data structure (data.frame)
   - **Dimensions**: 671 patients × 95 variables
   - **Key columns**: Survival, age, sex, ELN risk, cytogenetics

#### Mutation Data

6. **mutations_matrix.csv** (in `03_Results/10_Mutations/`)
   - **Format**: Binary matrix (CSV)
   - **Dimensions**: 478 samples × ~200 genes
   - **Values**: 0 (wild-type), 1 (mutant)
   - **Key genes**: NPM1, FLT3, DNMT3A, IDH1, IDH2, TET2, TP53, RUNX1, ASXL1, NRAS, KRAS

---

### **ANALYSIS-READY DATA (03_Results/05_Analysis_Ready_Data/)**

7. **gold_standard_cohort.rds**
   - **Purpose**: 478 samples with complete 4-way multi-omics integration
   - **Includes**: Expression + Drug + Clinical + Mutations
   - **Format**: R list containing all data types with aligned sample IDs

8. **survival_data_with_clusters.csv** (in `03_Results/08_Survival_Analysis/`)
   - **Dimensions**: 671 samples × ~15 columns
   - **Key columns**:
     - `sample_id`: RNA sample ID
     - `cluster`: Cluster assignment (1 or 2)
     - `OS_months`: Overall survival in months
     - `OS_event`: Event indicator (0=censored, 1=death)
     - `age`, `sex`: Demographics
     - Clinical variables

---

### **PHASE 5 DRUG VALIDATION FILES (03_Results/23_Drug_Validation/)**

#### Primary Results (Already described above)
- `all_drugs_differential_response.csv` (155 drugs)
- `drug_cluster_independence_SIMPLIFIED.csv` (20 drugs, R² analysis)

#### Biological Validation

9. **bcl2_pathway_expression_FIXED.csv** (1.4 KB)
   - **Dimensions**: 10 genes × 8 columns
   - **Genes**: BCL2, BCL2L1, MCL1, BCL2L11, BAX, BAK1, BID, BBC3, BAD, PMAIP1
   - **Columns**:
     - `gene`: Gene symbol
     - `mean_cluster1`, `mean_cluster2`: Mean log2 expression
     - `log2fc`: Log2 fold-change (C1 vs C2)
     - `pvalue`: Wilcoxon test p-value
     - `fdr`: FDR-corrected p-value
     - `direction`: Which cluster has higher expression
   - **Key Finding**: BCL2 higher in Cluster 1 (p=8.55×10⁻²⁵), correlates with Venetoclax sensitivity (ρ=-0.552)

10. **immune_checkpoint_expression_FIXED.csv** (737 bytes)
    - **Dimensions**: 5 genes × 8 columns
    - **Genes**: CD47, BTLA, CTLA4, HAVCR2 (TIM-3), LAG3
    - **Key Finding**: CD47 higher in Cluster 1 (p=5.65×10⁻²⁸); BTLA, CTLA4, TIM-3 higher in Cluster 2

11. **drug_class_enrichment_FIXED.csv** (925 bytes)
    - **Dimensions**: 9 drug classes × 7 columns
    - **Classes**: BCL-2 inhibitors, MEK inhibitors, HDAC inhibitors, mTOR inhibitors, etc.
    - **Columns**: Class name, drugs tested, drugs significant, enrichment p-value, FDR
    - **Key Finding**: BCL-2, MEK, HDAC, mTOR inhibitors show 100% differential response

#### Publication-Ready Tables (Supplementary)

12-16. **Supplementary_Table_S5-S9.csv**
    - Enhanced versions of primary results with additional annotations
    - Ready for journal submission

---

### **VALIDATION COHORT DATA**

#### TCGA-LAML (Adult Validation)

17. **tcga_cluster_assignments.csv** (in `03_Results/17_TCGA_Validation/`)
    - **Dimensions**: 151 samples × 3 columns
    - **Purpose**: Apply 50-gene classifier to TCGA cohort
    - **Columns**: Sample ID, cluster, consensus_score

18. **tcga_survival_results.csv**
    - **Contains**: Survival analysis results (HR=1.24, p=0.291)
    - **Interpretation**: Non-significant due to low power (36.8%)

#### TARGET-AML (Pediatric Validation)

19. **target_cluster_assignments.csv** (in `03_Results/18_TARGET_Validation/`)
    - **Dimensions**: 1,713 samples × 3 columns
    - **Purpose**: Test age-specificity of subtypes
    - **Key Finding**: OPPOSITE prognostic effect (HR=0.81 vs 1.35 in adults)

---

## 📋 DATA SCHEMAS

### Master Sample ID Mapping Schema

```
unified_sample_id | expression_id | drug_response_id | clinical_id | mutation_id | has_expression | has_drug_response | has_clinical | has_mutations | n_data_types | cohort_category
-----------------|---------------|------------------|-------------|-------------|----------------|-------------------|--------------|---------------|--------------|----------------
BA2000R          | BA2000R       | BA2000D          | BA2000R     | BA2000D     | TRUE           | TRUE              | TRUE         | TRUE          | 4            | complete_quad_omics
```

### Cluster Assignment Schema

```
sample_id | cluster | consensus_score
----------|---------|----------------
BA2000R   | 1       | 0.742
BA2003R   | 1       | 0.852
```

### Drug Response Schema

```
drug        | n_samples | n_cluster1 | n_cluster2 | mean_auc_cluster1 | mean_auc_cluster2 | wilcoxon_pvalue | cohens_d | fdr
------------|-----------|------------|------------|-------------------|-------------------|-----------------|----------|----
Venetoclax  | 367       | 184        | 183        | 107.35            | 192.00            | 2.78e-24        | -1.25    | 4.31e-22
```

### Cluster Independence Schema

```
drug       | n_samples | r2_mutations_only | r2_mutations_plus_cluster | r2_improvement | pct_improvement | p_cluster_adds_value | fdr_cluster
-----------|-----------|-------------------|---------------------------|----------------|-----------------|---------------------|------------
Venetoclax | 367       | 0.140             | 0.365                     | 0.225          | 161%            | 4.73e-24            | 9.46e-23
```

---

## 🔗 DATA INTEGRATION WORKFLOW

### Step 1: Load Raw Data
```r
# Expression
expr <- read.delim("01_Data/BeatAML_Downloaded_Data/beataml_expression.txt")
# Dimensions: 22,843 genes × 671 samples

# Drug Response
drug <- read.delim("01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt")
# Format: Long (sample × drug × AUC)

# Clinical
clinical <- readxl::read_excel("01_Data/BeatAML_Downloaded_Data/beataml_clinical.xlsx")
# Dimensions: 671 samples × 95 variables

# Mutations
mutations <- read.delim("01_Data/BeatAML_Downloaded_Data/beataml_mutations.txt")
# Format: Variant-level (11,720 variants)
```

### Step 2: Load Processed Data
```r
# Analysis-ready expression (batch-corrected)
expr_bc <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

# Drug response matrix
drug_matrix <- readRDS("03_Results/01_Processed_Data/drug_response_auc.rds")

# Sample mapping
mapping <- read.csv("03_Results/01_Processed_Data/master_sample_id_mapping.csv")
```

### Step 3: Load Cluster Assignments
```r
# Cluster labels for all 671 samples
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
```

### Step 4: Merge for Analysis
```r
# Example: Merge drug data with clusters
library(dplyr)

# Get drug response for samples with cluster assignments
drug_cluster <- drug_matrix %>%
  rownames_to_column("sample_id") %>%
  inner_join(clusters, by = "sample_id")

# Result: 520 samples × 166 drugs + cluster assignment
```

### Step 5: Load Phase 5 Results
```r
# All drug differential response results
drug_results <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")

# Cluster independence analysis (R² improvement)
independence <- read.csv("03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv")

# BCL-2 pathway validation
bcl2 <- read.csv("03_Results/23_Drug_Validation/bcl2_pathway_expression_FIXED.csv")
```

---

## 📊 DATA DIMENSIONS SUMMARY

| Data Type | Raw Dimensions | Processed Dimensions | Format | Size |
|-----------|----------------|----------------------|--------|------|
| **Expression** | 22,843 genes × 671 samples | Same, log2 + batch-corrected | RDS | 99 MB |
| **Drug Response** | 520 samples × 166 drugs | Wide matrix, AUC values | RDS | 4.8 MB |
| **Clinical** | 671 samples × 95 variables | Cleaned, survival formatted | RDS | 150 KB |
| **Mutations** | 11,720 variants | 478 samples × 200 genes (binary) | CSV | Variable |
| **Clusters** | 671 samples × 1 label | + consensus score | CSV | 73 KB |
| **Phase 5 Drug** | 155 drugs × 15 statistics | All differential tests | CSV | 33 KB |
| **Phase 5 Independence** | 20 drugs × 11 R² metrics | Mutation + cluster models | CSV | 4.5 KB |

---

## 🔑 CRITICAL VARIABLES

### Sample Identifiers
- **Primary Key**: `unified_sample_id` or `sample_id`
- **RNA samples**: `BA####R` format (e.g., BA2000R)
- **DNA samples**: `BA####D` format (e.g., BA2000D)
- **Conversion**: Replace 'R' with 'D' for RNA→DNA matching

### Survival Variables
- `overallSurvival`: Days from diagnosis to death/last contact
- `OS_months`: Overall survival in months (derived)
- `vitalStatus` or `OS_event`: 0=alive/censored, 1=deceased/event

### Cluster Variables
- `cluster` or `cluster_assignment`: 1 or 2
- `consensus_score`: 0-1 (clustering stability)

### Drug Response Variables
- `auc`: Area Under Curve (LOWER = MORE SENSITIVE)
- Range: typically 0-300
- Interpretation: <100 = sensitive, >200 = resistant

### Mutation Variables
- Binary: 0 (wild-type), 1 (mutant)
- Key mutations: NPM1, FLT3, DNMT3A, IDH1, IDH2, TET2, TP53, RUNX1, ASXL1, NRAS, KRAS

---

## 📁 FILE FORMATS

### R Data Structures (.rds)
- Binary format, R-specific
- Fast to load, preserves R objects
- Load with: `readRDS("file.rds")`
- Most expression and processed data

### CSV Files (.csv)
- Text format, universal
- Load with: `read.csv("file.csv")`
- Results tables, mappings, annotations

### Tab-Delimited (.txt)
- Text format, tab-separated
- Load with: `read.delim("file.txt")`
- Raw BeatAML data files

### Excel Files (.xlsx)
- Binary spreadsheet format
- Load with: `readxl::read_excel("file.xlsx")`
- Clinical data, drug families

---

## 🎯 QUICK START: LOADING KEY DATA

### Minimal Code to Get Started

```r
# Set working directory
setwd("D:/Projects/Project_AML")

# Load essential packages
library(dplyr)
library(survival)

# 1. Load cluster assignments (671 samples)
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
# Columns: sample_id, cluster (1 or 2), consensus_score

# 2. Load drug response data (520 samples, 166 drugs)
drug_response <- readRDS("03_Results/01_Processed_Data/drug_response_auc.rds")
# Matrix format: rows = samples, columns = drugs

# 3. Load Phase 5 drug validation results (155 drugs tested)
drug_results <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")
# Key columns: drug, wilcoxon_pvalue, cohens_d, fdr, cluster_more_sensitive

# 4. Load cluster independence results (20 drugs, R² analysis)
independence <- read.csv("03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv")
# Key columns: drug, r2_mutations_only, r2_mutations_plus_cluster, r2_improvement, p_cluster_adds_value

# 5. Load survival data with clusters (671 samples)
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
# Key columns: sample_id, cluster, OS_months, OS_event, age, sex

# 6. Load expression data (batch-corrected, 22,843 genes × 671 samples)
expression <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")
# Matrix format: rows = genes, columns = samples
```

### Example Analysis: Venetoclax by Cluster

```r
# Load necessary data
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")
drug_response <- readRDS("03_Results/01_Processed_Data/drug_response_auc.rds")

# Extract Venetoclax AUC values
venetoclax_auc <- drug_response[, "Venetoclax"]  # Adjust column name if needed

# Create data frame with cluster assignments
venetoclax_data <- data.frame(
  sample_id = rownames(drug_response),
  venetoclax_auc = venetoclax_auc
) %>%
  inner_join(clusters, by = "sample_id") %>%
  filter(!is.na(venetoclax_auc))

# Test differential response
wilcox.test(venetoclax_auc ~ cluster, data = venetoclax_data)

# Summary statistics
venetoclax_data %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_auc = mean(venetoclax_auc),
    sd_auc = sd(venetoclax_auc),
    median_auc = median(venetoclax_auc)
  )
```

---

## 📌 IMPORTANT NOTES

### Sample Size Variations
- **Complete multi-omics**: 478 samples (all 4 data types)
- **Drug response cohort**: 520 samples (expression + drug + clinical)
- **Survival cohort**: 671 samples (expression + survival)
- **Multivariate analysis**: 459 samples (after removing missing mutation data)

### Missing Data
- **Expression**: 0.53% missing (excellent)
- **Drug response**: ~6% missing per drug (varies)
- **Clinical**: 18.6% missing (expected for retrospective)
- **Mutations**: Varies by gene (187 samples excluded from multivariate due to incomplete mutation calls)

### Data Quality
- All raw data MD5-verified
- Expression data batch-corrected (ComBat)
- Drugs with <30 samples excluded from analysis
- Only samples passing QC included

### File Naming Conventions
- **Phase-specific**: Results organized by phase number (01-23)
- **FIXED suffix**: Indicates corrected version after debugging
- **SIMPLIFIED suffix**: Simplified analysis approach
- **V2, V3 suffix**: Updated versions of documents

---

## 🔄 DATA UPDATE HISTORY

| Date | Update | Files Affected |
|------|--------|----------------|
| Oct 2, 2025 | Initial data processing | 01_Processed_Data/ |
| Oct 15, 2025 | Phase 4 completion | 21_Manuscript_Prep/ |
| Oct 16, 2025 | Phase 5 drug validation | 23_Drug_Validation/ |
| Oct 25, 2025 | Phase 5 manuscript integration | Root-level summaries updated |

---

## 📞 DATA ACCESS CONTACTS

### Original Data Sources
- **BeatAML**: https://www.nature.com/articles/s41586-018-0623-z (Tyner et al. Nature 2018)
- **TCGA-LAML**: https://portal.gdc.cancer.gov/projects/TCGA-LAML
- **TARGET-AML**: https://portal.gdc.cancer.gov/projects/TARGET-AML

### Data Availability
- Raw BeatAML data: `01_Data/BeatAML_Downloaded_Data/` (338 MB, 6 files)
- Processed data ready for analysis: `03_Results/01_Processed_Data/`
- Publication materials: `03_Results/23_Drug_Validation/` + `04_Figures/22_Drug_Validation/`

---

## 🎓 CITATION

If using this data structure or analysis:

**Proposed Citation**:
> [Authors]. "Molecular Subtypes in Adult AML Predict Venetoclax Response Independent of Genomic Alterations with Validated BCL-2 Mechanism." [Journal], [Year].
>
> Data from: Tyner JW, et al. "Functional genomic landscape of acute myeloid leukaemia." *Nature* 562, 526-531 (2018). doi:10.1038/s41586-018-0623-z

---

**Document Version**: 1.0
**Last Updated**: October 25, 2025
**Contact**: See project README for maintainer information
**Status**: ✅ Complete - All 5 phases documented
