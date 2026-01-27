# Data Structure - Quick Reference Card

**Last Updated**: October 25, 2025

---

## 🎯 ESSENTIAL FILES (Top 10)

| # | File Path | Purpose | Dimensions | Format |
|---|-----------|---------|------------|--------|
| 1 | `03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv` | ⭐⭐⭐ Cluster labels | 671 samples | CSV |
| 2 | `03_Results/23_Drug_Validation/all_drugs_differential_response.csv` | ⭐⭐⭐ Drug sensitivity | 155 drugs | CSV |
| 3 | `03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv` | ⭐⭐⭐ R² improvement | 20 drugs | CSV |
| 4 | `03_Results/01_Processed_Data/drug_response_auc.rds` | Drug AUC matrix | 520×166 | RDS |
| 5 | `03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds` | Expression | 22,843×671 | RDS |
| 6 | `03_Results/01_Processed_Data/master_sample_id_mapping.csv` | Sample integration | 1,067 samples | CSV |
| 7 | `03_Results/08_Survival_Analysis/survival_data_with_clusters.csv` | Survival + clusters | 671 samples | CSV |
| 8 | `03_Results/23_Drug_Validation/bcl2_pathway_expression_FIXED.csv` | BCL-2 validation | 10 genes | CSV |
| 9 | `01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt` | Raw drug data | ~86K rows | TXT |
| 10 | `01_Data/BeatAML_Downloaded_Data/beataml_expression.txt` | Raw expression | 22,843×671 | TXT |

---

## 📊 SAMPLE COHORTS

```
┌─────────────────────────────────────────────────┐
│ TOTAL SAMPLES IN PROJECT: 1,067                │
├─────────────────────────────────────────────────┤
│ Complete Quad-Omics     478  (Expression + Drug + Clinical + Mutations)
│ Drug Response Cohort    520  (Expression + Drug + Clinical)
│ Survival Cohort         671  (Expression + Clinical + Survival)
│ Multivariate Analysis   459  (Complete mutation data)
└─────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────┐
│ VALIDATION COHORTS                              │
├─────────────────────────────────────────────────┤
│ TCGA-LAML (Adult)       151  (Survival validation)
│ TARGET-AML (Pediatric) 1,713 (Age-specificity test)
│ TOTAL ACROSS 3 COHORTS 2,535
└─────────────────────────────────────────────────┘
```

---

## 🔗 SAMPLE ID MAPPING

```
RNA Sample ID:  BA2000R  ─┐
                          ├─> unified_sample_id: BA2000R
DNA Sample ID:  BA2000D  ─┘

Conversion Rule: Replace 'R' with 'D' (RNA → DNA)
```

**Sample ID Components**:
- `BA`: BeatAML prefix
- `2000`: Patient number
- `R`: RNA-seq sample
- `D`: DNA-seq sample (mutations, drug response)

---

## 📋 DATA TYPES & DIMENSIONS

| Data Type | Genes/Drugs/Variables | Samples | Matrix Size | File Size |
|-----------|----------------------|---------|-------------|-----------|
| **Expression** | 22,843 genes | 671 | 22,843 × 671 | 99 MB |
| **Drug Response** | 166 drugs | 520 | 520 × 166 | 4.8 MB |
| **Clinical** | 95 variables | 671 | 671 × 95 | 150 KB |
| **Mutations** | ~200 genes | 478 | 478 × 200 | Variable |
| **Clusters** | 2 subtypes | 671 | 671 × 1 | 73 KB |

---

## 🗂️ DIRECTORY STRUCTURE (Simplified)

```
Project_AML/
│
├── 01_Data/
│   ├── BeatAML_Downloaded_Data/    ⭐ Raw data (6 files, 338 MB)
│   ├── TCGA_LAML/                  Validation cohort
│   └── Gene_Sets/                  Reference annotations
│
├── 03_Results/
│   ├── 01_Processed_Data/          ⭐ Clean data, master mapping
│   ├── 06_Molecular_Subtypes/      ⭐⭐⭐ CLUSTER ASSIGNMENTS
│   ├── 23_Drug_Validation/         ⭐⭐⭐ PHASE 5 DRUG RESULTS
│   ├── 21_Manuscript_Prep/         Phase 4 materials
│   └── [22 other result folders]
│
├── 04_Figures/
│   ├── 21_Main_Figures/            Figures 1-4
│   └── 22_Drug_Validation/         ⭐ Figure 5 (Drug response)
│
└── [Root summaries]                ⭐ Start here
    ├── COMPLETE_PROJECT_SUMMARY_ALL_PHASES_V3.md
    ├── DATA_STRUCTURE_DOCUMENTATION.md
    └── KEY_STATISTICS_PHASE5_QUICK_REFERENCE.md
```

---

## 🔑 KEY VARIABLES

### Cluster Assignment
- **Variable**: `cluster` or `cluster_assignment`
- **Values**: 1 (NPM1+, favorable) or 2 (TP53+/RUNX1+, adverse)
- **File**: `sample_cluster_assignments.csv`

### Drug Sensitivity
- **Variable**: `auc` (Area Under Curve)
- **Interpretation**: **LOWER = MORE SENSITIVE**
- **Range**: 0-300 (typical)
- **Example**: Venetoclax AUC 107 (sensitive) vs 192 (resistant)

### Survival
- **Time**: `OS_months` (overall survival in months)
- **Event**: `OS_event` (0=censored, 1=death)
- **Format**: Right-censored survival data

### Mutations
- **Format**: Binary (0=wild-type, 1=mutant)
- **Key genes**: NPM1, FLT3, DNMT3A, IDH1, IDH2, TET2, TP53, RUNX1, ASXL1, NRAS, KRAS

---

## 📊 PHASE 5 KEY RESULTS

### All Drugs Differential Response
**File**: `all_drugs_differential_response.csv`

| Column | Description | Example |
|--------|-------------|---------|
| `drug` | Drug name | "Venetoclax" |
| `n_samples` | Total samples | 367 |
| `mean_auc_cluster1` | Cluster 1 mean AUC | 107.35 |
| `mean_auc_cluster2` | Cluster 2 mean AUC | 192.00 |
| `wilcoxon_pvalue` | Differential test | 2.78×10⁻²⁴ |
| `cohens_d` | Effect size | -1.25 |
| `fdr` | FDR-corrected p | 4.31×10⁻²² |
| `cluster_more_sensitive` | Which cluster responds | "Cluster_1" |

### Cluster Independence Analysis
**File**: `drug_cluster_independence_SIMPLIFIED.csv`

| Column | Description | Venetoclax |
|--------|-------------|------------|
| `drug` | Drug name | "Venetoclax" |
| `n_samples` | Sample size | 367 |
| `r2_mutations_only` | R² from mutations | 0.140 (14%) |
| `r2_mutations_plus_cluster` | R² with cluster | 0.365 (36.5%) |
| `r2_improvement` | Absolute gain | +0.225 |
| `pct_improvement` | Percent gain | +161% |
| `p_cluster_adds_value` | ANOVA p-value | 4.73×10⁻²⁴ |
| `fdr_cluster` | FDR correction | 9.46×10⁻²³ |

**INTERPRETATION**: Adding cluster to mutation model increases R² from 14% to 36.5% (+161%)

---

## 💻 QUICK LOAD COMMANDS (R)

```r
# Set working directory
setwd("D:/Projects/Project_AML")

# ⭐ ESSENTIAL 3 FILES ⭐

# 1. Cluster assignments
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# 2. Drug differential response (Phase 5)
drugs <- read.csv("03_Results/23_Drug_Validation/all_drugs_differential_response.csv")

# 3. Cluster independence (Phase 5 - R² improvement)
independence <- read.csv("03_Results/23_Drug_Validation/drug_cluster_independence_SIMPLIFIED.csv")

# ADDITIONAL USEFUL FILES

# Drug response matrix (520 samples × 166 drugs)
drug_auc <- readRDS("03_Results/01_Processed_Data/drug_response_auc.rds")

# Expression matrix (batch-corrected)
expr <- readRDS("03_Results/04_Batch_Corrected/beataml_expression_batchcorrected.rds")

# Survival data with clusters
survival <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Sample mapping (integrates all data types)
mapping <- read.csv("03_Results/01_Processed_Data/master_sample_id_mapping.csv")
```

---

## 🎯 DATA ANALYSIS CHEAT SHEET

### Find Venetoclax-Sensitive Patients
```r
# Load data
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Venetoclax-sensitive = Cluster 1 (NPM1+)
venetoclax_sensitive <- clusters %>% filter(cluster == 1)
# Result: ~220 samples (42%)
```

### Calculate Drug Response by Cluster
```r
library(dplyr)

# Load drug AUC matrix
drug_auc <- readRDS("03_Results/01_Processed_Data/drug_response_auc.rds")

# Load clusters
clusters <- read.csv("03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv")

# Merge
drug_cluster <- drug_auc %>%
  rownames_to_column("sample_id") %>%
  inner_join(clusters, by = "sample_id")

# Test one drug
wilcox.test(Venetoclax ~ cluster, data = drug_cluster)
```

### Get Samples with Complete Data
```r
# Load mapping
mapping <- read.csv("03_Results/01_Processed_Data/master_sample_id_mapping.csv")

# Filter to complete quad-omics
complete <- mapping %>% filter(cohort_category == "complete_quad_omics")
# Result: 478 samples
```

---

## 📌 IMPORTANT REMINDERS

### Drug AUC Interpretation
- **Lower AUC** = **MORE sensitive** to drug
- Example: Venetoclax AUC 107 (Cluster 1) vs 192 (Cluster 2)
- Cluster 1 is MORE sensitive (79% lower AUC)

### Cluster Labels
- **Cluster 1**: NPM1+/DNMT3A+, favorable-risk, Venetoclax-sensitive
- **Cluster 2**: TP53+/RUNX1+/ASXL1+, adverse-risk, Venetoclax-resistant
- Never swap labels arbitrarily!

### Missing Data
- Drug response: ~6% missing (varies by drug)
- Not all samples have all data types
- Use `master_sample_id_mapping.csv` to check availability

### File Formats
- `.rds` = R-specific binary (fast, preserves structure)
- `.csv` = Universal text (slower but readable anywhere)
- `.txt` = Tab-delimited text (raw BeatAML data)

---

## 🔍 FINDING SPECIFIC DATA

**"I need cluster assignments"**
→ `03_Results/06_Molecular_Subtypes/sample_cluster_assignments.csv`

**"I need drug response data"**
→ Raw: `01_Data/BeatAML_Downloaded_Data/beataml_drug_auc.txt`
→ Processed: `03_Results/01_Processed_Data/drug_response_auc.rds`

**"I need Phase 5 drug validation results"**
→ `03_Results/23_Drug_Validation/` (all files)

**"I need to integrate expression + drug + cluster"**
→ Use `master_sample_id_mapping.csv` as bridge

**"I need Venetoclax-specific results"**
→ Row 2 in `all_drugs_differential_response.csv`
→ Row 2 in `drug_cluster_independence_SIMPLIFIED.csv`

**"I need survival data"**
→ `03_Results/08_Survival_Analysis/survival_data_with_clusters.csv`

**"I need mutation data"**
→ `03_Results/10_Mutations/mutation_matrix.csv`

---

## 📊 SAMPLE COUNTS BY ANALYSIS

| Analysis | N Samples | Data Required |
|----------|-----------|---------------|
| Clustering | 671 | Expression + Clinical |
| Drug Response (Phase 5) | 520 | Expression + Drug + Clinical |
| Multivariate (survival) | 459 | Expression + Clinical + Mutations (complete) |
| Complete Multi-Omics | 478 | Expression + Drug + Clinical + Mutations |
| TCGA Validation | 151 | External cohort |
| TARGET Validation | 1,713 | External cohort (pediatric) |

---

## ✅ DATA QUALITY CHECKS

**Expression Data**: ✅
- 0.53% missing (excellent)
- Batch-corrected (ComBat)
- Log2-transformed

**Drug Response Data**: ✅
- 6.04% missing (acceptable for screening data)
- 166 compounds tested
- Ex vivo assay (BeatAML protocol)

**Clinical Data**: ✅
- 18.6% missing (expected for retrospective)
- Survival data available for 671 samples
- ELN risk groups classified

**Mutation Data**: ✅
- Variant-level validated
- Binary matrix created
- Key driver genes covered

---

**Quick Reference Version**: 1.0
**Last Updated**: October 25, 2025
**For Full Details**: See `DATA_STRUCTURE_DOCUMENTATION.md`
