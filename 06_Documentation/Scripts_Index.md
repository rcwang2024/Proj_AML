# Analysis Scripts Index

**Project:** Beat AML Multi-Omics Integration

**Generated:** 2025-10-02

**Total Scripts:** 17

---

## 01_Data_Processing/

### 01_verify_files.py
**Location:** `02_Scripts/01_Data_Processing/01_verify_files.py`
**Purpose:** Verify presence and integrity of all downloaded BeatAML data files
**Inputs:** Raw data files in `01_Data/BeatAML_Downloaded_Data/`
**Outputs:** Console output with file verification results
**Runtime:** <1 minute

---

### 02_inspect_data.py
**Location:** `02_Scripts/01_Data_Processing/02_inspect_data.py`
**Purpose:** Initial exploration of all data types
**Inputs:** All 5 BeatAML data files
**Outputs:**
- Console output with data dimensions
- Basic statistics for each dataset
**Runtime:** ~2 minutes

---

### 03_sample_inventory.py
**Location:** `02_Scripts/01_Data_Processing/03_sample_inventory.py`
**Purpose:** Create unified sample inventory across all data types
**Inputs:** All BeatAML data files
**Outputs:**
- `03_Results/01_Processed_Data/sample_inventory.csv`
- Sample counts per data type
**Runtime:** ~2 minutes

---

### 04_drug_response_analysis.py
**Location:** `02_Scripts/01_Data_Processing/04_drug_response_analysis.py`
**Purpose:** Analyze drug response data and create summary statistics
**Inputs:** `beataml_drug_auc.txt`
**Outputs:**
- `drug_response_summary.csv`
- `samples_drug_counts.csv`
- `top_20_drugs.csv`
**Runtime:** ~2 minutes

---

### 05_clinical_data_analysis.py
**Location:** `02_Scripts/01_Data_Processing/05_clinical_data_analysis.py`
**Purpose:** Analyze clinical data completeness and distributions
**Inputs:** `beataml_clinical.xlsx`
**Outputs:**
- `clinical_data_summary.csv`
- Survival statistics
**Runtime:** ~2 minutes

---

### 06_expression_data_analysis.py
**Location:** `02_Scripts/01_Data_Processing/06_expression_data_analysis.py`
**Purpose:** Analyze gene expression data characteristics
**Inputs:** `beataml_expression.txt`
**Outputs:**
- Expression statistics
- `gene_detection_per_sample.csv`
**Runtime:** ~5 minutes
**Memory:** ~4 GB (large matrix)

---

### 07_mutation_data_analysis.py
**Location:** `02_Scripts/01_Data_Processing/07_mutation_data_analysis.py`
**Purpose:** Analyze mutation data and identify driver mutations
**Inputs:** `beataml_mutations.txt`
**Outputs:**
- `driver_mutation_frequencies.csv`
- Mutation burden statistics
**Runtime:** ~2 minutes

---

### 08_sample_overlap_and_mapping.py
**Location:** `02_Scripts/01_Data_Processing/08_sample_overlap_and_mapping.py`
**Purpose:** Create master sample mapping and analyze overlaps
**Inputs:** All data files
**Outputs:**
- `master_sample_id_mapping.csv`
- `sample_overlap_upset.png`
- Sample overlap statistics
**Dependencies:** upsetplot, matplotlib
**Runtime:** ~3 minutes

---

## 02_Quality_Control/

### 01_batch_effect_assessment.py
**Location:** `02_Scripts/02_Quality_Control/01_batch_effect_assessment.py`
**Purpose:** Detect batch effects in expression data via PCA and statistical testing
**Inputs:**
- `beataml_expression.txt`
- `beataml_clinical.xlsx`
**Outputs:**
- `batch_effect_assessment.txt`
- `pca_variance_explained.csv`
- `pca_biplot_pc1_pc2.png`
- `pca_by_batch_variable.png`
- `batch_effect_boxplots.png`
**Dependencies:** sklearn, scipy
**Runtime:** ~5 minutes
**Memory:** ~4 GB
**Key Finding:** Significant batch effects detected (centerID)

---

### 02_expression_quality_metrics.py
**Location:** `02_Scripts/02_Quality_Control/02_expression_quality_metrics.py`
**Purpose:** Comprehensive quality assessment of expression data
**Methods:**
- PCA-based outlier detection
- Sample correlation analysis
- Hierarchical clustering
**Inputs:** `beataml_expression.txt`
**Outputs:**
- `expression_outliers.csv`
- `pca_outliers.png`
- `sample_correlation_metrics.png`
- `sample_correlation_heatmap.png`
**Dependencies:** sklearn, scipy
**Runtime:** ~10 minutes
**Memory:** ~4 GB
**Key Finding:** 7 outliers identified (~1%)

---

### 03_comprehensive_qc_final.py
**Location:** `02_Scripts/02_Quality_Control/03_comprehensive_qc_final.py`
**Purpose:** Comprehensive QC for drug, clinical, mutation data, and missing data analysis
**Inputs:** All data files, master mapping
**Outputs:**
- `drug_response_qc.csv`
- `clinical_completeness.csv`
- `mutation_data_qc.csv`
- `missing_data_comprehensive_report.csv`
- Multiple QC figures
**Runtime:** ~5 minutes

---

## 03_Power_Analysis/

### 01_statistical_power_analysis.py
**Location:** `02_Scripts/03_Power_Analysis/01_statistical_power_analysis.py`
**Purpose:** Calculate statistical power for 6 major analysis types
**Inputs:**
- Master mapping
- Clinical data
- Driver frequencies
**Outputs:**
- `statistical_power_analysis.csv`
- Power estimates for all analyses
**Dependencies:** statsmodels
**Runtime:** ~2 minutes
**Key Finding:** All 6 analyses feasible (power ≥0.8)

---

### 02_comprehensive_analysis_roadmap.py
**Location:** `02_Scripts/03_Power_Analysis/02_comprehensive_analysis_roadmap.py`
**Purpose:** Generate initial comprehensive analysis roadmap (Tiers 1-3)
**Inputs:** Master mapping, driver frequencies
**Outputs:**
- `comprehensive_analysis_roadmap.csv`
- `Comprehensive_Analysis_Roadmap.md` (partial)
**Runtime:** ~3 minutes

---

### 03_clinical_integration_analyses.py
**Location:** `02_Scripts/03_Power_Analysis/03_clinical_integration_analyses.py`
**Purpose:** Detailed power analysis for clinical integration analyses (Tier 2 expansion)
**Inputs:** Master mapping, clinical data
**Outputs:**
- `clinical_integration_analyses.csv`
- `Clinical_Integration_Analyses.md`
**Runtime:** ~3 minutes

---

### 04_tier3_advanced_analyses.py
**Location:** `02_Scripts/03_Power_Analysis/04_tier3_advanced_analyses.py`
**Purpose:** Power analysis for Tier 3 advanced integrative analyses
**Inputs:** Master mapping
**Outputs:**
- `tier3_advanced_analyses.csv`
**Runtime:** ~2 minutes

---

### 05_generate_final_roadmap.py
**Location:** `02_Scripts/03_Power_Analysis/05_generate_final_roadmap.py`
**Purpose:** Consolidate all roadmap components into final comprehensive document
**Inputs:** All roadmap CSVs
**Outputs:**
- `Analysis_Roadmap.md` (final comprehensive)
- `consolidated_analysis_roadmap.csv`
**Runtime:** ~2 minutes

---

## 04_Documentation/

### 01_generate_inventory_report.py
**Location:** `02_Scripts/04_Documentation/01_generate_inventory_report.py`
**Purpose:** Generate comprehensive BeatAML Data Inventory Report
**Inputs:** All analysis results, QC outputs, power analysis
**Outputs:**
- `05_Reports/BeatAML_Data_Inventory_Report.md` (comprehensive)
**Runtime:** ~3 minutes
**Report Length:** ~400 lines, 10 sections

---

## Usage Instructions

### Running Scripts in Order

**Phase 1: Data Processing (Week 1)**
```bash
cd 02_Scripts/01_Data_Processing
python 01_verify_files.py
python 02_inspect_data.py
python 03_sample_inventory.py
python 04_drug_response_analysis.py
python 05_clinical_data_analysis.py
python 06_expression_data_analysis.py
python 07_mutation_data_analysis.py
python 08_sample_overlap_and_mapping.py
```

**Phase 2: Quality Control (Week 1-2)**
```bash
cd 02_Scripts/02_Quality_Control
python 01_batch_effect_assessment.py
python 02_expression_quality_metrics.py
python 03_comprehensive_qc_final.py
```

**Phase 3: Power Analysis (Week 2)**
```bash
cd 02_Scripts/03_Power_Analysis
python 01_statistical_power_analysis.py
python 02_comprehensive_analysis_roadmap.py
python 03_clinical_integration_analyses.py
python 04_tier3_advanced_analyses.py
python 05_generate_final_roadmap.py
```

**Phase 4: Documentation (Week 2)**
```bash
cd 02_Scripts/04_Documentation
python 01_generate_inventory_report.py
```

### Common Dependencies

All scripts require:
- Python 3.8+
- pandas
- numpy
- pathlib (standard library)

Additional packages by script:
- **Visualization:** matplotlib, seaborn
- **Statistics:** scipy, statsmodels
- **Machine Learning:** scikit-learn
- **Excel reading:** openpyxl

### Installing Dependencies

```bash
pip install pandas numpy matplotlib seaborn scipy scikit-learn statsmodels openpyxl upsetplot
```

Or use requirements file:
```bash
pip install -r requirements.txt
```

### Memory Requirements

- Most scripts: <2 GB RAM
- Expression analysis scripts: 4-8 GB RAM
- Network analysis (future): 8-16 GB RAM

### Runtime Estimates

- **Total Phase 1:** ~20 minutes
- **Total Phase 2:** ~20 minutes
- **Total Phase 3:** ~12 minutes
- **Total Phase 4:** ~5 minutes

**Overall:** ~1 hour for complete analysis pipeline

---

## Script Development Guidelines

### For Future Scripts

**Naming Convention:**
- `##_descriptive_name.py` (numbered by execution order)
- Use underscores, lowercase
- Be descriptive

**Required Elements:**
1. Docstring at top with purpose, author, date
2. UTF-8 encoding handling for Windows
3. Pathlib for paths (not strings)
4. Create output directories if needed
5. Print progress messages
6. Save outputs with descriptive names

**Example Template:**
```python
\"\"\"
Phase X: Description
Task X.X: Specific task name

Author: AML Multi-Omics Project
Date: YYYY-MM-DD
\"\"\"

import pandas as pd
from pathlib import Path
import sys

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Paths
DATA_DIR = Path("D:/Projects/Project_AML/01_Data/BeatAML_Downloaded_Data")
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/...")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TASK DESCRIPTION")
print("=" * 80)
print()

# Analysis code here...

print("=" * 80)
print("TASK COMPLETE")
print("=" * 80)
```

---

**For questions about scripts, refer to this index or contact project team.**
