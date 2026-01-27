# TARGET-AML Validation - Analysis Plan

## Status: IN PROGRESS (Background Download)

**Date**: 2025-10-11
**Script**: `02_Scripts/Phase3_CriticalValidation/10_target_aml_validation.R`
**Background Process ID**: 3dea27

---

## Overview

TARGET-AML provides an **independent pediatric AML cohort** (ages 0-20) for validating our molecular subtypes. This is a critical test of whether our adult-derived subtypes extend to pediatric populations.

### Why TARGET-AML is Important

1. **Age diversity**: Tests generalizability beyond adult populations
2. **Independent cohort**: Different institution, sequencing platform, treatment protocols
3. **Biological insight**: If subtypes work in pediatric AML, suggests age-independent mechanism
4. **Clinical relevance**: ~20% of AML cases are pediatric

---

## Data Download

### Dataset Characteristics
- **Project**: TARGET-AML (Therapeutically Applicable Research to Generate Effective Treatments)
- **Source**: NCI Genomic Data Commons (GDC)
- **Data type**: RNA-Seq (STAR - Counts)
- **Samples**: 3,227 files
- **Size**: 13.6 GB
- **Download time**: Estimated 30-60 minutes

### Technical Details
- Using `TCGAbiolinks` R package for automated download
- Data downloaded in 14 chunks to manage file size
- Workflow: STAR-aligned gene counts

---

## Analysis Pipeline

### 1. Data Processing
- Convert Ensembl IDs to gene symbols
- Match genes with BeatAML reference
- Normalize using log2(CPM+1) for consistency
- Filter to common genes across cohorts

### 2. Clinical Data Extraction
- **Outcome**: Overall survival (OS)
- **Variables**: Age, sex, vital status, follow-up time
- **Age range**: 0-20 years (pediatric)

### 3. Classifier Application
- Apply BeatAML-trained 50-gene Random Forest classifier
- Predict cluster assignments (Cluster 1 vs Cluster 2)
- Calculate prediction confidence scores
- Handle missing genes via mean imputation

### 4. Survival Validation
**Primary endpoint**: Log-rank test for survival difference
**Secondary analyses**:
- Cox proportional hazards regression
- Hazard ratio with 95% CI
- Median survival comparison
- Kaplan-Meier curves

### 5. Cross-Cohort Comparison
**Cohorts**:
- BeatAML (n=~700, adult, discovery)
- TCGA-AML (n=~150, adult, validation)
- TARGET-AML (n=?, pediatric, validation)

**Metrics**:
- Heterogeneity testing (Cochran's Q, I² statistic)
- Forest plot of hazard ratios
- Age-stratified analysis

---

## Expected Outcomes

### Scenario 1: Strong Validation ✓✓✓
- **Result**: Significant survival difference (p<0.05)
- **HR**: Consistent with adults (HR ~1.3-1.5)
- **Interpretation**: Subtypes show age-independent biology
- **Impact**: Supports universal applicability of classifier

### Scenario 2: Partial Validation ✓
- **Result**: Significant but different effect size
- **HR**: Outside adult range (e.g., HR<1.2 or >2.0)
- **Interpretation**: Age-specific effect modification
- **Impact**: Suggests pediatric-specific calibration needed

### Scenario 3: No Validation ✗
- **Result**: No significant survival difference (p>0.05)
- **Possible reasons**:
  1. Insufficient power (low sample size/events)
  2. Age-specific biology (subtypes don't exist in pediatric)
  3. Treatment differences mask prognostic effect
  4. Classifier doesn't transfer to pediatric population
- **Impact**: Limits generalizability claims

---

## Key Questions to Answer

1. **Do molecular subtypes exist in pediatric AML?**
   - Are cluster proportions similar to adults?
   - What is prediction confidence?

2. **Do subtypes predict survival in children?**
   - Is the effect significant?
   - Is the direction consistent with adults?

3. **Is the effect size similar to adults?**
   - Compare HR across cohorts
   - Test for heterogeneity

4. **What are age-specific considerations?**
   - Compare age distributions
   - Assess age × cluster interaction

5. **Does sample size provide adequate power?**
   - Calculate achieved power
   - Estimate required events for definitive test

---

## Output Files

### Data Files (when complete)
- `target_aml_counts_raw.rds` - Raw gene counts
- `target_aml_metadata.rds` - Clinical metadata
- `target_aml_expression_normalized.rds` - Normalized expression matrix
- `target_aml_clinical.csv` - Processed clinical data

### Results Files
- `target_cluster_assignments.csv` - Predicted clusters with confidence
- `target_survival_validation.csv` - Survival analysis results

### Figures
- `target_kaplan_meier.pdf` - Survival curves by cluster
- `forest_plot_all_cohorts.pdf` - Cross-cohort HR comparison

---

## Monitoring Progress

To check download progress:
```r
# In R console
library(BashOutput)
BashOutput("3dea27")  # Background process ID
```

Or check for output files:
```bash
ls -lh 03_Results/18_TARGET_Validation/
```

---

## Next Steps After Completion

### If Validation Successful
1. Update meta-analysis to include TARGET-AML
2. Add to manuscript as pediatric validation
3. Emphasize age-independent biology
4. Update composite forest plot

### If Validation Unsuccessful
1. Perform power analysis
2. Document as limitation
3. Discuss age-specific biology hypothesis
4. Recommend pediatric-specific classifier development

### Always
1. Update Phase3_CriticalValidation_Summary.md
2. Add TARGET results to cross-cohort comparison
3. Create comprehensive validation summary
4. Prepare supplementary figures

---

## Timeline

- **Download**: 30-60 minutes (in progress)
- **Processing**: 5-10 minutes
- **Analysis**: 2-5 minutes
- **Total**: ~1 hour

**Expected completion**: ~2025-10-11 23:00-23:30

---

## Contact/References

**TARGET-AML Project**: https://ocg.cancer.gov/programs/target/projects/acute-myeloid-leukemia
**GDC Portal**: https://portal.gdc.cancer.gov/projects/TARGET-AML
**Publication**: Bolouri et al. (2018) Nature Medicine - TARGET AML genomics

---

**Status**: ⏳ Downloading... (Started 2025-10-11 22:03)
**Last Updated**: 2025-10-11 22:15
**Estimated Completion**: 2025-10-11 23:00
