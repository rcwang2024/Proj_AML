# Analysis Cohort Definitions - AML Multi-Omics Project

**Document Version:** 1.0
**Date:** October 2, 2025
**Total Sample Universe:** 970 unique samples

---

## 📋 Table of Contents

1. [Cohort Overview](#cohort-overview)
2. [Primary Cohorts](#primary-cohorts)
3. [Secondary Cohorts](#secondary-cohorts)
4. [Exploratory Cohorts](#exploratory-cohorts)
5. [Statistical Power Analysis](#statistical-power-analysis)
6. [Cohort Selection Guide](#cohort-selection-guide)
7. [Quality Metrics](#quality-metrics)

---

## Cohort Overview

### Data Type Availability

| Data Type | Samples Available | Coverage |
|-----------|------------------|----------|
| **Expression** | 707 | 72.9% |
| **Drug Response** | 603 | 62.2% |
| **Clinical** | 934 | 96.3% |
| **Mutations** | 871 | 89.8% |

### Cohort Hierarchy

```
Total Samples (n=970)
    │
    ├── Gold Standard Cohort (n=478) ⭐ PRIMARY
    │   └── All 4 data types
    │
    ├── Triple Omics Cohorts (n=153)
    │   ├── Cohort A: E+D+C (n=16)
    │   └── Cohort B: E+M+C (n=137)
    │
    ├── Dual Omics Cohorts (n=195)
    │   ├── E+D (n=494 at-least, 0 exactly-2)
    │   └── E+M (n=615 at-least, 0 exactly-2)
    │
    └── Full Data Type Cohorts
        ├── Full Expression (n=707)
        ├── Full Clinical (n=934)
        └── Others...
```

---

## Primary Cohorts

### 1. Gold Standard Cohort ⭐

**Sample Size:** n = 478

**Data Types Included:**
- ✅ Gene Expression (22,843 genes)
- ✅ Drug Response (166 compounds)
- ✅ Clinical Data (95 variables)
- ✅ Somatic Mutations (WES, 3,333 genes)

**Intended Use:**
1. **Primary multi-omics integration analyses**
   - Comprehensive molecular characterization
   - Multi-omics clustering and subtyping
   - Pathway enrichment with complete context

2. **Drug response prediction**
   - Molecular features → drug sensitivity
   - Mutation + expression signatures → AUC/IC50
   - Machine learning models with full feature set

3. **Survival analysis with complete molecular context**
   - Cox regression with multi-omics features
   - Risk stratification models
   - Biomarker discovery across all data types

4. **Mutation-expression-drug correlation**
   - Driver mutations → expression changes → drug response
   - Pathway activity → drug sensitivity
   - Resistance mechanism identification

**Strengths:**
- ✅ **Complete molecular profiles** - no missing data types
- ✅ **Largest complete cohort** - 478 samples is excellent for multi-omics
- ✅ **Balanced representation** - 49% of total samples
- ✅ **Statistical power** - sufficient for complex models (10-15 variables)
- ✅ **Clinical outcomes** - 100% survival data available
- ✅ **Drug diversity** - mean 101 drugs tested per sample

**Limitations:**
- ⚠️ **Selection bias** - samples had to pass QC for all 4 assays
- ⚠️ **Cost constraints** - may exclude some patient populations
- ⚠️ **Reduced sample size** - 478 vs 970 total

**Statistical Power:**
- **Survival analysis:** ~287 events expected (60% mortality) - excellent
- **Drug response:** ~100 samples/drug average - good for correlation
- **Machine learning:** Train (287) / Validate (96) / Test (95) - adequate
- **Subgroup analysis:** 3-4 groups with ~120 samples each - good
- **Multivariate models:** 478/15 = 32 samples per variable (rule of thumb: 10-15) - good

**Recommended For:**
- 🌟 **All primary research questions**
- 🌟 **Publications requiring complete multi-omics**
- 🌟 **Drug response prediction models**
- 🌟 **Comprehensive survival models**

---

## Secondary Cohorts

### 2. Triple Omics Cohort A: Expression + Drug Response + Clinical

**Sample Size:** n = 16

**Data Types Included:**
- ✅ Gene Expression (22,843 genes)
- ✅ Drug Response (166 compounds)
- ✅ Clinical Data (95 variables)
- ❌ Somatic Mutations (not available)

**Intended Use:**
1. **Pharmacogenomics without mutation context**
   - Expression-drug correlation
   - Gene signature → drug response
   - Clinical + expression → drug sensitivity

2. **Validation cohort for expression-drug models**
   - Test models developed in gold standard cohort
   - Assess if mutations are required for prediction

**Strengths:**
- ✅ Complete pharmacogenomics data (expression + drug + clinical)
- ✅ Can assess expression-only drug prediction
- ✅ Clinical context for adjustment

**Limitations:**
- ⚠️ **Very small sample size** (n=16) - limited statistical power
- ⚠️ **No mutation data** - cannot assess mutation-drug relationships
- ⚠️ **Insufficient for primary analyses** - use for validation only

**Statistical Power:**
- **Limited** - n=16 only suitable for:
  - Validation of existing models
  - Descriptive analyses
  - Case studies
- **Not recommended for:**
  - De novo model building
  - Subgroup analyses
  - Survival analyses

**Recommended For:**
- ⚠️ **Validation only** - test expression-drug models
- ⚠️ **Exploratory** - assess need for mutation data
- ❌ **Not for primary analyses** - too small

---

### 3. Triple Omics Cohort B: Expression + Mutations + Clinical

**Sample Size:** n = 137

**Data Types Included:**
- ✅ Gene Expression (22,843 genes)
- ✅ Clinical Data (95 variables)
- ✅ Somatic Mutations (WES, 3,333 genes)
- ❌ Drug Response (not available)

**Intended Use:**
1. **Mutation-expression correlation**
   - Driver mutations → gene expression changes
   - Mutation → pathway activity
   - eQTL-like analysis (mutation as variant)

2. **Molecular subtyping without pharmacogenomics**
   - Expression + mutation clustering
   - Pathway-based classification
   - Molecular risk stratification

3. **Survival analysis with genomic + transcriptomic features**
   - Cox regression with expression + mutations
   - Prognostic signature development
   - Risk prediction models

**Strengths:**
- ✅ **Good sample size** (n=137) - adequate for many analyses
- ✅ **Complete genomic context** - mutations + expression
- ✅ **Clinical outcomes** - survival data available
- ✅ **Sufficient power** - can build multivariate models

**Limitations:**
- ⚠️ **No drug response data** - cannot predict drug sensitivity
- ⚠️ **Smaller than gold standard** - less power than n=478

**Statistical Power:**
- **Survival analysis:** ~82 events expected - good
- **Expression-mutation correlation:** n=137 adequate for discovery
- **Machine learning:** Train (82) / Validate (27) / Test (28) - borderline
- **Subgroup analysis:** 2-3 groups feasible
- **Multivariate models:** 137/15 = 9 samples per variable - borderline

**Recommended For:**
- ✅ **Mutation-expression correlation studies**
- ✅ **Molecular subtyping (non-drug)**
- ✅ **Prognostic model development**
- ✅ **Pathway analysis with genomic context**
- ❌ **Not for drug response prediction**

---

## Exploratory Cohorts

### 4. Expression-Drug Cohort: Expression + Drug Response (at least)

**Sample Size:** n = 494 (at least; includes gold standard)

**Data Types Included:**
- ✅ Gene Expression (22,843 genes)
- ✅ Drug Response (166 compounds)
- ✅ Clinical Data (494/494 have clinical - 100%)
- ⚠️ Mutations (478/494 have mutations - 97%)

**Composition:**
- Gold Standard (E+D+C+M): 478 samples
- Triple (E+D+C only): 16 samples
- Total: 494 samples

**Intended Use:**
1. **Maximized pharmacogenomics analyses**
   - Largest cohort for expression-drug correlation
   - Gene signature discovery for drug response
   - Expression-based drug sensitivity prediction

2. **Drug response analyses with maximum power**
   - ~83 samples per drug (494 samples / 166 drugs)
   - Better statistical power than gold standard alone

**Strengths:**
- ✅ **Largest pharmacogenomics cohort** (n=494)
- ✅ **Excellent power** for expression-drug correlation
- ✅ **100% clinical coverage** - all samples have clinical
- ✅ **97% mutation coverage** - most have mutations too

**Limitations:**
- ⚠️ **Heterogeneous** - 16 samples lack mutations
- ⚠️ **Need to account for missing data** in mutation analyses

**Statistical Power:**
- **Excellent** for expression-drug correlation
- ~99 samples per drug average (increased from 101 in gold standard)
- Multivariate models: 494/15 = 33 samples per variable - excellent

**Recommended For:**
- ✅ **Primary pharmacogenomics discovery**
- ✅ **Gene signature development for drug response**
- ✅ **Maximum power drug-expression correlation**
- ⚠️ **Be cautious** - account for 16 samples without mutations

---

### 5. Expression-Mutation Cohort: Expression + Mutations (at least)

**Sample Size:** n = 615 (at least; includes gold standard)

**Data Types Included:**
- ✅ Gene Expression (22,843 genes)
- ✅ Somatic Mutations (WES, 3,333 genes)
- ✅ Clinical Data (615/615 have clinical - 100%)
- ⚠️ Drug Response (478/615 have drug - 78%)

**Composition:**
- Gold Standard (E+D+C+M): 478 samples
- Triple (E+C+M only): 137 samples
- Total: 615 samples

**Intended Use:**
1. **Maximized mutation-expression correlation**
   - Largest cohort for genomic-transcriptomic integration
   - Driver mutation → expression changes
   - Pathway activity from mutations + expression

2. **Molecular classification with genomic context**
   - Clustering with expression + mutations
   - Subtype discovery
   - Pathway-based classification

3. **Survival analysis with complete molecular data**
   - Cox regression with expression + mutation features
   - Prognostic signatures

**Strengths:**
- ✅ **Largest genomic-transcriptomic cohort** (n=615)
- ✅ **Excellent power** for mutation-expression correlation
- ✅ **100% clinical coverage**
- ✅ **78% drug coverage** - reasonable for exploratory

**Limitations:**
- ⚠️ **Heterogeneous drug data** - 137 samples lack drug response
- ⚠️ **Cannot use for primary drug analyses** without imputation

**Statistical Power:**
- **Excellent** for mutation-expression correlation
- ~369 events expected for survival - excellent
- Multivariate models: 615/15 = 41 samples per variable - excellent

**Recommended For:**
- ✅ **Primary mutation-expression correlation**
- ✅ **Molecular subtyping with genomic context**
- ✅ **Large-scale survival modeling**
- ✅ **Pathway enrichment analysis**
- ⚠️ **Not primary for drug response** (use Cohort 4 instead)

---

### 6. Full Expression Cohort: All with Expression Data

**Sample Size:** n = 707

**Data Types Included:**
- ✅ Gene Expression (22,843 genes)
- ⚠️ Clinical Data (667/707 have clinical - 94%)
- ⚠️ Drug Response (494/707 have drug - 70%)
- ⚠️ Mutations (615/707 have mutations - 87%)

**Intended Use:**
1. **Expression-only analyses**
   - Gene expression patterns across all samples
   - Clustering and subtyping without other data
   - Batch effect assessment
   - Technical QC

2. **Maximum power expression studies**
   - Differential expression
   - Co-expression networks
   - Gene signature discovery

**Strengths:**
- ✅ **Largest expression cohort** (n=707)
- ✅ **Maximum statistical power** for expression analyses
- ✅ **94% clinical coverage** - most have outcomes

**Limitations:**
- ⚠️ **Incomplete other data types**
  - 30% lack drug response
  - 13% lack mutations
  - 6% lack clinical
- ⚠️ **Cannot integrate other omics** without subsetting

**Statistical Power:**
- **Excellent** for expression-only analyses
- Differential expression: very high power
- Clustering: large sample size enables fine subtype detection

**Recommended For:**
- ✅ **Expression-only discovery**
- ✅ **Technical QC and normalization**
- ✅ **Gene co-expression networks**
- ✅ **Initial clustering and exploration**
- ❌ **Not for multi-omics** (use appropriate subset)

---

### 7. Full Clinical Cohort: All with Clinical Data

**Sample Size:** n = 934

**Data Types Included:**
- ✅ Clinical Data (95 variables, 100% complete survival)
- ⚠️ Expression (667/934 have expression - 71%)
- ⚠️ Mutations (871/934 have mutations - 93%)
- ⚠️ Drug Response (599/934 have drug - 64%)

**Intended Use:**
1. **Clinical-only analyses**
   - Survival analysis with clinical variables only
   - Epidemiological analyses
   - Clinical risk models (without molecular data)

2. **Maximum power clinical studies**
   - Clinical factor associations
   - Outcome prediction from clinical features
   - Baseline characteristics

3. **Comparison cohort**
   - Compare molecular vs clinical-only models
   - Assess added value of omics data

**Strengths:**
- ✅ **Largest cohort** (n=934)
- ✅ **Maximum power** for clinical analyses
- ✅ **100% survival data** - ~560 events expected
- ✅ **93% mutation coverage** - most have genomic data

**Limitations:**
- ⚠️ **Incomplete molecular data**
  - 29% lack expression
  - 7% lack mutations
  - 36% lack drug response
- ⚠️ **Cannot integrate omics** without subsetting

**Statistical Power:**
- **Excellent** for clinical-only analyses
- Survival analysis: ~560 events - very high power
- Multivariate clinical models: 934 samples - excellent

**Recommended For:**
- ✅ **Clinical-only survival models**
- ✅ **Baseline cohort characterization**
- ✅ **Clinical risk factor studies**
- ✅ **Comparison with molecular models**
- ❌ **Not for omics integration** (use appropriate subset)

---

## Statistical Power Analysis

### Sample Size Requirements (General Guidelines)

| Analysis Type | Minimum n | Adequate n | Optimal n | Gold Standard | Other Cohorts |
|---------------|-----------|------------|-----------|---------------|---------------|
| **Cox Regression** | 100 | 200 | 500 | ✅ (478) | Cohort B (137) |
| **Multivariate (10-15 vars)** | 150 | 300 | 500 | ✅ (478) | Cohorts 4,5,6,7 |
| **Machine Learning** | 200 | 400 | 1000 | ✅ (478) | Cohorts 5,6,7 |
| **Subgroup Analysis** | 50/group | 100/group | 200/group | ✅ (120/group) | Cohort 5,6,7 |
| **Correlation Studies** | 30 | 100 | 300 | ✅ (478) | All |
| **Drug Response** | 50 | 100 | 200 | ✅ (101 avg) | Cohort 4 |

### Power by Cohort

| Cohort | Survival | Drug Response | Mutation-Expr | ML Models | Subgroups |
|--------|----------|---------------|---------------|-----------|-----------|
| **Gold Standard (478)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Triple A (16)** | ⭐ | ⭐ | N/A | ⭐ | ⭐ |
| **Triple B (137)** | ⭐⭐⭐ | N/A | ⭐⭐⭐ | ⭐⭐ | ⭐⭐ |
| **Expr-Drug (494)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| **Expr-Mut (615)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Full Expr (707)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Full Clinical (934)** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |

⭐ = Limited, ⭐⭐ = Borderline, ⭐⭐⭐ = Adequate, ⭐⭐⭐⭐ = Good, ⭐⭐⭐⭐⭐ = Excellent

---

## Cohort Selection Guide

### Decision Tree: Which Cohort to Use?

```
START: What is your research question?

Q1: Do you need drug response data?
├─ YES → Q2: Do you need mutation data too?
│   ├─ YES → Use Gold Standard (n=478) ⭐
│   └─ NO → Use Expression-Drug Cohort (n=494)
└─ NO → Q3: Do you need both expression and mutations?
    ├─ YES → Q4: Maximizing sample size important?
    │   ├─ YES → Use Expression-Mutation Cohort (n=615)
    │   └─ NO → Use Gold Standard (n=478) ⭐
    └─ NO → Q5: What's your primary data type?
        ├─ Expression → Use Full Expression Cohort (n=707)
        ├─ Mutations → Use Full Clinical Cohort with mutations (n=871)
        └─ Clinical → Use Full Clinical Cohort (n=934)
```

### By Analysis Type

| Analysis | Recommended Cohort | Alternative | Avoid |
|----------|-------------------|-------------|-------|
| **Drug Response Prediction** | Gold Standard (478) | Expr-Drug (494) | Triple A (16) |
| **Mutation-Expression** | Expr-Mut (615) | Gold Standard (478) | Triple A (16) |
| **Survival with Omics** | Gold Standard (478) | Expr-Mut (615) | Triple A (16) |
| **Molecular Subtyping** | Gold Standard (478) | Expr-Mut (615) | Triple A (16) |
| **Expression Signatures** | Full Expr (707) | Expr-Mut (615) | - |
| **Clinical Risk Models** | Full Clinical (934) | Gold Standard (478) | - |
| **Pathway Analysis** | Gold Standard (478) | Expr-Mut (615) | - |
| **Biomarker Discovery** | Gold Standard (478) | Expr-Drug (494) | Triple A (16) |

### By Research Goal

**Comprehensive Multi-Omics Integration:**
→ **Gold Standard Cohort (n=478)** ⭐

**Maximum Power Pharmacogenomics:**
→ **Expression-Drug Cohort (n=494)**

**Maximum Power Genomic-Transcriptomic:**
→ **Expression-Mutation Cohort (n=615)**

**Expression Discovery (no omics integration):**
→ **Full Expression Cohort (n=707)**

**Clinical-Only or Baseline Studies:**
→ **Full Clinical Cohort (n=934)**

**Validation Studies Only:**
→ **Triple Cohort A or B** (use with caution)

---

## Quality Metrics

### Data Completeness by Cohort

| Cohort | Expression | Drug | Clinical | Mutations | Overall |
|--------|-----------|------|----------|-----------|---------|
| **Gold Standard** | 100% | 100% | 100% | 100% | 100% ✅ |
| **Triple A** | 100% | 100% | 100% | 0% | 75% |
| **Triple B** | 100% | 0% | 100% | 100% | 75% |
| **Expr-Drug** | 100% | 100% | 100% | 97% | 99% ✅ |
| **Expr-Mut** | 100% | 78% | 100% | 100% | 95% |
| **Full Expr** | 100% | 70% | 94% | 87% | 88% |
| **Full Clinical** | 71% | 64% | 100% | 93% | 82% |

### Sample Overlap Summary

```
Total Samples: 970
│
├─ Complete (All 4): 478 (49.3%) ⭐
├─ Triple (3 types): 153 (15.8%)
├─ Dual (2 types): 195 (20.1%)
└─ Single (1 type): 39 (4.0%)

Missing: 105 samples (10.8%) with 0 data types
```

---

## Recommendations Summary

### Primary Analyses - Use These:

1. **🌟 Gold Standard Cohort (n=478)** - Default for all multi-omics
2. **Expression-Drug Cohort (n=494)** - Pharmacogenomics with max power
3. **Expression-Mutation Cohort (n=615)** - Genomic-transcriptomic with max power

### Secondary Analyses - Use These:

4. **Full Expression Cohort (n=707)** - Expression-only discovery
5. **Full Clinical Cohort (n=934)** - Clinical-only models

### Validation Only - Use Cautiously:

6. **Triple Cohort A (n=16)** - Too small for primary analyses
7. **Triple Cohort B (n=137)** - Borderline power, use for specific questions

---

## Implementation Notes

### Sample Selection

Use `master_sample_id_mapping.csv` to select cohorts:

```python
# Example: Gold Standard Cohort
gold = df[(df['has_expression']) &
          (df['has_drug_response']) &
          (df['has_clinical']) &
          (df['has_mutations'])]

# Example: Expression-Drug Cohort
expr_drug = df[(df['has_expression']) &
               (df['has_drug_response'])]
```

### Quality Control

For each cohort:
1. ✅ Verify sample counts match definitions
2. ✅ Check data completeness for required types
3. ✅ Assess missing data patterns
4. ✅ Evaluate potential selection bias
5. ✅ Document any deviations

### Reporting

When publishing, report:
- Cohort used and rationale
- Sample size and selection criteria
- Data completeness statistics
- Any exclusions or filtering applied
- Power calculations performed

---

## Cohort Summary Table

| # | Cohort Name | n | Data Types | Primary Use | Power | Priority |
|---|-------------|---|------------|-------------|-------|----------|
| 1 | **Gold Standard** | **478** | E+D+C+M | All multi-omics | ⭐⭐⭐⭐⭐ | **PRIMARY** |
| 2 | Triple A (E+D+C) | 16 | E+D+C | Validation only | ⭐ | Low |
| 3 | Triple B (E+M+C) | 137 | E+M+C | Mut-Expr correlation | ⭐⭐⭐ | Medium |
| 4 | Expression-Drug | 494 | E+D+(C+M) | Pharmacogenomics | ⭐⭐⭐⭐⭐ | **HIGH** |
| 5 | Expression-Mutation | 615 | E+M+C+(D) | Genomic-transcriptomic | ⭐⭐⭐⭐⭐ | **HIGH** |
| 6 | Full Expression | 707 | E+(others) | Expression discovery | ⭐⭐⭐⭐⭐ | Medium |
| 7 | Full Clinical | 934 | C+(others) | Clinical models | ⭐⭐⭐⭐⭐ | Medium |

**Legend:** E=Expression, D=Drug, C=Clinical, M=Mutations, +(X) = partially available

---

## Contact and Updates

**Document Owner:** AML Multi-Omics Project Team
**Last Updated:** October 2, 2025
**Next Review:** After Phase 4 completion

**Related Documents:**
- Sample overlap analysis: `sample_overlap_analysis.csv`
- Master ID mapping: `master_sample_id_mapping.csv`
- Gold standard list: `gold_standard_cohort_samples.txt`
- Phase 3 summary: `PHASE3_SAMPLE_INTEGRATION_SUMMARY.md`

---

**Document Version:** 1.0
**Status:** ✅ Complete and Ready for Use
