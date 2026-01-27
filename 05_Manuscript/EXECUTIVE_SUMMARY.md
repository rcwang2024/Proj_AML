# AML Molecular Subtyping Project - Executive Summary

## One-Sentence Summary
**Transcriptomic molecular subtypes in adult AML predict Venetoclax response independent of genomic mutations, providing a clinically actionable biomarker for precision therapy selection.**

---

## Project Overview

| Metric | Value |
|--------|-------|
| **Total Patients** | 2,535 |
| **Cohorts** | 3 (BeatAML, TCGA-LAML, TARGET-AML) |
| **Analysis Phases** | 6 complete |
| **Drugs Tested** | 155 |
| **Drugs Differential** | 72 (46.5%) |
| **Key Drug** | Venetoclax (p=2.78×10⁻²⁴) |

---

## The Two Molecular Subtypes

### Cluster 1 (45% of patients) - "NPM1-like / Venetoclax-Sensitive"
- **NPM1 mutations**: 42% (vs 10% in C2)
- **DNMT3A mutations**: 31% (vs 16% in C2)
- **Drug profile**: Hypersensitive to BCL-2 inhibitors
- **Venetoclax AUC**: 107 (highly sensitive)

### Cluster 2 (55% of patients) - "Adverse-Risk / Chemoresistant"
- **TP53 mutations**: 14% (vs 5% in C1)
- **RUNX1 mutations**: 18% (vs 5% in C1)
- **ASXL1 mutations**: 16% (vs 4% in C1)
- **Drug profile**: More sensitive to MEK/HDAC/mTOR inhibitors
- **Venetoclax AUC**: 192 (relatively resistant)

---

## Key Statistical Results

### Survival (Adult Meta-Analysis)
| Metric | Value |
|--------|-------|
| Pooled HR | 1.35 |
| 95% CI | 1.13-1.62 |
| p-value | 0.001 |
| I² heterogeneity | 0% |

### Multivariate Independence (Survival)
| Variable | HR | p-value |
|----------|-----|---------|
| **Cluster** | **1.06** | **0.649 (NOT independent)** |
| TP53 | 2.96 | 5.6×10⁻¹⁰ |
| TET2 | 1.42 | 0.031 |

### Drug Response Independence (Venetoclax)
| Model | R² | Improvement |
|-------|-----|-------------|
| Mutations only | 0.140 | — |
| Mutations + Cluster | 0.365 | **+161%** |
| p-value for cluster | 4.7×10⁻²⁴ | **INDEPENDENT** |

---

## Robustness Validation (Phase 6)

| Test | Result | Interpretation |
|------|--------|----------------|
| Bootstrap (10K) | 100% p<0.001 | Exceptional stability |
| LOOCV | 100% significant | No outlier-driven |
| Permutation (10K) | p<0.0001 | True signal |
| Sample-Split | p=3.2×10⁻¹² | Not circular |
| VIF | Max 1.59 | No multicollinearity |

**All 10 top drugs show "Exceptional" robustness**

---

## The Critical Distinction

### For PROGNOSIS (Survival Prediction):
❌ **Clusters are NOT independent** of mutations (p=0.649)
- TP53 dominates survival prediction
- Use standard ELN risk stratification

### For TREATMENT SELECTION (Drug Response):
✅ **Clusters ARE independent** of mutations
- 19/20 drugs: significant R² improvement (FDR<0.05)
- Mean +42% improvement over mutation-only models
- Venetoclax: +161% R² improvement

**This is NOT a paradox - different outcomes have different determinants**

---

## Venetoclax: The Headline Finding

| Metric | Value |
|--------|-------|
| **p-value** | 2.78×10⁻²⁴ |
| **Cohen's d** | -1.25 (very large) |
| **Fold difference** | 1.79× more sensitive (C1) |
| **R² improvement** | +161% beyond mutations |
| **BCL-2 mechanism** | 9/10 genes validated |
| **Bootstrap stability** | 100% (10,000 resamples) |
| **Held-out validation** | p=3.2×10⁻¹² |

### BCL-2 Pathway Validation
- BCL2 expression **higher in Cluster 1** (sensitive cluster)
- Higher target → greater dependency → better response
- 9/10 BCL-2 family genes differentially expressed (FDR<0.05)

---

## Clinical Translation Pathway

```
Current Status: RETROSPECTIVE VALIDATION COMPLETE
     ↓
Next Step: PROSPECTIVE CLINICAL TRIAL
     ↓
     Randomize patients by molecular subtype
     Test Venetoclax-based vs alternative regimens
     ↓
Future: COMPANION DIAGNOSTIC
     50-gene RT-qPCR panel
     Integrate with mutation testing
```

---

## Publication Strategy

### Recommended Title
**"Transcriptomic Molecular Subtypes Predict Venetoclax Response Independent of Genomic Alterations in Adult Acute Myeloid Leukemia"**

### Target Journals (Tier 1)
1. **Nature Medicine** (IF: 87) - Translational biomarker, FDA-approved drug
2. **Journal of Clinical Oncology** (IF: 51) - Clinical actionability
3. **Blood** (IF: 26) - Comprehensive AML study

### Key Selling Points
1. ✅ Large multi-cohort validation (2,535 patients)
2. ✅ Clinically actionable (Venetoclax is FDA-approved)
3. ✅ Mutation-INDEPENDENT for treatment (novel finding)
4. ✅ Mechanistically validated (BCL-2 pathway)
5. ✅ Exceptionally robust statistics (Phase 6)
6. ✅ Portable classifier (50 genes, 93% accuracy)

---

## Files Generated

### Phase 6 Robustness Validation
```
03_Results/24_Robustness_Validation/
├── Bootstrap_Top10_Drugs.csv
├── LOOCV_Top10_Drugs.csv
├── Permutation_Top10_Drugs.csv
├── Sample_Split_Validation.csv
├── VIF_Analysis_Summary.csv
├── Independence_Paradox_Explanation.txt
└── PHASE6_COMPREHENSIVE_REPORT.md
```

### Manuscript Materials
```
05_Manuscript/
├── MANUSCRIPT_DRAFT_COMPLETE.md (this document references)
└── EXECUTIVE_SUMMARY.md (this file)
```

---

## Bottom Line

**We discovered that transcriptomic molecular subtypes in adult AML:**

1. **Cannot replace** mutation testing for prognosis
2. **CAN augment** mutation testing for treatment selection
3. **Identify Venetoclax-responsive patients** with extraordinary accuracy
4. **Provide orthogonal information** to standard genomic profiling
5. **Are ready for prospective clinical validation**

**The 50-gene classifier is a potential companion diagnostic for Venetoclax-based therapy in adult AML.**

---

*Last Updated: November 2025*
*Project Status: PHASE 6 COMPLETE - PUBLICATION READY*
