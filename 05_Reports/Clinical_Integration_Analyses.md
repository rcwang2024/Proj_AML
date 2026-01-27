# Clinical Integration Analyses (TIER 2 Expansion)

**Generated:** 2025-10-02

This document provides detailed specifications for the clinical integration analyses, expanding TIER 2 of the comprehensive analysis roadmap.

---

## Analysis 2.4: Survival Analysis by Molecular Features

**Goal:** Identify prognostic molecular features in AML

**Methods:** Kaplan-Meier curves, Cox proportional hazards regression

**Data Requirements:** Survival + Expression + Mutations + Clinical

**Sample Size:**
- Required: n ≥ 100
- Available: n = 942
- With survival + expression: n = 671
- With survival + mutations: n = 615

**Statistical Power:** 0.85

**Feasibility:** YES
- *Justification:* n=942 total, 565 events (60.0%), 11/11 stratifications powered

**Stratifications:**
- Molecular Subtypes
- Mutation: FLT3
- Mutation: NPM1
- Mutation: TP53
- Mutation: DNMT3A
- Mutation: IDH1
- Mutation: IDH2
- Mutation: TET2
- Mutation: RUNX1
- Gene Expression Signatures
- Mutation Burden

**Timeline:** 2 weeks

**Scripts Location:** `02_Scripts/07_Survival_Analysis/`

**Deliverables:**
- Kaplan-Meier curves for all stratifications
- Hazard ratios with 95% confidence intervals
- Forest plots summarizing HR across features
- Univariate Cox regression results
- Multivariate Cox models (adjusted for clinical covariates)
- Log-rank test p-values
- Median survival times per group

---

## Analysis 2.5: Clinical-Molecular Correlation

**Goal:** Associate molecular features with clinical variables

**Methods:** Chi-square test, t-test, ANOVA, Fisher exact test

**Data Requirements:** Clinical + Expression + Mutations

**Sample Size:**
- Required: n ≥ 100
- Available: n = 615

**Statistical Power:** 0.72

**Feasibility:** YES
- *Justification:* n=615, 4/5 associations powered

**Associations to Test:**
- Mutations vs Age
- Mutations vs Gender
- Molecular Features vs De Novo/Secondary
- Gene Expression vs Blast %

**Timeline:** 1 week

**Scripts Location:** `02_Scripts/07_Survival_Analysis/`

**Deliverables:**
- Association tables with test statistics and p-values
- Visualization plots (boxplots, bar charts, heatmaps)
- Effect size calculations (Cohen's d, Cramer's V, odds ratios)
- Multiple testing correction (FDR, Bonferroni)
- Summary report of significant associations

---

## Analysis 2.6: Integrated Prognostic Model

**Goal:** Build comprehensive prognostic model combining clinical and molecular features

**Methods:** Multivariate Cox proportional hazards regression

**Data Requirements:** Survival + Clinical + Expression + Mutations

**Sample Size:**
- Required: n ≥ 100
- Available: n = 615
- Events (deaths): 368

**Statistical Power:** 0.85

**Feasibility:** YES
- *Justification:* n=615, 368 events, 26.3 events/predictor

**Feature Categories:**
- Clinical (5)
- Mutations (7)
- Expression (2)

**Model Recommendation:** Full model with all features feasible

**Timeline:** 2-3 weeks

**Scripts Location:** `02_Scripts/07_Survival_Analysis/`

**Deliverables:**
- Risk stratification model (low/intermediate/high risk)
- Nomogram for individual risk prediction
- C-index (concordance index) for model discrimination
- Calibration plots (predicted vs observed survival)
- Internal validation (bootstrap or cross-validation)
- Risk score formula and coefficients
- Comparison with existing prognostic models (e.g., ELN risk)

---

## Summary Table

| Analysis ID | Analysis Name | Sample Size | Power | Feasibility |
|-------------|---------------|-------------|-------|-------------|
| 2.4 | Survival Analysis by Molecular Features | 942 | 0.85 | YES |
| 2.5 | Clinical-Molecular Correlation | 615 | 0.72 | YES |
| 2.6 | Integrated Prognostic Model | 615 | 0.85 | YES |

**Summary:** 3/3 clinical integration analyses are fully feasible.
