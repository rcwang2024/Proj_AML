# Analysis Decisions Log

**Project:** Beat AML Multi-Omics Integration

**Version:** 1.0

**Last Updated:** 2025-10-02

---

## Purpose

This document records all major analysis decisions made during the BeatAML multi-omics integration project. Each decision is documented with:
- The decision made
- Rationale/justification
- Impact on analysis
- Date decided

---

## 1. Data Processing Decisions

### 1.1 Sample Inclusion Criteria

**Decision:** Include all samples with data available for each specific analysis

**Rationale:**
- Maximize sample size and statistical power
- Different analyses require different data type combinations
- Quality control will identify and handle outliers separately

**Impact:**
- Sample sizes vary by analysis type
- n=478 for gold standard quad-omics analyses
- n=707 for expression-only analyses
- Increases power for single data type analyses

**Date:** 2025-10-01

---

### 1.2 Outlier Handling Strategy

**Decision:**
- Identify outliers using multiple methods (PCA distance, correlation, clustering)
- Flag samples with ≥2 outlier flags for review
- Do NOT automatically exclude outliers
- Make case-by-case decisions before each major analysis

**Rationale:**
- Multiple outlier detection methods provide robust identification
- Biological vs technical outliers require different handling
- Some "outliers" may represent rare but valid biological subtypes
- Premature exclusion could eliminate interesting samples

**Impact:**
- 7 expression outliers identified (~1% of samples)
- Outliers flagged but retained for now
- Will revisit before molecular subtyping analysis

**Date:** 2025-10-01

---

### 1.3 Missing Data Strategy

**Decision:** Analysis-specific handling

**Approach:**
1. **Expression data:** Use complete cases (missing rate <0.1%)
2. **Drug response:** Accept incomplete drug panels (average 101/166 drugs)
3. **Clinical data:**
   - Use complete cases for survival analyses (100% complete)
   - Allow missingness for secondary variables
   - Do NOT impute clinical variables
4. **Mutations:** Absence of mutation = wild-type (not missing)

**Rationale:**
- Different data types have different missingness patterns
- Imputation not appropriate for categorical clinical data
- Drug data is intentionally incomplete (not all drugs on all samples)
- Missing mutation data means gene was wild-type in most cases

**Impact:**
- Maintains data integrity
- Sample size varies by analysis
- No artificial data creation through imputation

**Date:** 2025-10-01

---

### 1.4 Batch Correction Decision

**Decision:** APPLY batch correction to expression data before analysis

**Method:** ComBat (sva package) or limma::removeBatchEffect

**Batch Variable:** centerID (sequencing center)

**Rationale:**
- Significant batch effects detected:
  - PC1: F=24.5, p=9.07×10⁻³⁰
  - PC2: F=79.2, p=2.80×10⁻⁸³
- Batch effects will confound biological signal
- Critical for molecular subtyping and differential expression
- Standard practice in multi-center studies

**Impact:**
- Will create batch-corrected expression matrix
- Original data preserved for comparison
- Batch correction must be done BEFORE subtyping

**Implementation:**
- Correct for centerID while preserving biological variables
- Save corrected matrix as `beataml_expression_batchcorrected.txt`
- Validate correction via PCA

**Date:** 2025-10-02

**Status:** To be implemented in Phase 1

---

### 1.5 Gene Filtering

**Decision:** Retain all 22,843 genes for now; apply analysis-specific filters later

**Rationale:**
- Pre-filtering may remove biologically relevant genes
- Different analyses require different gene sets
- Will filter for:
  - Differential expression: genes with mean expression > threshold
  - Network analysis: high-variance genes
  - Drug prediction: feature selection during modeling

**Impact:**
- Full gene set available for all analyses
- Increases computational burden slightly
- Flexibility for different analyses

**Date:** 2025-10-01

---

### 1.6 Mutation Filtering

**Decision:** Apply minimal filtering

**Filters to Apply:**
- Keep: All variant types (missense, nonsense, frameshift, etc.)
- Keep: All VAF values (including low VAF <0.1 for subclonal mutations)
- Optional: Filter VAF <0.05 for driver-focused analyses

**Rationale:**
- Beat AML already applied quality filters
- Low VAF mutations may represent subclonal populations
- Different analyses have different VAF requirements

**Impact:**
- Comprehensive mutation data retained
- 11,721 total mutations across 871 samples
- Enables subclonal evolution analysis

**Date:** 2025-10-01

---

## 2. Cohort Definition Decisions

### 2.1 Primary Analysis Cohort

**Decision:** Define MULTIPLE cohorts for different analyses

**Cohorts Defined:**
1. **Gold Standard:** All 4 data types (n=478)
2. **Expression-Drug:** For drug prediction (n=494)
3. **Mutation-Expression:** For mut-expr integration (n=615)
4. **Mutation-Drug:** For mut-drug associations (n=583)
5. **Survival Cohort:** All with survival data (n=942)

**Rationale:**
- Maximizes sample size for each analysis type
- Gold standard cohort too restrictive for single-omics analyses
- Analysis-specific cohorts provide optimal power

**Impact:**
- Sample sizes vary by analysis (see cohort definitions document)
- More complex sample tracking required
- Increased statistical power

**Date:** 2025-10-01

---

### 2.2 Gold Standard Cohort Definition

**Decision:** Gold standard = all 4 data types (Expression + Mutations + Drug + Clinical)

**Sample Size:** n=478

**Use Cases:**
- Multi-omics integration
- Integrated prognostic models
- Personalized treatment frameworks

**Rationale:**
- Complete data enables full integration
- Sufficient sample size (power ≥0.85 for all major analyses)
- Represents ~49% of total cohort

**Impact:**
- Primary cohort for manuscript
- Enables comprehensive multi-omics analyses
- Trade-off between completeness and sample size

**Date:** 2025-10-01

---

## 3. Statistical Analysis Decisions

### 3.1 Significance Thresholds

**Decisions:**

**1. Differential Expression Analyses**
- FDR < 0.05 (Benjamini-Hochberg correction)
- |log2FC| > 1.0 for biological significance
- Mean expression > 1 for gene inclusion

**Rationale:**
- FDR controls false discovery rate in high-dimensional data
- log2FC threshold ensures biological relevance
- Standard thresholds in RNA-seq field

**2. Clinical Associations**
- p < 0.05 for univariate tests
- FDR < 0.10 for multiple comparisons (less stringent due to exploratory nature)

**Rationale:**
- Clinical associations are hypothesis-generating
- Balance Type I and Type II error

**3. Survival Analyses**
- p < 0.05 for log-rank tests
- Cox regression: p < 0.05 for individual features
- Multivariate models: likelihood ratio test p < 0.05

**4. Drug Prediction Models**
- Cross-validation R² > 0.3 for acceptable model
- R² > 0.5 for good model
- Test set validation required

**Date:** 2025-10-02

---

### 3.2 Multiple Testing Correction Methods

**Decision:** Use method appropriate for each analysis type

**Methods:**
1. **Benjamini-Hochberg FDR:** Default for most analyses
   - Differential expression
   - Gene-drug correlations
   - Mutation-expression associations

2. **Bonferroni:** For small number of pre-specified tests
   - Key driver mutation associations (n<10 tests)

3. **Permutation-based FDR:** For network analyses
   - When analytical FDR not appropriate

**Rationale:**
- FDR provides good balance of power and control
- Bonferroni too conservative for exploratory analyses
- Method should match analysis structure

**Date:** 2025-10-02

---

### 3.3 Sample Size Requirements

**Decision:** Minimum sample sizes for each analysis type

**Requirements:**
- Differential expression: ≥10 per group (prefer ≥20)
- Survival analysis: ≥50 events
- Machine learning: ≥100 total (prefer ≥150)
- Correlation analysis: ≥30 samples
- Network analysis: ≥50 samples

**Rationale:**
- Based on power calculations
- Ensures adequate statistical power (≥0.8)
- Standard recommendations from literature

**Impact:**
- Some subgroup analyses not feasible
- Rare mutations (<5%) have limited power

**Date:** 2025-10-02

---

## 4. Analysis-Specific Decisions

### 4.1 Molecular Subtyping

**Decision:** Use consensus clustering with k=3-5 clusters

**Method:**
- Consensus clustering (ConsensusClusterPlus)
- Features: Top 1000-5000 most variable genes
- Distance: 1 - Pearson correlation
- Clustering: Hierarchical (average linkage)
- Iterations: 1000

**Rationale:**
- Consensus clustering robust to sampling
- 3-5 clusters captures major AML subtypes
- High-variance genes drive biological subtypes

**Validation:**
- Silhouette width
- Gap statistic
- Clinical association (survival, ELN risk)

**Date:** 2025-10-02

---

### 4.2 Driver Mutation Definition

**Decision:** Use literature-based list of AML driver genes

**Driver Gene List:**
- Based on TCGA AML, Beat AML publications
- Genes with frequency ≥5% in cohort
- Known functional impact in AML

**Top Drivers:**
FLT3, NPM1, DNMT3A, TET2, IDH1, IDH2, TP53, NRAS, KRAS, RUNX1, ASXL1, SRSF2, U2AF1, CEBPA, WT1

**Rationale:**
- Evidence-based gene list
- Focus on recurrent, impactful mutations
- Balances comprehensiveness and focus

**Date:** 2025-10-01

---

### 4.3 Drug Response Thresholds

**Decision:** Define drug sensitivity based on AUC tertiles or quartiles

**Approach:**
- Sensitive: Bottom tertile (lowest AUC)
- Resistant: Top tertile (highest AUC)
- Or use continuous AUC values for correlation/regression

**Rationale:**
- No established AUC cutoffs in literature
- Tertile approach provides balanced groups
- Continuous approach preserves information

**Impact:**
- Flexible approach depending on analysis
- ~165 samples per group for tertile split

**Date:** 2025-10-02

---

### 4.4 Feature Selection for Predictive Modeling

**Decision:** Use multiple feature selection approaches

**Methods:**
1. Univariate filtering (top 1000 genes by variance)
2. Recursive feature elimination (RFE)
3. Elastic net automatic selection
4. Random forest importance

**Rationale:**
- No single method is universally best
- Ensemble of methods more robust
- Compare performance across methods

**Impact:**
- Multiple models to compare
- Computational burden increased
- Better model interpretability

**Date:** 2025-10-02

---

## 5. Software and Methods Decisions

### 5.1 Primary Analysis Software

**Decision:** R for statistical analyses, Python for machine learning

**R Packages:**
- Differential expression: DESeq2, limma
- Survival: survival, survminer
- Network: WGCNA
- Clustering: ConsensusClusterPlus
- Batch correction: sva (ComBat)

**Python Libraries:**
- Machine learning: scikit-learn
- Deep learning: TensorFlow/Keras (if needed)
- Data manipulation: pandas, numpy
- Visualization: seaborn, matplotlib

**Rationale:**
- R excellent for statistical genomics
- Python superior for ML/DL
- Leverage strengths of both languages

**Date:** 2025-10-01

---

### 5.2 Reproducibility Measures

**Decision:** Implement comprehensive reproducibility practices

**Practices:**
- Version control (Git)
- Set random seeds for all stochastic analyses
- Document package versions (sessionInfo())
- Script-based analysis (no manual steps)
- Docker containers for complex environments (if needed)

**Rationale:**
- Essential for scientific rigor
- Facilitates collaboration
- Enables replication by others

**Date:** 2025-10-01

---

## 6. Publication and Dissemination Decisions

### 6.1 Authorship

**Decision:** To be determined by project leadership

**Criteria:**
- Substantial intellectual contribution
- Writing/editing of manuscript
- Data analysis contribution
- Following ICMJE guidelines

**Date:** TBD

---

### 6.2 Data Sharing

**Decision:** Share processed data and code upon publication

**Plan:**
- Code: GitHub public repository
- Processed data: GEO or similar repository
- Raw data: Already available via dbGaP (phs001657)
- Supplementary tables: All significant associations

**Rationale:**
- Promotes reproducibility
- Benefits scientific community
- Increasingly required by journals

**Date:** 2025-10-02

---

## 7. Future Decisions

### Items Requiring Future Decision

1. **Exact clustering k value:** After initial analysis (Week 2-3)
2. **Outlier sample exclusions:** Before final analyses (Week 4)
3. **Final significance thresholds:** May adjust based on results
4. **External validation cohort:** Identify after initial results
5. **Specific drugs for detailed analysis:** Based on results

---

## Document History

**Version 1.0** - 2025-10-02 - Initial documentation of all major decisions

---

**For questions or to propose decision changes, contact project leadership.**
