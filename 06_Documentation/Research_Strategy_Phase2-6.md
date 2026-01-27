# BeatAML Multi-Omics Integration: Research Strategy for Phases 2-6

**Project Status:** Phase 1 COMPLETE ✅ | Ready for Phase 2 Molecular Subtyping
**Date:** 2025-10-04
**Target Timeline:** 12 weeks for complete analysis → manuscript submission

---

## I. RESEARCH CONTEXT & POSITIONING

### Current AML Research Landscape (2024-2025)

#### Major Advances:
- **12+ new approved drugs** since 2017 (venetoclax, FLT3i, IDH inhibitors, menin inhibitors)
- **Multi-omics classification** replacing genomic-only approaches
- **MRD detection** at 10⁻⁶ sensitivity predicting relapse
- **Precision medicine trials** (Beat AML Master Trial: genomic profiling in 7 days)

#### Key Research Gaps We're Addressing:

**Gap 1: Transcriptomic Subtyping Beyond Genetics**
- Current: ELN risk based on mutations/cytogenetics only
- Evidence: Transcriptomic subtypes reveal heterogeneity WITHIN genetic groups
- Recent: Studies show 8-17 transcriptomic clusters with distinct biology
- **Our Opportunity:** Define clinically actionable transcriptomic subtypes in Beat AML (n=707)

**Gap 2: Multi-Omics Integration for Drug Response**
- Current: Most models use genomics OR transcriptomics separately
- Evidence: Network/pathway approaches outperform gene-level models
- **Our Opportunity:** Integrate RNA-seq + mutations + clinical for robust prediction (n=478)

**Gap 3: Pathway-Level Analysis**
- Current: Gene-level analysis has high noise, limited interpretability
- Evidence: Pathway enrichment creates more stable, interpretable clusters
- **Our Opportunity:** Use pathway scores for molecular subtyping and drug response

**Gap 4: Clinical Validation of Molecular Subtypes**
- Current: Many studies identify subtypes without clinical utility validation
- Evidence: Need correlation with survival, treatment response, ELN risk
- **Our Opportunity:** Validate against survival and drug response in large cohort (n=942 survival, n=478 drugs)

**Gap 5: Limited Drug Response Integration**
- Current: Most AML studies lack ex vivo drug response data
- Evidence: Beat AML unique with 122 compounds tested
- **Our Opportunity:** Map drug response profiles to molecular subtypes

### Our Project's Strategic Advantages

**Unique Strengths:**
1. **Beat AML database access** (~970 AML samples, multi-omics)
2. **Comprehensive drug response** (122 compounds, n=478 with full data)
3. **Pure AML cohort** (no MPN contamination)
4. **Multi-omics integration potential** (expression + mutations + clinical + drugs)
5. **Large survival cohort** (n=942, 60% event rate)

**Realistic Constraints:**
1. **Drug response sample size** ~478 samples (not all 970)
2. **No proteomics** in our dataset (recent Beat AML proteomics not yet public)
3. **Batch effects** addressed ✅ (ComBat correction complete)
4. **Cross-sectional design** (no longitudinal data)
5. **Bulk tissue** (no single-cell resolution)

---

## II. COMPREHENSIVE ANALYSIS PLAN

### Phase 2: Molecular Subtyping (Weeks 3-4)

#### **Objective:**
Identify robust, biologically meaningful transcriptomic subtypes in AML using batch-corrected expression data

#### **Analysis 2.1: Gene Expression-Based Clustering**

**Method:** Consensus clustering with multiple algorithms

**Workflow:**
1. **Feature Selection:**
   - Top 5,000 most variable genes (by MAD or CV)
   - OR use published AML gene signatures (LSC17, immune signatures)
   - Filter: median expression > 1 (log2 scale)

2. **Consensus Clustering:**
   - Algorithms: Hierarchical, k-means, PAM, NMF
   - Distance: Euclidean, Pearson correlation
   - Test k=2 to k=10 clusters
   - Iterations: 1,000 with 80% subsampling
   - Tools: ConsensusClusterPlus R package

3. **Optimal k Selection:**
   - Consensus CDF (empirical cumulative distribution)
   - Delta Area plot (relative change in area under CDF)
   - Silhouette scores
   - Gap statistic
   - **Decision criteria:** Maximize stability + biological interpretability

4. **Cluster Stability:**
   - Silhouette width for each sample
   - Cluster purity across methods
   - Reproducibility in random subsets

**Deliverables:**
- Optimal number of clusters (k) with statistical justification
- Sample-to-cluster assignments for all 707 samples
- Consensus matrix heatmap
- PCA/UMAP colored by clusters
- Stability metrics (silhouette scores, consensus scores)

---

#### **Analysis 2.2: Pathway Enrichment-Based Clustering**

**Rationale:** More stable and interpretable than gene-level clustering

**Workflow:**
1. **Pathway Databases:**
   - KEGG pathways (metabolism, signaling, DNA repair)
   - Hallmark gene sets (MSigDB)
   - GO Biological Processes
   - **Custom AML pathways:**
     - Immune: NK cells, T-cell signaling, antigen presentation
     - Metabolism: glycolysis, OXPHOS, fatty acid metabolism
     - DNA repair: MMR, homologous recombination, NHEJ
     - Oncogenic: PI3K/AKT, MAPK, Wnt, Hedgehog, JAK-STAT
     - Stem cell: LSC signatures, hematopoietic stem cell markers

2. **Pathway Scoring:**
   - Method: ssGSEA (single-sample GSEA) or GSVA
   - Output: Pathway enrichment matrix (pathways × samples)
   - Normalization: Z-score across samples

3. **Clustering on Pathway Scores:**
   - Same consensus clustering as Analysis 2.1
   - Compare pathway-based vs gene-based clusters
   - Assess which provides better clinical separation

4. **Integration:**
   - Map pathway clusters to gene-based clusters
   - Identify concordance and discordance
   - Use both for comprehensive subtype definition

**Deliverables:**
- Pathway enrichment score matrix (pathways × 707 samples)
- Pathway-based cluster assignments
- Comparison heatmap (gene-based vs pathway-based)
- Pathway enrichment profiles per cluster

---

#### **Analysis 2.3: Subtype Characterization**

**Objective:** Deeply characterize each molecular subtype

**Tasks:**

1. **Differential Expression Analysis:**
   - For each cluster vs all others
   - Method: limma or DESeq2
   - Threshold: FDR < 0.05, |log2FC| > 1.0
   - Top 100 upregulated/downregulated genes per cluster
   - Functional annotation (GO, KEGG, Reactome)

2. **Mutation Enrichment:**
   - Fisher's exact test: which mutations enriched in which clusters?
   - **Expected patterns to validate:**
     - NPM1-mutant cluster
     - TP53-mutant cluster
     - FLT3-ITD enrichment
     - DNMT3A/TET2/IDH co-mutations
   - **Novel associations:** Any unexpected mutation-cluster pairs?

3. **Clinical Associations:**
   - Age distribution across clusters (ANOVA or Kruskal-Wallis)
   - ELN risk distribution (chi-square)
   - De novo vs secondary AML (chi-square)
   - FAB subtypes (chi-square)
   - WBC count, blast percentage (ANOVA)

4. **Biological Interpretation & Naming:**
   - Name clusters based on dominant biology:
     - "Immune-enriched" (high immune pathway scores)
     - "Proliferative/DNA repair-high"
     - "Stem-like/quiescent" (LSC signature high)
     - "Metabolically active" (glycolysis/OXPHOS high)
     - "NPM1-mutant" (if clearly defined)
     - "TP53-mutant/complex karyotype"

5. **Comparison to Published Subtypes:**
   - Cheng et al. 2022 (8 transcriptomic subgroups)
   - Song et al. 2025 (UAMOCS 3 subtypes)
   - Pathway enrichment studies
   - TCGA-LAML classifications

**Deliverables:**
- Differential expression tables (1 per cluster)
- Mutation enrichment bar plots
- Clinical association heatmaps
- Subtype naming and biological interpretation document
- Comparison table with published subtypes

**Timeline:** 2-3 weeks

---

### Phase 3: Survival & Clinical Outcome Analysis (Week 5)

#### **Analysis 3.1: Survival Analysis by Molecular Subtype**

**Objective:** Determine if molecular subtypes predict clinical outcome

**Workflow:**

1. **Kaplan-Meier Survival Curves:**
   - Overall survival (OS) by molecular subtype
   - Event-free survival (EFS) by molecular subtype
   - Log-rank test for significance
   - Median survival times with 95% CI
   - 1-year, 3-year, 5-year survival rates

2. **Cox Proportional Hazards Regression:**
   - **Univariate:** Subtype as predictor
     - Hazard ratios (HR) with 95% CI
     - P-values for each subtype vs reference

   - **Multivariate:** Adjust for established risk factors
     - Covariates: Age (<60 vs ≥60), ELN risk, treatment intensity
     - Assess if subtypes add value beyond ELN risk
     - Test for interactions (e.g., subtype × age)

3. **Stratified Analysis:**
   - **Within ELN risk groups:** Do subtypes further stratify survival?
     - Favorable risk: Can we identify high-risk patients?
     - Adverse risk: Can we identify lower-risk patients?

   - **Age-stratified:** (<60 vs ≥60 years)
   - **Treatment-stratified:** Intensive vs non-intensive therapy

4. **Prognostic Value Comparison:**
   - **C-index (concordance index):**
     - ELN risk alone
     - Molecular subtypes alone
     - ELN risk + subtypes combined
     - LSC17 score (if available)

   - **Decision curve analysis:** Clinical utility at different threshold probabilities

5. **Reclassification Analysis:**
   - Net reclassification improvement (NRI)
   - Integrated discrimination improvement (IDI)
   - Do subtypes correctly reclassify patients compared to ELN risk?

**Deliverables:**
- Kaplan-Meier curves (OS and EFS) with p-values
- Cox regression tables (univariate and multivariate)
- Forest plots showing hazard ratios
- C-index comparison table
- Reclassification tables
- Identify best/worst prognosis subtypes

---

#### **Analysis 3.2: Subtype Validation in TCGA-LAML**

**Objective:** Validate molecular subtypes in independent cohort

**Workflow:**

1. **Build Subtype Classifier:**
   - Training data: Beat AML subtypes
   - Features: Top 50-100 marker genes per cluster
   - Algorithm: Random Forest or SVM
   - Cross-validation: 5-fold CV within Beat AML (80/20 split)
   - Performance metrics: Accuracy, balanced accuracy, F1-score

2. **External Validation:**
   - Download TCGA-LAML expression data (n=151)
   - Normalize to match Beat AML (log2(TPM+1))
   - Apply trained classifier to TCGA samples
   - Assign TCGA samples to Beat AML subtypes

3. **Validation Checks:**
   - **Survival:** Do TCGA samples classified to same subtype show similar survival?
   - **Mutations:** Do TCGA samples show same mutation enrichments?
   - **Pathways:** Do TCGA samples show same pathway profiles?
   - **Clinical:** Do TCGA samples show same age/sex/FAB distributions?

4. **Reproducibility Assessment:**
   - Concordance between Beat AML and TCGA subtypes
   - Which subtypes replicate well? Which don't?
   - Robustness across cohorts

**Deliverables:**
- Subtype classifier model (saved for future use)
- TCGA validation results (survival curves, mutation enrichments)
- Reproducibility report
- Identification of robust vs cohort-specific subtypes

**Timeline:** 1 week

---

### Phase 4: Drug Response Integration (Weeks 6-7)

#### **Critical Consideration:**
**Sample size:** n=478 with all 4 data types (expression + mutations + clinical + drugs)
- This is adequate for hypothesis generation
- Findings will need prospective validation
- Report exact sample sizes for all comparisons
- Focus on effect sizes, not just p-values

---

#### **Analysis 4.1: Drug Response Data Preparation**

**Workflow:**

1. **Data Extraction:**
   - Load drug AUC/IC50 values for all 122 compounds
   - Identify samples with ≥50% of compounds tested
   - Create drug response matrix (samples × drugs)
   - Handle missing data: Complete case analysis per drug

2. **Drug Annotation:**
   - Drug families (FLT3i, BCL2i, IDHi, etc.)
   - Mechanisms of action
   - Clinical approval status (FDA-approved vs experimental)
   - **Prioritize clinically relevant drugs:**
     - FLT3 inhibitors: midostaurin, gilteritinib, quizartinib
     - BCL2 inhibitors: venetoclax
     - IDH inhibitors: ivosidenib, enasidenib
     - MEK inhibitors
     - mTOR inhibitors
     - HDAC inhibitors

3. **Response Classification:**
   - **Continuous:** Use AUC values as-is
   - **Binary:** Define responders vs non-responders
     - Method: Median split or clustering-based
     - Validate against known associations (FLT3-ITD + FLT3i)

4. **Quality Control:**
   - Check for outlier responses
   - Verify known associations:
     - FLT3-ITD mutations → sensitive to FLT3 inhibitors
     - IDH1/2 mutations → sensitive to IDH inhibitors
   - If expected associations missing → data quality issue

**Deliverables:**
- Drug response matrix (478 samples × 122 drugs)
- Drug annotation table with families/mechanisms
- List of prioritized drugs for analysis
- QC report validating known associations

---

#### **Analysis 4.2: Subtype-Drug Response Association**

**Objective:** Identify which drugs work best for which molecular subtypes

**Workflow:**

1. **Differential Drug Response:**
   - For each drug: Do subtypes differ in sensitivity?
   - **Statistical tests:**
     - Kruskal-Wallis test (overall difference)
     - Wilcoxon rank-sum (pairwise comparisons)
     - Linear model: AUC ~ subtype + age + sex
   - **Multiple testing correction:** FDR < 0.10 (more lenient given sample size)

2. **Drug Family Analysis:**
   - Aggregate drugs by mechanism/family
   - Test: Do certain subtypes respond better to drug classes?
   - **Hypotheses to test:**
     - "Immune-enriched" → better response to immunomodulators?
     - "DNA repair-high" → better response to DNA-damaging agents?
     - "Metabolically active" → better response to metabolic inhibitors?

3. **Heatmap Visualization:**
   - Rows = drugs (or drug families)
   - Columns = subtypes
   - Color = median AUC or z-score
   - Annotate with known mechanisms
   - Hierarchical clustering to identify drug-subtype patterns

4. **Case-Control Analysis for Key Drugs:**
   - **For each high-impact drug:**
     - Define responders (low AUC) vs non-responders (high AUC)
     - Compare subtype distribution
     - Compare mutation profiles
     - Compare pathway scores

   - **Focus drugs:**
     - Venetoclax (BCL2i): Which subtype most sensitive?
     - FLT3 inhibitors: Beyond FLT3-ITD, are there expression predictors?
     - IDH inhibitors: Beyond IDH mutations, expression biomarkers?

5. **Mechanistic Interpretation:**
   - For significant drug-subtype associations:
     - Check expression of drug targets in subtypes
     - Check pathway dependencies
     - Generate hypotheses for mechanism

**Statistical Considerations:**
- With n~50-100 per subtype, power is limited
- Report effect sizes (Cohen's d, median difference)
- Use FDR < 0.10 (exploratory threshold)
- Emphasize hypothesis generation over definitive conclusions
- Validate top findings with literature or cell line data

**Deliverables:**
- Subtype-drug response association table (122 drugs × k subtypes)
- Heatmap of drug sensitivities by subtype
- Prioritized drug-subtype combinations for follow-up
- Case-control analysis for top 10 drugs
- Mechanistic interpretation for significant associations

---

#### **Analysis 4.3: Predictive Modeling**

**Objective:** Build multi-omics models to predict drug response

**Approach:** Start simple, add complexity incrementally

**Workflow:**

1. **Feature Engineering:**
   - **Mutation features:** Binary matrix (has mutation yes/no)
   - **Expression features:**
     - Option A: Top 1000 most variable genes
     - Option B: Pathway enrichment scores (preferred - more stable)
   - **Clinical features:** Age, sex, WBC count, blast %
   - **Interaction terms:** Mutation × expression (e.g., NPM1-mutant × gene X)

2. **Model Training (per drug):**
   - **Target:** Drug response (AUC or binary responder)
   - **Algorithms:**
     - Baseline: Mutations-only (logistic/linear regression)
     - Elastic Net (L1+L2 regularization)
     - Random Forest (non-linear, feature importance)
     - XGBoost (gradient boosting)

   - **Cross-validation:** 5-fold CV (due to sample size ~100-400 per drug)
   - **Performance metrics:**
     - Continuous: R², RMSE, Spearman correlation
     - Binary: AUC-ROC, balanced accuracy, F1-score

3. **Feature Importance:**
   - Identify top predictors for each drug
   - Are they transcriptomic, mutational, or clinical?
   - Consistent predictors across drugs of same family?

4. **Model Comparison:**
   - **Baseline:** Mutations-only model
   - **Model 1:** Mutations + clinical
   - **Model 2:** Mutations + expression
   - **Model 3:** Mutations + expression + pathway scores
   - **Best:** Which features contribute most?

5. **Benchmark Against Published Models:**
   - Compare to MDREAM (npj Precision Oncology 2023)
   - Compare to NetAML (Advanced Science 2025)
   - Similar performance? Better? Novel biomarkers?

**Reality Check:**
- With n~100-400 samples per drug, expect R² ~0.3-0.5 (modest)
- Goal: Identify robust biomarkers, not perfect predictions
- Some drugs will be more predictable than others
- Compare predictability to published studies

**Deliverables:**
- Predictive models for top 20 drugs (clinically relevant)
- Feature importance rankings (top 20 features per drug)
- Model performance comparison table
- Identified biomarkers requiring validation
- Scatter plots: predicted vs observed response

**Timeline:** 2 weeks

---

### Phase 5: Integrative Analysis & Interpretation (Weeks 8-9)

#### **Analysis 5.1: Multi-Omics Integration Figure**

**Objective:** Create comprehensive "Figure 1" for manuscript

**Components:**
1. **Panel A:** Consensus clustering heatmap
   - Top variable genes × samples
   - Samples ordered by cluster
   - Color bars: cluster assignment, ELN risk, age, survival status

2. **Panel B:** Mutation oncoplot aligned by cluster
   - Top 20 mutated genes × samples
   - Samples ordered same as Panel A
   - Mutation enrichment bars on left

3. **Panel C:** Clinical annotations bar plot
   - Age, sex, ELN risk, FAB subtype, treatment
   - Samples ordered by cluster

4. **Panel D:** Survival curves by cluster
   - Kaplan-Meier with p-value
   - Color-coded by cluster

5. **Panel E:** Pathway enrichment heatmap
   - Key pathways × clusters
   - Z-scores of pathway enrichment

6. **Panel F:** Drug response profiles
   - Selected drugs × clusters
   - Median AUC with significance stars

**Tools:** ComplexHeatmap R package, ggplot2, survminer

**Deliverable:** Multi-panel integrative figure (publication-quality)

---

#### **Analysis 5.2: Clinical Utility Assessment**

**Objective:** Determine if your subtypes improve upon ELN risk stratification

**Workflow:**

1. **Reclassification Analysis:**
   - Within each ELN category, do subtypes further stratify survival?
   - Calculate C-index for:
     - ELN risk alone
     - Molecular subtypes alone
     - ELN risk + subtypes combined
   - Test improvement: C-index difference with 95% CI

2. **Net Reclassification Improvement (NRI):**
   - Do subtypes correctly reclassify patients into better risk categories?
   - Event NRI and non-event NRI

3. **Decision Curve Analysis:**
   - Clinical utility at different threshold probabilities
   - Does your classification change treatment decisions?
   - Net benefit compared to ELN risk

4. **Risk Score Development:**
   - Combine subtype + mutations + clinical into single prognostic score
   - Train Cox model with feature selection
   - Validate in TCGA cohort

**Deliverable:**
- Reclassification tables
- C-index comparison (with 95% CI)
- Decision curve plots
- Proposed risk score algorithm

---

#### **Analysis 5.3: Literature Contextualization**

**Objective:** Position findings within current AML research

**Tasks:**

1. **Compare to Published Subtypes:**
   - How do your clusters relate to:
     - Cheng et al. 2022: 8 transcriptomic subgroups
     - Song et al. 2025: UAMOCS 3 subtypes
     - TCGA-LAML clusters
   - Venn diagram of subtype overlaps
   - Concordance analysis

2. **Benchmark Drug Response Models:**
   - Compare your models to MDREAM, NetAML
   - Similar performance? Novel findings?

3. **Identify Novel Insights:**
   - What's NEW in your analysis that hasn't been reported?
   - Which subtype-drug associations are unreported?
   - Any unexpected biological findings?

**Deliverable:**
- Literature comparison table
- Venn diagrams
- List of novel findings (for discussion section)

**Timeline:** 2 weeks

---

### Phase 6: Manuscript Preparation (Weeks 10-12)

#### **Proposed Title:**
> "Multi-Omics Molecular Subtyping Reveals Clinically Actionable Drug Response Profiles in Acute Myeloid Leukemia"

#### **Proposed Structure:**

**Abstract (250 words):**
- Background: AML heterogeneity, limitations of genetic-only classification
- Methods: Multi-omics integration (n=478 gold standard, n=707 expression, n=942 survival)
- Results: X distinct molecular subtypes, survival stratification, drug response profiles
- Conclusions: Clinical utility, precision medicine implications

**Introduction:**
- AML heterogeneity and unmet clinical need
- Limitations of ELN risk classification
- Need for integrated multi-omics approach
- Study objectives

**Results:**
1. **Molecular Subtyping:**
   - X distinct transcriptomic subtypes identified
   - Biological characterization (mutations, pathways)
   - Comparison to published classifications

2. **Clinical Associations:**
   - Survival stratification beyond ELN risk
   - Age/treatment associations
   - Validation in TCGA-LAML

3. **Drug Response Profiles:**
   - Subtype-specific drug sensitivities
   - Predictive biomarkers for targeted therapies
   - Mechanistic interpretation

4. **Clinical Utility:**
   - Improvement over ELN risk (C-index, NRI)
   - Decision curve analysis
   - Proposed risk stratification algorithm

**Discussion:**
- Summary of key findings
- Comparison to recent literature (2024-2025)
- Clinical implications for treatment selection
- Limitations and future directions
- Conclusion

**Methods:**
- Data sources (Beat AML, TCGA)
- Quality control and batch correction
- Clustering methodology
- Statistical analyses
- Software and reproducibility

**Figures:**
- Figure 1: Multi-omics integration (6 panels)
- Figure 2: Survival analysis by subtype
- Figure 3: Drug response profiles
- Figure 4: Clinical utility (C-index, decision curves)
- Supplementary: 10-15 figures

**Tables:**
- Table 1: Patient characteristics (overall and by subtype)
- Table 2: Subtype characterization (mutations, pathways)
- Table 3: Cox regression results
- Supplementary: 5-10 tables

---

#### **Target Journals (in order):**

**Tier 1 (Aim High):**
- Nature Medicine (IF: 58.7)
- Blood (IF: 25.5)
- Leukemia (IF: 12.8)

**Tier 2 (Realistic):**
- Clinical Cancer Research (IF: 13.8)
- npj Precision Oncology (IF: 6.8)
- Molecular Cancer (IF: 27.7)

**Tier 3 (Backup):**
- Leukemia Research (IF: 3.0)
- Cancers (IF: 5.2)

**Timeline:** 3 weeks for complete draft

---

## III. STATISTICAL & METHODOLOGICAL RIGOR

### **Critical Success Factors:**

1. **Always Report Exact Sample Sizes**
   - Don't say "~100 samples" → say "103 samples"
   - Report overlap: "478 samples have all 4 data types"

2. **Multiple Testing Correction**
   - FDR < 0.05 for expression/pathway analyses
   - FDR < 0.10 for drug response (exploratory)
   - Bonferroni for focused hypothesis tests

3. **Effect Sizes, Not Just P-values**
   - Hazard ratios with 95% CI
   - Cohen's d for group comparisons
   - Median difference with IQR

4. **Confidence Intervals**
   - For all key estimates
   - Bootstrap for non-parametric metrics

5. **Biological Validation**
   - Check against known AML biology
   - Literature support for findings
   - Pathway interpretation

6. **Honest Limitations**
   - Limited drug response sample size
   - Cross-sectional design
   - Validation cohort differences

7. **Reproducibility**
   - Version control (Git)
   - Document all parameters
   - Provide processed data as supplements
   - Share code on GitHub

---

## IV. EXPECTED OUTCOMES & IMPACT

### **Primary Findings (Expected):**

1. **4-8 transcriptomic subtypes** with distinct biology
2. **Survival stratification** beyond ELN risk (HR range: 0.3-3.0)
3. **Subtype-specific drug sensitivities** for 10-20 compounds
4. **Predictive biomarkers** for FLT3i, BCL2i, IDH inhibitors
5. **Improved C-index:** ELN risk (0.65) → Combined (0.72-0.75)

### **Novel Contributions:**

1. **First comprehensive transcriptomic subtyping** of full Beat AML cohort
2. **Integration of subtypes with extensive drug response** data (n=478, 122 drugs)
3. **Pathway-based interpretation** of drug sensitivity mechanisms
4. **Clinical utility assessment** vs standard-of-care ELN risk
5. **Validation in TCGA-LAML** for reproducibility

### **Clinical Impact:**

- **Improved patient stratification** for clinical trial enrollment
- **Hypothesis-generating** for rational drug combinations
- **Foundation for prospective biomarker** validation studies
- **Contribution to precision medicine** in AML

### **Publication Potential:**

**Estimated Impact Factor:** 10-25 (Tier 1-2 journals)
**Expected Citations:** 50-100 in first 2 years
**Clinical Relevance:** High (addresses unmet need in AML)

---

## V. RISK MITIGATION & CONTINGENCY PLANS

| **Challenge** | **Solution** |
|---------------|--------------|
| Limited drug response samples | Focus on hypothesis generation; emphasize effect sizes |
| Batch effects in validation data | Rigorous normalization; batch correction |
| Subtypes don't validate in TCGA | Test biological consistency not just labels |
| Multiple testing burden | Pre-specify primary hypotheses; FDR correction |
| Clinical utility unclear | Compare to ELN risk; propose combined score |
| Weak drug response predictability | Focus on clinically relevant drugs; compare to benchmarks |

---

## VI. IMMEDIATE NEXT STEPS

### **This Week (Week 3):**

1. ✅ **Batch correction:** COMPLETE
2. **Review outlier samples:**
   - Load `03_Results/02_QC_Reports/expression_outliers.csv`
   - Investigate 7 outlier samples
   - Decision: Keep or exclude? (Document in Analysis_Decisions.md)

3. **Prepare gold standard cohort (n=478):**
   - Load master_sample_id_mapping.csv
   - Filter for `cohort_category == 'complete_quad_omics'`
   - Extract expression, mutations, clinical, drug data
   - Verify sample alignment across data types

4. **Begin Phase 2.1 - Consensus Clustering:**
   - Use batch-corrected expression: `beataml_expression_batchcorrected.rds`
   - Select top 5,000 variable genes
   - Run ConsensusClusterPlus (k=2-10)
   - Generate consensus matrix heatmaps

### **Resources Needed:**

**Computational:**
- R (v4.3+): ConsensusClusterPlus, GSVA, limma, survival, ComplexHeatmap
- Python (v3.9+): scikit-learn, pandas, seaborn
- RAM: 32GB+ for large matrix operations
- Storage: 500GB for intermediate files

**Biological Knowledge:**
- AML biology literature (2024-2025)
- Pathway databases: MSigDB, KEGG, Reactome
- Published AML gene signatures (LSC17, immune)

---

## VII. SUCCESS METRICS

### **Phase 2 (Molecular Subtyping) Success Criteria:**
✅ Identify 4-8 robust clusters with silhouette score > 0.3
✅ Clusters show significant mutation enrichment (p < 0.05)
✅ Pathway enrichment makes biological sense
✅ Clusters correlate with survival (log-rank p < 0.05)

### **Phase 3 (Survival) Success Criteria:**
✅ Subtypes stratify survival beyond ELN risk (multivariate p < 0.05)
✅ C-index improvement > 0.05 over ELN alone
✅ At least partial validation in TCGA-LAML

### **Phase 4 (Drug Response) Success Criteria:**
✅ Identify ≥5 significant subtype-drug associations (FDR < 0.10)
✅ Validate known associations (FLT3-ITD + FLT3i)
✅ Predictive models outperform mutation-only baseline

### **Phase 6 (Manuscript) Success Criteria:**
✅ Complete draft ready for submission
✅ All figures publication-quality
✅ Methods fully documented and reproducible
✅ Statistical review complete

---

## VIII. CONCLUSION

This research plan is designed to be:
- **Scientifically rigorous:** Follows 2024-2025 best practices
- **Clinically relevant:** Addresses AML precision medicine gaps
- **Realistic:** Acknowledges limitations, works within constraints
- **Publishable:** Structured for high-impact journal (IF: 10-25)
- **Reproducible:** Code, data, methods fully documented

**Current Status:** ✅ Phase 1 COMPLETE | Ready for Phase 2 Molecular Subtyping

**Timeline to First Results:** 2-3 weeks (molecular subtypes)
**Timeline to Manuscript Submission:** 12 weeks

**Science is iterative.** Be prepared to adjust based on what the data reveals.

---

**Document prepared:** 2025-10-04
**Next milestone:** Week 3 - Consensus clustering complete
**Project lead:** AML Multi-Omics Research Team
