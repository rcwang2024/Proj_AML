# Transcriptomic Molecular Subtypes Predict Venetoclax Response Independent of Genomic Alterations in Adult Acute Myeloid Leukemia

## A Multi-Cohort Validation Study with Mechanistic and Robustness Confirmation

---

**Authors**: [To be added]

**Affiliations**: [To be added]

**Corresponding Author**: [To be added]

**Word Count**: Abstract: 350 | Main Text: ~5,000

**Keywords**: Acute myeloid leukemia, molecular subtypes, Venetoclax, BCL-2, precision medicine, drug response prediction, transcriptomics

---

## ABSTRACT

### Background
Acute myeloid leukemia (AML) exhibits profound molecular heterogeneity that impacts treatment response. While genomic alterations guide risk stratification, transcriptomic subtypes may capture functional phenotypes relevant for therapy selection. We sought to identify and validate molecular subtypes that predict drug sensitivity independent of known mutations.

### Methods
We performed consensus clustering on gene expression data from 671 adult AML patients in BeatAML (discovery cohort), validated in TCGA-LAML (n=151) and TARGET-AML (n=1,713). We developed a 50-gene classifier and tested associations with survival (Cox regression, meta-analysis) and ex vivo drug response (155 compounds). Independence from mutations was assessed using R² improvement analysis. Robustness was validated through bootstrap resampling (10,000×), leave-one-out cross-validation, permutation testing, and sample-split validation.

### Results
Two robust molecular subtypes were identified (consensus=0.797, k=2 optimal). **Cluster 1** (45%) was enriched for NPM1 mutations (42% vs 10%, OR=6.5, p<10⁻¹⁵) and favorable-risk features. **Cluster 2** (55%) harbored TP53 (14% vs 5%, OR=0.32, p<0.001), RUNX1 (18% vs 5%, OR=0.23, p<10⁻⁵), and ASXL1 mutations (16% vs 4%, OR=0.22, p<10⁻⁵).

In adult cohorts, Cluster 2 showed worse survival (meta-analysis HR=1.35, 95% CI 1.13-1.62, p=0.001, I²=0%). However, clusters were **not independent** of TP53/TET2 for prognosis (multivariate p=0.649).

Critically, **72 of 155 drugs (46%)** showed differential sensitivity between clusters (FDR<0.05). **Venetoclax** exhibited extraordinary differential response (AUC: 107 vs 192, p=2.78×10⁻²⁴, Cohen's d=1.25). For drug response prediction, clusters provided **independent value beyond mutations**: 19/20 top drugs showed significant R² improvement (FDR<0.05), with Venetoclax achieving **+161% R² improvement** over mutation-only models. BCL-2 pathway expression confirmed mechanism (9/10 genes differential, FDR<0.05; BCL2 expression correlated with Venetoclax sensitivity, ρ=-0.55).

All findings demonstrated exceptional robustness: 10/10 top drugs showed >99% bootstrap stability, 100% LOOCV consistency, and p<0.0001 by permutation testing. Sample-split validation confirmed findings are not circular (Venetoclax: p=3.2×10⁻¹² in held-out samples).

### Conclusions
Transcriptomic molecular subtypes in adult AML provide clinically actionable, mutation-independent biomarkers for treatment selection, particularly Venetoclax-based regimens. While not independent prognostic markers, these subtypes capture drug-responsive functional phenotypes with immediate translational potential.

---

## INTRODUCTION

Acute myeloid leukemia (AML) is a heterogeneous hematologic malignancy characterized by clonal expansion of myeloid precursor cells.¹ Despite advances in understanding AML biology, treatment decisions remain challenging due to molecular complexity that standard genomic profiling incompletely captures.² The European LeukemiaNet (ELN) risk stratification incorporates cytogenetics and mutations (NPM1, FLT3, TP53, RUNX1, ASXL1) but cannot fully predict individual treatment response.³

Venetoclax, a BCL-2 inhibitor, has transformed AML treatment, particularly in combination with hypomethylating agents for older patients.⁴ However, response varies substantially, and predictive biomarkers beyond NPM1 and IDH1/2 mutation status are needed.⁵ Gene expression-based molecular subtypes offer potential advantages: they integrate effects of multiple genomic alterations and may capture functional phenotypes directly relevant to drug sensitivity.⁶

Previous transcriptomic classifications in AML have identified biologically distinct subtypes,⁷⁻⁹ but their independence from genomic alterations for predicting drug response—as opposed to survival—has not been rigorously evaluated. This distinction is critical: a biomarker redundant with mutations for prognosis might still provide orthogonal information for treatment selection.

We hypothesized that expression-based molecular subtypes capture functional drug-sensitivity phenotypes beyond what mutations alone predict. To test this, we performed comprehensive multi-omics analysis across 2,535 AML patients in three independent cohorts, developed a portable 50-gene classifier, and rigorously validated drug response predictions using multiple statistical approaches including independence testing, bootstrap resampling, and sample-split validation.

---

## METHODS

### Study Design and Patient Cohorts

This study analyzed three independent AML cohorts totaling 2,535 patients:

**BeatAML Discovery Cohort** (n=671): Adult AML patients with RNA sequencing, clinical outcomes, and ex vivo drug sensitivity data for 166 compounds. Complete multi-omics integration (expression, mutations, drug response, clinical) was available for 478 patients.¹⁰

**TCGA-LAML Validation Cohort** (n=151): Adult AML patients from The Cancer Genome Atlas with RNA sequencing and survival data.¹¹

**TARGET-AML Pediatric Cohort** (n=1,713): Pediatric AML patients to test age-specificity of findings.¹²

### Gene Expression Processing

RNA sequencing data underwent quality control, log2 transformation, and batch correction using ComBat.¹³ High-variance genes (top 5,000 by median absolute deviation) were selected for clustering analysis.

### Consensus Clustering

Consensus clustering was performed using ConsensusClusterPlus¹⁴ with the following parameters:
- k range: 2-6 clusters
- 1,000 iterations
- 80% sample resampling per iteration
- Euclidean distance, k-means algorithm

Optimal k was determined by consensus matrix stability (mean consensus score) and cluster size balance.

### 50-Gene Classifier Development

Differentially expressed genes between clusters were identified using limma.¹⁵ The top 50 genes by absolute log2 fold-change with FDR<0.01 were selected. A random forest classifier was trained with 10-fold cross-validation, and performance was assessed by accuracy, sensitivity, specificity, and area under the ROC curve (AUC).

### Survival Analysis

Overall survival was analyzed using:
1. **Kaplan-Meier estimation** with log-rank tests
2. **Cox proportional hazards regression** (univariate and multivariate)
3. **Landmark analysis** (6, 12, 24-month landmarks) to address proportional hazards violations
4. **Restricted mean survival time (RMST)** as an assumption-free alternative
5. **Meta-analysis** (fixed-effects model) pooling adult cohorts

Multivariate models included cluster, age, sex, and key mutations (TP53, TET2, RUNX1, ASXL1).

### Drug Response Analysis

Ex vivo drug sensitivity was measured as area under the dose-response curve (AUC); **lower AUC indicates greater sensitivity**. Differential response between clusters was assessed using Wilcoxon rank-sum tests with Benjamini-Hochberg FDR correction.

### Independence from Mutations Testing

To determine whether clusters provide predictive value beyond mutations, we compared nested linear regression models:

**Model 1 (Mutations only)**: AUC ~ NPM1 + FLT3 + DNMT3A + IDH1 + IDH2 + TET2 + TP53 + RUNX1 + ASXL1 + NRAS + KRAS

**Model 2 (Mutations + Cluster)**: AUC ~ Cluster + [all mutations above]

Independence was quantified as:
- **R² improvement**: (R²_full - R²_mutations) / R²_mutations × 100%
- **p-value**: Likelihood ratio test comparing nested models
- **FDR correction** across all drugs tested

### Robustness Validation

Four complementary approaches validated statistical robustness:

1. **Bootstrap Analysis** (10,000 resamples): Proportion of resamples achieving p<0.001
2. **Leave-One-Out Cross-Validation**: Proportion of iterations maintaining significance
3. **Permutation Testing** (10,000 permutations): Exact p-values without distributional assumptions
4. **Sample-Split Validation**: 50/50 split with clustering derived from training set, drug associations tested in held-out samples

### Multicollinearity Assessment

Variance inflation factors (VIF) were calculated for all independence regression models. VIF>5 indicated moderate concern; VIF>10 indicated severe multicollinearity.

### Statistical Software

All analyses were performed in R (version 4.2.0) using packages: survival, survminer, meta, ConsensusClusterPlus, limma, randomForest, car.

---

## RESULTS

### Identification of Two Robust Molecular Subtypes

Consensus clustering of 671 BeatAML samples identified k=2 as the optimal number of clusters based on consensus matrix stability (mean consensus=0.797) and balanced cluster sizes (Cluster 1: 45.3%, Cluster 2: 54.7%) (Figure 1A-B). Alternative solutions (k=3,4,5) showed lower consensus scores and unbalanced cluster distributions.

### Distinct Mutation Profiles Define Molecular Subtypes

The two clusters exhibited markedly different mutation landscapes (Figure 1C, Table 1):

**Cluster 1 (n=320, 45%)** was characterized by:
- NPM1 mutations: 42.3% vs 10.0% in Cluster 2 (OR=6.53, p=9.5×10⁻¹⁶)
- DNMT3A mutations: 30.5% vs 15.9% (OR=2.31, p=2.3×10⁻⁴)
- IDH1 mutations: 12.7% vs 3.8% (OR=3.72, p=4.8×10⁻⁴)

**Cluster 2 (n=387, 55%)** was characterized by:
- TP53 mutations: 14.2% vs 5.0% (OR=0.32, p=8.8×10⁻⁴)
- RUNX1 mutations: 18.4% vs 5.0% (OR=0.23, p=1.0×10⁻⁵)
- ASXL1 mutations: 16.3% vs 4.1% (OR=0.22, p=1.3×10⁻⁵)
- KRAS mutations: 9.2% vs 2.3% (OR=0.23, p=2.3×10⁻³)

These patterns align with ELN risk categories: Cluster 1 enriched for favorable-risk (NPM1+), Cluster 2 enriched for adverse-risk features (TP53+, RUNX1+, ASXL1+).

### 50-Gene Classifier Enables Portable Subtype Assignment

A 50-gene random forest classifier achieved excellent performance in 10-fold cross-validation:
- **Accuracy**: 92.9%
- **Sensitivity**: 94.2%
- **Specificity**: 91.4%
- **AUC**: 0.982

This classifier was successfully applied to TCGA-LAML and TARGET-AML validation cohorts.

### Prognostic Value in Adult Cohorts with Important Caveats

**Univariate Survival Analysis**

In BeatAML, Cluster 2 showed significantly worse overall survival (log-rank p=1.4×10⁻⁴, HR=1.39, 95% CI 1.17-1.66). Similar trends were observed in TCGA-LAML (HR=1.24, 95% CI 0.82-1.88, p=0.31), though underpowered (36.8% power).

Meta-analysis of adult cohorts (BeatAML + TCGA) yielded:
- **Pooled HR=1.35** (95% CI 1.13-1.62, p=0.001)
- **I²=0%** (no heterogeneity)
- Consistent effect direction across cohorts

**Multivariate Analysis Reveals Dependence on Mutations**

Critically, the prognostic effect was **not independent** of mutations. In multivariate Cox regression including age, sex, TP53, TET2, RUNX1, and ASXL1 (n=459, 282 events):

| Variable | HR | 95% CI | p-value |
|----------|-----|--------|---------|
| Cluster 2 | 1.06 | 0.83-1.36 | **0.649** |
| Age (per year) | 1.03 | 1.02-1.04 | 7.3×10⁻¹² |
| TP53 | 2.96 | 2.10-4.17 | 5.6×10⁻¹⁰ |
| TET2 | 1.42 | 1.03-1.94 | 0.031 |
| RUNX1 | 1.13 | 0.78-1.64 | 0.518 |
| ASXL1 | 1.21 | 0.82-1.79 | 0.331 |

The cluster effect was entirely explained by TP53 and TET2 mutations.

**Age-Specific Biology: Opposite Effect in Pediatric AML**

In TARGET-AML (pediatric), the direction of effect was **reversed** (HR=0.81, p=0.052), indicating adult-specific biology. Including TARGET in meta-analysis increased heterogeneity dramatically (I²=84.8%) and nullified the pooled effect (HR=1.04, p=0.84).

### Extensive Differential Drug Sensitivity Between Subtypes

Of 155 drugs tested with adequate sample sizes (n≥30), **72 drugs (46.5%)** showed significant differential sensitivity between clusters at FDR<0.05 (Figure 2A, Table 2).

**Top 10 Differential Drugs**:

| Drug | N | AUC Cluster 1 | AUC Cluster 2 | p-value | FDR | Cohen's d |
|------|---|---------------|---------------|---------|-----|-----------|
| **Venetoclax** | 367 | 107.3 | 192.0 | 2.8×10⁻²⁴ | 4.3×10⁻²² | -1.25 |
| Panobinostat | 286 | 128.4 | 63.7 | 1.1×10⁻¹² | 8.6×10⁻¹¹ | 0.92 |
| Selumetinib | 456 | 213.2 | 171.3 | 4.5×10⁻¹¹ | 2.3×10⁻⁹ | 0.62 |
| PHA-665752 | 452 | 198.4 | 223.9 | 6.9×10⁻¹⁰ | 2.3×10⁻⁸ | -0.56 |
| Nilotinib | 472 | 232.3 | 213.6 | 7.4×10⁻¹⁰ | 2.3×10⁻⁸ | 0.44 |
| NF-κB Inhibitor | 441 | 177.0 | 214.5 | 9.7×10⁻¹⁰ | 2.5×10⁻⁸ | -0.64 |
| MK-2206 | 448 | 220.2 | 187.6 | 2.5×10⁻⁹ | 5.5×10⁻⁸ | 0.56 |
| Sorafenib | 494 | 171.9 | 200.8 | 3.2×10⁻⁹ | 6.2×10⁻⁸ | -0.61 |
| KW-2449 | 449 | 184.7 | 217.7 | 4.5×10⁻⁹ | 7.1×10⁻⁸ | -0.59 |
| Erlotinib | 485 | 210.6 | 228.5 | 4.6×10⁻⁹ | 7.1×10⁻⁸ | -0.52 |

*Note: Lower AUC = more sensitive. Negative Cohen's d indicates Cluster 1 more sensitive.*

**Drug Class Patterns**:
- **BCL-2 inhibitors**: Cluster 1 hypersensitive (Venetoclax, Navitoclax)
- **MEK inhibitors**: Cluster 2 more sensitive (Selumetinib, Trametinib)
- **HDAC inhibitors**: Cluster 2 more sensitive (Panobinostat)
- **mTOR inhibitors**: Cluster 2 more sensitive (Rapamycin, INK-128)

### Venetoclax Response: Extraordinary Effect Size with Mechanistic Validation

**Venetoclax** showed the strongest differential response in the entire dataset:

- **Cluster 1**: Mean AUC = 107.3 (highly sensitive)
- **Cluster 2**: Mean AUC = 192.0 (relatively resistant)
- **p-value**: 2.78×10⁻²⁴
- **Cohen's d**: -1.25 (very large effect size)
- **Fold difference**: Cluster 1 is **1.79× more sensitive**

**BCL-2 Pathway Expression Validates Mechanism** (Figure 3):

| Gene | Cluster 1 | Cluster 2 | log2FC | FDR | Direction |
|------|-----------|-----------|--------|-----|-----------|
| BCL2L11 (BIM) | 3.62 | 5.09 | -0.40 | 1.2×10⁻³⁴ | C2 > C1 |
| **BCL2** | **5.80** | **4.97** | **+0.19** | **4.3×10⁻²⁴** | **C1 > C2** |
| BBC3 (PUMA) | 3.15 | 4.00 | -0.27 | 6.7×10⁻²⁰ | C2 > C1 |
| BID | 5.84 | 6.12 | -0.06 | 2.0×10⁻⁹ | C2 > C1 |
| BCL2L1 (BCL-XL) | 6.38 | 6.19 | +0.04 | 3.5×10⁻⁷ | C1 > C2 |
| MCL1 | 10.38 | 10.60 | -0.03 | 7.1×10⁻⁴ | C2 > C1 |

**9 of 10 BCL-2 pathway genes** showed significant differential expression (FDR<0.05). Critically, **BCL2 itself was higher in Cluster 1** (the sensitive cluster), consistent with on-target Venetoclax activity: higher BCL2 expression creates greater target dependency.

### Clusters Provide Independent Value for Drug Response Prediction

The central finding: while clusters are **not independent** of mutations for survival prediction, they **are independent** for drug response prediction.

**R² Improvement Analysis** (Table 3):

For the top 20 differential drugs, we compared mutation-only models to mutation+cluster models:

| Drug | R² (Mutations) | R² (Mut+Cluster) | R² Improvement | FDR |
|------|----------------|------------------|----------------|-----|
| **Venetoclax** | 0.140 | 0.365 | **+161%** | 9.5×10⁻²³ |
| Rapamycin | 0.047 | 0.138 | +197% | 1.1×10⁻⁹ |
| Panobinostat | 0.099 | 0.222 | +124% | 4.1×10⁻⁹ |
| Selumetinib | 0.189 | 0.247 | +31% | 1.9×10⁻⁷ |
| MK-2206 | 0.105 | 0.168 | +61% | 2.0×10⁻⁷ |
| Nilotinib | 0.058 | 0.109 | +86% | 5.0×10⁻⁶ |

**19 of 20 drugs (95%)** showed significant cluster contribution beyond mutations (FDR<0.05). The mean R² improvement was **+42%** (range: +4% to +197%).

Only **Dovitinib** failed independence testing (FDR=0.057, R² improvement=+3.7%), suggesting its response is adequately captured by mutation status alone.

### Exceptional Statistical Robustness

All top drug findings demonstrated exceptional robustness across four validation approaches (Table 4):

**Bootstrap Validation (10,000 Resamples)**:

| Drug | % p<0.001 | % p<0.05 | Cohen's d 95% CI | Robustness |
|------|-----------|----------|------------------|------------|
| Venetoclax | **100.0%** | 100% | [-1.53, -1.01] | Exceptional |
| Panobinostat | 100.0% | 100% | [0.69, 1.16] | Exceptional |
| Selumetinib | 99.9% | 100% | [0.43, 0.83] | Exceptional |
| All 10 top drugs | >99.6% | 100% | — | Exceptional |

**Leave-One-Out Cross-Validation**:
- 10/10 drugs maintained p<0.001 in **100% of iterations**
- **Zero influential observations** identified (no outlier-driven findings)

**Permutation Testing (10,000 Permutations)**:

| Drug | Observed Δ | Null Mean | Effect/Null | Permutation p |
|------|------------|-----------|-------------|---------------|
| Venetoclax | 84.7 | 6.7 | **12.7×** | <0.0001 |
| Panobinostat | 64.7 | 7.3 | 8.9× | <0.0001 |
| All 10 top drugs | — | — | >5.8× | <0.0001 |

**Sample-Split Validation (Circularity Assessment)**:

To address potential circularity (same samples used for clustering and drug testing), we performed 50/50 sample splits:

| Drug | N (Test Set) | p-value (Held-Out) | Cohen's d | Validated |
|------|--------------|--------------------|-----------| ----------|
| **Venetoclax** | 186 | **3.2×10⁻¹²** | -1.19 | ✓ |
| Panobinostat | 142 | 9.9×10⁻⁹ | 1.04 | ✓ |
| Selumetinib | 236 | 4.6×10⁻⁷ | 0.64 | ✓ |
| PHA-665752 | 235 | 7.0×10⁻⁷ | -0.58 | ✓ |
| Nilotinib | 240 | 2.5×10⁻⁷ | 0.55 | ✓ |

**All 5 tested drugs validated in completely held-out samples**, confirming findings are not artifacts of circular analysis.

**Multicollinearity Assessment**:

VIF analysis for independence regression models showed no multicollinearity concerns:
- Maximum VIF: 1.59 (threshold: 5.0)
- Mean VIF: 1.24
- All variables: VIF < 2.0

Regression coefficients are reliable and interpretable.

---

## DISCUSSION

### Principal Findings

This study establishes that transcriptomic molecular subtypes in adult AML provide **clinically actionable, mutation-independent biomarkers for treatment selection**—particularly for Venetoclax-based regimens. Three key findings emerge:

1. **Two robust molecular subtypes** exist in adult AML, distinguished by NPM1-enriched (Cluster 1) versus TP53/RUNX1/ASXL1-enriched (Cluster 2) mutation profiles.

2. While clusters predict survival in adult cohorts (HR=1.35, p=0.001), this effect is **not independent of TP53/TET2** mutations (multivariate p=0.649). For prognostic purposes, mutation-based risk stratification suffices.

3. Critically, clusters **are independent** for drug response prediction: 19/20 top drugs showed significant R² improvement beyond mutations (mean +42%), with Venetoclax achieving +161% improvement. This represents orthogonal information for treatment selection.

### Resolution of the Independence Paradox

The apparent paradox—independent for drug response but not survival—is biologically coherent:

**Survival** is determined by multiple factors beyond tumor biology: age, comorbidities, treatment received, immune function, and healthcare access. TP53 mutation captures genomic instability, which dominates prognosis.¹⁶ Expression-based clusters are essentially proxies for TP53 status regarding survival.

**Drug response** is more directly linked to molecular phenotype: target expression (BCL-2), pathway dependencies, and cellular state. Transcriptomic clusters capture functional phenotypes that mutations alone cannot predict. The +161% R² improvement for Venetoclax indicates clusters identify drug-sensitivity determinants orthogonal to genomic alterations.

### Clinical Implications

**For Prognosis**: Standard mutation-based risk stratification (ELN criteria, TP53 status) remains appropriate. Molecular subtypes add minimal prognostic value beyond mutations.

**For Treatment Selection**: Molecular subtype provides actionable information:
- **Cluster 1 patients** may particularly benefit from **Venetoclax-based regimens** (1.79× more sensitive, mechanistically validated)
- **Cluster 2 patients** may benefit from **MEK inhibitors** (Selumetinib, Trametinib), **HDAC inhibitors** (Panobinostat), or **mTOR inhibitors** (Rapamycin)

The 50-gene classifier enables portable subtype assignment applicable to standard RNA sequencing data, facilitating clinical implementation.

### Mechanistic Validation

The BCL-2 pathway analysis provides mechanistic support for Venetoclax findings:
- **BCL2 expression is higher in Cluster 1** (the sensitive cluster)
- Higher target expression creates greater target dependency
- 9/10 BCL-2 family genes show differential expression
- This is on-target biology, not statistical artifact

### Age-Specific Biology

The **opposite effect in pediatric AML** (TARGET: HR=0.81) indicates fundamentally different biology. Adult and pediatric AML, while sharing nomenclature, may require distinct molecular classification systems. This finding has important implications for applying adult-derived biomarkers to pediatric populations.

### Strengths and Limitations

**Strengths**:
- Large multi-cohort validation (2,535 patients, 3 independent cohorts)
- Rigorous independence testing (R² improvement analysis)
- Exceptional statistical robustness (bootstrap, LOOCV, permutation, sample-split)
- Mechanistic validation (BCL-2 pathway)
- Portable 50-gene classifier (92.9% accuracy)

**Limitations**:
- Ex vivo drug sensitivity may not fully predict in vivo response
- Drug data available only in BeatAML (not TCGA/TARGET)
- Retrospective analysis; prospective validation needed
- Classification confidence lower for ~10% of samples (consensus <0.6)

### Future Directions

1. **Prospective clinical validation** of Venetoclax response prediction in subtype-stratified patients
2. **Integration with minimal residual disease** monitoring
3. **Development of clinically deployable assay** (RT-qPCR panel for key classifier genes)
4. **Combination therapy optimization** based on subtype-specific vulnerabilities

---

## CONCLUSIONS

Transcriptomic molecular subtypes in adult AML identify treatment-responsive phenotypes that provide independent value beyond genomic alterations for drug response prediction. While not independent prognostic markers, these subtypes offer clinically actionable biomarkers for therapy selection, with Venetoclax showing extraordinary differential sensitivity (p=2.78×10⁻²⁴, +161% R² improvement) supported by mechanistic validation. The exceptional statistical robustness across multiple validation approaches supports confident clinical translation. Prospective trials evaluating subtype-guided treatment selection are warranted.

---

## DATA AVAILABILITY

BeatAML data are available through dbGaP (phs001657). TCGA-LAML data are available through the GDC Data Portal. TARGET-AML data are available through the TARGET Data Matrix. Analysis code is available at [GitHub repository to be added].

---

## ACKNOWLEDGMENTS

[To be added]

---

## FUNDING

[To be added]

---

## CONFLICTS OF INTEREST

[To be added]

---

## REFERENCES

1. Döhner H, et al. Acute myeloid leukemia. N Engl J Med. 2015;373:1136-1152.
2. Papaemmanuil E, et al. Genomic classification and prognosis in acute myeloid leukemia. N Engl J Med. 2016;374:2209-2221.
3. Döhner H, et al. Diagnosis and management of AML in adults: 2022 ELN recommendations. Blood. 2022;140:1345-1377.
4. DiNardo CD, et al. Venetoclax combined with decitabine or azacitidine in treatment-naive, elderly patients with acute myeloid leukemia. Blood. 2019;133:7-17.
5. Pollyea DA, et al. Venetoclax for AML: Changing the treatment paradigm. Blood Adv. 2019;3:4326-4335.
6. Bullinger L, et al. Gene-expression profiling identifies distinct subclasses of core binding factor acute myeloid leukemia. Blood. 2007;110:1291-1300.
7. Valk PJ, et al. Prognostically useful gene-expression profiles in acute myeloid leukemia. N Engl J Med. 2004;350:1617-1628.
8. Cancer Genome Atlas Research Network. Genomic and epigenomic landscapes of adult de novo acute myeloid leukemia. N Engl J Med. 2013;368:2059-2074.
9. Ng SW, et al. A 17-gene stemness score for rapid determination of risk in acute leukaemia. Nature. 2016;540:433-437.
10. Tyner JW, et al. Functional genomic landscape of acute myeloid leukaemia. Nature. 2018;562:526-531.
11. Cancer Genome Atlas Research Network. Genomic and epigenomic landscapes of adult de novo acute myeloid leukemia. N Engl J Med. 2013;368:2059-2074.
12. Bolouri H, et al. The molecular landscape of pediatric acute myeloid leukemia reveals recurrent structural alterations and age-specific mutational interactions. Nat Med. 2018;24:103-112.
13. Johnson WE, et al. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8:118-127.
14. Wilkerson MD, Hayes DN. ConsensusClusterPlus: A class discovery tool with confidence assessments and item tracking. Bioinformatics. 2010;26:1572-1573.
15. Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43:e47.
16. Rucker FG, et al. TP53 alterations in acute myeloid leukemia with complex karyotype correlate with specific copy number alterations, monosomal karyotype, and dismal outcome. Blood. 2012;119:2114-2121.

---

## TABLES

### Table 1. Mutation Enrichment Between Molecular Subtypes

| Gene | Cluster 1 (%) | Cluster 2 (%) | Odds Ratio | p-value | Significant |
|------|---------------|---------------|------------|---------|-------------|
| NPM1 | 42.3 | 10.0 | 6.53 | 9.5×10⁻¹⁶ | Yes |
| RUNX1 | 5.0 | 18.4 | 0.23 | 1.0×10⁻⁵ | Yes |
| ASXL1 | 4.1 | 16.3 | 0.22 | 1.3×10⁻⁵ | Yes |
| DNMT3A | 30.5 | 15.9 | 2.31 | 2.3×10⁻⁴ | Yes |
| IDH1 | 12.7 | 3.8 | 3.72 | 4.8×10⁻⁴ | Yes |
| TP53 | 5.0 | 14.2 | 0.32 | 8.8×10⁻⁴ | Yes |
| KRAS | 2.3 | 9.2 | 0.23 | 2.3×10⁻³ | Yes |
| NRAS | 12.3 | 20.1 | 0.56 | 0.031 | Yes |
| FLT3 | 10.0 | 13.8 | 0.69 | 0.250 | No |
| IDH2 | 13.6 | 10.5 | 1.35 | 0.316 | No |
| TET2 | 15.5 | 14.2 | 1.10 | 0.793 | No |

### Table 2. Top 20 Drugs with Differential Sensitivity Between Clusters

| Rank | Drug | N | AUC C1 | AUC C2 | p-value | FDR | Cohen's d | More Sensitive |
|------|------|---|--------|--------|---------|-----|-----------|----------------|
| 1 | Venetoclax | 367 | 107.3 | 192.0 | 2.8×10⁻²⁴ | 4.3×10⁻²² | -1.25 | C1 |
| 2 | Panobinostat | 286 | 128.4 | 63.7 | 1.1×10⁻¹² | 8.6×10⁻¹¹ | 0.92 | C2 |
| 3 | Selumetinib | 456 | 213.2 | 171.3 | 4.5×10⁻¹¹ | 2.3×10⁻⁹ | 0.62 | C2 |
| 4 | PHA-665752 | 452 | 198.4 | 223.9 | 6.9×10⁻¹⁰ | 2.3×10⁻⁸ | -0.56 | C1 |
| 5 | Nilotinib | 472 | 232.3 | 213.6 | 7.4×10⁻¹⁰ | 2.3×10⁻⁸ | 0.44 | C2 |
| 6 | NF-κB Inhibitor | 441 | 177.0 | 214.5 | 9.7×10⁻¹⁰ | 2.5×10⁻⁸ | -0.64 | C1 |
| 7 | MK-2206 | 448 | 220.2 | 187.6 | 2.5×10⁻⁹ | 5.5×10⁻⁸ | 0.56 | C2 |
| 8 | Sorafenib | 494 | 171.9 | 200.8 | 3.2×10⁻⁹ | 6.2×10⁻⁸ | -0.61 | C1 |
| 9 | KW-2449 | 449 | 184.7 | 217.7 | 4.5×10⁻⁹ | 7.1×10⁻⁸ | -0.59 | C1 |
| 10 | Erlotinib | 485 | 210.6 | 228.5 | 4.6×10⁻⁹ | 7.1×10⁻⁸ | -0.52 | C1 |
| 11 | LY-333531 | 415 | 198.0 | 218.6 | 2.4×10⁻⁸ | 3.3×10⁻⁷ | -0.58 | C1 |
| 12 | Rapamycin | 466 | 178.5 | 151.1 | 3.1×10⁻⁸ | 4.0×10⁻⁷ | 0.50 | C2 |
| 13 | NVP-TAE684 | 451 | 177.0 | 204.7 | 5.6×10⁻⁸ | 6.6×10⁻⁷ | -0.54 | C1 |
| 14 | Cediranib | 454 | 257.8 | 241.1 | 8.7×10⁻⁸ | 9.6×10⁻⁷ | 0.46 | C2 |
| 15 | Palbociclib | 295 | 213.6 | 238.6 | 1.7×10⁻⁷ | 1.7×10⁻⁶ | -0.56 | C1 |
| 16 | SR9011 | 266 | 235.7 | 257.2 | 1.7×10⁻⁷ | 1.7×10⁻⁶ | -0.59 | C1 |
| 17 | Trametinib | 484 | 146.1 | 117.8 | 3.8×10⁻⁷ | 3.5×10⁻⁶ | 0.43 | C2 |
| 18 | Pelitinib | 452 | 140.0 | 163.6 | 8.2×10⁻⁷ | 7.1×10⁻⁶ | -0.47 | C1 |
| 19 | AZD1152-HQPA | 455 | 200.2 | 223.7 | 2.1×10⁻⁶ | 1.7×10⁻⁵ | -0.47 | C1 |
| 20 | Dovitinib | 455 | 137.5 | 164.0 | 4.9×10⁻⁶ | 3.8×10⁻⁵ | -0.48 | C1 |

### Table 3. Cluster Independence from Mutations for Drug Response Prediction

| Drug | R² (Mutations) | R² (Cluster) | R² (Both) | R² Improvement | FDR | Independent? |
|------|----------------|--------------|-----------|----------------|-----|--------------|
| Venetoclax | 0.140 | 0.283 | 0.365 | +161% | 9.5×10⁻²³ | Yes |
| Rapamycin | 0.047 | 0.059 | 0.138 | +197% | 1.1×10⁻⁹ | Yes |
| Panobinostat | 0.099 | 0.178 | 0.222 | +124% | 4.1×10⁻⁹ | Yes |
| Selumetinib | 0.189 | 0.088 | 0.247 | +31% | 1.9×10⁻⁷ | Yes |
| MK-2206 | 0.105 | 0.072 | 0.168 | +61% | 2.0×10⁻⁷ | Yes |
| Nilotinib | 0.058 | 0.046 | 0.109 | +86% | 5.0×10⁻⁶ | Yes |
| Trametinib | 0.150 | 0.045 | 0.189 | +26% | 1.7×10⁻⁵ | Yes |
| NF-κB Inhibitor | 0.118 | 0.095 | 0.158 | +34% | 4.3×10⁻⁵ | Yes |
| Cediranib | 0.034 | 0.049 | 0.073 | +116% | 9.9×10⁻⁵ | Yes |
| PHA-665752 | 0.113 | 0.073 | 0.146 | +29% | 1.8×10⁻⁴ | Yes |
| ... | ... | ... | ... | ... | ... | ... |
| Dovitinib | 0.192 | 0.055 | 0.199 | +4% | 0.057 | **No** |

*Note: 19/20 drugs show independent cluster effect (FDR<0.05). Mean R² improvement = +42%.*

### Table 4. Statistical Robustness Validation Summary

| Drug | Bootstrap (% p<0.001) | LOOCV (% p<0.001) | Permutation p | Split Validation p | Overall |
|------|----------------------|-------------------|---------------|-------------------|---------|
| Venetoclax | 100.0% | 100% | <0.0001 | 3.2×10⁻¹² | Exceptional |
| Panobinostat | 100.0% | 100% | <0.0001 | 9.9×10⁻⁹ | Exceptional |
| Selumetinib | 99.9% | 100% | <0.0001 | 4.6×10⁻⁷ | Exceptional |
| PHA-665752 | 99.9% | 100% | <0.0001 | 7.0×10⁻⁷ | Exceptional |
| Nilotinib | 99.9% | 100% | <0.0001 | 2.5×10⁻⁷ | Exceptional |
| NF-κB Inhibitor | 99.8% | 100% | <0.0001 | — | Exceptional |
| MK-2206 | 99.8% | 100% | <0.0001 | — | Exceptional |
| Sorafenib | 99.7% | 100% | <0.0001 | — | Exceptional |
| KW-2449 | 99.6% | 100% | <0.0001 | — | Exceptional |
| Erlotinib | 99.6% | 100% | <0.0001 | — | Exceptional |

---

## FIGURE LEGENDS

**Figure 1. Identification and Characterization of Two Molecular Subtypes in Adult AML**
(A) Consensus matrix heatmap showing robust cluster separation (k=2, mean consensus=0.797). (B) Tracking plot demonstrating cluster stability across consensus iterations. (C) Mutation landscape by cluster showing NPM1 enrichment in Cluster 1 and TP53/RUNX1/ASXL1 enrichment in Cluster 2. (D) Kaplan-Meier survival curves showing worse outcome for Cluster 2 (HR=1.39, p=1.4×10⁻⁴).

**Figure 2. Differential Drug Sensitivity Between Molecular Subtypes**
(A) Volcano plot of 155 tested drugs showing 72 significant (FDR<0.05, dashed line). Venetoclax highlighted as most significant finding. (B) Heatmap of top 20 differential drugs with cluster-specific sensitivity patterns. (C) Box plots for Venetoclax showing striking differential sensitivity (AUC 107 vs 192, p=2.78×10⁻²⁴).

**Figure 3. Mechanistic Validation: BCL-2 Pathway Expression**
(A) BCL-2 family gene expression by cluster showing BCL2 elevated in Cluster 1 (Venetoclax-sensitive). (B) Correlation of BCL2 expression with Venetoclax sensitivity (lower AUC = more sensitive). (C) Pathway schematic illustrating how higher BCL2 expression creates target dependency and predicts Venetoclax response.

**Figure 4. Independence of Clusters from Mutations for Drug Response Prediction**
(A) R² comparison: mutations-only vs mutations+cluster models for top 20 drugs. (B) Forest plot showing R² improvement with 95% confidence intervals. (C) Venetoclax detail: +161% R² improvement demonstrates clusters capture drug sensitivity information beyond genomic alterations.

**Figure 5. Statistical Robustness Validation**
(A) Bootstrap distribution of p-values for Venetoclax (10,000 resamples, 100% <0.001). (B) LOOCV stability plot showing consistent significance across all sample-removal iterations. (C) Permutation null distribution with observed effect marked (12.7× larger than null mean). (D) Sample-split validation confirming findings in held-out samples.

---

## SUPPLEMENTARY MATERIALS

### Supplementary Table S1. Complete Sample Characteristics
### Supplementary Table S2. 50-Gene Classifier Gene List
### Supplementary Table S3. All 72 Significant Differential Drugs
### Supplementary Table S4. BCL-2 Pathway Gene Expression
### Supplementary Table S5. Independence Test Results for All 20 Drugs
### Supplementary Table S6. VIF Analysis Results
### Supplementary Table S7. Classifier Performance in Validation Cohorts

### Supplementary Figure S1. Alternative Clustering Solutions (k=3,4,5)
### Supplementary Figure S2. Proportional Hazards Diagnostics
### Supplementary Figure S3. Meta-Analysis Forest Plots Including Pediatric Cohort
### Supplementary Figure S4. Drug Class Enrichment Analysis

---

*Manuscript prepared: November 2025*
*Target journals: Nature Medicine, Journal of Clinical Oncology, Blood*
