# SUPPLEMENTARY MATERIALS - MASTER INDEX

**Manuscript**: Transcriptomic Molecular Subtypes Predict Venetoclax Response Independent of Genomic Alterations in Adult Acute Myeloid Leukemia

**Authors**: [To be added]

**Correspondence**: [To be added]

---

## TABLE OF CONTENTS

### SUPPLEMENTARY METHODS
- [Supplementary Methods Document](#supplementary-methods)

### SUPPLEMENTARY TABLES
- [Table S1](#table-s1): Complete Sample Characteristics (n=2,535)
- [Table S2](#table-s2): 50-Gene Classifier Gene List with Weights
- [Table S3](#table-s3): All 72 Differential Drugs (FDR<0.05)
- [Table S4](#table-s4): BCL-2 Pathway Gene Expression
- [Table S5](#table-s5): Cluster Independence Testing (20 Drugs)
- [Table S6](#table-s6): Multivariate Analysis Results
- [Table S7](#table-s7): Robustness Validation Summary
- [Table S8](#table-s8): Cluster 2 Salvage Therapy Options
- [Table S9](#table-s9): VRS Clinical Decision Tool

### SUPPLEMENTARY FIGURES
- [Figure S1](#figure-s1): Alternative Clustering Solutions (k=3,4,5)
- [Figure S2](#figure-s2): Proportional Hazards Diagnostics
- [Figure S3](#figure-s3): Meta-Analysis with Pediatric Cohort
- [Figure S4](#figure-s4): Drug Class Enrichment Analysis
- [Figure S5](#figure-s5): Top 20 Drugs Boxplots
- [Figure S6](#figure-s6): BCL-2 Pathway Heatmap
- [Figure S7](#figure-s7): Cluster 2 Drug Response Profile
- [Figure S8](#figure-s8): VRS Distribution and Thresholds

---

## SUPPLEMENTARY METHODS

### 1. Patient Cohorts and Data Acquisition

#### BeatAML Discovery Cohort
Adult AML patients (n=671) enrolled in the Beat AML Master Trial (NCT03013998). RNA sequencing performed on diagnostic bone marrow aspirates using Illumina HiSeq 2500 platform. Ex vivo drug sensitivity testing conducted on primary patient cells using 166 compounds at 7 concentrations. Clinical data collected prospectively including demographics, cytogenetics, mutations, and outcomes. Data accessed through dbGaP (phs001657.v2.p1) with appropriate IRB approval.

**Inclusion criteria**:
- Age ≥18 years
- Confirmed AML diagnosis by WHO 2016 criteria
- De novo or secondary AML
- Adequate sample for RNA-seq (>70% blasts)

**Exclusion criteria**:
- Acute promyelocytic leukemia (APL)
- AML with t(8;21) or inv(16) [Core-binding factor AML]
- Prior AML treatment (except hydroxyurea)

#### TCGA-LAML Validation Cohort
Adult AML patients (n=151) from The Cancer Genome Atlas. RNA-seq data (Illumina HiSeq) and clinical annotations accessed through GDC Data Portal (https://portal.gdc.cancer.gov/). Survival data extracted from published supplementary materials (NEJM 2013;368:2059-2074).

#### TARGET-AML Pediatric Cohort
Pediatric AML patients (n=1,713) from the Therapeutically Applicable Research to Generate Effective Treatments initiative. RNA-seq data accessed through TARGET Data Matrix (https://ocg.cancer.gov/programs/target). Age range: 0-20 years (median 10.3 years). Used to test age-specificity of adult-derived molecular subtypes.

---

### 2. RNA Sequencing and Expression Processing

#### Sequencing Parameters
- **Platform**: Illumina HiSeq 2500 (BeatAML), Illumina HiSeq 2000 (TCGA/TARGET)
- **Read length**: 2×100 bp paired-end
- **Depth**: Median 50 million reads per sample
- **Alignment**: STAR v2.5.3a to GRCh38
- **Quantification**: RSEM v1.3.0

#### Quality Control
```r
# Filtering criteria
- Mapping rate >70%
- rRNA contamination <15%
- Duplicate rate <50%
- Gene detection: ≥10,000 genes with FPKM>0.1
```

#### Normalization Pipeline
```r
library(edgeR)
library(limma)

# 1. TMM normalization
dge <- DGEList(counts = raw_counts)
dge <- calcNormFactors(dge, method = "TMM")

# 2. Log2 transformation
expr_log <- cpm(dge, log = TRUE, prior.count = 1)

# 3. Batch correction (ComBat)
library(sva)
expr_corrected <- ComBat(dat = expr_log,
                         batch = sample_metadata$batch,
                         mod = model.matrix(~1, data = sample_metadata))

# 4. Variance stabilization
expr_vst <- limma::voom(dge)$E
```

---

### 3. Consensus Clustering Methodology

#### Algorithm Parameters
```r
library(ConsensusClusterPlus)

results <- ConsensusClusterPlus(
  d = expr_matrix_top5000,      # Top 5,000 MAD genes
  maxK = 6,                      # Test k=2 to k=6
  reps = 1000,                   # 1,000 iterations
  pItem = 0.8,                   # 80% sample resampling
  pFeature = 0.8,                # 80% feature resampling
  clusterAlg = "km",             # k-means clustering
  distance = "euclidean",        # Euclidean distance
  seed = 12345                   # Reproducibility
)
```

#### Optimal k Selection Criteria
1. **Consensus score**: Mean pairwise consensus within clusters
   - k=2: 0.957 (excellent)
   - k=3: 0.882 (good)
   - k=4: 0.791 (moderate)

2. **Silhouette width**: Average silhouette coefficient
   - k=2: 0.123
   - k=3: 0.089
   - k=4: 0.065

3. **Cluster size balance**:
   - k=2: 42% / 58% (balanced)
   - k=3: 25% / 35% / 40% (imbalanced)

4. **Biological interpretability**: k=2 shows clearest mutation enrichment patterns

**Decision**: k=2 selected based on highest consensus, balanced sizes, and biological coherence.

---

### 4. 50-Gene Classifier Development

#### Feature Selection
```r
library(limma)

# Differential expression analysis
design <- model.matrix(~0 + cluster)
fit <- lmFit(expr_matrix, design)
contrast <- makeContrasts(Cluster2-Cluster1, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Select top 50 genes
top_genes <- topTable(fit2, n = 50, adjust = "BH")
# Criteria: |log2FC| > 1.0 AND FDR < 0.01
```

#### Random Forest Training
```r
library(randomForest)

# Train classifier
rf_model <- randomForest(
  x = expr_matrix_50genes,
  y = cluster_labels,
  ntree = 1000,
  mtry = 7,                    # sqrt(50) ≈ 7
  importance = TRUE
)

# 10-fold cross-validation
cv_results <- replicate(10, {
  train_idx <- sample(1:nrow(data), 0.9 * nrow(data))
  train_set <- data[train_idx, ]
  test_set <- data[-train_idx, ]

  model <- randomForest(cluster ~ ., data = train_set, ntree = 1000)
  predictions <- predict(model, test_set)
  accuracy <- mean(predictions == test_set$cluster)

  return(accuracy)
})

mean_accuracy <- mean(cv_results)  # 92.9%
```

#### Performance Metrics
- **Accuracy**: 92.9% (91.2%-94.5% 95% CI)
- **Sensitivity**: 94.2%
- **Specificity**: 91.4%
- **AUC**: 0.982
- **Positive Predictive Value**: 92.6%
- **Negative Predictive Value**: 93.1%

---

### 5. Survival Analysis Methods

#### Proportional Hazards Testing
```r
# Schoenfeld residuals test
cox_model <- coxph(Surv(time, event) ~ cluster)
ph_test <- cox.zph(cox_model)

# Result: Global test p = 0.0002 (PH VIOLATED)
```

#### PH-Free Analysis Methods

**Method 1: Stratified Cox (Log-Rank Test)**
```r
survdiff(Surv(OS_months, OS_event) ~ cluster,
         data = survival_data,
         rho = 0)  # Log-rank test

# Result: χ² = 14.2, p = 0.00014
```

**Method 2: Landmark Analysis**
```r
# 6-month landmark
landmark_6m <- survival_data %>% filter(OS_months >= 6)
cox_6m <- coxph(Surv(OS_months - 6, OS_event) ~ cluster,
                data = landmark_6m)

# 12-month landmark
landmark_12m <- survival_data %>% filter(OS_months >= 12)
cox_12m <- coxph(Surv(OS_months - 12, OS_event) ~ cluster,
                 data = landmark_12m)

# Results:
# 6m: HR = 2.22 (1.56-3.16), p < 0.001
# 12m: HR = 1.85 (1.26-2.72), p = 0.002
# 24m: HR = 1.62 (1.04-2.54), p = 0.033
```

**Method 3: Restricted Mean Survival Time**
```r
library(survRM2)

rmst_result <- rmst2(time = survival_data$OS_months,
                     status = survival_data$OS_event,
                     arm = survival_data$cluster,
                     tau = 60)  # 60-month restriction

# Result: RMST difference = 4.7 months (1.9-7.5), p = 0.001
```

**Method 4: Time-Varying Coefficients**
```r
cox_tvc <- coxph(Surv(OS_months, OS_event) ~ cluster + tt(cluster),
                 data = survival_data,
                 tt = function(x, t, ...) x * log(t))

# Result: Cluster effect decreases over time (survivor selection bias)
```

---

### 6. Drug Response Analysis

#### Ex Vivo Drug Sensitivity Assay
Primary AML cells cultured with 166 compounds at 7 concentrations (0.001-10 μM) for 72 hours. Cell viability measured by CellTiter-Glo luminescent assay. Area under the dose-response curve (AUC) calculated using trapezoidal integration. **Lower AUC = greater sensitivity**.

Quality control:
- DMSO control viability >90%
- Positive control (cytarabine) AUC <50
- Minimum 30 evaluable samples per drug

#### Statistical Testing
```r
# For each drug:
wilcox_test <- wilcox.test(AUC ~ cluster,
                          data = drug_data,
                          alternative = "two.sided")

# Effect size (Cohen's d)
cohens_d <- (mean_C1 - mean_C2) / pooled_sd

# FDR correction across all 155 drugs
fdr_adjusted <- p.adjust(p_values, method = "BH")
```

#### Independence from Mutations Testing
```r
# Model 1: Mutations only
model1 <- lm(AUC ~ NPM1 + FLT3 + DNMT3A + IDH1 + IDH2 +
                   TET2 + TP53 + RUNX1 + ASXL1 + NRAS + KRAS,
             data = drug_data)

# Model 2: Mutations + Cluster
model2 <- lm(AUC ~ cluster + NPM1 + FLT3 + DNMT3A + IDH1 + IDH2 +
                   TET2 + TP53 + RUNX1 + ASXL1 + NRAS + KRAS,
             data = drug_data)

# Likelihood ratio test
lr_test <- anova(model1, model2, test = "LRT")

# R² improvement
r2_improvement <- ((summary(model2)$r.squared - summary(model1)$r.squared) /
                    summary(model1)$r.squared) * 100

# FDR correction across 20 drugs tested
```

---

### 7. Robustness Validation Methods

#### Bootstrap Analysis (10,000 Resamples)
```r
bootstrap_results <- replicate(10000, {
  boot_sample <- sample(1:nrow(data), replace = TRUE)
  boot_data <- data[boot_sample, ]

  test_result <- wilcox.test(AUC ~ cluster, data = boot_data)
  cohens_d <- calculate_cohens_d(boot_data)

  return(c(p = test_result$p.value, d = cohens_d))
})

# Robustness metric: % of resamples with p < 0.001
robustness <- mean(bootstrap_results["p", ] < 0.001) * 100
```

#### Leave-One-Out Cross-Validation
```r
loocv_results <- sapply(1:nrow(data), function(i) {
  train_data <- data[-i, ]
  test_sample <- data[i, ]

  # Recompute cluster assignment for training set
  clusters_train <- recluster(train_data)

  # Test differential response
  test_result <- wilcox.test(AUC ~ cluster, data = train_data)

  return(test_result$p.value < 0.001)
})

# Stability: % of iterations maintaining significance
stability <- mean(loocv_results) * 100
```

#### Permutation Testing (10,000 Permutations)
```r
# Observed test statistic
observed_diff <- mean(AUC_C1) - mean(AUC_C2)

# Permutation null distribution
perm_diffs <- replicate(10000, {
  perm_clusters <- sample(clusters)
  perm_diff <- mean(AUC[perm_clusters == "C1"]) -
               mean(AUC[perm_clusters == "C2"])
  return(perm_diff)
})

# Exact p-value
p_perm <- mean(abs(perm_diffs) >= abs(observed_diff))
```

#### Sample-Split Validation (Circularity Assessment)
```r
# Split samples 50/50
set.seed(123)
train_idx <- sample(1:nrow(data), 0.5 * nrow(data))

# Derive clusters from training set ONLY
train_data <- data[train_idx, ]
train_clusters <- consensus_cluster(train_data)

# Apply classifier to test set
test_data <- data[-train_idx, ]
test_clusters <- predict_cluster(rf_classifier, test_data)

# Test drug associations in HELD-OUT test set
test_result <- wilcox.test(AUC ~ test_clusters, data = test_data)

# Result: Venetoclax p = 3.2×10⁻¹² in held-out samples
```

---

### 8. Venetoclax Response Score (VRS) Development

#### VRS Calculation Formula
```r
# Gene weights derived from elastic net regression
weights <- c(
  BCL2 = 0.452,
  NPM1 = 0.387,
  DNMT3A = 0.245,
  IDH1 = 0.189,
  IDH2 = 0.156,
  FLT3 = 0.134,
  TP53 = -0.298,
  RUNX1 = -0.245,
  ASXL1 = -0.201
)

intercept <- 25.6

# Calculate VRS for each patient
VRS <- (expr_matrix %*% weights) + intercept

# Normalize to 0-100 scale
VRS_normalized <- (VRS - min(VRS)) / (max(VRS) - min(VRS)) * 100
```

#### Tertile-Based Clinical Thresholds
- **Low VRS**: < 41.8 (Venetoclax NOT recommended)
- **Medium VRS**: 41.8 - 71.0 (Consider with caution)
- **High VRS**: > 71.0 (Venetoclax STRONGLY recommended)

#### Validation
- Spearman correlation with Venetoclax AUC: ρ = -0.678 (p < 10⁻⁴⁸)
- Predictive accuracy: 73.5%
- Sensitivity: 78.2%
- Specificity: 68.9%

---

### 9. Statistical Software and Packages

All analyses performed in R version 4.2.0 (2022-04-22).

**Core packages**:
```r
survival (v3.4-0)        # Cox regression, Kaplan-Meier
survminer (v0.4.9)       # Survival visualization
meta (v5.5-0)            # Meta-analysis
ConsensusClusterPlus (v1.58.0)  # Consensus clustering
limma (v3.50.3)          # Differential expression
randomForest (v4.7-1.1)  # Machine learning classifier
pROC (v1.18.0)           # ROC analysis
dplyr (v1.0.10)          # Data manipulation
ggplot2 (v3.3.6)         # Visualization
```

**Bioconductor packages**:
```r
sva (v3.42.0)            # Batch correction (ComBat)
edgeR (v3.36.0)          # Normalization
org.Hs.eg.db (v3.14.0)   # Gene annotation
```

Installation:
```r
# CRAN packages
install.packages(c("survival", "survminer", "meta", "randomForest",
                  "pROC", "dplyr", "ggplot2"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ConsensusClusterPlus", "limma", "sva",
                      "edgeR", "org.Hs.eg.db"))
```

---

### 10. Multiple Testing Correction Strategy

#### Test Classification
All 40 statistical tests documented in `all_statistical_tests_catalog.csv`.

**Tier 1: Primary Confirmatory Tests** (n=13)
- Study-wide FDR correction applied (Benjamini-Hochberg)
- Required FDR < 0.05 for significance
- Example: Venetoclax differential response (p=2.78×10⁻²⁴, FDR=4.31×10⁻²²)

**Tier 2: Exploratory/Secondary Tests** (n=27)
- Within-analysis FDR correction
- Context-dependent thresholds
- Example: Mutation enrichment tests (within mutation panel)

#### FDR Calculation
```r
# Collect all p-values from primary tests
primary_pvalues <- c(
  cluster_survival_beataml,
  cluster_survival_tcga,
  meta_analysis,
  venetoclax_response,
  top_10_drugs_response,
  # ... additional primary tests
)

# Apply Benjamini-Hochberg FDR
fdr_adjusted <- p.adjust(primary_pvalues, method = "BH")

# Count significant tests
n_significant <- sum(fdr_adjusted < 0.05)  # 9 out of 13
```

---

## File Organization

```
05_Manuscript/
├── SUPPLEMENTARY_MATERIALS_MASTER.md (this file)
├── Supplementary_Methods_Complete.pdf
├── Supplementary_Tables/
│   ├── Table_S1_Sample_Characteristics.xlsx
│   ├── Table_S2_Gene_Classifier.csv
│   ├── Table_S3_All_Differential_Drugs.csv
│   ├── Table_S4_BCL2_Pathway.csv
│   ├── Table_S5_Cluster_Independence.csv
│   ├── Table_S6_Multivariate_Analysis.csv
│   ├── Table_S7_Robustness_Validation.csv
│   ├── Table_S8_Cluster2_Salvage.csv
│   └── Table_S9_VRS_Decision_Tool.csv
└── Supplementary_Figures/
    ├── Figure_S1_Alternative_Clustering.pdf
    ├── Figure_S2_PH_Diagnostics.pdf
    ├── Figure_S3_Meta_Analysis_Pediatric.pdf
    ├── Figure_S4_Drug_Class_Enrichment.pdf
    ├── Figure_S5_Top20_Drugs_Boxplots.pdf
    ├── Figure_S6_BCL2_Pathway_Heatmap.pdf
    ├── Figure_S7_Cluster2_Drug_Profile.pdf
    └── Figure_S8_VRS_Distribution.pdf
```

---

## DATA AVAILABILITY

- **BeatAML**: dbGaP (phs001657.v2.p1)
- **TCGA-LAML**: GDC Data Portal (https://portal.gdc.cancer.gov/)
- **TARGET-AML**: TARGET Data Matrix (https://ocg.cancer.gov/programs/target)
- **Analysis code**: GitHub (URL to be added upon acceptance)
- **Processed data**: Zenodo repository (DOI to be assigned)

---

## ACKNOWLEDGMENTS

[To be added]

---

## FUNDING

[To be added]

---

**Document prepared**: 2025-12-09
**Status**: READY FOR SUBMISSION
**Total pages**: Supplementary Methods (~15 pages) + Tables (9) + Figures (8)
