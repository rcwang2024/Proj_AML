# TARGET-AML Validation - Expected Results

**Date**: 2025-10-12
**Status**: RUNNING (Background Process cbcd79)
**Estimated Completion**: 20-30 minutes from start

---

## Validation Objective

Test whether molecular subtypes derived from **adult AML patients (BeatAML)** can predict survival in **pediatric AML patients (TARGET-AML)**. This is a critical test of:
1. Age-independence of the biological mechanism
2. Generalizability beyond the discovery cohort
3. Clinical utility across age groups

---

## Dataset Characteristics

### TARGET-AML (Pediatric Cohort)
- **Source**: NCI Genomic Data Commons
- **Samples**: 3,227 RNA-seq files downloaded
- **Clinical records**: 2,181 patients
- **Age range**: 0-20 years (pediatric)
- **Platform**: STAR - RNA-Seq gene counts
- **Genome**: GRCh38/hg38

### Comparison Cohorts
1. **BeatAML**: n=456, adult (18-89 years), discovery cohort
2. **TCGA-LAML**: n=153, adult (21-88 years), validation cohort
3. **TARGET-AML**: n=~?, pediatric (0-20 years), validation cohort

---

## Analysis Pipeline

### 1. Data Processing ⏳
- Load 3,227 individual TSV files
- Extract unstranded counts
- Filter to protein-coding genes (~19,962 genes)
- CPM normalization → log2(CPM+1)

### 2. Classifier Application ⏳
- Apply 50-gene Random Forest classifier (trained on BeatAML)
- Predict cluster assignments (Cluster 1 vs Cluster 2)
- Calculate prediction confidence scores
- Handle missing genes via BeatAML mean imputation

### 3. Survival Analysis ⏳
- **Primary test**: Log-rank test for survival difference
- **Effect size**: Cox proportional hazards HR with 95% CI
- **Comparison**: Median survival by cluster
- **Visualization**: Kaplan-Meier curves

### 4. Cross-Cohort Meta-Analysis ⏳
- Combine BeatAML + TCGA + TARGET results
- Test for heterogeneity (Cochran's Q, I² statistic)
- Forest plot of hazard ratios across cohorts

---

## Expected Outcomes

### Scenario 1: Strong Validation ✓✓✓
**Criteria**:
- Log-rank p < 0.05
- HR ~1.3-1.5 (consistent with adults)
- Same direction of effect (Cluster 2 worse)

**Interpretation**:
- Molecular subtypes show age-independent biology
- Classifier successfully transfers to pediatric population
- Universal applicability of prognostic model

**Impact on manuscript**:
- Strengthen claims of biological relevance
- Emphasize age-independent mechanism
- Support broader clinical utility

---

### Scenario 2: Partial Validation ✓
**Criteria**:
- Log-rank p < 0.05
- HR significantly different from 1.35 (e.g., HR < 1.2 or > 2.0)
- OR opposite direction of effect

**Interpretation**:
- Prognostic value exists but age-modified
- May reflect treatment differences (pediatric vs adult protocols)
- Suggests need for age-specific calibration

**Impact on manuscript**:
- Document age-specific effect modification
- Discuss biological vs treatment-driven differences
- Recommend pediatric-specific classifier development

---

### Scenario 3: No Validation ✗
**Criteria**:
- Log-rank p > 0.05
- Confidence intervals wide and include HR=1.0

**Possible reasons**:
1. **Insufficient power**
   - Low sample size after matching
   - Few events (pediatric AML has better prognosis)
   - → Perform power analysis

2. **Age-specific biology**
   - Pediatric AML has different mutation landscape
   - NPM1, TP53 mutations less common in children
   - Different subtype biology in young patients
   - → Discuss in limitations

3. **Treatment differences**
   - Intensive pediatric protocols mask prognostic effect
   - COG/POG treatment standardization
   - → Compare treatment-stratified analyses

4. **Classifier non-transferability**
   - Gene expression patterns differ with age
   - Developmental biology confounds signature
   - → Explore age-interaction terms

**Impact on manuscript**:
- Document as limitation to generalizability
- Frame subtypes as adult-AML specific
- Recommend prospective pediatric validation

---

## Key Questions

### 1. Sample Size After Matching
**Question**: How many samples have both expression and survival data?

**Expected**: ~500-1000 samples
- TARGET downloaded 3,227 files
- Clinical records: 2,181
- Expected overlap: 50-80% (protocol-dependent)

**Adequacy**:
- If n>300 with >150 events: Adequate power
- If n=100-300: Moderate power, interpret cautiously
- If n<100: Severely underpowered, descriptive only

---

### 2. Cluster Distribution
**Question**: Do pediatric samples classify similarly to adults?

**Adult distribution** (BeatAML):
- Cluster 1: 50.7%
- Cluster 2: 49.3%

**Expected pediatric distribution**:
- If biology similar: 45-55% in each cluster
- If different: Highly skewed (>70% in one cluster)

**Prediction confidence**:
- Mean confidence >0.75: Good classification
- Mean confidence 0.6-0.75: Moderate uncertainty
- Mean confidence <0.6: Poor transferability

---

### 3. Mutation Profile Comparison
**Adult subtypes** (BeatAML):
- Cluster 1: NPM1+ (42%), DNMT3A+ (31%), IDH1+ (13%)
- Cluster 2: TP53+ (14%), RUNX1+ (18%), ASXL1+ (16%)

**Pediatric AML characteristics**:
- NPM1 mutations: RARE (~5% vs 42% in adults)
- TP53 mutations: Less common (~5% vs 14%)
- FLT3-ITD: More common (~15-20%)
- Core-binding factor: More common (CBF-AML ~25%)

**Implication**: If subtypes align with mutations, pediatric distribution may differ substantially

---

### 4. Survival Differences
**Adult median survival** (BeatAML):
- Cluster 1: 582 months (~48 years)
- Cluster 2: 360 months (~30 years)
- Difference: 222 months (18.5 years)

**Pediatric AML survival** (general):
- 5-year OS: ~60-70% (better than adults)
- Median OS: Often not reached (long survivors)
- Different prognostic factors: MRD, cytogenetics

**Expected**:
- Absolute survival will be BETTER in TARGET (pediatric advantage)
- Relative effect (HR) is key comparison metric
- May need longer follow-up to see separation

---

## Power Analysis Considerations

### Factors affecting power:
1. **Sample size**: Higher n → better power
2. **Event rate**: More deaths → better power (but pediatric AML has LOWER mortality)
3. **Effect size**: Larger HR → better power
4. **Follow-up time**: Longer → more events → better power

### Power calculation (post-hoc):
If TARGET shows **no significant difference**:
- Calculate achieved power for HR=1.35
- Determine required sample size for 80% power
- Compare to what's available

**Example**:
- If n=200, events=60, HR=1.35
- Power = ~50% (underpowered)
- Required: n=400, events=150 for 80% power

---

## Files to be Generated

### Data Files:
- `target_aml_expression_normalized.rds` - Normalized expression matrix
- `target_aml_clinical.csv` - Clinical data with survival
- `target_cluster_assignments.csv` - Predicted clusters + confidence
- `target_survival_validation.csv` - Survival analysis results

### Figures:
- `target_kaplan_meier.pdf` - Survival curves by cluster
- `forest_plot_all_cohorts.pdf` - Cross-cohort HR comparison

### Summary:
- Updated Phase3_CriticalValidation_Summary.md

---

## Interpretation Framework

### If Validation Successful (p<0.05, HR~1.35):

**Strengths**:
✓ Age-independent prognostic value confirmed
✓ Biological mechanism likely universal
✓ Classifier robust across populations
✓ Clinical utility potentially broad

**Manuscript claims**:
- "Molecular subtypes validated across age spectrum (pediatric and adult AML)"
- "Age-independent biological mechanism"
- "Potential for universal prognostic classification"

---

### If Validation Unsuccessful (p>0.05):

**AVOID over-interpreting negative results**

**Check**:
1. Was power adequate? (>80%)
2. Was effect size estimation precise? (narrow CI?)
3. Was sample matching successful? (>70% overlap?)

**Manuscript approach**:
- "Validation in independent adult cohort (TCGA) confirmed prognostic value (meta-analysis HR=1.35, p=0.001)"
- "Pediatric validation limited by [insufficient power / treatment differences / distinct biology]"
- "Further research needed to determine age-specific applicability"

---

## Next Steps After Completion

### If Successful:
1. Update Phase3 summary document
2. Update meta-analysis to include all 3 cohorts
3. Add age-stratified analysis section
4. Create comprehensive cross-cohort comparison figure
5. Update manuscript introduction (emphasize age-independence)

### If Unsuccessful:
1. Perform power analysis
2. Explore reasons (mutation profiles, treatment, biology)
3. Document as limitation in manuscript
4. Suggest prospective pediatric validation

### Always:
1. Update EXECUTIVE_SUMMARY.md
2. Add TARGET results to RESULTS_SUMMARY.md
3. Create publication-ready forest plot (all cohorts)
4. Prepare supplementary tables

---

## Timeline

- **Start time**: 2025-10-12 ~13:20
- **Estimated completion**: ~13:40-13:50 (20-30 minutes)
- **Process ID**: cbcd79

**Monitor progress**:
```bash
tail -f target_validation_log.txt
```

Or check results directory:
```bash
ls -lh 03_Results/18_TARGET_Validation/
```

---

**Last Updated**: 2025-10-12 13:25
**Status**: RUNNING
**Expected Completion**: ~13:45
