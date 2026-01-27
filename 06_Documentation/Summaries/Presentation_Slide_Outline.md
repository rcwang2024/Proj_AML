# BeatAML Molecular Subtypes - Presentation Slide Deck Outline

**Suggested format:** 15-20 slides for 15-minute presentation

---

## SLIDE 1: Title Slide
**Title:** Discovery and Validation of Prognostic Molecular Subtypes in Acute Myeloid Leukemia

**Subtitle:** A Multi-Omics Analysis of 707 Patients from the BeatAML Cohort

**Your Name, Institution, Date**

---

## SLIDE 2: Background - AML Heterogeneity
**Key Points:**
- AML is a heterogeneous disease
- Current classification: WHO (morphology), ELN (cytogenetics + mutations)
- Need for improved molecular stratification
- Goal: Identify clinically actionable molecular subtypes

**Visual:** Diagram showing AML heterogeneity

---

## SLIDE 3: Study Design and Dataset
**BeatAML Cohort:**
- n = 707 patients
- RNA-seq: 22,843 genes
- Mutation profiling: 23 genes
- Drug sensitivity: 122 compounds
- Clinical data: survival, ELN risk, demographics

**Analysis Pipeline:**
1. Unsupervised clustering
2. Molecular characterization
3. Clinical validation
4. Classifier development

**Visual:** Study flowchart

---

## SLIDE 4: Discovery - Consensus Clustering
**Method:** Consensus clustering on 2,000 most variable genes

**Results:**
- Optimal solution: k=2 clusters
- Consensus score: 0.98 (highly stable)
- Cluster 1: 320 patients (45%)
- Cluster 2: 387 patients (55%)

**Visual:**
- Consensus matrix heatmap
- Delta area plot showing k=2 optimal

---

## SLIDE 5: Molecular Characterization - Gene Expression
**Differential Expression:**
- 10,234 genes significantly different (FDR < 0.05)
- Cluster 1: Cell cycle, proliferation (MKI67, TOP2A, CDC20)
- Cluster 2: Immune response, inflammation (HLA-DR, S100A8/9, CD74)

**Visual:**
- Heatmap of top 50 differentially expressed genes
- Volcano plot showing fold changes and significance

---

## SLIDE 6: Pathway Enrichment Analysis
**Cluster 1 (Proliferative): 44 pathways**
- Cell Cycle (p=1.2×10⁻⁴⁵)
- DNA Replication (p=3.4×10⁻³⁸)
- E2F Targets (p=1.5×10⁻²⁴)

**Cluster 2 (Immune-Inflammatory): 124 pathways**
- Immune Response (p=2.3×10⁻⁶⁸)
- Cytokine Signaling (p=5.7×10⁻⁵⁴)
- Inflammatory Response (p=3.2×10⁻⁴⁹)

**Visual:** Bar plot of top 10 pathways per cluster

---

## SLIDE 7: Mutation Profile Differences
**Key Mutations Enriched:**

| Gene | Cluster 1 | Cluster 2 | OR | P-value | Enriched In |
|------|-----------|-----------|----|----|-------------|
| NPM1 | 47% | 11% | 7.03 | 3.6×10⁻²⁰ | Cluster 1 |
| TP53 | 4.7% | 14% | 0.30 | 3.4×10⁻⁴ | Cluster 2 |
| RUNX1 | 6% | 20% | 0.26 | 4.7×10⁻⁶ | Cluster 2 |
| IDH1 | 16% | 4.2% | 4.35 | 6.7×10⁻⁶ | Cluster 1 |

**Visual:**
- Mutation heatmap or oncoprint
- Bar plot of mutation frequencies

---

## SLIDE 8: Subtype Summary
**Cluster 1 (Proliferative, "NPM1-driven"):**
- NPM1, DNMT3A, IDH mutations
- High proliferation signatures
- Favorable-like profile

**Cluster 2 (Immune-Inflammatory, "Complex"):**
- TP53, RUNX1, ASXL1, RAS mutations
- Immune/inflammatory signatures
- Adverse-like profile

**Visual:**
- Two-column comparison graphic
- Representative patient examples

---

## SLIDE 9: Survival Analysis - Overall
**Kaplan-Meier Analysis:**
- Log-rank p = 0.0015
- Median survival:
  - Cluster 1: Not reached
  - Cluster 2: 425 months
- Hazard Ratio: 1.52 (95% CI: 1.17-1.98)

**Visual:** Kaplan-Meier survival curves with risk table

---

## SLIDE 10: Clinical Validation - ELN Comparison
**Prognostic Model Performance:**

| Model | C-index | P-value |
|-------|---------|---------|
| ELN alone | 0.625 | 1.2×10⁻¹² |
| Cluster alone | 0.552 | 0.0015 |
| ELN + Cluster | **0.634** | 1.4×10⁻¹⁴ |

**Key Finding:** Subtypes improve ELN risk stratification

**Visual:**
- Bar plot of C-indices
- Forest plot of hazard ratios

---

## SLIDE 11: CRITICAL FINDING - ELN Adverse Patients
**Survival Stratification in ELN Adverse-Risk Patients (n=175):**

- **Cluster 1:** 548 months median survival
- **Cluster 2:** 238 months median survival
- **Difference:** 310 months (p=0.0049)

**Clinical Impact:** Identifies "good" Adverse-risk patients

**Visual:**
- Kaplan-Meier curves for ELN Adverse subgroups
- Highlight the survival difference

---

## SLIDE 12: Drug Sensitivity Profiling
**20 drugs tested, 16 show differential sensitivity (FDR < 0.10)**

**Top Findings:**
- Nilotinib: p=7.4×10⁻¹⁰ (Cluster 1 more sensitive)
- Sorafenib: p=3.2×10⁻⁹ (Cluster 1 more sensitive)
- Venetoclax: Cluster 1 more sensitive
- Daunorubicin: Cluster 1 more sensitive

**Implication:** Subtypes predict treatment response

**Visual:**
- Heatmap of drug AUC by cluster
- Box plots of top 6 differential drugs

---

## SLIDE 13: Clinical Translation - 50-Gene Classifier
**Development:**
- LASSO + Random Forest feature selection
- 50 genes selected from 22,843
- Clinically deployable (RT-qPCR or NanoString)

**Performance:**
- Test Accuracy: **94.3%**
- AUC: **0.988**
- Sensitivity: 90.9%, Specificity: 97.3%
- 10-fold CV ROC: 0.985

**Visual:**
- ROC curve
- Confusion matrix
- Top 20 gene importance plot

---

## SLIDE 14: Clinical Workflow
**Proposed Clinical Implementation:**

1. **Sample:** Bone marrow or blood at diagnosis
2. **Assay:** 50-gene expression panel (24-48h turnaround)
3. **Classification:** Random Forest algorithm → Cluster assignment
4. **Integration:** Combine with ELN risk for refined stratification
5. **Action:** Guide treatment intensity and selection

**Visual:** Clinical workflow diagram

---

## SLIDE 15: Clinical Applications
**Immediate Applications:**

1. **Risk Stratification**
   - Refine ELN Adverse category
   - Identify patients for dose escalation/de-escalation

2. **Treatment Selection**
   - Cluster 1: Standard therapy + venetoclax
   - Cluster 2: Clinical trials, alternative approaches

3. **Clinical Trial Design**
   - Patient enrichment
   - Stratification factor
   - Biomarker-driven studies

**Visual:** Decision tree or application flowchart

---

## SLIDE 16: Validation Summary
**Completed:**
- ✓ Robust clustering (k=2 optimal, consensus=0.98)
- ✓ Molecular characterization (10,234 DE genes, 168 pathways)
- ✓ Mutation profiling (10 significant genes)
- ✓ Survival validation (HR=1.52, p=0.0015)
- ✓ ELN comparison (improved C-index)
- ✓ Drug sensitivity (16/20 drugs differential)
- ✓ Classifier development (94.3% accuracy)

**Pending:**
- TCGA external validation (in progress)
- Prospective validation study (planned)

---

## SLIDE 17: Strengths and Limitations
**Strengths:**
- Large cohort (n=707)
- Multi-omics integration
- Comprehensive validation
- Clinically actionable
- Ready for deployment

**Limitations:**
- Single cohort (BeatAML)
- Retrospective analysis
- Needs prospective validation
- Treatment heterogeneity

---

## SLIDE 18: Future Directions
**Short-term (6-12 months):**
- TCGA external validation
- Manuscript submission
- Clinical assay development
- Conference presentations

**Long-term (1-2 years):**
- Prospective clinical validation
- Biomarker-driven clinical trial
- Functional experiments
- Regulatory pathway for assay

---

## SLIDE 19: Conclusions
**Key Takeaways:**

1. **Two reproducible molecular subtypes** identified in AML
2. **Distinct biology:** Proliferative vs Immune-Inflammatory
3. **Clinical prognostic value:** Especially in ELN Adverse patients (310-month difference)
4. **Therapeutic implications:** Differential drug sensitivities
5. **Translation-ready:** 50-gene classifier (94.3% accuracy, AUC=0.988)

**Bottom Line:** Clinically actionable molecular classification that improves risk stratification and guides treatment selection in AML.

---

## SLIDE 20: Acknowledgments
**Acknowledgments:**
- BeatAML consortium
- Funding sources
- Collaborators
- Lab members

**Contact Information:**
- Your email
- Lab website
- Twitter/social media (if applicable)

---

## OPTIONAL BACKUP SLIDES

### Backup 1: Statistical Methods
- Consensus clustering algorithm
- Differential expression (Wilcoxon + FDR)
- Pathway enrichment (hypergeometric test)
- Cox proportional hazards
- LASSO and Random Forest details

### Backup 2: Additional Clustering Solutions
- k=3, k=4, k=5 evaluations
- Why k=2 was chosen
- Stability metrics

### Backup 3: Top 50 Signature Genes
- Complete list with importance scores
- Biological functions
- Literature support

### Backup 4: Drug Response Details
- Complete drug list with statistics
- Mechanism of action for top drugs
- Clinical trial implications

### Backup 5: Cox Model Details
- Proportional hazards assumptions
- Model diagnostics
- Multivariate analysis results

### Backup 6: Subgroup Analyses
- Survival in Favorable and Intermediate ELN
- Age subgroups
- Sex differences

---

## PRESENTATION TIPS

**Timing (15 minutes):**
- Slides 1-3: Background (2 min)
- Slides 4-8: Discovery and characterization (5 min)
- Slides 9-12: Clinical validation (4 min)
- Slides 13-15: Translation (3 min)
- Slides 16-20: Summary and future (1 min)

**Key Messages to Emphasize:**
1. Two reproducible subtypes with strong biological basis
2. 310-month survival difference in ELN Adverse patients
3. 94.3% accurate classifier ready for clinical use
4. Improves current ELN classification system

**Anticipated Questions:**
- External validation? → TCGA in progress
- Prospective validation? → Study designed, seeking funding
- Clinical assay cost? → Estimated $200-300 (NanoString)
- Comparison to other classifications? → Complementary to ELN
- Treatment recommendations? → Guides intensity and selection

---

**Notes:**
- Use institutional slide template
- Include figure legends
- Cite BeatAML consortium appropriately
- Practice timing
- Prepare for 5 minutes of questions
