# CLINICAL TRIAL PROTOCOL

## Molecular Cluster-Guided Treatment Selection in Adult Acute Myeloid Leukemia: A Phase II Randomized Trial

---

**Protocol ID**: AML-CLUSTER-001
**Protocol Version**: 1.0
**Protocol Date**: 2025-12-09
**Trial Phase**: Phase II
**Study Design**: Randomized, Open-Label, Multi-Center

---

## PROTOCOL SYNOPSIS

| Element | Description |
|---------|-------------|
| **Study Title** | Molecular Cluster-Guided vs Standard Treatment Selection in Adult AML |
| **Short Title** | CLUSTER-AML Trial |
| **Phase** | Phase II (Efficacy/Safety) |
| **Study Design** | Randomized, open-label, parallel-group, multi-center |
| **Primary Objective** | Compare complete remission (CR/CRi) rate between cluster-guided vs standard treatment |
| **Secondary Objectives** | Overall survival, progression-free survival, minimal residual disease, safety |
| **Study Population** | Newly diagnosed adult AML patients (age ≥60 years OR unfit for intensive therapy) |
| **Sample Size** | 200 patients (100 per arm) |
| **Treatment Arms** | **ARM A**: Cluster-guided treatment<br>**ARM B**: Standard of care |
| **Primary Endpoint** | CR/CRi rate after 2 cycles |
| **Key Secondary Endpoints** | OS at 12 months, PFS, MRD negativity rate |
| **Study Duration** | Enrollment: 24 months<br>Follow-up: 24 months<br>Total: 48 months |
| **Sponsor** | [Institution name] |
| **Principal Investigator** | [Name, MD] |

---

## 1. BACKGROUND AND RATIONALE

### 1.1 Disease Background

Acute myeloid leukemia (AML) is a heterogeneous hematologic malignancy with poor outcomes in older adults. Standard treatment includes intensive chemotherapy (7+3) for fit patients or hypomethylating agents (HMA) with venetoclax for unfit patients. However, response rates vary substantially, and predictive biomarkers beyond NPM1 mutation status are urgently needed.

### 1.2 Preliminary Data

Our comprehensive multi-omics analysis of 2,535 AML patients across 3 independent cohorts identified two molecular subtypes with distinct treatment sensitivities:

**Cluster 1 (45% of patients)**:
- NPM1-enriched (42% vs 10%, OR=6.5, p<10⁻¹⁵)
- **Extraordinary Venetoclax sensitivity** (p=2.78×10⁻²⁴, Cohen's d=1.25)
- Mean AUC: 107 (highly sensitive)
- Venetoclax Response Score (VRS) >71 (High)

**Cluster 2 (55% of patients)**:
- TP53/RUNX1/ASXL1-enriched (adverse-risk mutations)
- **Venetoclax resistance** (AUC: 192 vs 107, 1.79× higher)
- **Alternative drug sensitivities**: Panobinostat (p=1.1×10⁻¹²), Selumetinib (p=4.5×10⁻¹¹)
- VRS <42 (Low)

**Key Finding**: Clusters provide **independent predictive value** beyond mutations for drug response (19/20 drugs, mean +42% R² improvement, FDR<0.05).

### 1.3 Rationale for This Trial

Current standard of care does NOT account for molecular cluster assignment. This trial tests whether **cluster-guided treatment selection** improves outcomes compared to mutation-based standard of care.

**Hypothesis**: Patients assigned to cluster-appropriate therapy (Venetoclax for Cluster 1, alternatives for Cluster 2) will achieve higher CR/CRi rates than standard treatment.

### 1.4 Clinical Equipoise

**Current practice**: Venetoclax + HMA given to most unfit patients regardless of cluster
**Proposed practice**: Cluster 1 → Venetoclax + HMA; Cluster 2 → Panobinostat-based therapy

**Equipoise exists because**:
- Panobinostat shows strong efficacy in Cluster 2 (ex vivo data)
- Cluster 2 patients have poor Venetoclax response (AUC 192)
- No prospective data comparing cluster-guided vs standard treatment

---

## 2. STUDY OBJECTIVES

### 2.1 Primary Objective

To compare the **complete remission rate (CR + CRi)** after 2 cycles of treatment between:
- **ARM A**: Molecular cluster-guided treatment selection
- **ARM B**: Standard of care treatment selection

### 2.2 Secondary Objectives

1. Compare **overall survival (OS)** at 12 and 24 months
2. Compare **progression-free survival (PFS)**
3. Compare **minimal residual disease (MRD) negativity rate** after 2 cycles
4. Evaluate **safety and tolerability** of cluster-guided treatment
5. Assess **quality of life** (EORTC QLQ-C30)
6. Validate **Venetoclax Response Score (VRS)** prospectively
7. Assess **cost-effectiveness** of cluster-guided approach

### 2.3 Exploratory Objectives

1. Correlate cluster assignment with treatment-emergent mutations
2. Assess cluster stability at relapse
3. Identify predictors of response within each cluster
4. Evaluate drug combinations (Panobinostat + Selumetinib) in Cluster 2

---

## 3. STUDY DESIGN

### 3.1 Study Schema

```
                      SCREENING
                  (Bone marrow biopsy
                   RNA-seq, mutations)
                          |
                    ENROLLMENT
                    (N = 200)
                          |
                   RANDOMIZATION
                    (1:1 ratio)
              ____________|____________
             |                         |
         ARM A                     ARM B
   CLUSTER-GUIDED              STANDARD OF CARE
       (n=100)                     (n=100)
             |                         |
     Cluster Assignment         ELN Risk/NPM1
             |                         |
    ┌────────┴────────┐               |
    |                 |               |
CLUSTER 1       CLUSTER 2       NPM1+/-
Venetoclax      Panobinostat    Standard
+ HMA           + Selumetinib   Treatment
    |                 |               |
    └────────┬────────┴───────────────┘
             |
      TREATMENT CYCLES
      (Up to 6 cycles)
             |
      RESPONSE ASSESSMENT
      (After cycle 2, 4, 6)
             |
      FOLLOW-UP
      (Every 3 months × 24 months)
```

### 3.2 Randomization

- **Timing**: Within 7 days of enrollment
- **Method**: Centralized web-based randomization (REDCap)
- **Ratio**: 1:1 (cluster-guided : standard of care)
- **Stratification factors**:
  1. Age (<75 vs ≥75 years)
  2. ELN 2022 Risk (Favorable vs Intermediate vs Adverse)
  3. De novo vs secondary AML

### 3.3 Blinding

**Open-label** design:
- Patients and investigators aware of treatment assignment
- Response assessments by blinded central review
- MRD assessments by blinded core laboratory

---

## 4. STUDY POPULATION

### 4.1 Inclusion Criteria

1. **Age**: ≥60 years
2. **Diagnosis**: Newly diagnosed AML by WHO 2016 criteria
   - De novo AML OR
   - Secondary AML (MDS-related, therapy-related)
3. **Performance Status**: ECOG 0-2
4. **Fitness**: Unfit for intensive 7+3 chemotherapy (≥1 of):
   - Age ≥75 years
   - ECOG PS 2
   - Cardiac comorbidity (EF <50%)
   - Pulmonary comorbidity (DLCO <65%)
   - Hepatic dysfunction (bilirubin >1.5× ULN)
   - Renal dysfunction (CrCl 30-60 mL/min)
5. **Organ Function**:
   - Total bilirubin ≤2.0× ULN (unless Gilbert's)
   - AST/ALT ≤3.0× ULN
   - Creatinine clearance ≥30 mL/min
   - Absolute neutrophil count ≥0.5 × 10⁹/L (unless due to AML)
6. **Sample Availability**: Adequate bone marrow for RNA-seq (≥20% blasts)
7. **Informed Consent**: Able and willing to provide written informed consent

### 4.2 Exclusion Criteria

1. **AML Subtype**: Acute promyelocytic leukemia (APL, t(15;17))
2. **Prior Treatment**: Any prior AML-directed therapy (except hydroxyurea)
3. **CNS Involvement**: Active CNS leukemia
4. **HIV/Hepatitis**: Known HIV, active HBV, or HCV
5. **Other Malignancy**: Active malignancy requiring treatment
6. **Pregnancy**: Pregnant or breastfeeding women
7. **Contraindications**: Known hypersensitivity to study drugs
8. **Psychiatric**: Unable to comply with study procedures

### 4.3 Sample Size Justification

**Primary endpoint**: CR/CRi rate after 2 cycles

**Assumptions**:
- Standard arm CR/CRi rate: 50% (based on VIALE-A trial)
- Cluster-guided arm CR/CRi rate: 65% (15% absolute improvement)
- Alpha: 0.05 (two-sided)
- Power: 80%
- Dropout rate: 10%

**Calculation**:
- Required per arm: 91 patients
- With 10% dropout: 100 patients per arm
- **Total sample size: 200 patients**

**Expected cluster distribution** (Arm A only):
- Cluster 1 (Venetoclax): ~45 patients
- Cluster 2 (Panobinostat): ~55 patients

---

## 5. TREATMENT PLAN

### 5.1 ARM A: Cluster-Guided Treatment

#### ARM A - Cluster 1 (Venetoclax-Sensitive)

**Treatment**: Venetoclax + Azacitidine

**Dosing**:
- **Venetoclax**: Days 1-28
  - Ramp-up: 100 mg (D1), 200 mg (D2), 400 mg (D3-28)
  - Cycles 2+: 400 mg daily × 28 days
- **Azacitidine**: 75 mg/m² SC/IV days 1-7
- **Cycle length**: 28 days

**Dose Modifications**:
- Hematologic toxicity: Hold until ANC >500 and platelets >50K
- Non-hematologic grade ≥3: Hold until recovery to grade ≤1

#### ARM A - Cluster 2 (Venetoclax-Resistant)

**Treatment**: Panobinostat + Selumetinib + Azacitidine

**Dosing**:
- **Panobinostat**: 20 mg PO three times weekly (M/W/F), weeks 1-2 of each 28-day cycle
- **Selumetinib**: 75 mg PO BID, days 1-14 (may escalate to 100 mg BID if tolerated)
- **Azacitidine**: 75 mg/m² SC/IV days 1-7
- **Cycle length**: 28 days

**Rationale**: Combination targets complementary pathways (HDAC + MEK)

**Dose Modifications**:
- Diarrhea: Aggressive anti-diarrheal support
- Cardiac: Monitor QTc weekly, hold if QTc >480 ms
- Hematologic: Same as Cluster 1

### 5.2 ARM B: Standard of Care

**Treatment assignment based on NPM1 status**:

#### NPM1-Mutated Patients

**Treatment**: Venetoclax + Azacitidine (same as Cluster 1 regimen)

#### NPM1-Wild-Type Patients

**Treatment**: Standard HMA-based therapy per investigator choice:
- **Option 1**: Azacitidine 75 mg/m² days 1-7 × 28-day cycles
- **Option 2**: Decitabine 20 mg/m² days 1-5 × 28-day cycles
- **Option 3**: Venetoclax + HMA (if investigator discretion)

### 5.3 Treatment Duration

**Both Arms**:
- **Induction**: Cycles 1-2 (minimum)
- **Consolidation**: Cycles 3-6 (if CR/CRi achieved)
- **Maintenance**: Continue until relapse, unacceptable toxicity, or patient/physician decision
- **Maximum**: 24 cycles (2 years)

### 5.4 Supportive Care

**Required for all patients**:
- **Tumor lysis prophylaxis**: Allopurinol 300 mg daily (or rasburicase if high risk)
- **Antimicrobial prophylaxis**:
  - Antibacterial: Levofloxacin 500 mg daily
  - Antifungal: Fluconazole 400 mg daily (or posaconazole if neutropenic)
  - Antiviral: Acyclovir 400 mg BID
- **Growth factors**: G-CSF allowed per ASCO guidelines
- **Transfusions**: Maintain Hgb >8 g/dL, platelets >10K (prophylactic)

### 5.5 Prohibited Medications

- Strong CYP3A4 inhibitors (increase Venetoclax levels)
- Strong CYP3A4 inducers (decrease Venetoclax levels)
- QTc-prolonging agents (for Panobinostat arm)
- Other investigational agents

---

## 6. STUDY ASSESSMENTS

### 6.1 Screening (Days -14 to -1)

- Informed consent
- Medical history, demographics
- Physical examination, vital signs, ECOG PS
- **Bone marrow biopsy/aspirate**:
  - Morphology, flow cytometry, cytogenetics
  - **RNA extraction for RNA-seq** (send to central lab)
  - Mutation panel (NGS): NPM1, FLT3, IDH1/2, TP53, RUNX1, ASXL1, etc.
- Complete blood count with differential
- Comprehensive metabolic panel
- Coagulation studies (PT/PTT/INR)
- Pregnancy test (women of childbearing potential)
- ECG, echocardiogram (baseline EF)
- Quality of life: EORTC QLQ-C30

**Turnaround time for RNA-seq**: 7-10 days

### 6.2 Baseline (Day 1, Cycle 1)

- Physical examination, vital signs, ECOG PS
- CBC with differential
- CMP
- Randomization and treatment initiation

### 6.3 During Treatment

**Weekly (Cycles 1-2)**:
- CBC with differential
- CMP (weekly × 2, then every 2 weeks)
- Adverse events assessment

**Every Cycle (Days 1)**:
- Physical examination, ECOG PS
- CBC, CMP
- Adverse events
- Quality of life (Q2 cycles)

### 6.4 Response Assessments

**Timing**: After Cycle 2, 4, 6, then every 3 cycles

**Bone Marrow Assessment**:
- Morphology (blast percentage)
- Flow cytometry (MRD by multi-parameter flow, sensitivity 10⁻⁴)
- Cytogenetics (if abnormal at baseline)
- **Molecular MRD** (if mutation present at diagnosis)

**Response Criteria**: Modified ELN 2022
- **CR**: <5% blasts, ANC ≥1.0, platelets ≥100K, no extramedullary disease
- **CRi**: CR but ANC <1.0 or platelets <100K
- **MLFS**: <5% blasts, no count recovery
- **PR**: 50% decrease in blasts (but >5%)
- **Resistant disease**: Failure to achieve CR/CRi/MLFS after 2 cycles

**MRD Definition**:
- **MRD-negative**: <0.1% blasts by flow cytometry AND no detectable mutations by NGS (VAF <0.01%)
- **MRD-positive**: ≥0.1% blasts or detectable mutations

### 6.5 End of Treatment

- Physical examination
- CBC, CMP
- Bone marrow assessment (if not done recently)
- Quality of life
- Reason for treatment discontinuation

### 6.6 Follow-Up

**Frequency**: Every 3 months × 24 months, then every 6 months × 2 years

**Assessments**:
- Physical examination
- CBC with differential
- AML relapse assessment (bone marrow if clinical suspicion)
- Survival status
- Subsequent therapies

---

## 7. ENDPOINTS

### 7.1 Primary Endpoint

**Complete Remission Rate (CR + CRi) after 2 cycles**

**Definition**: Proportion of patients achieving CR or CRi per ELN 2022 criteria after 2 treatment cycles

**Analysis**: Chi-square test comparing ARM A vs ARM B

### 7.2 Secondary Endpoints

1. **Overall Survival (OS)**
   - Time from randomization to death from any cause
   - Kaplan-Meier curves, log-rank test, Cox regression

2. **Progression-Free Survival (PFS)**
   - Time from randomization to relapse or death
   - Kaplan-Meier curves, log-rank test

3. **MRD Negativity Rate**
   - Proportion achieving MRD-negative status after 2 cycles
   - Chi-square test

4. **Duration of Response (DOR)**
   - Time from CR/CRi to relapse or death
   - Median DOR by Kaplan-Meier

5. **Safety**
   - Adverse events (CTCAE v5.0)
   - Treatment-related mortality (TRM) at Day 30 and Day 60
   - Dose modifications and discontinuations

6. **Quality of Life**
   - Change from baseline in EORTC QLQ-C30 scores
   - Mixed models for repeated measures

### 7.3 Exploratory Endpoints

1. **VRS Validation**
   - Correlation of VRS with Venetoclax response in ARM A Cluster 1 patients
   - Spearman correlation, ROC analysis

2. **Cost-Effectiveness**
   - Incremental cost-effectiveness ratio (ICER)
   - Cost per quality-adjusted life year (QALY)

3. **Biomarker Discovery**
   - Identify predictors of response within each cluster
   - Machine learning models

---

## 8. STATISTICAL ANALYSIS PLAN

### 8.1 Analysis Populations

- **Intent-to-Treat (ITT)**: All randomized patients (primary analysis)
- **Per-Protocol (PP)**: Patients completing ≥2 cycles per protocol
- **Safety**: All patients receiving ≥1 dose of study drug

### 8.2 Primary Analysis

**Null Hypothesis**: CR/CRi rate in ARM A = CR/CRi rate in ARM B

**Statistical Test**: Chi-square test (or Fisher's exact if cell counts <5)

**Significance Level**: α = 0.05 (two-sided)

**Sample Size**: 200 patients provides 80% power to detect 15% improvement (50% → 65%)

### 8.3 Secondary Analyses

**Overall Survival**:
- Kaplan-Meier curves with 95% CI
- Log-rank test for comparison
- Cox proportional hazards model adjusting for stratification factors
- HR with 95% CI

**MRD Negativity**:
- Chi-square test
- Logistic regression adjusting for baseline factors

**Safety**:
- Descriptive statistics (frequencies, percentages)
- Comparison of Grade 3-4 AE rates (Chi-square)

### 8.4 Subgroup Analyses

**Pre-specified subgroups** (interaction tests):
1. Age (<75 vs ≥75 years)
2. ELN Risk (Favorable/Intermediate vs Adverse)
3. De novo vs secondary AML
4. TP53 mutation status
5. **Cluster assignment (ARM A only)**
6. **VRS tertile (ARM A only)**

**Cluster-specific analysis (ARM A)**:
- Compare CR/CRi rate: Cluster 1 (Venetoclax) vs Cluster 2 (Panobinostat)
- Expected Cluster 1 CR/CRi: 75% (based on VIALE-A NPM1+ subgroup)
- Expected Cluster 2 CR/CRi: 55% (estimated from ex vivo data)

### 8.5 Interim Analysis

**Timing**: After 100 patients complete 2 cycles (~12 months)

**Purpose**: Futility and safety monitoring

**Futility Rule**: Stop if conditional power <20% for primary endpoint

**Safety Rule**: Stop if treatment-related mortality >15% in either arm

**Statistical Adjustment**: O'Brien-Fleming spending function (α=0.003 for interim)

**DSMB**: Independent Data Safety Monitoring Board reviews interim data

### 8.6 Missing Data

**Primary Endpoint (CR/CRi)**:
- Missing considered as non-responders (conservative)
- Sensitivity analysis: Multiple imputation

**Time-to-Event Endpoints (OS, PFS)**:
- Censored at last known alive date
- Sensitivity analysis: Competing risks (death without relapse)

---

## 9. SAFETY MONITORING

### 9.1 Adverse Event Reporting

**Grading**: CTCAE v5.0

**Expected Toxicities**:
- **Venetoclax + HMA**: Cytopenias (90%), GI symptoms (40%), tumor lysis (5%)
- **Panobinostat + Selumetinib**: Diarrhea (60%), nausea/vomiting (50%), cytopenias (70%), QTc prolongation (15%)

**Serious Adverse Events (SAEs)**:
- Report within 24 hours to sponsor
- Report to IRB per institutional policy
- Report to FDA within 7 days (unexpected SAEs)

### 9.2 Dose-Limiting Toxicities (DLTs)

**Definition** (first 2 cycles):
- Grade 4 neutropenia >14 days
- Grade 4 thrombocytopenia >14 days with bleeding
- Any Grade ≥3 non-hematologic toxicity (except controlled nausea/vomiting)
- QTc prolongation >500 ms or >60 ms increase from baseline

**DLT Management**:
- Hold study drug until recovery
- Resume at reduced dose per protocol
- If recurrent DLT, consider permanent discontinuation

### 9.3 Data Safety Monitoring Board (DSMB)

**Composition**:
- 3 independent members (2 oncologists, 1 biostatistician)
- No conflicts of interest with study sponsor/investigators

**Meetings**:
- Every 6 months during enrollment
- Ad hoc if safety concerns arise

**Review**:
- Enrollment, protocol deviations
- SAEs, DLTs, treatment-related deaths
- Interim efficacy analysis (blinded)
- Recommendations: Continue, modify, or terminate trial

---

## 10. ETHICAL CONSIDERATIONS

### 10.1 Institutional Review Board (IRB)

- Protocol approved by IRB at each participating site
- Annual continuing review
- Protocol amendments require IRB re-approval

### 10.2 Informed Consent

- Written informed consent obtained before any study procedures
- Consent form approved by IRB
- Patient receives copy of signed consent
- Re-consent if protocol amended with significant changes

### 10.3 Risks and Benefits

**Risks**:
- Standard chemotherapy risks (cytopenias, infections)
- Venetoclax: Tumor lysis syndrome (mitigated with prophylaxis)
- Panobinostat: Diarrhea, QTc prolongation (mitigated with monitoring)
- RNA-seq: Minimal additional risk (already standard at many centers)

**Benefits**:
- Potential for improved treatment outcomes
- Access to personalized therapy
- Close monitoring and supportive care

**Risk-Benefit Ratio**: Acceptable given potential for improved outcomes and manageable toxicity profile

### 10.4 Data Privacy

- All data de-identified (coded with study ID)
- HIPAA-compliant data management
- Data stored on secure, encrypted servers
- Access limited to authorized study personnel

---

## 11. STUDY ORGANIZATION

### 11.1 Coordinating Center

**[Institution Name]**
- Protocol development
- Centralized randomization
- Data management (REDCap)
- Statistical analysis

### 11.2 Central Laboratories

**RNA-Seq Laboratory**: [Lab name]
- RNA extraction, quality control
- RNA sequencing (Illumina platform)
- Cluster assignment (50-gene classifier)
- VRS calculation
- Turnaround time: 7-10 days

**MRD Core Laboratory**: [Lab name]
- Multiparameter flow cytometry (10-color panel)
- Molecular MRD (NGS, sensitivity 10⁻⁴)

**Pathology Core**: [Lab name]
- Central review of bone marrow morphology
- Cytogenetics interpretation

### 11.3 Participating Sites

**Estimated**: 15-20 sites (US multi-center)

**Site Selection Criteria**:
- Academic medical center with AML expertise
- Enrollment capacity: ≥10 patients per year
- IRB approval and research infrastructure
- Ability to ship samples to central labs

---

## 12. REGULATORY CONSIDERATIONS

### 12.1 IND Application

- **Venetoclax**: FDA-approved for AML (no IND required)
- **Panobinostat + Selumetinib**: Combination requires IND
- **Sponsor**: [Institution] submits IND to FDA
- **IND Number**: [To be assigned]

### 12.2 Protocol Amendments

- Major amendments: FDA and IRB approval required before implementation
- Minor amendments: Notify FDA/IRB within 30 days

### 12.3 Data Monitoring

- Annual IND safety reports to FDA
- Immediate reporting of unexpected SAEs

---

## 13. PUBLICATION PLAN

### 13.1 Primary Manuscript

**Target Journal**: *Blood* or *Journal of Clinical Oncology*

**Primary Outcome**: CR/CRi rate at 2 cycles

**Expected Timeline**: Submit within 6 months of trial completion

### 13.2 Secondary Manuscripts

1. **Cluster-specific efficacy analysis** (by cluster)
2. **VRS validation study**
3. **Cost-effectiveness analysis**
4. **Long-term follow-up** (OS at 2 years)

### 13.3 Presentations

- Abstract submission to ASH Annual Meeting
- Oral presentation at major conferences (ASCO, EHA)

### 13.4 Data Sharing

- De-identified data deposited in public repository (dbGaP) after publication
- Analysis code shared on GitHub

---

## 14. FUNDING AND BUDGET

### 14.1 Estimated Budget

**Total Study Cost**: $4.5 - 6.0 million

**Major Cost Components**:
- **RNA-seq**: $800 × 200 = $160,000
- **MRD testing**: $500 × 400 = $200,000
- **Study drug costs**:
  - Venetoclax: Covered by standard of care
  - Panobinostat: $5,000/cycle × 6 cycles × 55 patients = $1,650,000
  - Selumetinib: $3,000/cycle × 6 cycles × 55 patients = $990,000
- **Coordinating center**: $500,000
- **Central labs**: $300,000
- **Site costs**: $5,000/patient × 200 = $1,000,000
- **Data management**: $200,000
- **Monitoring**: $150,000

### 14.2 Funding Sources

**Target Funders**:
- NIH/NCI (R01 or U01 mechanism)
- Leukemia & Lymphoma Society (TRP or SCOR)
- Industry partnership (Novartis - Panobinostat; AstraZeneca - Selumetinib)

**Grant Submission**: Prepare R01 application, target February 2026 cycle

---

## 15. STUDY TIMELINE

### 15.1 Pre-Enrollment Phase (Months 1-6)

- Month 1-2: IND preparation and submission
- Month 2-3: Site selection and contracting
- Month 3-4: IRB approvals at all sites
- Month 4-5: Central lab setup and validation
- Month 5-6: Investigator meeting, site training

### 15.2 Enrollment Phase (Months 7-30)

- Month 7: First patient enrolled
- Month 12: Interim analysis (n=100)
- Month 30: Last patient enrolled (target: 200 patients)

**Enrollment Rate**: ~8 patients/month across 15 sites

### 15.3 Treatment Phase (Months 7-36)

- Patients treated for up to 6 cycles (24 weeks)
- Primary endpoint assessed after cycle 2 (~8 weeks)

### 15.4 Follow-Up Phase (Months 30-54)

- Minimum 24 months follow-up for all patients
- OS and PFS assessment

### 15.5 Analysis and Publication (Months 54-60)

- Data lock, statistical analysis
- Manuscript preparation and submission
- Conference presentations

**Total Study Duration**: 60 months (5 years)

---

## 16. CONTINGENCY PLANS

### 16.1 Slow Enrollment

**Trigger**: <6 patients/month for 3 consecutive months

**Actions**:
- Add additional sites (target: 5 more sites)
- Expand eligibility (consider age ≥18 if safe)
- Increase site activation and recruitment efforts

### 16.2 High Toxicity

**Trigger**: Treatment-related mortality >15% in either arm

**Actions**:
- DSMB review and recommendation
- Consider dose reductions
- Modify eligibility criteria (exclude highest-risk patients)
- Possible trial suspension

### 16.3 Futility

**Trigger**: Conditional power <20% at interim analysis

**Actions**:
- DSMB recommendation to stop for futility
- Analyze accumulated data
- Publish negative results

---

## 17. CONCLUSIONS

This Phase II randomized trial will prospectively test whether **molecular cluster-guided treatment selection** improves outcomes in adult AML patients. The trial leverages robust preclinical data showing cluster-specific drug sensitivities and provides a pathway for precision medicine implementation.

**Key Innovations**:
1. **First prospective trial** of transcriptomic subtype-guided therapy in AML
2. **Cluster 2 salvage strategy** (Panobinostat + Selumetinib) offers alternative to ineffective Venetoclax
3. **VRS validation** enables rapid clinical implementation
4. **Cost-effectiveness analysis** informs healthcare policy

**Expected Impact**:
- If positive: Establish cluster assignment as standard biomarker for AML treatment selection
- Pathway to Phase III registration trial
- Foundation for expanded precision medicine in AML

---

## PROTOCOL APPROVAL

**Protocol Version**: 1.0
**Date**: 2025-12-09

**Principal Investigator**: _________________ Date: _______

**Sponsor Representative**: _________________ Date: _______

**IRB Chairperson**: _________________ Date: _______

---

## REFERENCES

1. DiNardo CD, et al. Venetoclax combined with decitabine or azacitidine in treatment-naive, elderly patients with acute myeloid leukemia. *Blood*. 2019;133:7-17.
2. Döhner H, et al. Diagnosis and management of AML in adults: 2022 ELN recommendations. *Blood*. 2022;140:1345-1377.
3. Tyner JW, et al. Functional genomic landscape of acute myeloid leukaemia. *Nature*. 2018;562:526-531.
4. [This manuscript - to be added upon publication]

---

**Protocol END**

**Total Pages**: 25
**Document Status**: READY FOR REGULATORY REVIEW
**Next Steps**: IND preparation, grant submission, site selection
