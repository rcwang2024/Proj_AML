# BeatAML Multi-Omics Data Dictionary

**Version:** 1.0

**Last Updated:** 2025-10-02

---

## 1. Gene Expression Data

### File Information
- **File Name:** `beataml_expression.txt`
- **Location:** `01_Data/BeatAML_Downloaded_Data/`
- **File Size:** ~269 MB
- **Format:** Tab-delimited text file

### Dimensions
- **Samples (columns):** 707 patients
- **Genes (rows):** 22,843 protein-coding genes

### Gene Identifiers
- **Type:** HGNC gene symbols
- **Column Name:** First column labeled as gene identifier
- **Examples:** TP53, FLT3, NPM1, DNMT3A
- **Unique:** Yes, one row per gene

### Sample Identifiers
- **Format:** BA####R (where # = digits, R = RNA)
- **Examples:** BA0001R, BA0534R, BA1234R
- **Unique Identifiers:** Yes, one column per sample
- **Mapping:** Can be linked to DNA samples (BA####D) via clinical file

### Data Values
- **Type:** Continuous numeric
- **Transformation:** Log2-transformed normalized expression
- **Units:** Log2(FPKM+1) or Log2(TPM+1)
- **Typical Range:** 0-20
- **Interpretation:**
  - 0: Gene not expressed
  - <5: Low expression
  - 5-10: Moderate expression
  - >10: High expression
- **Normalization:** Pre-normalized by Beat AML consortium
- **Platform:** Illumina RNA-seq

### Processing Pipeline
1. RNA extraction from patient samples
2. Library preparation and sequencing
3. Read alignment to human genome (GRCh38/hg38)
4. Quantification (FPKM or TPM)
5. Log2 transformation
6. Normalization across samples

### Missing Data
- **Encoding:** NA, blank, or 0
- **Frequency:** Rare (< 0.1% of values)
- **Handling:** Genes with excessive missing data excluded

### Quality Control Applied
- **By Beat AML:** Yes (samples with low quality excluded)
- **Outlier Detection:** Performed in this project
- **Batch Effects:** Detected (centerID variable), correction recommended

### Source and Version
- **Source:** Beat AML consortium
- **Data Version:** Public release
- **dbGaP Accession:** phs001657
- **Citation:** Tyner et al., Nature 2018; Bottomly et al., Cancer Cell 2022

---

## 2. Drug Response Data

### File Information
- **File Name:** `beataml_drug_auc.txt`
- **Location:** `01_Data/BeatAML_Downloaded_Data/`
- **File Size:** ~3.5 MB
- **Format:** Tab-delimited text (long format)

### File Structure
- **Format:** Long format (one row per sample-drug pair)
- **Total Rows:** 63,395 sample-drug combinations
- **Columns:** 3 primary columns

### Column Definitions

**1. lab_id / sample_id**
- **Type:** Character string
- **Format:** Sample identifier
- **Examples:** Patient/sample codes
- **Links to:** Clinical and expression data

**2. inhibitor**
- **Type:** Character string
- **Format:** Drug/compound name
- **Total Unique Drugs:** 166
- **Examples:** Sorafenib, Venetoclax, Gilteritinib, Azacitidine
- **Drug Classes:**
  - Kinase inhibitors (FLT3, JAK, BTK, etc.)
  - Epigenetic modifiers (IDH, DNMT, BET)
  - BCL2 inhibitors
  - Chemotherapy agents
  - Other targeted therapies

**3. auc (Area Under Curve)**
- **Type:** Continuous numeric
- **Definition:** Area under the dose-response curve
- **Units:** Arbitrary units
- **Typical Range:** 0-1000+
- **Interpretation:**
  - Lower AUC = Greater drug sensitivity
  - Higher AUC = Drug resistance
  - AUC < 100: Highly sensitive
  - AUC 100-500: Moderate sensitivity
  - AUC > 500: Resistant
- **Calculation:** Integrated from dose-response experiments

### Drug Information

**Drug Name Type:** Generic or compound names

**Drug Families/Classes:**
- FLT3 inhibitors: Sorafenib, Gilteritinib, Quizartinib
- IDH inhibitors: Ivosidenib (IDH1), Enasidenib (IDH2)
- BCL2 inhibitors: Venetoclax
- Hypomethylating agents: Azacitidine, Decitabine
- JAK inhibitors: Ruxolitinib
- MEK inhibitors: Trametinib
- Others: Various targeted and chemotherapy agents

**Target Information:**
- Refer to DrugBank or DGIdb for detailed target information
- Most drugs target kinases, epigenetic modifiers, or apoptosis pathways

### Missing Data
- **Encoding:** Absent rows (not all drugs tested on all samples)
- **Average Coverage:** ~101 drugs per sample
- **Range:** 50-166 drugs per sample
- **Reason:** Not all samples screened against all drugs

### Quality Metrics
- **Quality Flags:** Not explicitly provided
- **Extreme Values:** Checked (no invalid values detected)
- **Replicate Information:** Not provided in public data

### Units and Ranges
- **AUC Units:** Arbitrary (relative)
- **Valid Range:** 0 to ~3000 (99% < 1000)
- **Outliers:** Rare

### Source and Version
- **Source:** Beat AML high-throughput drug screening
- **Platform:** Ex vivo drug sensitivity assay
- **Citation:** Tyner et al., Nature 2018

---

## 3. Clinical Data

### File Information
- **File Name:** `beataml_clinical.xlsx`
- **Location:** `01_Data/BeatAML_Downloaded_Data/`
- **File Size:** ~1.2 MB
- **Format:** Excel spreadsheet (.xlsx)

### Overview
- **Samples:** 942 patients
- **Variables:** 95 clinical features
- **Data Types:** Mixed (numeric, categorical, dates)

### Key Variables with Definitions

#### Sample Identifiers
**dbgap_rnaseq_sample**
- RNA sample ID (BA####R format)
- Links to expression data

**dbgap_dnaseq_sample**
- DNA sample ID (BA####D format)
- Links to mutation data

#### Demographics

**ageAtDiagnosis**
- Age at AML diagnosis
- Type: Numeric
- Units: Years
- Range: 18-90+
- Completeness: 97.2%

**consensus_sex**
- Biological sex
- Type: Categorical
- Values: Male, Female, Unknown
- Completeness: 99.9%

#### Survival Variables

**overallSurvival**
- Overall survival time
- Type: Numeric
- Units: Days
- Range: 0-3000+
- Completeness: 100%
- Note: Time from diagnosis to death or last follow-up

**vitalStatus**
- Vital status at last follow-up
- Type: Categorical
- Values: Dead, Alive
- Completeness: 100%

#### Disease Characteristics

**isDenovo**
- De novo vs secondary AML
- Type: Boolean
- Values: TRUE (de novo), FALSE (secondary/therapy-related)
- Completeness: 94.7%

**%.Blasts.in.BM**
- Blast percentage in bone marrow
- Type: Numeric
- Units: Percentage (0-100)
- Completeness: 74.6%

#### Laboratory Values

**WBC**
- White blood cell count
- Type: Numeric
- Units: × 10⁹/L
- Completeness: Variable

#### Treatment Information

**priorTreatment**
- Prior therapy status
- Type: Categorical/text
- Completeness: Variable

### Missing Data Encoding
- **NA:** Standard R/Excel missing value
- **Empty cells:** Also indicate missing
- **-1 or -999:** Sometimes used for missing numeric values
- **"Unknown" or "Not Available":** Text missing indicators

### Data Types and Valid Ranges

**Numeric Variables:**
- Age: 0-120 years
- Survival: 0-5000 days
- Lab values: Biologically plausible ranges

**Categorical Variables:**
- Predefined levels (Male/Female, Dead/Alive, etc.)
- Check for typos or inconsistent coding

**Date Variables:**
- Format: Various (YYYY-MM-DD or Excel date format)
- May require parsing

### ELN Risk Definitions

**2017_ELN_risk** (if available)
- European LeukemiaNet risk stratification
- Values: Favorable, Intermediate, Adverse
- Based on cytogenetics and molecular markers

### Survival Time Units
- **Overall Survival:** Days from diagnosis
- **Event:** Death (1 = dead, 0 = alive/censored)

### Vital Status Encoding
- **Dead:** Patient died
- **Alive:** Patient alive at last follow-up (censored)

### Source
- **Source:** Beat AML clinical database
- **Citation:** Tyner et al., Nature 2018

---

## 4. Mutation Data

### File Information
- **File Name:** `beataml_mutations.txt`
- **Location:** `01_Data/BeatAML_Downloaded_Data/`
- **File Size:** ~5.8 MB
- **Format:** Tab-delimited text (MAF-like format)

### File Structure
- **Format:** Long format (one row per mutation)
- **Total Rows:** 11,721 somatic mutations
- **Samples:** 871 patients
- **Genes:** 615 unique genes mutated

### Key Column Definitions

**dbgap_sample_id**
- Sample identifier
- Format: BA####D (D = DNA)
- Links to clinical data

**symbol**
- Gene symbol (HGNC)
- Examples: TP53, FLT3, NPM1, DNMT3A
- Type: Character

**variant_classification**
- Mutation consequence type
- Values:
  - **Missense_Mutation:** Amino acid change
  - **Nonsense_Mutation:** Premature stop codon
  - **Frame_Shift_Del:** Frameshift deletion
  - **Frame_Shift_Ins:** Frameshift insertion
  - **In_Frame_Del:** In-frame deletion
  - **In_Frame_Ins:** In-frame insertion
  - **Splice_Site:** Affects splicing
  - **Silent:** Synonymous (no AA change)

**chromosome**
- Chromosome location
- Values: 1-22, X, Y, MT

**start_position / end_position**
- Genomic coordinates (GRCh37 or GRCh38)
- Type: Integer
- Units: Base pairs

**reference_allele**
- Reference base(s) from genome
- Type: Character (A, C, G, T, or longer)

**tumor_seq_allele1 / tumor_seq_allele2**
- Observed alleles in tumor
- Usually tumor_seq_allele2 is the variant

**t_vaf (Tumor Variant Allele Frequency)**
- Proportion of sequencing reads with variant
- Type: Numeric (0.0 to 1.0)
- Definition: (variant reads) / (total reads)
- **Interpretation:**
  - VAF ~0.5: Heterozygous mutation in majority of cells
  - VAF ~1.0: Homozygous or loss of heterozygosity
  - VAF <0.1: Subclonal mutation (present in minority)
  - VAF 0.4-0.6: Typical for clonal heterozygous mutations

### Mutation Type Definitions

**Missense Mutation**
- Single nucleotide change resulting in different amino acid
- Example: TP53 p.R273H

**Nonsense Mutation**
- Change to stop codon, truncating protein
- Example: TP53 p.Q192*

**Frameshift**
- Insertion/deletion not multiple of 3
- Shifts reading frame, usually truncates protein

**In-frame Indel**
- Insertion/deletion multiple of 3
- Changes amino acids but preserves reading frame

**Splice Site**
- Affects canonical splice sites (GT-AG)
- May cause exon skipping or intron retention

### VAF (Variant Allele Frequency) Definition

**Definition:** Fraction of sequencing reads supporting the variant allele

**Range:** 0.0 (no variant reads) to 1.0 (all reads have variant)

**Typical Distribution:**
- Peak at 0.5 for heterozygous clonal mutations
- Lower VAF for subclonal mutations
- Higher VAF (>0.5) for copy-neutral LOH or subclonal expansion

**Quality Considerations:**
- Low VAF (<0.05) may be sequencing errors
- Very high VAF (>0.95) may indicate homozygous variant

### Quality Metrics

**Sequencing Coverage:** Variable across samples and positions

**Filtering Applied:**
- Germline variants removed (tumor-normal comparison or population filter)
- Low-quality calls filtered by Beat AML
- Likely passenger mutations may be included

### Consequence Predictions

**SIFT, PolyPhen:** May be available in extended annotations

**ClinVar:** Pathogenicity annotations (if available)

### Driver Mutation Definitions

**Driver Genes:** Genes recurrently mutated in AML and known to contribute to disease

**Examples:**
- FLT3, NPM1, DNMT3A, IDH1, IDH2, TET2, TP53, NRAS, KRAS

**Identification:**
- Based on published AML driver gene lists
- High frequency in cohort
- Known functional impact

### Filtering Applied

**By Beat AML:**
- Germline variants removed
- Low-quality variants filtered
- Likely artifacts excluded

**Recommended Additional Filtering:**
- VAF > 0.05 (remove low-frequency variants)
- Exclude synonymous mutations for driver analyses

### Source and Version
- **Source:** Beat AML targeted or whole-exome sequencing
- **Genome Build:** GRCh37/hg19 or GRCh38/hg38
- **Pipeline:** GATK MuTect2 or similar
- **Citation:** Tyner et al., Nature 2018

---

## General Notes

### Cross-Dataset Linking

**Sample ID Mapping:**
- Expression: BA####R
- Mutations: BA####D
- Clinical: Contains both RNA and DNA IDs
- Use clinical file as bridge to link datasets

**Unified Sample IDs:**
- See `master_sample_id_mapping.csv` for unified identifiers

### Data Access

**Source:** dbGaP controlled access

**Accession:** phs001657

**Access Requirements:** Data Use Agreement, IRB approval

### Data Processing

**Scripts Location:** `02_Scripts/`

**Processed Data:** `03_Results/01_Processed_Data/`

**Quality Control:** `03_Results/02_QC_Reports/`

### Questions or Issues

**Contact:** AML Multi-Omics Project Team

**Documentation:** See `05_Reports/` for comprehensive reports

---

**END OF DATA DICTIONARY**
