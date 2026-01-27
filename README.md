# Molecular Subtyping and Prognostic Validation in Acute Myeloid Leukemia

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Language-Python-green.svg)](https://www.python.org/)

This repository contains the analysis scripts and documentation for a comprehensive multi-omics study of adult Acute Myeloid Leukemia (AML). The project identifies robust molecular subtypes and establishes their clinical utility for treatment selection, specifically predicting response to **Venetoclax** and identifying salvage therapies for resistant cases.

## 🚀 Key Findings

* **Two Robust Subtypes**: Identified two molecular subtypes (k=2, consensus=0.957) across 2,535 patients in 3 independent cohorts.
* **Venetoclax Sensitivity**: Cluster 1 is extraordinarily sensitive to Venetoclax (p = 2.78×10⁻²⁴), with a **161% improvement** in predictive accuracy over mutation-only models.
* **Independent Utility**: While subtypes are not independent prognostic markers for survival, they provide **independent predictive value for 19/20 tested drugs** (mean +42% R² improvement).
* **Clinical Utility Tools**: Includes a 50-gene classifier for subtype assignment and the **Venetoclax Response Score (VRS)** with defined clinical thresholds.
* **Salvage Therapies**: Identified **Panobinostat** and **Selumetinib** as optimal treatments for Cluster 2 (Venetoclax-resistant) patients.

## 📁 Repository Structure

* `02_Scripts/`: Core analysis pipeline.
  * `Phase1-5/`: Primary subtyping and validation.
  * `Phase7_Enhancements/`: Clinical utility tools, VRS development, and salvage therapy analysis.
* `04_Figures/`: Publication-ready visualizations.
* `05_Manuscript/`: LaTeX source and drafts, including a full clinical trial protocol (**CLUSTER-AML Trial**).
* `06_Documentation/`: Detailed methodology and data inventory.

## 🛠️ Setup & Requirements

### Dependencies

- **R**: `dplyr`, `ggplot2`, `readr`, `tidyr`, `immunedeconv`, `data.table`
* **Python**: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`

### Data Availability

This project utilizes data from the **Beat AML Database**, **TCGA-LAML**, and **TARGET-AML**. Users must download raw data from the respective portals (e.g., GDC) to reproduce the analysis.

## 📊 Analysis Pipeline

1. **Data Integration**: Merges RNA-seq, somatic mutations, and drug response AUC data.
2. **Molecular Subtyping**: Consensus clustering to identify k=2 optimal clusters.
3. **Validation**: PH-free survival analysis and meta-analysis across adult and pediatric cohorts.
4. **Drug Screening**: Hierarchical regression to test independent predictive value of subtypes.
5. **Clinical Translation**: Generation of VRS thresholds and salvage therapy pathways.

## 📄 Documentation

For a detailed walkthrough of the implementation and results, see:
* [Complete Project Summary](06_Documentation/COMPLETE_PROJECT_SUMMARY_ALL_PHASES_V3.md)
* [Clinical Enhancement Summary](06_Documentation/COMPLETE_ENHANCEMENT_SUMMARY.md)
* [Clinical Trial Protocol](05_Manuscript/CLINICAL_TRIAL_PROTOCOL.md)

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

For questions or collaborations, please contact the project lead via the institutional details provided in the manuscript.
