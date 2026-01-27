# Molecular Subtyping and Prognostic Validation in Acute Myeloid Leukemia

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Language-Python-green.svg)](https://www.python.org/)

This repository contains the analysis scripts for a multi-omics study of adult Acute Myeloid Leukemia (AML). The project identifies robust molecular subtypes and establishes their clinical utility for treatment selection, specifically predicting response to **Venetoclax**.

## 🚀 Key Findings

* **Two Robust Subtypes**: Identified two molecular subtypes (k=2) cross-validated in independent cohorts.
* **Venetoclax Sensitivity**: Cluster 1 shows high sensitivity to Venetoclax (p = 2.78×10⁻²⁴).
* **Independent Utility**: Subtypes provide independent predictive value for drug response (+42% R² improvement over mutation-only models).
* **Clinical Tools**: Includes a 50-gene classifier and the **Venetoclax Response Score (VRS)** thresholds.

## 📁 Repository Structure

* `02_Scripts/`: Analysis pipeline.
  * `01_Data_Processing/`: Initial data cleaning and multi-omics integration.
  * `Phase5_DrugValidation/`: Primary drug response analysis and BCL-2 pathway validation.
  * `Phase7_Enhancements/`: Clinical utility tools, VRS development, and salvage therapy analysis.

## 🛠️ Setup & Requirements

### Dependencies

* **R**: `dplyr`, `ggplot2`, `readr`, `tidyr`, `immunedeconv`, `data.table`
* **Python**: `pandas`, `numpy`, `matplotlib`, `seaborn`, `scikit-learn`

### Data Availability

This project utilizes data from the **Beat AML Database**, **TCGA-LAML**, and **TARGET-AML**. Users must download raw data from the respective portals to reproduce the analysis.

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
