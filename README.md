# Molecular Subtyping and Prognostic Validation in Acute Myeloid Leukemia

[![Language-R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![Language-Python](https://img.shields.io/badge/Language-Python-green.svg)](https://www.python.org/)

This repository contains the analysis suite for a multi-omics study of adult Acute Myeloid Leukemia (AML). The project identifies robust molecular subtypes and establishes their clinical utility for treatment selection, specifically predicting response to **Venetoclax**.

## ?? Key Findings
* **Two Robust Subtypes**: Identified two molecular subtypes (k=2) cross-validated in independent cohorts.
* **Venetoclax Predictive Value**: Cluster 1 shows profound hypersensitivity to Venetoclax (p = 2.78×10?²4).
* **Independent Utility**: Subtypes provide independent predictive value (+42% R² improvement) beyond genomic risk classification.
* **Clinical Decision Tool**: Includes a 50-gene classifier and the **Venetoclax Response Score (VRS)**.

## ?? Repository Structure
The project is organized into a clean 01--06 structured hierarchy:

1. **01_Data/**: Raw multi-omics datasets (BeatAML, TCGA, TARGET).
2. **02_Scripts/**: Full analysis pipeline organized by Phase (1--9).
3. **03_Results/**: Processed tables, manifestations, and QC reports.
4. **04_Figures/**: Publication-ready figures and visualizations.
5. **05_Submission/**: Final deliverables hub including Manuscript, Figures, and Supplementary Tables.
6. **06_Revision/**: Tracking for peer review responses and revisions.

## ??? Getting Started
The entire pipeline can be executed via the master orchestrator:
``bash
Rscript RUN_MASTER_PIPELINE.R
``

## ?? License
This project is licensed under the MIT License.
