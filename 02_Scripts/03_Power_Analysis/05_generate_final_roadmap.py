"""
Phase 5: Statistical Power and Feasibility Assessment
Task 5.2 Final: Generate Comprehensive Analysis Roadmap

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
from pathlib import Path
import sys

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Paths
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/03_Power_Analysis")
REPORT_DIR = Path("D:/Projects/Project_AML/05_Reports")

print("=" * 80)
print("GENERATING FINAL COMPREHENSIVE ANALYSIS ROADMAP")
print("=" * 80)
print()

# Load all analysis CSVs
print("Loading analysis specifications...")

# Original roadmap (Tier 1 + partial Tier 2 & 3)
df_original = pd.read_csv(OUTPUT_DIR / "comprehensive_analysis_roadmap.csv")

# Clinical integration (Tier 2 expansion)
df_clinical = pd.read_csv(OUTPUT_DIR / "clinical_integration_analyses.csv")

# Tier 3 advanced analyses
df_tier3 = pd.read_csv(OUTPUT_DIR / "tier3_advanced_analyses.csv")

# Combine all
df_all = pd.concat([df_original, df_clinical, df_tier3], ignore_index=True)

# Sort by analysis_id
df_all = df_all.sort_values('analysis_id')

print(f"✓ Loaded {len(df_all)} analyses")
print()

# ============================================================================
# GENERATE COMPREHENSIVE MARKDOWN ROADMAP
# ============================================================================

print("Generating comprehensive roadmap document...")

report_lines = []

# Header
report_lines.append("# AML Multi-Omics Integration: Comprehensive Analysis Roadmap")
report_lines.append("")
report_lines.append("**Project:** Beat AML Multi-Omics Integration Study")
report_lines.append("")
report_lines.append("**Generated:** 2025-10-02")
report_lines.append("")
report_lines.append("**Version:** 1.0")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# Table of Contents
report_lines.append("## Table of Contents")
report_lines.append("")
report_lines.append("1. [Executive Summary](#executive-summary)")
report_lines.append("2. [Cohort Overview](#cohort-overview)")
report_lines.append("3. [TIER 1: Core Multi-Omics Analyses](#tier-1-core-multi-omics-analyses-highest-priority)")
report_lines.append("4. [TIER 2: Clinical Integration Analyses](#tier-2-clinical-integration-analyses-high-priority)")
report_lines.append("5. [TIER 3: Advanced Integrative Analyses](#tier-3-advanced-integrative-analyses-exploratory)")
report_lines.append("6. [Implementation Timeline](#implementation-timeline)")
report_lines.append("7. [Feasibility Summary](#feasibility-summary)")
report_lines.append("8. [Resource Requirements](#resource-requirements)")
report_lines.append("9. [Success Metrics](#success-metrics)")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# Executive Summary
report_lines.append("## Executive Summary")
report_lines.append("")

tier1_analyses = df_all[df_all['tier'] == 1]
tier2_analyses = df_all[df_all['tier'] == 2]
tier3_analyses = df_all[df_all['tier'] == 3]

total_feasible = (df_all['feasibility'] == 'YES').sum()

report_lines.append(f"This roadmap outlines **{len(df_all)} comprehensive analyses** for the Beat AML multi-omics integration project, organized into 3 priority tiers:")
report_lines.append("")
report_lines.append(f"- **TIER 1 (Highest Priority):** {len(tier1_analyses)} core multi-omics analyses")
report_lines.append(f"- **TIER 2 (High Priority):** {len(tier2_analyses)} clinical integration analyses")
report_lines.append(f"- **TIER 3 (Exploratory):** {len(tier3_analyses)} advanced integrative analyses")
report_lines.append("")
report_lines.append(f"### Overall Feasibility")
report_lines.append("")
report_lines.append(f"**{total_feasible}/{len(df_all)} analyses are fully feasible** with the available Beat AML cohort.")
report_lines.append("")
report_lines.append("All analyses have been power-calculated and are supported by robust sample sizes across multiple data type combinations.")
report_lines.append("")

# Key statistics
mean_power_all = df_all['statistical_power'].mean()
report_lines.append(f"### Key Statistics")
report_lines.append("")
report_lines.append(f"- **Mean Statistical Power:** {mean_power_all:.2f}")
report_lines.append(f"- **Gold Standard Cohort (all 4 data types):** n = 478")
report_lines.append(f"- **Expression + Drug:** n = 494")
report_lines.append(f"- **Survival Data:** n = 942 (565 events, 60%)")
report_lines.append("")

# Timeline estimate
report_lines.append("### Estimated Timeline")
report_lines.append("")
report_lines.append("| Tier | Sequential | Parallel |")
report_lines.append("|------|------------|----------|")
report_lines.append("| **Tier 1** | 8-12 weeks | 6-8 weeks |")
report_lines.append("| **Tier 2** | 7-10 weeks | 5-7 weeks |")
report_lines.append("| **Tier 3** | 11-16 weeks | 8-10 weeks |")
report_lines.append("| **TOTAL** | **26-38 weeks** | **19-25 weeks** |")
report_lines.append("")
report_lines.append("With parallelization and adequate resources, the project can be completed in **~5-6 months**.")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# Cohort Overview
report_lines.append("## Cohort Overview")
report_lines.append("")
report_lines.append("### Available Sample Sizes by Data Type Combination")
report_lines.append("")
report_lines.append("| Data Type Combination | Sample Size | Key Analyses |")
report_lines.append("|----------------------|-------------|--------------|")
report_lines.append("| Expression only | n = 707 | Molecular subtyping |")
report_lines.append("| Mutations only | n = 871 | Mutation landscape |")
report_lines.append("| Drug Response only | n = 603 | Drug screening |")
report_lines.append("| Clinical only | n = 934 | Demographics |")
report_lines.append("| Expression + Mutations | n = 615 | Mutation-expression integration |")
report_lines.append("| Expression + Drug | n = 494 | Drug prediction |")
report_lines.append("| Mutations + Drug | n = 583 | Mutation-drug associations |")
report_lines.append("| Expression + Clinical | n = 671 | Clinical-molecular correlation |")
report_lines.append("| **Gold Standard (all 4)** | **n = 478** | **Multi-omics integration** |")
report_lines.append("| Survival data | n = 942 (565 events) | Prognostic analyses |")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# TIER 1
report_lines.append("## TIER 1: Core Multi-Omics Analyses (HIGHEST PRIORITY)")
report_lines.append("")
report_lines.append("These foundational analyses establish the molecular landscape and provide essential insights for all downstream analyses.")
report_lines.append("")

for _, analysis in tier1_analyses.iterrows():
    report_lines.append(f"### Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
    report_lines.append("")

    # Status badge
    if analysis['feasibility'] == 'YES':
        status = "✓ **FEASIBLE**"
    else:
        status = "⚠ **LIMITED**"
    report_lines.append(f"**Status:** {status} | **Power:** {analysis['statistical_power']:.2f} | **Timeline:** {analysis['timeline_weeks']}")
    report_lines.append("")

    report_lines.append(f"**Goal:** {analysis['goal']}")
    report_lines.append("")
    report_lines.append(f"**Methods:** {analysis['methods']}")
    report_lines.append("")

    report_lines.append(f"**Sample Size:**")
    report_lines.append(f"- Required: n ≥ {analysis['required_n']}")
    report_lines.append(f"- Available: n = {analysis['available_n']}")
    report_lines.append("")

    report_lines.append(f"**Data Requirements:** {analysis['data_requirements']}")
    report_lines.append("")

    report_lines.append(f"**Justification:** {analysis['justification']}")
    report_lines.append("")

    report_lines.append(f"**Scripts Location:** `{analysis['scripts_location']}`")
    report_lines.append("")

    report_lines.append("**Deliverables:**")
    for deliverable in analysis['deliverables'].split('; '):
        report_lines.append(f"- {deliverable}")
    report_lines.append("")

    # Add key comparisons if available
    if pd.notna(analysis.get('key_comparisons')):
        report_lines.append(f"**Key Comparisons:** {analysis['key_comparisons']}")
        report_lines.append("")

    report_lines.append("---")
    report_lines.append("")

# TIER 2
report_lines.append("## TIER 2: Clinical Integration Analyses (HIGH PRIORITY)")
report_lines.append("")
report_lines.append("These analyses integrate molecular features with clinical outcomes to identify prognostic biomarkers and validate clinical utility.")
report_lines.append("")

for _, analysis in tier2_analyses.iterrows():
    report_lines.append(f"### Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
    report_lines.append("")

    # Status badge
    if analysis['feasibility'] == 'YES':
        status = "✓ **FEASIBLE**"
    else:
        status = "⚠ **LIMITED**"
    report_lines.append(f"**Status:** {status} | **Power:** {analysis['statistical_power']:.2f} | **Timeline:** {analysis['timeline_weeks']}")
    report_lines.append("")

    report_lines.append(f"**Goal:** {analysis['goal']}")
    report_lines.append("")
    report_lines.append(f"**Methods:** {analysis['methods']}")
    report_lines.append("")

    report_lines.append(f"**Sample Size:**")
    report_lines.append(f"- Required: n ≥ {analysis['required_n']}")
    report_lines.append(f"- Available: n = {analysis['available_n']}")
    report_lines.append("")

    report_lines.append(f"**Data Requirements:** {analysis['data_requirements']}")
    report_lines.append("")

    report_lines.append(f"**Justification:** {analysis['justification']}")
    report_lines.append("")

    # Add stratifications if available
    if pd.notna(analysis.get('stratifications')):
        report_lines.append(f"**Stratifications:**")
        strats = analysis['stratifications'].split('; ')
        for strat in strats[:5]:  # Limit to first 5
            report_lines.append(f"- {strat}")
        if len(strats) > 5:
            report_lines.append(f"- ... and {len(strats)-5} more")
        report_lines.append("")

    # Add associations if available
    if pd.notna(analysis.get('associations')):
        report_lines.append(f"**Associations to Test:**")
        for assoc in analysis['associations'].split('; '):
            report_lines.append(f"- {assoc}")
        report_lines.append("")

    # Add feature categories if available
    if pd.notna(analysis.get('feature_categories')):
        report_lines.append(f"**Feature Categories:**")
        for cat in analysis['feature_categories'].split('; '):
            report_lines.append(f"- {cat}")
        report_lines.append("")

    # Add key pairs if available
    if pd.notna(analysis.get('key_pairs')):
        report_lines.append(f"**Key Associations:**")
        pairs = analysis['key_pairs'].split('; ')
        for pair in pairs[:5]:
            report_lines.append(f"- {pair}")
        if len(pairs) > 5:
            report_lines.append(f"- ... and {len(pairs)-5} more")
        report_lines.append("")

    report_lines.append(f"**Scripts Location:** `{analysis['scripts_location']}`")
    report_lines.append("")

    report_lines.append("**Deliverables:**")
    for deliverable in analysis['deliverables'].split('; '):
        report_lines.append(f"- {deliverable}")
    report_lines.append("")

    report_lines.append("---")
    report_lines.append("")

# TIER 3
report_lines.append("## TIER 3: Advanced Integrative Analyses (EXPLORATORY)")
report_lines.append("")
report_lines.append("These exploratory analyses provide advanced insights and develop translational tools for precision medicine.")
report_lines.append("")

for _, analysis in tier3_analyses.iterrows():
    report_lines.append(f"### Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
    report_lines.append("")

    # Status badge
    if analysis['feasibility'] == 'YES':
        status = "✓ **FEASIBLE**"
    else:
        status = "⚠ **LIMITED**"
    report_lines.append(f"**Status:** {status} | **Power:** {analysis['statistical_power']:.2f} | **Timeline:** {analysis['timeline_weeks']}")
    report_lines.append("")

    report_lines.append(f"**Goal:** {analysis['goal']}")
    report_lines.append("")
    report_lines.append(f"**Methods:** {analysis['methods']}")
    report_lines.append("")

    report_lines.append(f"**Sample Size:**")
    report_lines.append(f"- Required: n ≥ {analysis['required_n']}")
    report_lines.append(f"- Available: n = {analysis['available_n']}")
    report_lines.append("")

    report_lines.append(f"**Data Requirements:** {analysis['data_requirements']}")
    report_lines.append("")

    report_lines.append(f"**Justification:** {analysis['justification']}")
    report_lines.append("")

    # Add network approaches if available
    if pd.notna(analysis.get('network_approaches')):
        report_lines.append(f"**Network Approaches:**")
        for approach in analysis['network_approaches'].split('; '):
            report_lines.append(f"- {approach}")
        report_lines.append("")

    # Add analysis components if available
    if pd.notna(analysis.get('analysis_components')):
        report_lines.append(f"**Analysis Components:**")
        for comp in analysis['analysis_components'].split('; '):
            report_lines.append(f"- {comp}")
        report_lines.append("")

    # Add framework components if available
    if pd.notna(analysis.get('framework_components')):
        report_lines.append(f"**Framework Components:**")
        for comp in analysis['framework_components'].split('; '):
            report_lines.append(f"- {comp}")
        report_lines.append("")

    # Add dependencies if available
    if pd.notna(analysis.get('dependencies')):
        report_lines.append(f"**Dependencies:** {analysis['dependencies']}")
        report_lines.append("")

    # Add caveats if available
    if pd.notna(analysis.get('caveat')):
        report_lines.append(f"⚠ **Caveat:** {analysis['caveat']}")
        report_lines.append("")

    report_lines.append(f"**Scripts Location:** `{analysis['scripts_location']}`")
    report_lines.append("")

    report_lines.append("**Deliverables:**")
    for deliverable in analysis['deliverables'].split('; '):
        report_lines.append(f"- {deliverable}")
    report_lines.append("")

    report_lines.append("---")
    report_lines.append("")

# Implementation Timeline
report_lines.append("## Implementation Timeline")
report_lines.append("")
report_lines.append("### Recommended Execution Order")
report_lines.append("")
report_lines.append("#### Phase 1: Foundation (Weeks 1-8)")
report_lines.append("**Objective:** Establish molecular landscape")
report_lines.append("")
report_lines.append("- Analysis 1.1: Molecular Subtyping (parallel)")
report_lines.append("- Analysis 1.2: Mutation Landscape (parallel)")
report_lines.append("- Analysis 1.3: Mutation-Expression Integration (after 1.1, 1.2)")
report_lines.append("")
report_lines.append("#### Phase 2: Prediction & Clinical Integration (Weeks 9-16)")
report_lines.append("**Objective:** Build predictive models and clinical correlations")
report_lines.append("")
report_lines.append("- Analysis 1.4: Drug Response Prediction (after Phase 1)")
report_lines.append("- Analysis 2.1: Survival by Molecular Features (parallel)")
report_lines.append("- Analysis 2.4: Survival Analysis (after 1.1)")
report_lines.append("- Analysis 2.5: Clinical-Molecular Correlation (parallel)")
report_lines.append("")
report_lines.append("#### Phase 3: Advanced Integration (Weeks 17-24)")
report_lines.append("**Objective:** Advanced analyses and prognostic modeling")
report_lines.append("")
report_lines.append("- Analysis 2.2: Mutation-Drug Associations (parallel)")
report_lines.append("- Analysis 2.3: Network Analysis (parallel)")
report_lines.append("- Analysis 2.6: Integrated Prognostic Model (after 2.4, 2.5)")
report_lines.append("- Analysis 3.1: Subtype-Specific Drug Sensitivities (after 1.1)")
report_lines.append("")
report_lines.append("#### Phase 4: Exploratory & Translational (Weeks 17-30)")
report_lines.append("**Objective:** Advanced discoveries and clinical tools")
report_lines.append("")
report_lines.append("- Analysis 3.2: Clinical-Molecular Correlations (parallel)")
report_lines.append("- Analysis 3.3: Drug Combination Synergy (exploratory)")
report_lines.append("- Analysis 3.4: Multi-Omics Network Analysis (parallel)")
report_lines.append("- Analysis 3.5: Drug Mechanism Discovery (parallel)")
report_lines.append("- Analysis 3.6: Personalized Treatment Framework (final)")
report_lines.append("")

# Feasibility Summary Table
report_lines.append("## Feasibility Summary")
report_lines.append("")
report_lines.append("### All Analyses")
report_lines.append("")
report_lines.append("| Analysis ID | Analysis Name | Sample Size | Power | Feasibility |")
report_lines.append("|-------------|---------------|-------------|-------|-------------|")
for _, analysis in df_all.iterrows():
    feas_icon = '✓' if analysis['feasibility'] == 'YES' else '⚠'
    report_lines.append(f"| {analysis['analysis_id']} | {analysis['analysis_name']} | {analysis['available_n']} | {analysis['statistical_power']:.2f} | {feas_icon} {analysis['feasibility']} |")
report_lines.append("")

# Summary by tier
report_lines.append("### Summary by Tier")
report_lines.append("")
for tier_num in [1, 2, 3]:
    tier_data = df_all[df_all['tier'] == tier_num]
    tier_feasible = (tier_data['feasibility'] == 'YES').sum()
    tier_mean_power = tier_data['statistical_power'].mean()
    tier_name = ['', 'TIER 1', 'TIER 2', 'TIER 3'][tier_num]
    report_lines.append(f"**{tier_name}:** {tier_feasible}/{len(tier_data)} feasible (mean power: {tier_mean_power:.2f})")
report_lines.append("")
report_lines.append(f"**OVERALL:** {total_feasible}/{len(df_all)} analyses are fully feasible")
report_lines.append("")

# Resource Requirements
report_lines.append("## Resource Requirements")
report_lines.append("")
report_lines.append("### Personnel")
report_lines.append("")
report_lines.append("- **Bioinformatician (Lead):** 1 FTE for entire project")
report_lines.append("- **Computational Biologist:** 1 FTE for advanced analyses (Tier 2-3)")
report_lines.append("- **Statistician:** 0.5 FTE for power analyses and validation")
report_lines.append("- **Clinical Collaborator:** 0.25 FTE for clinical interpretation")
report_lines.append("")
report_lines.append("### Computational Resources")
report_lines.append("")
report_lines.append("- **High-performance computing:** Required for network analysis and ML models")
report_lines.append("- **Storage:** ~500 GB for processed data and results")
report_lines.append("- **RAM:** ≥64 GB recommended for expression matrix operations")
report_lines.append("")
report_lines.append("### Software & Tools")
report_lines.append("")
report_lines.append("- **R packages:** DESeq2, limma, WGCNA, survival, caret")
report_lines.append("- **Python libraries:** scikit-learn, pandas, seaborn, lifelines")
report_lines.append("- **Pathway databases:** MSigDB, KEGG, Reactome")
report_lines.append("- **Visualization:** Cytoscape, R Shiny (for personalized framework)")
report_lines.append("")

# Success Metrics
report_lines.append("## Success Metrics")
report_lines.append("")
report_lines.append("### Tier 1 Completion Criteria")
report_lines.append("")
report_lines.append("- [ ] Molecular subtypes identified and validated")
report_lines.append("- [ ] Comprehensive mutation landscape characterized")
report_lines.append("- [ ] Mutation-expression associations identified (FDR < 0.05)")
report_lines.append("- [ ] Drug prediction models achieve R² > 0.3")
report_lines.append("")
report_lines.append("### Tier 2 Completion Criteria")
report_lines.append("")
report_lines.append("- [ ] Prognostic molecular features identified (p < 0.05)")
report_lines.append("- [ ] Clinical-molecular associations validated")
report_lines.append("- [ ] Integrated prognostic model C-index > 0.7")
report_lines.append("- [ ] Mutation-drug associations confirmed")
report_lines.append("")
report_lines.append("### Tier 3 Completion Criteria")
report_lines.append("")
report_lines.append("- [ ] Multi-omics networks constructed and validated")
report_lines.append("- [ ] Drug mechanisms elucidated")
report_lines.append("- [ ] Personalized treatment framework deployed")
report_lines.append("")
report_lines.append("### Publication Targets")
report_lines.append("")
report_lines.append("- **High-impact journals:** Nature Communications, Cell Reports, Blood")
report_lines.append("- **Minimum publications:** 3-5 primary research articles")
report_lines.append("- **Data release:** GEO/dbGaP deposition for reproducibility")
report_lines.append("")

# Footer
report_lines.append("---")
report_lines.append("")
report_lines.append("## Document Information")
report_lines.append("")
report_lines.append("**Version:** 1.0")
report_lines.append("")
report_lines.append("**Last Updated:** 2025-10-02")
report_lines.append("")
report_lines.append("**Contact:** AML Multi-Omics Project Team")
report_lines.append("")
report_lines.append("**Related Documents:**")
report_lines.append("- Statistical Power Analysis Report")
report_lines.append("- Clinical Integration Analyses Specification")
report_lines.append("- Analysis Cohort Definitions")
report_lines.append("")

# Save markdown
output_file = REPORT_DIR / "Analysis_Roadmap.md"
with open(output_file, 'w', encoding='utf-8') as f:
    f.write('\n'.join(report_lines))

print(f"✓ Saved: {output_file}")
print()

# Also save consolidated CSV
consolidated_file = OUTPUT_DIR / "consolidated_analysis_roadmap.csv"
df_all.to_csv(consolidated_file, index=False)
print(f"✓ Saved: {consolidated_file}")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("COMPREHENSIVE ANALYSIS ROADMAP COMPLETE")
print("=" * 80)
print()
print(f"TOTAL ANALYSES: {len(df_all)}")
print(f"  Tier 1: {len(tier1_analyses)}")
print(f"  Tier 2: {len(tier2_analyses)}")
print(f"  Tier 3: {len(tier3_analyses)}")
print()
print(f"FEASIBILITY: {total_feasible}/{len(df_all)} analyses fully feasible")
print(f"MEAN POWER: {mean_power_all:.2f}")
print()
print("OUTPUTS:")
print(f"  1. {output_file.name} (in 05_Reports/)")
print(f"  2. {consolidated_file.name}")
print()
print("=" * 80)
