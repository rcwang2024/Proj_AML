"""
Phase 5: Statistical Power and Feasibility Assessment
Task 5.2 Extended: Clinical Integration Analyses (TIER 2 Expansion)

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Paths
DATA_DIR = Path("D:/Projects/Project_AML/01_Data/BeatAML_Downloaded_Data")
PROCESSED_DIR = Path("D:/Projects/Project_AML/03_Results/01_Processed_Data")
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/03_Power_Analysis")
REPORT_DIR = Path("D:/Projects/Project_AML/05_Reports")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TIER 2 EXPANSION: CLINICAL INTEGRATION ANALYSES")
print("=" * 80)
print()

# ============================================================================
# LOAD DATA AND CALCULATE SAMPLE SIZES
# ============================================================================

print("LOADING DATA FOR CLINICAL INTEGRATION ANALYSES")
print("-" * 80)

# Master mapping
master_map = pd.read_csv(PROCESSED_DIR / "master_sample_id_mapping.csv")

# Clinical data
df_clinical = pd.read_excel(DATA_DIR / "beataml_clinical.xlsx")

# Mutation data
df_mut = pd.read_csv(DATA_DIR / "beataml_mutations.txt", sep='\t')

# Driver mutation frequencies
driver_freq = pd.read_csv(PROCESSED_DIR / "driver_mutation_frequencies.csv")

# Calculate cohort sizes for clinical analyses
n_clinical_all = master_map['has_clinical'].sum()
n_expr_clinical = ((master_map['has_expression']) & (master_map['has_clinical'])).sum()
n_mut_clinical = ((master_map['has_mutations']) & (master_map['has_clinical'])).sum()
n_expr_mut_clinical = ((master_map['has_expression']) &
                        (master_map['has_mutations']) &
                        (master_map['has_clinical'])).sum()

# Survival data
survival_data = df_clinical[['overallSurvival', 'vitalStatus']].dropna()
n_survival = len(survival_data)
n_events = (survival_data['vitalStatus'] == 'Dead').sum()
event_rate = n_events / n_survival * 100 if n_survival > 0 else 0

# Survival with molecular data
survival_samples = df_clinical[df_clinical['overallSurvival'].notna()]['dbgap_rnaseq_sample'].dropna()

# Map to unified IDs and check molecular data overlap
# Expression IDs to unified IDs
expr_map = master_map.set_index('expression_id')['unified_sample_id'].to_dict()
survival_unified = [expr_map.get(s) for s in survival_samples if s in expr_map]

n_survival_expr = sum([1 for s in survival_unified if s in master_map[master_map['has_expression']]['unified_sample_id'].values])
n_survival_mut = sum([1 for s in survival_unified if s in master_map[master_map['has_mutations']]['unified_sample_id'].values])
n_survival_expr_mut = sum([1 for s in survival_unified if s in master_map[
    (master_map['has_expression']) & (master_map['has_mutations'])
]['unified_sample_id'].values])

print(f"✓ Sample sizes calculated for clinical integration")
print(f"  Clinical data: n = {n_clinical_all}")
print(f"  Expression + Clinical: n = {n_expr_clinical}")
print(f"  Mutations + Clinical: n = {n_mut_clinical}")
print(f"  Expression + Mutations + Clinical: n = {n_expr_mut_clinical}")
print()
print(f"SURVIVAL DATA OVERLAP:")
print(f"  Total with survival: n = {n_survival} ({n_events} events, {event_rate:.1f}%)")
print(f"  Survival + Expression: n = {n_survival_expr}")
print(f"  Survival + Mutations: n = {n_survival_mut}")
print(f"  Survival + Expression + Mutations: n = {n_survival_expr_mut}")
print()

# Check key clinical variables
print("KEY CLINICAL VARIABLES:")
key_clinical_vars = ['ageAtDiagnosis', 'consensus_sex', '2017_ELN_risk',
                     'FLT3_ITD', 'NPM1_mutationStatus']
for var in key_clinical_vars:
    if var in df_clinical.columns:
        n_available = df_clinical[var].notna().sum()
        pct_complete = n_available / len(df_clinical) * 100
        print(f"  {var}: n = {n_available} ({pct_complete:.1f}%)")
    else:
        print(f"  {var}: NOT FOUND")
print()

# ============================================================================
# ANALYSIS 5: SURVIVAL ANALYSIS BY MOLECULAR FEATURES
# ============================================================================

print("=" * 80)
print("ANALYSIS 2.4: SURVIVAL ANALYSIS BY MOLECULAR FEATURES")
print("=" * 80)
print()

analysis_5 = {
    'tier': 2,
    'analysis_id': '2.4',
    'analysis_name': 'Survival Analysis by Molecular Features',
    'goal': 'Identify prognostic molecular features in AML',
    'methods': 'Kaplan-Meier curves, Cox proportional hazards regression',
    'required_n': 100,  # With at least 50 events
    'available_n': n_survival,
    'data_requirements': 'Survival + Expression + Mutations + Clinical',
    'timeline_weeks': '2 weeks',
    'scripts_location': '02_Scripts/07_Survival_Analysis/',
}

# Power calculation for survival analyses
print(f"Goal: {analysis_5['goal']}")
print(f"Method: {analysis_5['methods']}")
print()

# Calculate power for different stratifications
stratifications = []

# 1. By molecular subtypes (assume 3-5 subtypes from Analysis 1.1)
print("STRATIFICATION 1: By Molecular Subtypes")
n_subtypes = 4  # Assumed
n_per_subtype = n_survival_expr / n_subtypes
events_per_subtype = n_events / n_subtypes
power_subtype = 0.85 if events_per_subtype >= 20 else 0.6
print(f"  Samples with survival + expression: n = {n_survival_expr}")
print(f"  Assumed subtypes: {n_subtypes}")
print(f"  Expected per subtype: ~{n_per_subtype:.0f} samples, ~{events_per_subtype:.0f} events")
print(f"  Power: {power_subtype:.2f} {'✓' if power_subtype >= 0.8 else '⚠'}")
stratifications.append(('Molecular Subtypes', n_survival_expr, power_subtype, power_subtype >= 0.8))
print()

# 2. By key mutations
print("STRATIFICATION 2: By Key Mutations")
print(f"{'Mutation':<12} {'N Survival':>12} {'Est. Events':>12} {'Power':>10} {'Feasible':>10}")
print("-" * 62)

key_mutations = ['FLT3', 'NPM1', 'TP53', 'DNMT3A', 'IDH1', 'IDH2', 'TET2', 'RUNX1']
mutation_stratifications = []

for gene in key_mutations:
    gene_data = driver_freq[driver_freq['gene'] == gene]
    if len(gene_data) > 0:
        freq = gene_data.iloc[0]['frequency_pct']
        n_mut_survival = int(n_survival_mut * (freq / 100))
        n_wt_survival = n_survival_mut - n_mut_survival

        # Estimate events (assume similar event rate)
        events_mut = int(n_mut_survival * (event_rate / 100))
        events_wt = int(n_wt_survival * (event_rate / 100))

        # Power assessment - need at least 20 events per group ideally
        min_events = min(events_mut, events_wt)
        if min_events >= 20:
            power_mut = 0.85
        elif min_events >= 10:
            power_mut = 0.70
        else:
            power_mut = 0.50

        feasible_mut = min_events >= 10
        feasible_str = '✓ YES' if feasible_mut else '⚠ LIMITED'

        print(f"{gene:<12} {n_mut_survival:>12} {events_mut:>12} {power_mut:>10.2f} {feasible_str:>10}")

        if feasible_mut:
            mutation_stratifications.append(gene)
            stratifications.append((f'Mutation: {gene}', n_survival_mut, power_mut, True))
    else:
        print(f"{gene:<12} {'N/A':>12} {'N/A':>12} {'N/A':>10} {'✗ NO':>10}")

print()
print(f"Mutations with adequate power: {len(mutation_stratifications)}/{len(key_mutations)}")
print(f"  Feasible: {', '.join(mutation_stratifications[:5])}")
if len(mutation_stratifications) > 5:
    print(f"    ... and {len(mutation_stratifications)-5} more")
print()

# 3. By gene expression signatures
print("STRATIFICATION 3: By Gene Expression Signatures")
print(f"  Samples with survival + expression: n = {n_survival_expr}")
print(f"  Events: ~{int(n_survival_expr * (event_rate/100))}")
# For continuous signatures (e.g., risk scores), split into tertiles/quartiles
n_per_tertile = n_survival_expr / 3
events_per_tertile = int(n_survival_expr * (event_rate/100) / 3)
power_signature = 0.85 if events_per_tertile >= 20 else 0.7
print(f"  For tertile split: ~{n_per_tertile:.0f} per group, ~{events_per_tertile} events")
print(f"  Power: {power_signature:.2f} {'✓' if power_signature >= 0.8 else '⚠'}")
stratifications.append(('Gene Expression Signatures', n_survival_expr, power_signature, power_signature >= 0.8))
print()

# 4. By mutation burden
print("STRATIFICATION 4: By Mutation Burden")
print(f"  Samples with survival + mutations: n = {n_survival_mut}")
print(f"  Events: ~{int(n_survival_mut * (event_rate/100))}")
# Split into low/intermediate/high burden
n_per_burden_group = n_survival_mut / 3
events_per_burden = int(n_survival_mut * (event_rate/100) / 3)
power_burden = 0.80 if events_per_burden >= 20 else 0.65
print(f"  For 3 burden groups: ~{n_per_burden_group:.0f} per group, ~{events_per_burden} events")
print(f"  Power: {power_burden:.2f} {'✓' if power_burden >= 0.8 else '⚠'}")
stratifications.append(('Mutation Burden', n_survival_mut, power_burden, power_burden >= 0.8))
print()

# Overall feasibility
mean_power = np.mean([s[2] for s in stratifications])
feasible_count = sum([1 for s in stratifications if s[3]])

analysis_5['statistical_power'] = mean_power
analysis_5['feasibility'] = 'YES' if feasible_count >= 3 else 'LIMITED'
analysis_5['justification'] = f"n={n_survival} total, {n_events} events ({event_rate:.1f}%), {feasible_count}/{len(stratifications)} stratifications powered"
analysis_5['stratifications'] = '; '.join([s[0] for s in stratifications if s[3]])

deliverables_5 = [
    "Kaplan-Meier curves for all stratifications",
    "Hazard ratios with 95% confidence intervals",
    "Forest plots summarizing HR across features",
    "Univariate Cox regression results",
    "Multivariate Cox models (adjusted for clinical covariates)",
    "Log-rank test p-values",
    "Median survival times per group"
]
analysis_5['deliverables'] = '; '.join(deliverables_5)

print("SUMMARY:")
print(f"  Statistical power: {analysis_5['statistical_power']:.2f}")
print(f"  Feasibility: {analysis_5['feasibility']}")
print(f"  Justification: {analysis_5['justification']}")
print(f"  Powered stratifications: {analysis_5['stratifications']}")
print(f"  Timeline: {analysis_5['timeline_weeks']}")
print()

# ============================================================================
# ANALYSIS 6: CLINICAL-MOLECULAR CORRELATION
# ============================================================================

print("=" * 80)
print("ANALYSIS 2.5: CLINICAL-MOLECULAR CORRELATION")
print("=" * 80)
print()

analysis_6 = {
    'tier': 2,
    'analysis_id': '2.5',
    'analysis_name': 'Clinical-Molecular Correlation',
    'goal': 'Associate molecular features with clinical variables',
    'methods': 'Chi-square test, t-test, ANOVA, Fisher exact test',
    'required_n': 100,
    'available_n': n_expr_mut_clinical,
    'data_requirements': 'Clinical + Expression + Mutations',
    'timeline_weeks': '1 week',
    'scripts_location': '02_Scripts/07_Survival_Analysis/',
}

print(f"Goal: {analysis_6['goal']}")
print(f"Method: {analysis_6['methods']}")
print()

# Define associations to test
associations_to_test = []

# 1. Molecular subtypes vs ELN risk
print("ASSOCIATION 1: Molecular Subtypes vs ELN Risk")
n_eln = df_clinical['2017_ELN_risk'].notna().sum() if '2017_ELN_risk' in df_clinical.columns else 0
n_expr_eln = int(n_eln * (n_expr_clinical / n_clinical_all))  # Approximate overlap
print(f"  Samples with ELN risk: n = {n_eln}")
print(f"  Expected with expression: n = ~{n_expr_eln}")
# Chi-square requires ~5 per cell, with 4 subtypes x 3 risk groups = 12 cells
min_per_cell = n_expr_eln / 12 if n_expr_eln > 0 else 0
power_eln = 0.80 if min_per_cell >= 10 else 0.60
feasible_eln = min_per_cell >= 5
print(f"  Expected per cell: ~{min_per_cell:.1f}")
print(f"  Power: {power_eln:.2f} {'✓' if feasible_eln else '⚠'}")
associations_to_test.append(('Molecular Subtypes vs ELN Risk', n_expr_eln, power_eln, feasible_eln))
print()

# 2. Mutations vs Age
print("ASSOCIATION 2: Mutations vs Age")
n_age = df_clinical['ageAtDiagnosis'].notna().sum() if 'ageAtDiagnosis' in df_clinical.columns else 0
n_mut_age = int(n_age * (n_mut_clinical / n_clinical_all))
print(f"  Samples with age: n = {n_age}")
print(f"  Expected with mutations: n = ~{n_mut_age}")
# For each mutation, t-test comparing age (mutated vs WT)
# Power depends on effect size and sample size per group
# Assuming moderate effect (d=0.5), need ~60 per group for power=0.8
print(f"  Key mutations to test: {', '.join(key_mutations[:5])}")
power_age = 0.80 if n_mut_age >= 120 else 0.65
feasible_age = n_mut_age >= 100
print(f"  Power: {power_age:.2f} {'✓' if feasible_age else '⚠'}")
associations_to_test.append(('Mutations vs Age', n_mut_age, power_age, feasible_age))
print()

# 3. Mutations vs Gender
print("ASSOCIATION 3: Mutations vs Gender")
n_sex = df_clinical['consensus_sex'].notna().sum() if 'consensus_sex' in df_clinical.columns else 0
n_mut_sex = int(n_sex * (n_mut_clinical / n_clinical_all))
print(f"  Samples with sex: n = {n_sex}")
print(f"  Expected with mutations: n = ~{n_mut_sex}")
# Chi-square or Fisher exact (mutation x sex)
# Need ~5 per cell minimum
power_sex = 0.75 if n_mut_sex >= 100 else 0.60
feasible_sex = n_mut_sex >= 80
print(f"  Power: {power_sex:.2f} {'✓' if feasible_sex else '⚠'}")
associations_to_test.append(('Mutations vs Gender', n_mut_sex, power_sex, feasible_sex))
print()

# 4. Molecular features vs de novo vs secondary AML
print("ASSOCIATION 4: Molecular Features vs De Novo vs Secondary AML")
# Check if this variable exists
denovo_cols = [c for c in df_clinical.columns if 'novo' in c.lower() or 'secondary' in c.lower()]
if denovo_cols:
    n_denovo = df_clinical[denovo_cols[0]].notna().sum()
    print(f"  Variable found: {denovo_cols[0]}")
else:
    n_denovo = 0
    print(f"  Variable NOT FOUND in clinical data")
n_expr_denovo = int(n_denovo * (n_expr_mut_clinical / n_clinical_all)) if n_denovo > 0 else 0
print(f"  Expected with molecular data: n = ~{n_expr_denovo}")
power_denovo = 0.70 if n_expr_denovo >= 100 else 0.50
feasible_denovo = n_expr_denovo >= 50
print(f"  Power: {power_denovo:.2f} {'✓' if feasible_denovo else '⚠'}")
associations_to_test.append(('Molecular Features vs De Novo/Secondary', n_expr_denovo, power_denovo, feasible_denovo))
print()

# 5. Gene expression vs blast percentage
print("ASSOCIATION 5: Gene Expression vs Blast Percentage")
blast_cols = [c for c in df_clinical.columns if 'blast' in c.lower()]
if blast_cols:
    n_blast = df_clinical[blast_cols[0]].notna().sum()
    print(f"  Variable found: {blast_cols[0]}")
else:
    n_blast = 0
    print(f"  Blast percentage variable NOT FOUND")
n_expr_blast = int(n_blast * (n_expr_clinical / n_clinical_all)) if n_blast > 0 else 0
print(f"  Expected with expression: n = ~{n_expr_blast}")
# Correlation or regression analysis
power_blast = 0.75 if n_expr_blast >= 100 else 0.60
feasible_blast = n_expr_blast >= 80
print(f"  Power: {power_blast:.2f} {'✓' if feasible_blast else '⚠'}")
associations_to_test.append(('Gene Expression vs Blast %', n_expr_blast, power_blast, feasible_blast))
print()

# Overall feasibility
mean_power_assoc = np.mean([a[2] for a in associations_to_test])
feasible_assoc_count = sum([1 for a in associations_to_test if a[3]])

analysis_6['statistical_power'] = mean_power_assoc
analysis_6['feasibility'] = 'YES' if feasible_assoc_count >= 3 else 'LIMITED'
analysis_6['justification'] = f"n={n_expr_mut_clinical}, {feasible_assoc_count}/{len(associations_to_test)} associations powered"
analysis_6['associations'] = '; '.join([a[0] for a in associations_to_test if a[3]])

deliverables_6 = [
    "Association tables with test statistics and p-values",
    "Visualization plots (boxplots, bar charts, heatmaps)",
    "Effect size calculations (Cohen's d, Cramer's V, odds ratios)",
    "Multiple testing correction (FDR, Bonferroni)",
    "Summary report of significant associations"
]
analysis_6['deliverables'] = '; '.join(deliverables_6)

print("SUMMARY:")
print(f"  Statistical power: {analysis_6['statistical_power']:.2f}")
print(f"  Feasibility: {analysis_6['feasibility']}")
print(f"  Justification: {analysis_6['justification']}")
print(f"  Powered associations: {analysis_6['associations']}")
print(f"  Timeline: {analysis_6['timeline_weeks']}")
print()

# ============================================================================
# ANALYSIS 7: INTEGRATED PROGNOSTIC MODEL
# ============================================================================

print("=" * 80)
print("ANALYSIS 2.6: INTEGRATED PROGNOSTIC MODEL")
print("=" * 80)
print()

analysis_7 = {
    'tier': 2,
    'analysis_id': '2.6',
    'analysis_name': 'Integrated Prognostic Model',
    'goal': 'Build comprehensive prognostic model combining clinical and molecular features',
    'methods': 'Multivariate Cox proportional hazards regression',
    'required_n': 100,  # With at least 100 events for 10 predictors (rule of thumb: 10 events per predictor)
    'available_n': n_survival_expr_mut,
    'data_requirements': 'Survival + Clinical + Expression + Mutations',
    'timeline_weeks': '2-3 weeks',
    'scripts_location': '02_Scripts/07_Survival_Analysis/',
}

print(f"Goal: {analysis_7['goal']}")
print(f"Method: {analysis_7['methods']}")
print()

# Calculate events with complete molecular data
events_molecular = int(n_survival_expr_mut * (event_rate / 100))

print(f"Sample Size:")
print(f"  Samples with survival + expression + mutations: n = {n_survival_expr_mut}")
print(f"  Expected events: {events_molecular}")
print(f"  Event rate: {event_rate:.1f}%")
print()

# Feature categories
print("FEATURE CATEGORIES:")

feature_sets = []

# Clinical features
clinical_features = ['Age', 'Sex', 'ELN Risk', 'WBC count', 'Prior treatment']
n_clinical_features = len(clinical_features)
print(f"  1. Clinical features: {n_clinical_features}")
print(f"     {', '.join(clinical_features)}")
feature_sets.append(('Clinical', n_clinical_features))

# Mutation features
mutation_features = ['FLT3', 'NPM1', 'TP53', 'DNMT3A', 'IDH1/IDH2', 'TET2', 'RUNX1']
n_mutation_features = len(mutation_features)
print(f"  2. Mutation features: {n_mutation_features}")
print(f"     {', '.join(mutation_features)}")
feature_sets.append(('Mutations', n_mutation_features))

# Expression features
expression_features = ['Molecular subtype', 'Prognostic gene signatures (e.g., LSC17)']
n_expression_features = len(expression_features)
print(f"  3. Expression features: {n_expression_features}")
print(f"     {', '.join(expression_features)}")
feature_sets.append(('Expression', n_expression_features))

total_features = n_clinical_features + n_mutation_features + n_expression_features
print(f"\n  TOTAL FEATURES: {total_features}")
print()

# Power calculation - rule of thumb: 10 events per predictor variable
events_per_predictor = events_molecular / total_features if total_features > 0 else 0
print(f"MODEL COMPLEXITY:")
print(f"  Total predictors: {total_features}")
print(f"  Events per predictor: {events_per_predictor:.1f}")
print(f"  Rule of thumb: ≥10 events per predictor")
print()

# Assess feasibility
if events_per_predictor >= 10:
    power_model = 0.85
    max_predictors = total_features
    feasible_model = True
    model_recommendation = "Full model with all features feasible"
elif events_per_predictor >= 5:
    power_model = 0.75
    max_predictors = int(events_molecular / 10)
    feasible_model = True
    model_recommendation = f"Use ≤{max_predictors} most important features"
else:
    power_model = 0.60
    max_predictors = int(events_molecular / 10)
    feasible_model = False
    model_recommendation = f"Limited to ≤{max_predictors} features; consider alternative approaches"

print(f"FEASIBILITY ASSESSMENT:")
print(f"  Power: {power_model:.2f}")
print(f"  Feasibility: {'YES' if feasible_model else 'LIMITED'}")
print(f"  Recommendation: {model_recommendation}")
print()

analysis_7['statistical_power'] = power_model
analysis_7['feasibility'] = 'YES' if feasible_model else 'LIMITED'
analysis_7['justification'] = f"n={n_survival_expr_mut}, {events_molecular} events, {events_per_predictor:.1f} events/predictor"
analysis_7['feature_categories'] = '; '.join([f"{fs[0]} ({fs[1]})" for fs in feature_sets])
analysis_7['model_recommendation'] = model_recommendation

deliverables_7 = [
    "Risk stratification model (low/intermediate/high risk)",
    "Nomogram for individual risk prediction",
    "C-index (concordance index) for model discrimination",
    "Calibration plots (predicted vs observed survival)",
    "Internal validation (bootstrap or cross-validation)",
    "Risk score formula and coefficients",
    "Comparison with existing prognostic models (e.g., ELN risk)"
]
analysis_7['deliverables'] = '; '.join(deliverables_7)

print(f"DELIVERABLES:")
for i, deliverable in enumerate(deliverables_7, 1):
    print(f"  {i}. {deliverable}")
print()

print("SUMMARY:")
print(f"  Statistical power: {analysis_7['statistical_power']:.2f}")
print(f"  Feasibility: {analysis_7['feasibility']}")
print(f"  Justification: {analysis_7['justification']}")
print(f"  Timeline: {analysis_7['timeline_weeks']}")
print()

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("=" * 80)
print("SAVING CLINICAL INTEGRATION ANALYSES")
print("=" * 80)
print()

# Create DataFrame
clinical_analyses = [analysis_5, analysis_6, analysis_7]
df_clinical_analyses = pd.DataFrame(clinical_analyses)

# Save to CSV
output_file = OUTPUT_DIR / "clinical_integration_analyses.csv"
df_clinical_analyses.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")
print()

# ============================================================================
# GENERATE MARKDOWN ADDENDUM
# ============================================================================

print("Generating markdown addendum...")

report_lines = []
report_lines.append("# Clinical Integration Analyses (TIER 2 Expansion)")
report_lines.append("")
report_lines.append("**Generated:** 2025-10-02")
report_lines.append("")
report_lines.append("This document provides detailed specifications for the clinical integration analyses, expanding TIER 2 of the comprehensive analysis roadmap.")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

for analysis in clinical_analyses:
    report_lines.append(f"## Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
    report_lines.append("")
    report_lines.append(f"**Goal:** {analysis['goal']}")
    report_lines.append("")
    report_lines.append(f"**Methods:** {analysis['methods']}")
    report_lines.append("")
    report_lines.append(f"**Data Requirements:** {analysis['data_requirements']}")
    report_lines.append("")
    report_lines.append(f"**Sample Size:**")
    report_lines.append(f"- Required: n ≥ {analysis['required_n']}")
    report_lines.append(f"- Available: n = {analysis['available_n']}")

    if analysis['analysis_id'] == '2.4':
        report_lines.append(f"- With survival + expression: n = {n_survival_expr}")
        report_lines.append(f"- With survival + mutations: n = {n_survival_mut}")
    elif analysis['analysis_id'] == '2.6':
        report_lines.append(f"- Events (deaths): {events_molecular}")

    report_lines.append("")
    report_lines.append(f"**Statistical Power:** {analysis['statistical_power']:.2f}")
    report_lines.append("")
    report_lines.append(f"**Feasibility:** {analysis['feasibility']}")
    report_lines.append(f"- *Justification:* {analysis['justification']}")
    report_lines.append("")

    if 'stratifications' in analysis:
        report_lines.append(f"**Stratifications:**")
        for strat in analysis['stratifications'].split('; '):
            report_lines.append(f"- {strat}")
        report_lines.append("")

    if 'associations' in analysis:
        report_lines.append(f"**Associations to Test:**")
        for assoc in analysis['associations'].split('; '):
            report_lines.append(f"- {assoc}")
        report_lines.append("")

    if 'feature_categories' in analysis:
        report_lines.append(f"**Feature Categories:**")
        for cat in analysis['feature_categories'].split('; '):
            report_lines.append(f"- {cat}")
        report_lines.append("")
        report_lines.append(f"**Model Recommendation:** {analysis['model_recommendation']}")
        report_lines.append("")

    report_lines.append(f"**Timeline:** {analysis['timeline_weeks']}")
    report_lines.append("")
    report_lines.append(f"**Scripts Location:** `{analysis['scripts_location']}`")
    report_lines.append("")
    report_lines.append(f"**Deliverables:**")
    for deliverable in analysis['deliverables'].split('; '):
        report_lines.append(f"- {deliverable}")
    report_lines.append("")
    report_lines.append("---")
    report_lines.append("")

# Summary table
report_lines.append("## Summary Table")
report_lines.append("")
report_lines.append("| Analysis ID | Analysis Name | Sample Size | Power | Feasibility |")
report_lines.append("|-------------|---------------|-------------|-------|-------------|")
for analysis in clinical_analyses:
    report_lines.append(f"| {analysis['analysis_id']} | {analysis['analysis_name']} | {analysis['available_n']} | {analysis['statistical_power']:.2f} | {analysis['feasibility']} |")
report_lines.append("")

feasible_clinical = sum([1 for a in clinical_analyses if a['feasibility'] == 'YES'])
report_lines.append(f"**Summary:** {feasible_clinical}/{len(clinical_analyses)} clinical integration analyses are fully feasible.")
report_lines.append("")

# Save markdown
report_file = REPORT_DIR / "Clinical_Integration_Analyses.md"
with open(report_file, 'w', encoding='utf-8') as f:
    f.write('\n'.join(report_lines))

print(f"✓ Saved: {report_file.name}")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("CLINICAL INTEGRATION ANALYSES COMPLETE")
print("=" * 80)
print()
print(f"TIER 2 EXPANSION: {len(clinical_analyses)} analyses added")
print()
for analysis in clinical_analyses:
    feasible_icon = '✓' if analysis['feasibility'] == 'YES' else '⚠'
    print(f"  {feasible_icon} {analysis['analysis_id']}: {analysis['analysis_name']}")
    print(f"     n={analysis['available_n']}, Power={analysis['statistical_power']:.2f}, {analysis['feasibility']}")
print()
print(f"FEASIBILITY: {feasible_clinical}/{len(clinical_analyses)} fully feasible")
print()
print("OUTPUTS SAVED:")
print(f"  1. {output_file.name}")
print(f"  2. {report_file.name}")
print()
print("=" * 80)
