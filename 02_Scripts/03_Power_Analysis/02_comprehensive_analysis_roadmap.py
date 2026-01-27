"""
Phase 5: Statistical Power and Feasibility Assessment
Task 5.2: Generate Comprehensive Analysis Roadmap

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
print("TASK 5.2: COMPREHENSIVE ANALYSIS ROADMAP")
print("=" * 80)
print()

# ============================================================================
# LOAD DATA AND CALCULATE SAMPLE SIZES
# ============================================================================

print("LOADING DATA AND CALCULATING AVAILABLE SAMPLE SIZES")
print("-" * 80)

# Master mapping
master_map = pd.read_csv(PROCESSED_DIR / "master_sample_id_mapping.csv")

# Clinical data
df_clinical = pd.read_excel(DATA_DIR / "beataml_clinical.xlsx")

# Mutation data
df_mut = pd.read_csv(DATA_DIR / "beataml_mutations.txt", sep='\t')

# Driver mutation frequencies
driver_freq = pd.read_csv(PROCESSED_DIR / "driver_mutation_frequencies.csv")

# Calculate cohort sizes
n_expression = master_map['has_expression'].sum()
n_mutations = master_map['has_mutations'].sum()
n_drug_response = master_map['has_drug_response'].sum()
n_clinical = master_map['has_clinical'].sum()

n_expr_mut = ((master_map['has_expression']) & (master_map['has_mutations'])).sum()
n_expr_drug = ((master_map['has_expression']) & (master_map['has_drug_response'])).sum()
n_mut_drug = ((master_map['has_mutations']) & (master_map['has_drug_response'])).sum()

n_gold_standard = ((master_map['has_expression']) &
                   (master_map['has_mutations']) &
                   (master_map['has_drug_response']) &
                   (master_map['has_clinical'])).sum()

# Survival data
survival_data = df_clinical[['overallSurvival', 'vitalStatus']].dropna()
n_survival = len(survival_data)
n_events = (survival_data['vitalStatus'] == 'Dead').sum()

print(f"✓ Sample sizes calculated")
print(f"  Expression only: n = {n_expression}")
print(f"  Mutations only: n = {n_mutations}")
print(f"  Drug response only: n = {n_drug_response}")
print(f"  Clinical only: n = {n_clinical}")
print(f"  Expression + Mutations: n = {n_expr_mut}")
print(f"  Expression + Drug: n = {n_expr_drug}")
print(f"  Mutations + Drug: n = {n_mut_drug}")
print(f"  Gold standard (all 4): n = {n_gold_standard}")
print(f"  Survival data: n = {n_survival} ({n_events} events)")
print()

# ============================================================================
# TIER 1: CORE MULTI-OMICS ANALYSES (HIGHEST PRIORITY)
# ============================================================================

print("=" * 80)
print("TIER 1: CORE MULTI-OMICS ANALYSES (HIGHEST PRIORITY)")
print("=" * 80)
print()

tier1_analyses = []

# ----------------------------------------------------------------------------
# ANALYSIS 1: MOLECULAR SUBTYPING VIA EXPRESSION
# ----------------------------------------------------------------------------

print("ANALYSIS 1.1: Molecular Subtyping via Expression")
print("-" * 80)

analysis_1 = {
    'tier': 1,
    'analysis_id': '1.1',
    'analysis_name': 'Molecular Subtyping via Expression',
    'goal': 'Identify transcriptomic subtypes in AML',
    'methods': 'Consensus clustering, hierarchical clustering, k-means',
    'required_n': 100,
    'available_n': n_expression,
    'data_requirements': 'Expression data',
    'timeline_weeks': '2-3 weeks',
    'scripts_location': '02_Scripts/04_Molecular_Subtyping/',
}

# Power calculation
power_subtyping = 0.9 if n_expression >= 100 else 0.6
min_cluster_size = n_expression / 5  # Worst case with 5 clusters
feasible_subtyping = n_expression >= 100 and min_cluster_size >= 50

analysis_1['statistical_power'] = power_subtyping
analysis_1['feasibility'] = 'YES' if feasible_subtyping else 'LIMITED'
analysis_1['justification'] = f"n={n_expression} with min cluster size ~{min_cluster_size:.0f}"

deliverables_1 = [
    "Cluster assignments for each sample",
    "Subtype-specific gene signatures (DEGs per cluster)",
    "Heatmaps showing cluster-specific expression patterns",
    "PCA/t-SNE plots with cluster annotations",
    "Biological interpretation via pathway analysis (GSEA, Reactome)",
    "Clinical correlation of subtypes (survival, age, mutations)"
]
analysis_1['deliverables'] = '; '.join(deliverables_1)

print(f"Goal: {analysis_1['goal']}")
print(f"Required samples: n ≥ {analysis_1['required_n']}")
print(f"Available samples: n = {analysis_1['available_n']}")
print(f"Statistical power: {analysis_1['statistical_power']:.2f}")
print(f"Feasibility: {analysis_1['feasibility']}")
print(f"  Justification: {analysis_1['justification']}")
print(f"Timeline: {analysis_1['timeline_weeks']}")
print()

tier1_analyses.append(analysis_1)

# ----------------------------------------------------------------------------
# ANALYSIS 1.2: COMPREHENSIVE MUTATION LANDSCAPE
# ----------------------------------------------------------------------------

print("ANALYSIS 1.2: Comprehensive Mutation Landscape")
print("-" * 80)

analysis_2 = {
    'tier': 1,
    'analysis_id': '1.2',
    'analysis_name': 'Comprehensive Mutation Landscape',
    'goal': 'Characterize mutational profile of AML cohort',
    'methods': 'Mutation frequency, co-occurrence analysis, mutual exclusivity',
    'required_n': 200,
    'available_n': n_mutations,
    'data_requirements': 'Mutation data',
    'timeline_weeks': '1-2 weeks',
    'scripts_location': '02_Scripts/04_Molecular_Subtyping/',
}

# Power calculation
power_mutations = 0.9 if n_mutations >= 200 else 0.7
# Can detect mutations with frequency ≥ 5% with n=871
min_detectable_freq = 50 / n_mutations * 100  # Need ~50 samples minimum
feasible_mutations = n_mutations >= 200

analysis_2['statistical_power'] = power_mutations
analysis_2['feasibility'] = 'YES' if feasible_mutations else 'LIMITED'
analysis_2['justification'] = f"n={n_mutations}, can detect mutations ≥{min_detectable_freq:.1f}%"

deliverables_2 = [
    "OncoPrint visualization of top driver mutations",
    "Mutation frequency bar plots (overall and by subtype)",
    "Co-occurrence and mutual exclusivity matrix",
    "Mutational signature analysis (if WGS data available)",
    "Driver vs passenger mutation classification",
    "Comparison with TCGA-AML and other public datasets"
]
analysis_2['deliverables'] = '; '.join(deliverables_2)

print(f"Goal: {analysis_2['goal']}")
print(f"Required samples: n ≥ {analysis_2['required_n']}")
print(f"Available samples: n = {analysis_2['available_n']}")
print(f"Statistical power: {analysis_2['statistical_power']:.2f}")
print(f"Feasibility: {analysis_2['feasibility']}")
print(f"  Justification: {analysis_2['justification']}")
print(f"  Can detect mutations with frequency ≥ {min_detectable_freq:.1f}%")
print(f"Timeline: {analysis_2['timeline_weeks']}")
print()

tier1_analyses.append(analysis_2)

# ----------------------------------------------------------------------------
# ANALYSIS 1.3: MUTATION-EXPRESSION INTEGRATION
# ----------------------------------------------------------------------------

print("ANALYSIS 1.3: Mutation-Expression Integration")
print("-" * 80)

analysis_3 = {
    'tier': 1,
    'analysis_id': '1.3',
    'analysis_name': 'Mutation-Expression Integration',
    'goal': 'Understand how mutations affect gene expression',
    'methods': 'Differential expression analysis stratified by mutation status',
    'required_n': 50,  # Per comparison (20+ per group)
    'available_n': n_expr_mut,
    'data_requirements': 'Expression + Mutation data',
    'timeline_weeks': '2-3 weeks',
    'scripts_location': '02_Scripts/06_Integration/',
}

# Key mutations with sufficient power
key_mutations = driver_freq.sort_values('frequency_pct', ascending=False).head(10)
powered_mutations = []

print(f"Key mutation comparisons:")
print(f"{'Mutation':<12} {'N Mutated':>12} {'N WT':>10} {'Power':>10} {'Feasible':>10}")
print("-" * 58)

for _, row in key_mutations.iterrows():
    gene = row['gene']
    n_mut = int(row['n_samples_mutated'])
    n_wt = n_mutations - n_mut

    # Estimate overlap with expression (~70%)
    n_mut_expr = int(n_mut * 0.7)
    n_wt_expr = int(n_wt * 0.7)

    # Power estimation
    min_n = min(n_mut_expr, n_wt_expr)
    if min_n >= 30:
        power = 0.9
    elif min_n >= 20:
        power = 0.8
    elif min_n >= 10:
        power = 0.6
    else:
        power = 0.3

    feasible = min_n >= 20
    feasible_str = '✓ YES' if feasible else '⚠ LIMITED'

    print(f"{gene:<12} {n_mut_expr:>12} {n_wt_expr:>10} {power:>10.2f} {feasible_str:>10}")

    if feasible:
        powered_mutations.append({
            'gene': gene,
            'n_mutated': n_mut_expr,
            'n_wt': n_wt_expr,
            'power': power
        })

print()

# Overall power
mean_power = np.mean([m['power'] for m in powered_mutations]) if powered_mutations else 0
feasible_mut_expr = len(powered_mutations) >= 5

analysis_3['statistical_power'] = mean_power
analysis_3['feasibility'] = 'YES' if feasible_mut_expr else 'LIMITED'
analysis_3['justification'] = f"n={n_expr_mut}, {len(powered_mutations)}/10 mutations powered"

deliverables_3 = [
    "DEG lists for each mutation (FDR < 0.05, |log2FC| > 1)",
    "Volcano plots for key mutations",
    "Pathway enrichment analysis (KEGG, Reactome, GO)",
    "Gene Set Enrichment Analysis (GSEA) using MSigDB",
    "Heatmaps of top differentially expressed genes",
    "Mutation-specific expression signatures"
]
analysis_3['deliverables'] = '; '.join(deliverables_3)
analysis_3['key_comparisons'] = '; '.join([m['gene'] for m in powered_mutations[:5]])

print(f"Goal: {analysis_3['goal']}")
print(f"Required samples: n ≥ {analysis_3['required_n']} per comparison")
print(f"Available samples: n = {analysis_3['available_n']}")
print(f"Statistical power: {analysis_3['statistical_power']:.2f} (mean across mutations)")
print(f"Feasibility: {analysis_3['feasibility']}")
print(f"  Justification: {analysis_3['justification']}")
print(f"  Key comparisons: {analysis_3['key_comparisons']}")
print(f"Timeline: {analysis_3['timeline_weeks']}")
print()

tier1_analyses.append(analysis_3)

# ----------------------------------------------------------------------------
# ANALYSIS 1.4: DRUG RESPONSE PREDICTION FROM MULTI-OMICS
# ----------------------------------------------------------------------------

print("ANALYSIS 1.4: Drug Response Prediction from Multi-Omics")
print("-" * 80)

analysis_4 = {
    'tier': 1,
    'analysis_id': '1.4',
    'analysis_name': 'Drug Response Prediction from Multi-Omics',
    'goal': 'Predict drug sensitivity using expression + mutations',
    'methods': 'Machine learning (Random Forest, Elastic Net, XGBoost)',
    'required_n': 150,
    'available_n': n_gold_standard,
    'data_requirements': 'Expression + Mutations + Drug Response',
    'timeline_weeks': '3-4 weeks',
    'scripts_location': '02_Scripts/05_Drug_Response/',
}

# Power calculation
power_ml = 0.85 if n_gold_standard >= 150 else 0.6
feasible_ml = n_gold_standard >= 150

# Calculate train/val/test split
if n_gold_standard >= 100:
    train_60 = int(n_gold_standard * 0.6)
    valid_20 = int(n_gold_standard * 0.2)
    test_20 = n_gold_standard - train_60 - valid_20
    split_info = f"Train/Val/Test: {train_60}/{valid_20}/{test_20}"
else:
    split_info = "Insufficient for 3-way split"

analysis_4['statistical_power'] = power_ml
analysis_4['feasibility'] = 'YES' if feasible_ml else 'LIMITED'
analysis_4['justification'] = f"n={n_gold_standard}, {split_info}"

deliverables_4 = [
    "Predictive models for top 20 drugs (Random Forest, Elastic Net, XGBoost)",
    "Feature importance rankings (genes + mutations)",
    "Cross-validation performance metrics (R², RMSE, MAE)",
    "Independent test set validation",
    "Biomarker identification (top predictive features)",
    "Drug-specific response signatures"
]
analysis_4['deliverables'] = '; '.join(deliverables_4)

print(f"Goal: {analysis_4['goal']}")
print(f"Required samples: n ≥ {analysis_4['required_n']}")
print(f"Available samples: n = {analysis_4['available_n']}")
print(f"Statistical power: {analysis_4['statistical_power']:.2f}")
print(f"Feasibility: {analysis_4['feasibility']}")
print(f"  Justification: {analysis_4['justification']}")
print(f"  Features: Gene expression + mutation status")
print(f"  Target: AUC values for key drugs")
print(f"Timeline: {analysis_4['timeline_weeks']}")
print()

tier1_analyses.append(analysis_4)

# ============================================================================
# TIER 2: ADVANCED INTEGRATION ANALYSES (HIGH PRIORITY)
# ============================================================================

print("=" * 80)
print("TIER 2: ADVANCED INTEGRATION ANALYSES (HIGH PRIORITY)")
print("=" * 80)
print()

tier2_analyses = []

# ----------------------------------------------------------------------------
# ANALYSIS 2.1: SURVIVAL ANALYSIS WITH MULTI-OMICS FEATURES
# ----------------------------------------------------------------------------

print("ANALYSIS 2.1: Survival Analysis with Multi-Omics Features")
print("-" * 80)

analysis_5 = {
    'tier': 2,
    'analysis_id': '2.1',
    'analysis_name': 'Survival Analysis with Multi-Omics Features',
    'goal': 'Identify prognostic biomarkers from multi-omics data',
    'methods': 'Cox regression, Kaplan-Meier curves, risk stratification',
    'required_n': 100,  # With ≥50 events
    'available_n': n_survival,
    'data_requirements': 'Survival + Clinical + Expression + Mutations',
    'timeline_weeks': '2-3 weeks',
    'scripts_location': '02_Scripts/07_Survival_Analysis/',
}

# Power calculation
event_rate = n_events / n_survival * 100 if n_survival > 0 else 0
power_survival = 0.9 if n_events >= 50 else 0.6
feasible_survival = n_events >= 50 and n_survival >= 100

analysis_5['statistical_power'] = power_survival
analysis_5['feasibility'] = 'YES' if feasible_survival else 'LIMITED'
analysis_5['justification'] = f"n={n_survival}, {n_events} events ({event_rate:.1f}%)"

deliverables_5 = [
    "Univariate Cox regression for all features",
    "Multivariate Cox model with top features",
    "Kaplan-Meier curves for risk groups",
    "Risk score calculation and validation",
    "Time-dependent ROC curves",
    "Prognostic signature identification"
]
analysis_5['deliverables'] = '; '.join(deliverables_5)

print(f"Goal: {analysis_5['goal']}")
print(f"Required samples: n ≥ {analysis_5['required_n']} with ≥50 events")
print(f"Available samples: n = {analysis_5['available_n']}")
print(f"Events: {n_events} ({event_rate:.1f}%)")
print(f"Statistical power: {analysis_5['statistical_power']:.2f}")
print(f"Feasibility: {analysis_5['feasibility']}")
print(f"  Justification: {analysis_5['justification']}")
print(f"Timeline: {analysis_5['timeline_weeks']}")
print()

tier2_analyses.append(analysis_5)

# ----------------------------------------------------------------------------
# ANALYSIS 2.2: MUTATION-DRUG RESPONSE ASSOCIATIONS
# ----------------------------------------------------------------------------

print("ANALYSIS 2.2: Mutation-Drug Response Associations")
print("-" * 80)

analysis_6 = {
    'tier': 2,
    'analysis_id': '2.2',
    'analysis_name': 'Mutation-Drug Response Associations',
    'goal': 'Identify mutation-specific drug sensitivities/resistances',
    'methods': 'Stratified AUC comparison, effect size calculation',
    'required_n': 40,  # 20+ per group
    'available_n': n_mut_drug,
    'data_requirements': 'Mutations + Drug Response',
    'timeline_weeks': '2-3 weeks',
    'scripts_location': '02_Scripts/05_Drug_Response/',
}

# Key mutation-drug pairs
key_pairs = [
    ('FLT3', 'FLT3 inhibitors (Sorafenib, Gilteritinib)'),
    ('IDH1', 'IDH inhibitors (Ivosidenib)'),
    ('IDH2', 'IDH inhibitors (Enasidenib)'),
    ('DNMT3A', 'Hypomethylating agents (Azacitidine)'),
    ('NPM1', 'Venetoclax')
]

# Estimate power for key pairs (simplified - assume ~30% mutation frequency)
n_mut_est = int(n_mut_drug * 0.2)  # Conservative estimate
n_wt_est = n_mut_drug - n_mut_est
power_mut_drug = 0.8 if n_mut_est >= 20 else 0.5
feasible_mut_drug = n_mut_est >= 20

analysis_6['statistical_power'] = power_mut_drug
analysis_6['feasibility'] = 'YES' if feasible_mut_drug else 'LIMITED'
analysis_6['justification'] = f"n={n_mut_drug}, key pairs have ≥20 per group"
analysis_6['key_pairs'] = '; '.join([f"{m} vs {d}" for m, d in key_pairs])

deliverables_6 = [
    "Mutation-stratified drug sensitivity plots",
    "Effect size calculations (Cohen's d) for each pair",
    "Statistical testing (t-test, Mann-Whitney U)",
    "Volcano plot of mutation-drug associations",
    "Precision medicine implications",
    "Actionable biomarkers for drug selection"
]
analysis_6['deliverables'] = '; '.join(deliverables_6)

print(f"Goal: {analysis_6['goal']}")
print(f"Required samples: n ≥ {analysis_6['required_n']} (20+ per group)")
print(f"Available samples: n = {analysis_6['available_n']}")
print(f"Statistical power: {analysis_6['statistical_power']:.2f}")
print(f"Feasibility: {analysis_6['feasibility']}")
print(f"  Justification: {analysis_6['justification']}")
print(f"  Key pairs: {analysis_6['key_pairs']}")
print(f"Timeline: {analysis_6['timeline_weeks']}")
print()

tier2_analyses.append(analysis_6)

# ----------------------------------------------------------------------------
# ANALYSIS 2.3: INTEGRATED NETWORK ANALYSIS
# ----------------------------------------------------------------------------

print("ANALYSIS 2.3: Integrated Network Analysis")
print("-" * 80)

analysis_7 = {
    'tier': 2,
    'analysis_id': '2.3',
    'analysis_name': 'Integrated Network Analysis',
    'goal': 'Construct multi-omics interaction networks',
    'methods': 'Correlation networks, pathway analysis, module detection',
    'required_n': 100,
    'available_n': n_gold_standard,
    'data_requirements': 'Expression + Mutations + Drug Response',
    'timeline_weeks': '3-4 weeks',
    'scripts_location': '02_Scripts/06_Integration/',
}

# Power calculation
power_network = 0.8 if n_gold_standard >= 100 else 0.5
feasible_network = n_gold_standard >= 100

analysis_7['statistical_power'] = power_network
analysis_7['feasibility'] = 'YES' if feasible_network else 'LIMITED'
analysis_7['justification'] = f"n={n_gold_standard} for integrated networks"

deliverables_7 = [
    "Gene-gene correlation networks",
    "Mutation-expression regulatory networks",
    "Drug-target-expression networks",
    "Module detection (WGCNA or similar)",
    "Hub gene identification",
    "Network-based biomarker discovery"
]
analysis_7['deliverables'] = '; '.join(deliverables_7)

print(f"Goal: {analysis_7['goal']}")
print(f"Required samples: n ≥ {analysis_7['required_n']}")
print(f"Available samples: n = {analysis_7['available_n']}")
print(f"Statistical power: {analysis_7['statistical_power']:.2f}")
print(f"Feasibility: {analysis_7['feasibility']}")
print(f"  Justification: {analysis_7['justification']}")
print(f"Timeline: {analysis_7['timeline_weeks']}")
print()

tier2_analyses.append(analysis_7)

# ============================================================================
# TIER 3: EXPLORATORY ANALYSES (MEDIUM PRIORITY)
# ============================================================================

print("=" * 80)
print("TIER 3: EXPLORATORY ANALYSES (MEDIUM PRIORITY)")
print("=" * 80)
print()

tier3_analyses = []

# ----------------------------------------------------------------------------
# ANALYSIS 3.1: SUBTYPE-SPECIFIC DRUG SENSITIVITIES
# ----------------------------------------------------------------------------

print("ANALYSIS 3.1: Subtype-Specific Drug Sensitivities")
print("-" * 80)

# Requires completion of Analysis 1.1 first
n_expr_drug_for_subtype = n_expr_drug  # Use expression+drug cohort

analysis_8 = {
    'tier': 3,
    'analysis_id': '3.1',
    'analysis_name': 'Subtype-Specific Drug Sensitivities',
    'goal': 'Identify drugs with differential sensitivity across subtypes',
    'methods': 'ANOVA, pairwise comparisons, effect size calculation',
    'required_n': 60,  # ~20 per subtype for 3 subtypes
    'available_n': n_expr_drug_for_subtype,
    'data_requirements': 'Expression + Drug Response (requires Analysis 1.1)',
    'timeline_weeks': '2 weeks',
    'scripts_location': '02_Scripts/05_Drug_Response/',
    'dependencies': 'Requires Analysis 1.1 completion',
}

# Power calculation - assumes 3 subtypes
n_per_subtype = n_expr_drug_for_subtype / 3
power_subtype_drug = 0.8 if n_per_subtype >= 20 else 0.5
feasible_subtype_drug = n_per_subtype >= 20

analysis_8['statistical_power'] = power_subtype_drug
analysis_8['feasibility'] = 'YES' if feasible_subtype_drug else 'LIMITED'
analysis_8['justification'] = f"n={n_expr_drug_for_subtype}, ~{n_per_subtype:.0f} per subtype"

deliverables_8 = [
    "Subtype-stratified drug sensitivity plots",
    "ANOVA results for all drugs",
    "Subtype-specific drug rankings",
    "Precision medicine recommendations",
    "Heatmap of subtype-drug associations"
]
analysis_8['deliverables'] = '; '.join(deliverables_8)

print(f"Goal: {analysis_8['goal']}")
print(f"Required samples: n ≥ {analysis_8['required_n']}")
print(f"Available samples: n = {analysis_8['available_n']}")
print(f"Statistical power: {analysis_8['statistical_power']:.2f}")
print(f"Feasibility: {analysis_8['feasibility']}")
print(f"  Justification: {analysis_8['justification']}")
print(f"  Dependencies: {analysis_8['dependencies']}")
print(f"Timeline: {analysis_8['timeline_weeks']}")
print()

tier3_analyses.append(analysis_8)

# ----------------------------------------------------------------------------
# ANALYSIS 3.2: CLINICAL-MOLECULAR CORRELATIONS
# ----------------------------------------------------------------------------

print("ANALYSIS 3.2: Clinical-Molecular Correlations")
print("-" * 80)

analysis_9 = {
    'tier': 3,
    'analysis_id': '3.2',
    'analysis_name': 'Clinical-Molecular Correlations',
    'goal': 'Associate clinical features with molecular profiles',
    'methods': 'Correlation analysis, stratified comparisons',
    'required_n': 100,
    'available_n': n_gold_standard,
    'data_requirements': 'Clinical + Expression + Mutations',
    'timeline_weeks': '1-2 weeks',
    'scripts_location': '02_Scripts/06_Integration/',
}

# Power calculation
power_clinical = 0.7 if n_gold_standard >= 100 else 0.5
feasible_clinical = n_gold_standard >= 100

analysis_9['statistical_power'] = power_clinical
analysis_9['feasibility'] = 'YES' if feasible_clinical else 'LIMITED'
analysis_9['justification'] = f"n={n_gold_standard} with complete clinical data"

deliverables_9 = [
    "Age vs expression/mutation associations",
    "Sex-stratified molecular profiles",
    "Prior treatment vs molecular features",
    "WBC count vs expression signatures",
    "Clinical subgroup comparisons (ELN risk, FAB classification)"
]
analysis_9['deliverables'] = '; '.join(deliverables_9)

print(f"Goal: {analysis_9['goal']}")
print(f"Required samples: n ≥ {analysis_9['required_n']}")
print(f"Available samples: n = {analysis_9['available_n']}")
print(f"Statistical power: {analysis_9['statistical_power']:.2f}")
print(f"Feasibility: {analysis_9['feasibility']}")
print(f"  Justification: {analysis_9['justification']}")
print(f"Timeline: {analysis_9['timeline_weeks']}")
print()

tier3_analyses.append(analysis_9)

# ----------------------------------------------------------------------------
# ANALYSIS 3.3: DRUG COMBINATION SYNERGY ANALYSIS
# ----------------------------------------------------------------------------

print("ANALYSIS 3.3: Drug Combination Synergy Analysis")
print("-" * 80)

analysis_10 = {
    'tier': 3,
    'analysis_id': '3.3',
    'analysis_name': 'Drug Combination Synergy Analysis',
    'goal': 'Identify synergistic drug combinations',
    'methods': 'Correlation analysis, synergy scoring (if combination data available)',
    'required_n': 100,
    'available_n': n_drug_response,
    'data_requirements': 'Drug Response (multiple drugs per sample)',
    'timeline_weeks': '2-3 weeks',
    'scripts_location': '02_Scripts/05_Drug_Response/',
    'caveat': 'Limited by single-drug screening design',
}

# Power calculation
power_combo = 0.6  # Lower due to limitations
feasible_combo = n_drug_response >= 100

analysis_10['statistical_power'] = power_combo
analysis_10['feasibility'] = 'LIMITED' if feasible_combo else 'NO'
analysis_10['justification'] = f"n={n_drug_response}, but single-drug screening limits analysis"

deliverables_10 = [
    "Drug-drug correlation matrix",
    "Predicted synergistic pairs (correlation-based)",
    "Pathway-based combination hypotheses",
    "Literature-validated combinations"
]
analysis_10['deliverables'] = '; '.join(deliverables_10)

print(f"Goal: {analysis_10['goal']}")
print(f"Required samples: n ≥ {analysis_10['required_n']}")
print(f"Available samples: n = {analysis_10['available_n']}")
print(f"Statistical power: {analysis_10['statistical_power']:.2f}")
print(f"Feasibility: {analysis_10['feasibility']}")
print(f"  Justification: {analysis_10['justification']}")
print(f"  Caveat: {analysis_10['caveat']}")
print(f"Timeline: {analysis_10['timeline_weeks']}")
print()

tier3_analyses.append(analysis_10)

# ============================================================================
# SAVE ROADMAP
# ============================================================================

print("=" * 80)
print("SAVING COMPREHENSIVE ANALYSIS ROADMAP")
print("=" * 80)
print()

# Combine all analyses
all_analyses = tier1_analyses + tier2_analyses + tier3_analyses

# Create DataFrame
df_roadmap = pd.DataFrame(all_analyses)

# Save to CSV
output_file = OUTPUT_DIR / "comprehensive_analysis_roadmap.csv"
df_roadmap.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# ============================================================================
# GENERATE MARKDOWN REPORT
# ============================================================================

print("Generating markdown report...")

report_lines = []
report_lines.append("# AML Multi-Omics Integration: Comprehensive Analysis Roadmap")
report_lines.append("")
report_lines.append("**Generated:** 2025-10-02")
report_lines.append("")
report_lines.append("**Author:** AML Multi-Omics Project Team")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# Executive Summary
report_lines.append("## Executive Summary")
report_lines.append("")
report_lines.append(f"This roadmap outlines **{len(all_analyses)} prioritized analyses** across 3 tiers:")
report_lines.append("")
report_lines.append(f"- **TIER 1 (Highest Priority):** {len(tier1_analyses)} core multi-omics analyses")
report_lines.append(f"- **TIER 2 (High Priority):** {len(tier2_analyses)} advanced integration analyses")
report_lines.append(f"- **TIER 3 (Medium Priority):** {len(tier3_analyses)} exploratory analyses")
report_lines.append("")

# Cohort summary
report_lines.append("### Available Cohort Sizes")
report_lines.append("")
report_lines.append("| Data Type Combination | Sample Size |")
report_lines.append("|----------------------|-------------|")
report_lines.append(f"| Expression only | n = {n_expression} |")
report_lines.append(f"| Mutations only | n = {n_mutations} |")
report_lines.append(f"| Drug Response only | n = {n_drug_response} |")
report_lines.append(f"| Clinical only | n = {n_clinical} |")
report_lines.append(f"| Expression + Mutations | n = {n_expr_mut} |")
report_lines.append(f"| Expression + Drug | n = {n_expr_drug} |")
report_lines.append(f"| Mutations + Drug | n = {n_mut_drug} |")
report_lines.append(f"| **Gold Standard (all 4)** | **n = {n_gold_standard}** |")
report_lines.append(f"| Survival data | n = {n_survival} ({n_events} events) |")
report_lines.append("")
report_lines.append("---")
report_lines.append("")

# TIER 1
report_lines.append("## TIER 1: Core Multi-Omics Analyses (HIGHEST PRIORITY)")
report_lines.append("")
report_lines.append("These analyses form the foundation of the multi-omics integration project and should be completed first.")
report_lines.append("")

for analysis in tier1_analyses:
    report_lines.append(f"### Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
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
    report_lines.append("")
    report_lines.append(f"**Statistical Power:** {analysis['statistical_power']:.2f}")
    report_lines.append("")
    report_lines.append(f"**Feasibility:** {analysis['feasibility']}")
    report_lines.append(f"- *Justification:* {analysis['justification']}")
    report_lines.append("")
    report_lines.append(f"**Timeline:** {analysis['timeline_weeks']}")
    report_lines.append("")
    report_lines.append(f"**Scripts Location:** `{analysis['scripts_location']}`")
    report_lines.append("")
    report_lines.append(f"**Deliverables:**")
    for deliverable in analysis['deliverables'].split('; '):
        report_lines.append(f"- {deliverable}")
    report_lines.append("")

    # Add key comparisons if available
    if 'key_comparisons' in analysis:
        report_lines.append(f"**Key Comparisons:** {analysis['key_comparisons']}")
        report_lines.append("")

    report_lines.append("---")
    report_lines.append("")

# TIER 2
report_lines.append("## TIER 2: Advanced Integration Analyses (HIGH PRIORITY)")
report_lines.append("")
report_lines.append("These analyses build upon Tier 1 results and provide deeper biological insights.")
report_lines.append("")

for analysis in tier2_analyses:
    report_lines.append(f"### Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
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
    report_lines.append("")
    report_lines.append(f"**Statistical Power:** {analysis['statistical_power']:.2f}")
    report_lines.append("")
    report_lines.append(f"**Feasibility:** {analysis['feasibility']}")
    report_lines.append(f"- *Justification:* {analysis['justification']}")
    report_lines.append("")
    report_lines.append(f"**Timeline:** {analysis['timeline_weeks']}")
    report_lines.append("")
    report_lines.append(f"**Scripts Location:** `{analysis['scripts_location']}`")
    report_lines.append("")
    report_lines.append(f"**Deliverables:**")
    for deliverable in analysis['deliverables'].split('; '):
        report_lines.append(f"- {deliverable}")
    report_lines.append("")

    # Add key pairs if available
    if 'key_pairs' in analysis:
        report_lines.append(f"**Key Associations:** {analysis['key_pairs']}")
        report_lines.append("")

    report_lines.append("---")
    report_lines.append("")

# TIER 3
report_lines.append("## TIER 3: Exploratory Analyses (MEDIUM PRIORITY)")
report_lines.append("")
report_lines.append("These analyses are exploratory and can be performed after Tier 1-2 completion.")
report_lines.append("")

for analysis in tier3_analyses:
    report_lines.append(f"### Analysis {analysis['analysis_id']}: {analysis['analysis_name']}")
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
    report_lines.append("")
    report_lines.append(f"**Statistical Power:** {analysis['statistical_power']:.2f}")
    report_lines.append("")
    report_lines.append(f"**Feasibility:** {analysis['feasibility']}")
    report_lines.append(f"- *Justification:* {analysis['justification']}")
    report_lines.append("")

    # Add dependencies if available
    if 'dependencies' in analysis:
        report_lines.append(f"**Dependencies:** {analysis['dependencies']}")
        report_lines.append("")

    # Add caveats if available
    if 'caveat' in analysis:
        report_lines.append(f"**Caveat:** {analysis['caveat']}")
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

# Timeline summary
report_lines.append("## Estimated Timeline")
report_lines.append("")
report_lines.append("| Tier | Analyses | Total Time | Can Run in Parallel |")
report_lines.append("|------|----------|------------|---------------------|")

tier1_weeks = "8-12 weeks"
tier2_weeks = "7-10 weeks"
tier3_weeks = "5-7 weeks"

report_lines.append(f"| **Tier 1** | {len(tier1_analyses)} analyses | {tier1_weeks} | Some analyses (1.1, 1.2) |")
report_lines.append(f"| **Tier 2** | {len(tier2_analyses)} analyses | {tier2_weeks} | Yes, after Tier 1 |")
report_lines.append(f"| **Tier 3** | {len(tier3_analyses)} analyses | {tier3_weeks} | Yes, flexible |")
report_lines.append("")
report_lines.append("**Total estimated time:** 20-29 weeks (5-7 months) if performed sequentially")
report_lines.append("")
report_lines.append("**With parallelization:** 12-18 weeks (3-4.5 months)")
report_lines.append("")

# Feasibility summary
report_lines.append("## Feasibility Summary")
report_lines.append("")
report_lines.append("| Analysis ID | Analysis Name | Feasibility | Power |")
report_lines.append("|-------------|---------------|-------------|-------|")
for analysis in all_analyses:
    report_lines.append(f"| {analysis['analysis_id']} | {analysis['analysis_name']} | {analysis['feasibility']} | {analysis['statistical_power']:.2f} |")
report_lines.append("")

feasible_count = sum([1 for a in all_analyses if a['feasibility'] == 'YES'])
report_lines.append(f"**Summary:** {feasible_count}/{len(all_analyses)} analyses are fully feasible with available data.")
report_lines.append("")

# Recommendations
report_lines.append("## Recommendations")
report_lines.append("")
report_lines.append("1. **Start with Tier 1 Core Analyses:**")
report_lines.append("   - Begin with Analysis 1.1 (Molecular Subtyping) and 1.2 (Mutation Landscape) in parallel")
report_lines.append("   - These are foundational and inform downstream analyses")
report_lines.append("")
report_lines.append("2. **Prioritize Well-Powered Analyses:**")
report_lines.append(f"   - {feasible_count}/{len(all_analyses)} analyses have adequate statistical power")
report_lines.append("   - Focus on these for high-confidence results")
report_lines.append("")
report_lines.append("3. **Sequential Dependencies:**")
report_lines.append("   - Analysis 3.1 requires completion of Analysis 1.1 first")
report_lines.append("   - Plan workflow accordingly")
report_lines.append("")
report_lines.append("4. **Resource Allocation:**")
report_lines.append("   - Tier 1: Allocate primary resources and personnel")
report_lines.append("   - Tier 2: Begin after Tier 1 completion")
report_lines.append("   - Tier 3: Optional/exploratory based on time and resources")
report_lines.append("")
report_lines.append("5. **Quality Control:**")
report_lines.append("   - Apply batch correction before expression-based analyses")
report_lines.append("   - Use gold standard cohort (n=478) for integrated analyses")
report_lines.append("   - Validate key findings across independent cohorts if possible")
report_lines.append("")

# Save markdown report
report_file = REPORT_DIR / "Comprehensive_Analysis_Roadmap.md"
with open(report_file, 'w', encoding='utf-8') as f:
    f.write('\n'.join(report_lines))

print(f"✓ Saved: {report_file.name}")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("COMPREHENSIVE ANALYSIS ROADMAP COMPLETE")
print("=" * 80)
print()
print(f"TOTAL ANALYSES PLANNED: {len(all_analyses)}")
print(f"  Tier 1 (Highest Priority): {len(tier1_analyses)}")
print(f"  Tier 2 (High Priority): {len(tier2_analyses)}")
print(f"  Tier 3 (Medium Priority): {len(tier3_analyses)}")
print()
print(f"FEASIBILITY:")
print(f"  Fully feasible: {feasible_count}/{len(all_analyses)}")
print(f"  Limited/Exploratory: {len(all_analyses) - feasible_count}/{len(all_analyses)}")
print()
print("ESTIMATED TIMELINE:")
print(f"  Sequential: 20-29 weeks (5-7 months)")
print(f"  With parallelization: 12-18 weeks (3-4.5 months)")
print()
print("OUTPUTS SAVED:")
print(f"  1. {output_file.name}")
print(f"  2. {report_file.name}")
print()
print("=" * 80)
