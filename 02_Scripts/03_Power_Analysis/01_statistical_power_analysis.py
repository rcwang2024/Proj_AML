"""
Phase 5: Statistical Power and Feasibility Assessment
Task 5.1: Statistical Power Analysis

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
import numpy as np
from scipy import stats
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
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/02_QC_Reports")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TASK 5.1: STATISTICAL POWER ANALYSIS")
print("=" * 80)
print()

# ============================================================================
# LOAD DATA
# ============================================================================

print("LOADING DATA FOR POWER ANALYSIS")
print("-" * 80)

# Master mapping
master_map = pd.read_csv(PROCESSED_DIR / "master_sample_id_mapping.csv")

# Clinical data
df_clinical = pd.read_excel(DATA_DIR / "beataml_clinical.xlsx")

# Mutation data
df_mut = pd.read_csv(DATA_DIR / "beataml_mutations.txt", sep='\t')

# Drug response
df_drug = pd.read_csv(DATA_DIR / "beataml_drug_auc.txt", sep='\t')

# Driver mutation frequencies
driver_freq = pd.read_csv(PROCESSED_DIR / "driver_mutation_frequencies.csv")

print(f"✓ Data loaded")
print()

# ============================================================================
# ANALYSIS 1: MULTI-OMICS INTEGRATION
# ============================================================================

print("1. MULTI-OMICS INTEGRATION ANALYSIS")
print("-" * 80)

# Gold standard cohort (all 4 data types)
gold_standard = master_map[
    (master_map['has_expression']) &
    (master_map['has_drug_response']) &
    (master_map['has_clinical']) &
    (master_map['has_mutations'])
]

n_multi_omics = len(gold_standard)

print(f"Cohort: Samples with ALL FOUR data types")
print(f"Sample size: n = {n_multi_omics}")
print()

# Power assessment
min_exploratory = n_multi_omics >= 50
good_robust = n_multi_omics >= 100
excellent_pub = n_multi_omics >= 200

print(f"Assessment:")
print(f"  Is n ≥ 50 (minimum for exploratory)? {'✓ YES' if min_exploratory else '✗ NO'}")
print(f"  Is n ≥ 100 (good for robust analysis)? {'✓ YES' if good_robust else '✗ NO'}")
print(f"  Is n ≥ 200 (excellent for publication)? {'✓ YES' if excellent_pub else '✗ NO'}")
print()

# Calculate power for moderate effect (d=0.5, two-sample t-test)
from statsmodels.stats.power import TTestIndPower
power_analysis = TTestIndPower()
# Assume equal groups
n_per_group = n_multi_omics / 2
power_moderate = power_analysis.solve_power(effect_size=0.5, nobs1=n_per_group,
                                            ratio=1.0, alpha=0.05, alternative='two-sided')

print(f"Power for detecting moderate effects (d=0.5):")
print(f"  Power (two-sample t-test): {power_moderate:.3f} ({power_moderate*100:.1f}%)")
print(f"  Interpretation: {'Adequate (≥0.8)' if power_moderate >= 0.8 else 'Limited (<0.8)'}")
print()

feasible_multi = n_multi_omics >= 100 and power_moderate >= 0.8
print(f"Recommendation: {'✓ FEASIBLE' if feasible_multi else '⚠ LIMITED POWER'}")
if feasible_multi:
    print(f"  Justification: n={n_multi_omics} provides adequate power for multi-omics integration")
else:
    print(f"  Justification: n={n_multi_omics} may limit detection of subtle effects")
print()

# ============================================================================
# ANALYSIS 2: MOLECULAR SUBTYPING
# ============================================================================

print("2. MOLECULAR SUBTYPING (EXPRESSION-BASED)")
print("-" * 80)

# Samples with expression
n_expression = (master_map['has_expression']).sum()

print(f"Method: Unsupervised clustering (hierarchical, k-means, consensus)")
print(f"Required: n ≥ 100 samples with expression")
print(f"Available: n = {n_expression}")
print()

sufficient_clustering = n_expression >= 100

print(f"Assessment:")
print(f"  Can identify 3-5 stable clusters? {'✓ YES' if sufficient_clustering else '⚠ BORDERLINE'}")

# Estimate cluster sizes for k=3,4,5
for k in [3, 4, 5]:
    avg_cluster_size = n_expression / k
    print(f"    For {k} clusters: ~{avg_cluster_size:.0f} samples/cluster")

print()

# Cluster-specific signatures
min_cluster_size = n_expression / 5  # Worst case with 5 clusters
signature_power = min_cluster_size >= 50

print(f"  Power to detect cluster-specific signatures?")
print(f"    Minimum cluster size (k=5): ~{min_cluster_size:.0f}")
print(f"    {'✓ Adequate (≥50)' if signature_power else '⚠ Limited (<50)'}")
print()

validation_ok = n_expression >= 200
print(f"  Sufficient for validation? {'✓ YES (n≥200)' if validation_ok else '⚠ BORDERLINE'}")
print()

feasible_subtype = sufficient_clustering and signature_power
print(f"Recommendation: {'✓ FEASIBLE' if feasible_subtype else '⚠ LIMITED'}")
print()

# ============================================================================
# ANALYSIS 3: MUTATION-EXPRESSION CORRELATION
# ============================================================================

print("3. MUTATION-EXPRESSION CORRELATION")
print("-" * 80)

print(f"Method: Differential expression by mutation status")
print()

# Get samples with both expression and mutations
expr_mut_samples = master_map[
    (master_map['has_expression']) &
    (master_map['has_mutations'])
]['unified_sample_id']

n_expr_mut = len(expr_mut_samples)
print(f"Samples with expression + mutations: {n_expr_mut}")
print()

# Analyze key mutations
key_mutations = driver_freq.sort_values('frequency_pct', ascending=False).head(10)

print(f"Power assessment for key mutations:")
print(f"{'Gene':<12} {'Mutated':>10} {'Wild-type':>10} {'Power':>10} {'Feasible':>10}")
print("-" * 56)

mut_power_results = []

for _, row in key_mutations.iterrows():
    gene = row['gene']
    n_mutated = int(row['n_samples_mutated'])
    n_wt = 871 - n_mutated  # Total mutation samples - mutated

    # Adjust for expression overlap (assume ~70% have expression)
    n_mut_expr = int(n_mutated * 0.7)
    n_wt_expr = int(n_wt * 0.7)

    # Power for 2-fold change (effect size ~1.0 for log2FC=1)
    min_n = min(n_mut_expr, n_wt_expr)

    # Simple power estimation
    if min_n >= 30:
        power_est = 0.9
    elif min_n >= 20:
        power_est = 0.8
    elif min_n >= 10:
        power_est = 0.6
    else:
        power_est = 0.3

    feasible = min_n >= 10

    print(f"{gene:<12} {n_mut_expr:>10} {n_wt_expr:>10} {power_est:>10.2f} {'✓ YES' if feasible else '✗ NO':>10}")

    mut_power_results.append({
        'gene': gene,
        'n_mutated_expr': n_mut_expr,
        'n_wt_expr': n_wt_expr,
        'power_estimate': power_est,
        'feasible': feasible
    })

print()
feasible_genes = [r['gene'] for r in mut_power_results if r['feasible']]
print(f"Recommendation: {len(feasible_genes)}/{len(key_mutations)} mutations have sufficient power")
print(f"  Feasible: {', '.join(feasible_genes[:5])}")
if len(feasible_genes) > 5:
    print(f"    ... and {len(feasible_genes)-5} more")
print()

# ============================================================================
# ANALYSIS 4: MUTATION-DRUG RESPONSE ASSOCIATION
# ============================================================================

print("4. MUTATION-DRUG RESPONSE ASSOCIATION")
print("-" * 80)

# Samples with both mutations and drug data
mut_drug_samples = master_map[
    (master_map['has_mutations']) &
    (master_map['has_drug_response'])
]

n_mut_drug = len(mut_drug_samples)

print(f"Method: Compare drug sensitivity by mutation status")
print(f"Required: ≥20 samples per group")
print(f"Available samples with mutation + drug data: n = {n_mut_drug}")
print()

# Key mutation-drug pairs
key_pairs = [
    ('FLT3', 'FLT3 inhibitors', ['Sorafenib', 'Sunitinib', 'Gilteritinib']),
    ('IDH1', 'IDH inhibitors', ['Ivosidenib', 'AG-221']),
    ('IDH2', 'IDH inhibitors', ['Enasidenib', 'AG-221']),
    ('DNMT3A', 'Hypomethylating agents', ['Azacitidine', 'Decitabine']),
    ('NPM1', 'Targeted agents', ['Venetoclax'])
]

print(f"Key mutation-drug associations:")
print(f"{'Mutation':<12} {'Drug Class':<25} {'Est. N Mut':>12} {'Est. N WT':>10} {'Feasible':>10}")
print("-" * 75)

drug_assoc_results = []

for gene, drug_class, _ in key_pairs:
    # Get mutation frequency
    gene_freq = driver_freq[driver_freq['gene'] == gene]
    if len(gene_freq) > 0:
        n_mut = int(gene_freq.iloc[0]['n_samples_mutated'])

        # Estimate overlap with drug data (~80%)
        n_mut_drug_est = int(n_mut * 0.8)
        n_wt_drug_est = n_mut_drug - n_mut_drug_est

        feasible = n_mut_drug_est >= 20 and n_wt_drug_est >= 20

        print(f"{gene:<12} {drug_class:<25} {n_mut_drug_est:>12} {n_wt_drug_est:>10} {'✓ YES' if feasible else '⚠ LIMITED':>10}")

        drug_assoc_results.append({
            'mutation': gene,
            'drug_class': drug_class,
            'n_mut_drug': n_mut_drug_est,
            'n_wt_drug': n_wt_drug_est,
            'feasible': feasible
        })
    else:
        print(f"{gene:<12} {drug_class:<25} {'N/A':>12} {'N/A':>10} {'✗ NO':>10}")

print()
feasible_pairs = [r for r in drug_assoc_results if r['feasible']]
print(f"Recommendation: {len(feasible_pairs)}/{len(key_pairs)} associations have sufficient power")
if feasible_pairs:
    print(f"  Feasible pairs:")
    for pair in feasible_pairs:
        print(f"    • {pair['mutation']} vs {pair['drug_class']}")
print()

# ============================================================================
# ANALYSIS 5: SURVIVAL ANALYSIS
# ============================================================================

print("5. SURVIVAL ANALYSIS")
print("-" * 80)

# Get survival data
survival_data = df_clinical[['overallSurvival', 'vitalStatus']].dropna()
n_survival = len(survival_data)

# Count events (deaths)
# vitalStatus: dead=1, alive=0 (usually)
events = (survival_data['vitalStatus'] == 'Dead').sum()
event_rate = events / n_survival * 100 if n_survival > 0 else 0

print(f"Samples with survival data: n = {n_survival}")
print(f"Events (deaths): n = {events}")
print(f"Event rate: {event_rate:.1f}%")
print()

print(f"Assessment:")

# Cox regression with 5-10 variables
min_events_cox = 50
cox_5var = events >= min_events_cox
cox_10var = events >= min_events_cox * 2  # More conservative for 10 vars

print(f"  Cox regression with 5-10 variables:")
print(f"    Need n ≥ 50 events: {'✓ YES' if cox_5var else '✗ NO'} (have {events})")
print(f"    5 variables: {'✓ Feasible' if cox_5var else '✗ Underpowered'}")
print(f"    10 variables: {'✓ Feasible' if cox_10var else '⚠ Limited'}")
print()

# Kaplan-Meier by groups
n_per_km_group = n_survival / 3  # Assume 3 risk groups
km_feasible = n_per_km_group >= 20

print(f"  Kaplan-Meier by groups:")
print(f"    Need n ≥ 20 per group: {'✓ YES' if km_feasible else '✗ NO'}")
print(f"    For 3 groups: ~{n_per_km_group:.0f} samples/group")
print(f"    {'✓ Feasible' if km_feasible else '✗ Underpowered'}")
print()

survival_feasible = cox_5var and km_feasible
print(f"Recommendation: {'✓ FEASIBLE' if survival_feasible else '⚠ LIMITED'}")
if survival_feasible:
    print(f"  • Cox regression with ≤5-7 variables")
    print(f"  • KM curves for 2-3 risk groups")
    print(f"  • Multivariate models with clinical + molecular features")
print()

# ============================================================================
# ANALYSIS 6: PREDICTIVE MODELING
# ============================================================================

print("6. PREDICTIVE MODELING (DRUG RESPONSE)")
print("-" * 80)

# Samples with expression + drug response
expr_drug = master_map[
    (master_map['has_expression']) &
    (master_map['has_drug_response'])
]

n_expr_drug = len(expr_drug)

print(f"Machine learning approach: Random Forest, Elastic Net")
print(f"Training samples needed: ≥100 (preferably ≥200)")
print(f"Available: n = {n_expr_drug} with expression + drug response")
print()

# Also count with mutations
expr_drug_mut = master_map[
    (master_map['has_expression']) &
    (master_map['has_drug_response']) &
    (master_map['has_mutations'])
]
n_full = len(expr_drug_mut)

print(f"With all features (expression + mutations): n = {n_full}")
print()

# Power assessment
min_ml = n_expr_drug >= 100
good_ml = n_expr_drug >= 200
cv_feasible = n_expr_drug >= 150  # For 5-fold CV

print(f"Assessment:")
print(f"  Minimum for ML (n≥100): {'✓ YES' if min_ml else '✗ NO'}")
print(f"  Good for ML (n≥200): {'✓ YES' if good_ml else '⚠ BORDERLINE'}")
print(f"  5-fold cross-validation: {'✓ Feasible' if cv_feasible else '⚠ Limited'}")
print()

# Training/validation/test split
if n_expr_drug >= 100:
    train_60 = int(n_expr_drug * 0.6)
    valid_20 = int(n_expr_drug * 0.2)
    test_20 = n_expr_drug - train_60 - valid_20

    print(f"  Recommended split (60/20/20):")
    print(f"    Training: {train_60}")
    print(f"    Validation: {valid_20}")
    print(f"    Test: {test_20}")
    print()

ml_feasible = n_expr_drug >= 150
print(f"Recommendation: {'✓ FEASIBLE' if ml_feasible else '⚠ LIMITED'}")
if ml_feasible:
    print(f"  • Cross-validated models feasible")
    print(f"  • Can use regularization (Elastic Net, LASSO)")
    print(f"  • Consider ensemble methods (Random Forest, XGBoost)")
else:
    print(f"  • Sample size limits model complexity")
    print(f"  • Use simple models with strong regularization")
    print(f"  • External validation recommended")
print()

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("SAVING POWER ANALYSIS RESULTS")
print("-" * 80)

power_results = []

# 1. Multi-omics integration
power_results.append({
    'analysis_type': 'Multi-omics Integration',
    'method': 'Integrated analysis (all 4 data types)',
    'sample_size': n_multi_omics,
    'min_required': 100,
    'power_estimate': power_moderate,
    'feasible': feasible_multi,
    'notes': f'Power={power_moderate:.2f} for d=0.5'
})

# 2. Molecular subtyping
power_results.append({
    'analysis_type': 'Molecular Subtyping',
    'method': 'Unsupervised clustering',
    'sample_size': n_expression,
    'min_required': 100,
    'power_estimate': 0.9 if feasible_subtype else 0.6,
    'feasible': feasible_subtype,
    'notes': f'Can identify 3-5 clusters'
})

# 3. Mutation-expression correlation
power_results.append({
    'analysis_type': 'Mutation-Expression Correlation',
    'method': 'Differential expression',
    'sample_size': n_expr_mut,
    'min_required': 20,  # 10 per group
    'power_estimate': np.mean([r['power_estimate'] for r in mut_power_results]),
    'feasible': len(feasible_genes) >= 5,
    'notes': f'{len(feasible_genes)} genes with power≥0.6'
})

# 4. Mutation-drug association
power_results.append({
    'analysis_type': 'Mutation-Drug Association',
    'method': 'Drug sensitivity by mutation',
    'sample_size': n_mut_drug,
    'min_required': 40,  # 20 per group
    'power_estimate': 0.8 if len(feasible_pairs) >= 3 else 0.5,
    'feasible': len(feasible_pairs) >= 3,
    'notes': f'{len(feasible_pairs)} associations feasible'
})

# 5. Survival analysis
power_results.append({
    'analysis_type': 'Survival Analysis',
    'method': 'Cox regression / Kaplan-Meier',
    'sample_size': n_survival,
    'min_required': 50,  # events
    'power_estimate': 0.9 if survival_feasible else 0.6,
    'feasible': survival_feasible,
    'notes': f'{events} events ({event_rate:.1f}%)'
})

# 6. Predictive modeling
power_results.append({
    'analysis_type': 'Predictive Modeling (Drug Response)',
    'method': 'Machine learning (RF, Elastic Net)',
    'sample_size': n_expr_drug,
    'min_required': 150,
    'power_estimate': 0.85 if ml_feasible else 0.6,
    'feasible': ml_feasible,
    'notes': f'Train/Val/Test: {train_60}/{valid_20}/{test_20}' if n_expr_drug >= 100 else 'Limited'
})

# Save to CSV
df_power = pd.DataFrame(power_results)
output_file = OUTPUT_DIR / "statistical_power_analysis.csv"
df_power.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")
print()

# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 80)
print("POWER ANALYSIS SUMMARY")
print("=" * 80)
print()

print(f"{'Analysis Type':<40} {'N':<8} {'Power':<10} {'Feasible':<10}")
print("-" * 72)
for _, row in df_power.iterrows():
    feasible_str = '✓ YES' if row['feasible'] else '⚠ LIMITED'
    print(f"{row['analysis_type']:<40} {row['sample_size']:<8} {row['power_estimate']:<10.2f} {feasible_str:<10}")

print()
print(f"FEASIBLE ANALYSES: {df_power['feasible'].sum()}/{len(df_power)}")
print()

feasible_list = df_power[df_power['feasible']]['analysis_type'].tolist()
if feasible_list:
    print("✓ Well-powered analyses:")
    for analysis in feasible_list:
        print(f"  • {analysis}")

print()

limited_list = df_power[~df_power['feasible']]['analysis_type'].tolist()
if limited_list:
    print("⚠ Limited power (proceed with caution):")
    for analysis in limited_list:
        print(f"  • {analysis}")

print()
print("=" * 80)
print("POWER ANALYSIS COMPLETE")
print("=" * 80)
