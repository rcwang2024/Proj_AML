"""
Phase 4: Data Quality Assessment
Tasks 4.3-4.6: Comprehensive QC for all data types

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/02_QC_Reports")
FIGURE_DIR = Path("D:/Projects/Project_AML/04_Figures/01_QC_Figures")

print("=" * 80)
print("COMPREHENSIVE DATA QUALITY ASSESSMENT (Tasks 4.3-4.6)")
print("=" * 80)
print()

# ============================================================================
# TASK 4.3: DRUG RESPONSE QC
# ============================================================================

print("TASK 4.3: DRUG RESPONSE DATA QUALITY")
print("-" * 80)

# Load previously generated summaries
drug_summary = pd.read_csv(OUTPUT_DIR.parent / "01_Processed_Data" / "drug_response_summary.csv")
drug_counts = pd.read_csv(OUTPUT_DIR.parent / "01_Processed_Data" / "samples_drug_counts.csv")

# Load raw drug data
df_drug = pd.read_csv(DATA_DIR / "beataml_drug_auc.txt", sep='\t')

# Check for extreme AUC values
auc_col = 'auc'
extreme_low = (df_drug[auc_col] < 0).sum()
extreme_high = (df_drug[auc_col] > 1000).sum()

print(f"Extreme AUC values:")
print(f"  AUC < 0: {extreme_low} ({extreme_low/len(df_drug)*100:.2f}%)")
print(f"  AUC > 1000: {extreme_high} ({extreme_high/len(df_drug)*100:.2f}%)")

# Drugs tested per sample
drug_col = 'n_drugs_tested' if 'n_drugs_tested' in drug_counts.columns else 'n_drugs'
print(f"\nDrugs tested per sample:")
print(f"  Mean: {drug_counts[drug_col].mean():.1f}")
print(f"  Median: {drug_counts[drug_col].median():.1f}")
print(f"  Min: {drug_counts[drug_col].min()}")
print(f"  Max: {drug_counts[drug_col].max()}")

# Samples with <50% drugs
threshold_50pct = 166 * 0.5
low_coverage = (drug_counts[drug_col] < threshold_50pct).sum()
print(f"  Samples with <50% drugs (<83): {low_coverage}")

# Drug completeness
drug_sample_counts = df_drug.groupby('inhibitor').size()
low_tested_drugs = (drug_sample_counts < 60).sum()  # <10% of 603 samples
print(f"\nDrug completeness:")
print(f"  Drugs tested in <10% samples (<60): {low_tested_drugs}")

# Save QC report
drug_qc = pd.DataFrame({
    'metric': ['Extreme AUC < 0', 'Extreme AUC > 1000', 'Mean drugs/sample',
               'Samples <50% drugs', 'Drugs in <10% samples'],
    'value': [extreme_low, extreme_high, drug_counts[drug_col].mean(),
              low_coverage, low_tested_drugs]
})
drug_qc.to_csv(OUTPUT_DIR / "drug_response_qc.csv", index=False)

# Visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.hist(drug_counts[drug_col], bins=30, edgecolor='black', alpha=0.7, color='coral')
ax.axvline(threshold_50pct, color='red', linestyle='--', label=f'50% threshold ({threshold_50pct:.0f})')
ax.set_xlabel('Number of Drugs Tested')
ax.set_ylabel('Number of Samples')
ax.set_title('Drugs Tested per Sample')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
top_20_drugs = drug_sample_counts.sort_values(ascending=False).head(20)
ax.barh(range(len(top_20_drugs)), top_20_drugs.values, color='steelblue')
ax.set_yticks(range(len(top_20_drugs)))
ax.set_yticklabels(top_20_drugs.index, fontsize=8)
ax.set_xlabel('Number of Samples')
ax.set_title('Top 20 Most Tested Drugs')
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURE_DIR / "drugs_per_sample_histogram.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"\n✓ Task 4.3 complete")
print()

# ============================================================================
# TASK 4.4: CLINICAL DATA COMPLETENESS
# ============================================================================

print("TASK 4.4: CLINICAL DATA COMPLETENESS")
print("-" * 80)

df_clinical = pd.read_excel(DATA_DIR / "beataml_clinical.xlsx")

# Calculate missingness
missing_pct = (df_clinical.isna().sum() / len(df_clinical) * 100).sort_values(ascending=False)

# Key variables - use actual column names
key_vars = {}
if 'overallSurvival' in df_clinical.columns:
    key_vars['overallSurvival'] = df_clinical['overallSurvival'].notna().sum() / len(df_clinical) * 100
if 'vitalStatus' in df_clinical.columns:
    key_vars['vitalStatus'] = df_clinical['vitalStatus'].notna().sum() / len(df_clinical) * 100
if 'ageAtDiagnosis' in df_clinical.columns:
    key_vars['ageAtDiagnosis'] = df_clinical['ageAtDiagnosis'].notna().sum() / len(df_clinical) * 100
if 'consensus_sex' in df_clinical.columns:
    key_vars['consensus_sex'] = df_clinical['consensus_sex'].notna().sum() / len(df_clinical) * 100
if 'FLT3_ITD' in df_clinical.columns:
    key_vars['FLT3_ITD'] = df_clinical['FLT3_ITD'].notna().sum() / len(df_clinical) * 100
if 'NPM1_mutationStatus' in df_clinical.columns:
    key_vars['NPM1_mutationStatus'] = df_clinical['NPM1_mutationStatus'].notna().sum() / len(df_clinical) * 100
if '2017_ELN_risk' in df_clinical.columns:
    key_vars['2017_ELN_risk'] = df_clinical['2017_ELN_risk'].notna().sum() / len(df_clinical) * 100

print("Key variables completeness:")
for var, pct in key_vars.items():
    print(f"  {var:20s}: {pct:5.1f}% complete")

# Samples with minimal clinical data
n_present = df_clinical.notna().sum(axis=1)
minimal_data = (n_present < 3).sum()
print(f"\nSamples with <3 variables: {minimal_data}")

# Save completeness
completeness_df = pd.DataFrame({
    'variable': missing_pct.index,
    'missing_pct': missing_pct.values,
    'complete_pct': 100 - missing_pct.values
})
completeness_df.to_csv(OUTPUT_DIR / "clinical_completeness.csv", index=False)

# Heatmap
fig, ax = plt.subplots(figsize=(12, 8))
# Select key variables for heatmap
key_cols = ['overallSurvival', 'vitalStatus', 'ageAtDiagnosis', 'consensus_sex', 'FLT3_ITD', 'NPM1_mutationStatus', '2017_ELN_risk']
key_cols = [c for c in key_cols if c in df_clinical.columns]
heatmap_data = df_clinical[key_cols].head(100).notna().astype(int)  # First 100 samples
sns.heatmap(heatmap_data.T, cmap=['red', 'green'], cbar_kws={'label': 'Data Available'},
           yticklabels=key_cols, xticklabels=False)
ax.set_title('Clinical Data Completeness (First 100 Samples)', fontweight='bold')
ax.set_xlabel('Samples')
plt.tight_layout()
plt.savefig(FIGURE_DIR / "clinical_completeness_heatmap.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Task 4.4 complete")
print()

# ============================================================================
# TASK 4.5: MUTATION DATA QC
# ============================================================================

print("TASK 4.5: MUTATION DATA QUALITY ASSESSMENT")
print("-" * 80)

df_mut = pd.read_csv(DATA_DIR / "beataml_mutations.txt", sep='\t')

# VAF distribution
vaf_col = 't_vaf'
vaf_values = df_mut[vaf_col].dropna()

print(f"VAF distribution:")
print(f"  Mean: {vaf_values.mean():.3f}")
print(f"  Median: {vaf_values.median():.3f}")
print(f"  Peak around 0.5 (heterozygous): {((vaf_values > 0.4) & (vaf_values < 0.6)).sum()} ({((vaf_values > 0.4) & (vaf_values < 0.6)).sum()/len(vaf_values)*100:.1f}%)")
print(f"  Low VAF (<0.1): {(vaf_values < 0.1).sum()} ({(vaf_values < 0.1).sum()/len(vaf_values)*100:.1f}%)")

# Mutation burden
burden = df_mut.groupby('dbgap_sample_id').size()
print(f"\nMutation burden:")
print(f"  Median: {burden.median():.0f}")
print(f"  Samples <5 mutations: {(burden < 5).sum()}")
print(f"  Samples >100 mutations: {(burden > 100).sum()}")

# Save QC
mut_qc = pd.DataFrame({
    'metric': ['Total mutations', 'Mean VAF', 'Median VAF', 'Low VAF (<0.1)',
               'Median burden', 'Samples <5 mut', 'Samples >100 mut'],
    'value': [len(df_mut), vaf_values.mean(), vaf_values.median(),
              (vaf_values < 0.1).sum(), burden.median(),
              (burden < 5).sum(), (burden > 100).sum()]
})
mut_qc.to_csv(OUTPUT_DIR / "mutation_data_qc.csv", index=False)

# Visualizations
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.hist(vaf_values, bins=50, edgecolor='black', alpha=0.7, color='green')
ax.axvline(0.5, color='red', linestyle='--', label='0.5 (het)')
ax.set_xlabel('Variant Allele Frequency')
ax.set_ylabel('Count')
ax.set_title('VAF Distribution')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.hist(burden, bins=50, edgecolor='black', alpha=0.7, color='purple')
ax.set_xlabel('Mutations per Sample')
ax.set_ylabel('Count')
ax.set_title('Mutation Burden Distribution')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(FIGURE_DIR / "mutation_burden_histogram.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Task 4.5 complete")
print()

# ============================================================================
# TASK 4.6: COMPREHENSIVE MISSING DATA ANALYSIS
# ============================================================================

print("TASK 4.6: COMPREHENSIVE MISSING DATA ANALYSIS")
print("-" * 80)

# Load master mapping
master_map = pd.read_csv(OUTPUT_DIR.parent / "01_Processed_Data" / "master_sample_id_mapping.csv")

# Calculate completeness
completeness_summary = []
for _, row in master_map.iterrows():
    n_types = row['n_data_types']
    has_all = (n_types == 4)

    # Check if any data type >30% missing (N/A for now as we track presence/absence)
    fail_threshold = False  # Placeholder

    completeness_summary.append({
        'sample_id': row['unified_sample_id'],
        'n_data_types': n_types,
        'has_expression': row['has_expression'],
        'has_drug': row['has_drug_response'],
        'has_clinical': row['has_clinical'],
        'has_mutations': row['has_mutations'],
        'fail_qc': fail_threshold or (n_types < 2)
    })

df_completeness = pd.DataFrame(completeness_summary)

# Summary
print(f"Sample completeness:")
print(f"  Complete (4 types): {(df_completeness['n_data_types'] == 4).sum()}")
print(f"  Partial (2-3 types): {((df_completeness['n_data_types'] >= 2) & (df_completeness['n_data_types'] < 4)).sum()}")
print(f"  Minimal (<2 types): {(df_completeness['n_data_types'] < 2).sum()}")
print(f"  Fail QC threshold: {df_completeness['fail_qc'].sum()}")

# Save
df_completeness.to_csv(OUTPUT_DIR / "missing_data_comprehensive_report.csv", index=False)

# Visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Data type availability
ax = axes[0, 0]
data_avail = df_completeness[['has_expression', 'has_drug', 'has_clinical', 'has_mutations']].sum()
ax.bar(range(len(data_avail)), data_avail.values, color=['steelblue', 'coral', 'green', 'purple'])
ax.set_xticks(range(len(data_avail)))
ax.set_xticklabels(['Expression', 'Drug', 'Clinical', 'Mutations'])
ax.set_ylabel('Number of Samples')
ax.set_title('Data Type Availability')
ax.grid(axis='y', alpha=0.3)

# N data types distribution
ax = axes[0, 1]
type_counts = df_completeness['n_data_types'].value_counts().sort_index()
ax.bar(type_counts.index, type_counts.values, color='skyblue', edgecolor='black')
ax.set_xlabel('Number of Data Types')
ax.set_ylabel('Number of Samples')
ax.set_title('Sample Distribution by Data Types')
ax.grid(axis='y', alpha=0.3)

# Completeness heatmap
ax = axes[1, 0]
heatmap_data = df_completeness[['has_expression', 'has_drug', 'has_clinical', 'has_mutations']].head(100).astype(int)
sns.heatmap(heatmap_data.T, cmap=['red', 'green'], cbar=False, ax=ax,
           yticklabels=['Expression', 'Drug', 'Clinical', 'Mutations'], xticklabels=False)
ax.set_title('Sample Completeness (First 100)')

# QC status
ax = axes[1, 1]
qc_status = df_completeness['fail_qc'].value_counts()
ax.pie(qc_status.values, labels=['Pass', 'Fail'], autopct='%1.1f%%', colors=['green', 'red'])
ax.set_title('QC Status')

plt.tight_layout()
plt.savefig(FIGURE_DIR / "sample_data_completeness.png", dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Task 4.6 complete")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("COMPREHENSIVE QC COMPLETE (Tasks 4.3-4.6)")
print("=" * 80)
print()
print("OUTPUTS SAVED:")
print("  Task 4.3: drug_response_qc.csv, drugs_per_sample_histogram.png")
print("  Task 4.4: clinical_completeness.csv, clinical_completeness_heatmap.png")
print("  Task 4.5: mutation_data_qc.csv, mutation_burden_histogram.png")
print("  Task 4.6: missing_data_comprehensive_report.csv, sample_data_completeness.png")
print()
print("=" * 80)
