"""
Phase 3: Sample Overlap and Integration Strategy
Tasks 3.1 & 3.2: Sample overlap analysis and master ID mapping

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
from itertools import combinations
from matplotlib_venn import venn3, venn3_circles
warnings.filterwarnings('ignore')

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Paths
DATA_DIR = Path("D:/Projects/Project_AML/01_Data/BeatAML_Downloaded_Data")
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/01_Processed_Data")
FIGURE_DIR = Path("D:/Projects/Project_AML/04_Figures/01_QC_Figures")

# Create directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("PHASE 3: SAMPLE OVERLAP AND INTEGRATION STRATEGY")
print("=" * 80)
print()

# ============================================================================
# STEP 1: LOAD ALL DATASETS AND EXTRACT SAMPLE IDs
# ============================================================================

print("STEP 1: LOADING ALL DATASETS")
print("-" * 80)

# 1. Expression data
print("Loading expression data...")
expr_file = DATA_DIR / "beataml_expression.txt"
df_expr = pd.read_csv(expr_file, sep='\t', nrows=1)  # Just load header
expr_samples = [col for col in df_expr.columns if col.startswith('BA')]
print(f"  ✓ Expression: {len(expr_samples)} samples")
print(f"    ID format: RNA samples (BA####R)")
print(f"    Example: {expr_samples[:3]}")

# 2. Drug response data
print("\nLoading drug response data...")
drug_file = DATA_DIR / "beataml_drug_auc.txt"
df_drug = pd.read_csv(drug_file, sep='\t')
# Get DNA sample IDs
drug_samples_dna = df_drug['dbgap_dnaseq_sample'].dropna().unique().tolist()
# Get RNA sample IDs (if available)
drug_samples_rna = df_drug['dbgap_rnaseq_sample'].dropna().unique().tolist()
print(f"  ✓ Drug response (DNA IDs): {len(drug_samples_dna)} samples")
print(f"    ID format: DNA samples (BA####D)")
print(f"    Example: {drug_samples_dna[:3]}")

# 3. Clinical data
print("\nLoading clinical data...")
clinical_file = DATA_DIR / "beataml_clinical.xlsx"
df_clinical = pd.read_excel(clinical_file)
# Get both RNA and DNA IDs
clinical_samples_rna = df_clinical['dbgap_rnaseq_sample'].dropna().unique().tolist()
clinical_samples_dna = df_clinical['dbgap_dnaseq_sample'].dropna().unique().tolist()
print(f"  ✓ Clinical (RNA IDs): {len(clinical_samples_rna)} samples")
print(f"  ✓ Clinical (DNA IDs): {len(clinical_samples_dna)} samples")
print(f"    Example RNA: {clinical_samples_rna[:3]}")
print(f"    Example DNA: {clinical_samples_dna[:3]}")

# 4. Mutation data
print("\nLoading mutation data...")
mut_file = DATA_DIR / "beataml_mutations.txt"
df_mut = pd.read_csv(mut_file, sep='\t')
mut_samples = df_mut['dbgap_sample_id'].unique().tolist()
print(f"  ✓ Mutations: {len(mut_samples)} samples")
print(f"    ID format: DNA samples (BA####D)")
print(f"    Example: {mut_samples[:3]}")

print()

# ============================================================================
# STEP 2: CREATE RNA-DNA MAPPING
# ============================================================================

print("STEP 2: CREATING RNA-DNA ID MAPPING")
print("-" * 80)

def rna_to_dna(rna_id):
    """Convert RNA ID (BA####R) to DNA ID (BA####D)"""
    if isinstance(rna_id, str) and rna_id.endswith('R'):
        return rna_id[:-1] + 'D'
    return rna_id

def dna_to_rna(dna_id):
    """Convert DNA ID (BA####D) to RNA ID (BA####R)"""
    if isinstance(dna_id, str) and dna_id.endswith('D'):
        return dna_id[:-1] + 'R'
    return dna_id

# Create mapping from clinical file (has both RNA and DNA)
rna_dna_map = {}
dna_rna_map = {}

for _, row in df_clinical.iterrows():
    rna_id = row['dbgap_rnaseq_sample']
    dna_id = row['dbgap_dnaseq_sample']

    if pd.notna(rna_id) and pd.notna(dna_id):
        rna_dna_map[rna_id] = dna_id
        dna_rna_map[dna_id] = rna_id

print(f"RNA-DNA mapping from clinical file:")
print(f"  {len(rna_dna_map)} RNA → DNA mappings")
print(f"  {len(dna_rna_map)} DNA → RNA mappings")
print()

# ============================================================================
# STEP 3: CREATE UNIFIED SAMPLE SET
# ============================================================================

print("STEP 3: CREATING UNIFIED SAMPLE SET")
print("-" * 80)

# Collect all unique sample IDs (using RNA as primary)
all_samples = set()

# Add expression samples (already RNA format)
all_samples.update(expr_samples)

# Add drug response samples (convert DNA to RNA)
for dna_id in drug_samples_dna:
    rna_id = dna_rna_map.get(dna_id, dna_to_rna(dna_id))
    all_samples.add(rna_id)

# Add clinical samples (already have RNA)
all_samples.update(clinical_samples_rna)

# Add mutation samples (convert DNA to RNA)
for dna_id in mut_samples:
    rna_id = dna_rna_map.get(dna_id, dna_to_rna(dna_id))
    all_samples.add(rna_id)

all_samples = sorted(list(all_samples))
print(f"Total unique samples (RNA ID format): {len(all_samples)}")
print()

# ============================================================================
# STEP 4: CREATE MASTER SAMPLE MAPPING TABLE
# ============================================================================

print("STEP 4: CREATING MASTER SAMPLE MAPPING TABLE")
print("-" * 80)

# Convert to sets for fast lookup
expr_set = set(expr_samples)
drug_dna_set = set(drug_samples_dna)
clinical_rna_set = set(clinical_samples_rna)
clinical_dna_set = set(clinical_samples_dna)
mut_set = set(mut_samples)

master_mapping = []

for rna_id in all_samples:
    # Get corresponding DNA ID
    dna_id = rna_dna_map.get(rna_id, rna_to_dna(rna_id))

    # Check presence in each dataset
    has_expr = rna_id in expr_set
    has_drug = dna_id in drug_dna_set
    has_clinical = rna_id in clinical_rna_set or dna_id in clinical_dna_set
    has_mut = dna_id in mut_set

    # Count data types
    n_data_types = sum([has_expr, has_drug, has_clinical, has_mut])

    # Determine cohort category
    if n_data_types == 4:
        cohort = "complete_quad_omics"
    elif n_data_types == 3:
        if has_expr and has_drug and has_clinical:
            cohort = "triple_expr_drug_clin"
        elif has_expr and has_drug and has_mut:
            cohort = "triple_expr_drug_mut"
        elif has_expr and has_clinical and has_mut:
            cohort = "triple_expr_clin_mut"
        elif has_drug and has_clinical and has_mut:
            cohort = "triple_drug_clin_mut"
        else:
            cohort = "triple_omics"
    elif n_data_types == 2:
        cohort = "dual_omics"
    elif n_data_types == 1:
        cohort = "single_omics"
    else:
        cohort = "no_data"

    master_mapping.append({
        'unified_sample_id': rna_id,
        'expression_id': rna_id if has_expr else np.nan,
        'drug_response_id': dna_id if has_drug else np.nan,
        'clinical_id': rna_id if has_clinical else np.nan,
        'mutation_id': dna_id if has_mut else np.nan,
        'has_expression': has_expr,
        'has_drug_response': has_drug,
        'has_clinical': has_clinical,
        'has_mutations': has_mut,
        'n_data_types': n_data_types,
        'cohort_category': cohort
    })

df_master = pd.DataFrame(master_mapping)

print(f"Master mapping table created:")
print(f"  Total samples: {len(df_master)}")
print(f"  Columns: {len(df_master.columns)}")
print()

# Summary by cohort
print("COHORT SUMMARY:")
cohort_counts = df_master['cohort_category'].value_counts()
for cohort, count in cohort_counts.items():
    pct = count / len(df_master) * 100
    print(f"  {cohort:30s}: {count:4d} ({pct:5.1f}%)")
print()

# ============================================================================
# STEP 5: DETAILED OVERLAP ANALYSIS (15 COMBINATIONS)
# ============================================================================

print("STEP 5: CALCULATING ALL 15 OVERLAP COMBINATIONS")
print("-" * 80)

overlap_results = []

# Single data types (4 combinations)
single_combos = [
    ('Expression only', [True, False, False, False]),
    ('Drug Response only', [False, True, False, False]),
    ('Clinical only', [False, False, True, False]),
    ('Mutations only', [False, False, False, True])
]

for name, (e, d, c, m) in single_combos:
    count = len(df_master[
        (df_master['has_expression'] == e) &
        (df_master['has_drug_response'] == d) &
        (df_master['has_clinical'] == c) &
        (df_master['has_mutations'] == m)
    ])
    overlap_results.append({'combination': name, 'n_samples': count, 'type': 'single'})
    print(f"  {name:40s}: {count:4d}")

print()

# Dual data types (6 combinations)
dual_combos = [
    ('Expression + Drug Response', [True, True, False, False]),
    ('Expression + Clinical', [True, False, True, False]),
    ('Expression + Mutations', [True, False, False, True]),
    ('Drug Response + Clinical', [False, True, True, False]),
    ('Drug Response + Mutations', [False, True, False, True]),
    ('Clinical + Mutations', [False, False, True, True])
]

for name, (e, d, c, m) in dual_combos:
    # At least these two, regardless of others
    count = len(df_master[
        (df_master['has_expression'] == e) &
        (df_master['has_drug_response'] == d) &
        (df_master['has_clinical'] == c) &
        (df_master['has_mutations'] == m)
    ])
    overlap_results.append({'combination': name, 'n_samples': count, 'type': 'dual'})
    print(f"  {name:40s}: {count:4d}")

print()

# Triple data types (4 combinations)
triple_combos = [
    ('Expression + Drug Response + Clinical', [True, True, True, False]),
    ('Expression + Drug Response + Mutations', [True, True, False, True]),
    ('Expression + Clinical + Mutations', [True, False, True, True]),
    ('Drug Response + Clinical + Mutations', [False, True, True, True])
]

for name, (e, d, c, m) in triple_combos:
    count = len(df_master[
        (df_master['has_expression'] == e) &
        (df_master['has_drug_response'] == d) &
        (df_master['has_clinical'] == c) &
        (df_master['has_mutations'] == m)
    ])
    overlap_results.append({'combination': name, 'n_samples': count, 'type': 'triple'})
    print(f"  {name:40s}: {count:4d}")

print()

# All four (1 combination) - GOLD STANDARD
quad_count = len(df_master[
    (df_master['has_expression'] == True) &
    (df_master['has_drug_response'] == True) &
    (df_master['has_clinical'] == True) &
    (df_master['has_mutations'] == True)
])
overlap_results.append({
    'combination': 'ALL FOUR (Gold Standard)',
    'n_samples': quad_count,
    'type': 'quad'
})
print(f"  {'ALL FOUR (Gold Standard)':40s}: {quad_count:4d} ⭐")

print()

df_overlap = pd.DataFrame(overlap_results)

# ============================================================================
# STEP 6: CALCULATE AT-LEAST COMBINATIONS
# ============================================================================

print("STEP 6: CALCULATING 'AT LEAST' COMBINATIONS")
print("-" * 80)

# Samples with at least each data type
at_least = {
    'Expression': (df_master['has_expression'] == True).sum(),
    'Drug Response': (df_master['has_drug_response'] == True).sum(),
    'Clinical': (df_master['has_clinical'] == True).sum(),
    'Mutations': (df_master['has_mutations'] == True).sum()
}

for data_type, count in at_least.items():
    pct = count / len(df_master) * 100
    print(f"  At least {data_type:15s}: {count:4d} ({pct:5.1f}%)")

print()

# ============================================================================
# STEP 7: CREATE VISUALIZATIONS
# ============================================================================

print("STEP 7: CREATING VISUALIZATIONS")
print("-" * 80)

# Figure 1: UpSet plot (best for 4-way)
try:
    from upsetplot import from_memberships, UpSet

    # Prepare data for UpSet plot
    memberships = []
    for _, row in df_master.iterrows():
        member = []
        if row['has_expression']:
            member.append('Expression')
        if row['has_drug_response']:
            member.append('Drug Response')
        if row['has_clinical']:
            member.append('Clinical')
        if row['has_mutations']:
            member.append('Mutations')
        if member:  # Only add if at least one data type
            memberships.append(tuple(member))

    upset_data = from_memberships(memberships)

    upset = UpSet(upset_data,
                  subset_size='count',
                  show_counts=True,
                  sort_by='cardinality',
                  sort_categories_by=None)

    upset.plot()
    plt.suptitle('Sample Overlap Across Four Data Types\n(UpSet Plot)',
                 fontsize=14, fontweight='bold', y=1.02)

    output_file = FIGURE_DIR / "sample_overlap_upset.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file.name}")
    plt.close()

except ImportError:
    print("  ⚠ upsetplot not installed - skipping UpSet plot")
    print("    Install with: pip install upsetplot")

# Figure 2: Venn diagram (3-way: Expression, Drug, Clinical + Mutations as note)
print("\nCreating Venn diagram (3-way)...")

# Get sample sets
expr_samples_set = set(df_master[df_master['has_expression']]['unified_sample_id'])
drug_samples_set = set(df_master[df_master['has_drug_response']]['unified_sample_id'])
clinical_samples_set = set(df_master[df_master['has_clinical']]['unified_sample_id'])
mut_samples_set = set(df_master[df_master['has_mutations']]['unified_sample_id'])

plt.figure(figsize=(12, 8))
venn = venn3([expr_samples_set, drug_samples_set, clinical_samples_set],
             set_labels=('Expression\n(n={})'.format(len(expr_samples_set)),
                        'Drug Response\n(n={})'.format(len(drug_samples_set)),
                        'Clinical\n(n={})'.format(len(clinical_samples_set))))

# Add circles
venn3_circles([expr_samples_set, drug_samples_set, clinical_samples_set],
              linewidth=1.5)

plt.title('Sample Overlap: Expression, Drug Response, and Clinical Data\n' +
          f'Mutations: {len(mut_samples_set)} samples (see UpSet plot for full 4-way overlap)',
          fontsize=14, fontweight='bold', pad=20)

# Add note about gold standard cohort
plt.text(0.5, -0.15, f'Gold Standard Cohort (All 4 data types): {quad_count} samples',
         ha='center', transform=plt.gca().transAxes,
         fontsize=12, fontweight='bold',
         bbox=dict(boxstyle='round', facecolor='gold', alpha=0.5))

output_file = FIGURE_DIR / "sample_overlap_venn.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_file.name}")
plt.close()

# Figure 3: Bar plot of all combinations
plt.figure(figsize=(14, 10))

# Sort by type then by count
df_overlap_sorted = df_overlap.sort_values(['type', 'n_samples'], ascending=[True, False])

colors = {'single': 'lightblue', 'dual': 'lightgreen', 'triple': 'orange', 'quad': 'gold'}
bar_colors = [colors[t] for t in df_overlap_sorted['type']]

plt.barh(range(len(df_overlap_sorted)), df_overlap_sorted['n_samples'],
         color=bar_colors, edgecolor='black', linewidth=1)
plt.yticks(range(len(df_overlap_sorted)), df_overlap_sorted['combination'])
plt.xlabel('Number of Samples', fontsize=12)
plt.title('Sample Counts for All Data Type Combinations\n(15 combinations total)',
          fontsize=14, fontweight='bold', pad=20)
plt.grid(axis='x', alpha=0.3)

# Add value labels
for i, (idx, row) in enumerate(df_overlap_sorted.iterrows()):
    plt.text(row['n_samples'] + 5, i, str(row['n_samples']),
             va='center', fontweight='bold')

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='lightblue', edgecolor='black', label='Single data type'),
    Patch(facecolor='lightgreen', edgecolor='black', label='Dual data types'),
    Patch(facecolor='orange', edgecolor='black', label='Triple data types'),
    Patch(facecolor='gold', edgecolor='black', label='Quad (Gold Standard)')
]
plt.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
output_file = FIGURE_DIR / "sample_overlap_barplot.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_file.name}")
plt.close()

# Figure 4: Summary statistics
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Pie chart: n_data_types distribution
ax = axes[0, 0]
data_type_counts = df_master['n_data_types'].value_counts().sort_index()
ax.pie(data_type_counts.values, labels=[f'{i} data types' for i in data_type_counts.index],
       autopct='%1.1f%%', startangle=90, colors=['lightgray', 'lightblue', 'lightgreen', 'orange', 'gold'])
ax.set_title('Distribution by Number of Data Types', fontweight='bold')

# Bar chart: Individual data type coverage
ax = axes[0, 1]
coverage_data = pd.DataFrame({
    'Data Type': ['Expression', 'Drug\nResponse', 'Clinical', 'Mutations'],
    'Samples': [at_least['Expression'], at_least['Drug Response'],
                at_least['Clinical'], at_least['Mutations']],
    'Percentage': [at_least['Expression']/len(df_master)*100,
                   at_least['Drug Response']/len(df_master)*100,
                   at_least['Clinical']/len(df_master)*100,
                   at_least['Mutations']/len(df_master)*100]
})
bars = ax.bar(coverage_data['Data Type'], coverage_data['Samples'],
              color=['steelblue', 'coral', 'green', 'purple'], alpha=0.7, edgecolor='black')
ax.set_ylabel('Number of Samples', fontweight='bold')
ax.set_title('Coverage by Data Type', fontweight='bold')
ax.grid(axis='y', alpha=0.3)
# Add percentage labels
for bar, pct in zip(bars, coverage_data['Percentage']):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{pct:.1f}%', ha='center', va='bottom', fontweight='bold')

# Heatmap: Cohort categories
ax = axes[1, 0]
cohort_summary = df_master['cohort_category'].value_counts().sort_values(ascending=True)
ax.barh(range(len(cohort_summary)), cohort_summary.values, color='skyblue', edgecolor='black')
ax.set_yticks(range(len(cohort_summary)))
ax.set_yticklabels(cohort_summary.index)
ax.set_xlabel('Number of Samples', fontweight='bold')
ax.set_title('Sample Counts by Cohort Category', fontweight='bold')
ax.grid(axis='x', alpha=0.3)
for i, val in enumerate(cohort_summary.values):
    ax.text(val + 2, i, str(val), va='center', fontweight='bold')

# Table: Key statistics
ax = axes[1, 1]
ax.axis('tight')
ax.axis('off')
table_data = [
    ['Metric', 'Value'],
    ['Total unique samples', f'{len(df_master):,}'],
    ['Gold standard (4 data types)', f'{quad_count:,}'],
    ['At least 3 data types', f'{(df_master["n_data_types"] >= 3).sum():,}'],
    ['At least 2 data types', f'{(df_master["n_data_types"] >= 2).sum():,}'],
    ['Expression coverage', f'{at_least["Expression"]:,} ({at_least["Expression"]/len(df_master)*100:.1f}%)'],
    ['Drug response coverage', f'{at_least["Drug Response"]:,} ({at_least["Drug Response"]/len(df_master)*100:.1f}%)'],
    ['Clinical coverage', f'{at_least["Clinical"]:,} ({at_least["Clinical"]/len(df_master)*100:.1f}%)'],
    ['Mutations coverage', f'{at_least["Mutations"]:,} ({at_least["Mutations"]/len(df_master)*100:.1f}%)']
]
table = ax.table(cellText=table_data, cellLoc='left', loc='center',
                colWidths=[0.6, 0.4])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)
# Style header row
for i in range(2):
    table[(0, i)].set_facecolor('#4CAF50')
    table[(0, i)].set_text_props(weight='bold', color='white')
ax.set_title('Key Statistics', fontweight='bold', pad=20)

plt.suptitle('Sample Overlap Analysis Summary', fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()
output_file = FIGURE_DIR / "sample_overlap_summary.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_file.name}")
plt.close()

print()

# ============================================================================
# STEP 8: SAVE OUTPUT FILES
# ============================================================================

print("STEP 8: SAVING OUTPUT FILES")
print("-" * 80)

# 1. Master sample ID mapping
output_file = OUTPUT_DIR / "master_sample_id_mapping.csv"
df_master.to_csv(output_file, index=False)
print(f"  ✓ Saved: {output_file.name}")
print(f"    Rows: {len(df_master):,}")
print(f"    Columns: {len(df_master.columns)}")

# 2. Sample overlap analysis table
output_file = OUTPUT_DIR / "sample_overlap_analysis.csv"
df_overlap.to_csv(output_file, index=False)
print(f"  ✓ Saved: {output_file.name}")
print(f"    Combinations: {len(df_overlap)}")

# 3. Summary statistics
summary_stats = {
    'metric': [
        'Total unique samples',
        'Samples with expression',
        'Samples with drug response',
        'Samples with clinical',
        'Samples with mutations',
        'Gold standard (all 4)',
        'At least 3 data types',
        'At least 2 data types',
        'Expression-only samples',
        'Drug-only samples',
        'Clinical-only samples',
        'Mutations-only samples'
    ],
    'value': [
        len(df_master),
        at_least['Expression'],
        at_least['Drug Response'],
        at_least['Clinical'],
        at_least['Mutations'],
        quad_count,
        (df_master['n_data_types'] >= 3).sum(),
        (df_master['n_data_types'] >= 2).sum(),
        len(df_master[(df_master['has_expression']) & (df_master['n_data_types'] == 1)]),
        len(df_master[(df_master['has_drug_response']) & (df_master['n_data_types'] == 1)]),
        len(df_master[(df_master['has_clinical']) & (df_master['n_data_types'] == 1)]),
        len(df_master[(df_master['has_mutations']) & (df_master['n_data_types'] == 1)])
    ]
}
df_summary = pd.DataFrame(summary_stats)
output_file = OUTPUT_DIR / "sample_overlap_summary_stats.csv"
df_summary.to_csv(output_file, index=False)
print(f"  ✓ Saved: {output_file.name}")

# 4. Gold standard cohort sample list
gold_standard_samples = df_master[
    (df_master['has_expression']) &
    (df_master['has_drug_response']) &
    (df_master['has_clinical']) &
    (df_master['has_mutations'])
]['unified_sample_id'].tolist()

output_file = OUTPUT_DIR / "gold_standard_cohort_samples.txt"
with open(output_file, 'w') as f:
    for sample in gold_standard_samples:
        f.write(f"{sample}\n")
print(f"  ✓ Saved: {output_file.name}")
print(f"    Samples: {len(gold_standard_samples)}")

print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("PHASE 3 COMPLETE: SAMPLE OVERLAP AND INTEGRATION STRATEGY")
print("=" * 80)
print()
print("TASK 3.1: SAMPLE OVERLAP ANALYSIS ✓")
print(f"  • Total unique samples: {len(df_master):,}")
print(f"  • All 15 combinations calculated")
print(f"  • Gold standard cohort: {quad_count:,} samples")
print()
print("TASK 3.2: MASTER SAMPLE ID MAPPING ✓")
print(f"  • Unified sample IDs: {len(df_master):,}")
print(f"  • RNA-DNA mapping: {len(rna_dna_map):,} pairs")
print(f"  • Cohort categories: {len(cohort_counts)} types")
print()
print("DATA FILES SAVED (4):")
print("  • master_sample_id_mapping.csv")
print("  • sample_overlap_analysis.csv")
print("  • sample_overlap_summary_stats.csv")
print("  • gold_standard_cohort_samples.txt")
print()
print("FIGURES SAVED (4):")
print("  • sample_overlap_upset.png (if upsetplot available)")
print("  • sample_overlap_venn.png")
print("  • sample_overlap_barplot.png")
print("  • sample_overlap_summary.png")
print()
print("GOLD STANDARD COHORT:")
print(f"  🌟 {quad_count} samples with all 4 data types")
print(f"     Expression + Drug Response + Clinical + Mutations")
print()
print("=" * 80)
