"""
Task 2.4: Mutation Data Analysis
Comprehensive analysis of somatic mutation data

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
from collections import Counter
warnings.filterwarnings('ignore')

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Paths
DATA_DIR = Path("D:/Projects/Project_AML/01_Data/BeatAML_Downloaded_Data")
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/01_Processed_Data")
FIGURE_DIR = OUTPUT_DIR

# Create output directory if needed
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input file
mut_file = DATA_DIR / "beataml_mutations.txt"

print("=" * 80)
print("TASK 2.4: MUTATION DATA ANALYSIS")
print("=" * 80)
print(f"\nAnalyzing: {mut_file.name}")
print(f"File size: {mut_file.stat().st_size / (1024**2):.2f} MB")
print()

# ============================================================================
# 1. LOAD DATA AND BASIC STRUCTURE
# ============================================================================

print("1. LOADING MUTATION DATA...")
print("-" * 80)

# Load mutation data
df_mut = pd.read_csv(mut_file, sep='\t', low_memory=False)

print(f"✓ Data loaded successfully")
print(f"  Shape: {df_mut.shape[0]:,} rows × {df_mut.shape[1]:,} columns")
print()

# Display columns
print(f"Columns present ({len(df_mut.columns)}):")
for i, col in enumerate(df_mut.columns, 1):
    print(f"  {i:2d}. {col}")
print()

# Determine file format
print("FILE FORMAT ANALYSIS:")
if 'Hugo_Symbol' in df_mut.columns and 'Variant_Classification' in df_mut.columns:
    file_format = "MAF-like format"
elif 'gene' in df_mut.columns or 'gene_symbol' in df_mut.columns:
    file_format = "Custom mutation format"
elif 'CHROM' in df_mut.columns and 'POS' in df_df.columns:
    file_format = "VCF-like format"
else:
    file_format = "Custom format"
print(f"  Inferred format: {file_format}")
print()

# ============================================================================
# 2. SAMPLE AND GENE COVERAGE
# ============================================================================

print("2. SAMPLE AND GENE COVERAGE")
print("-" * 80)

# Identify sample ID column
sample_cols = [col for col in df_mut.columns if 'sample' in col.lower() or
               'dbgap' in col.lower() or 'tumor' in col.lower()]
if sample_cols:
    sample_col = sample_cols[0]
elif 'Tumor_Sample_Barcode' in df_mut.columns:
    sample_col = 'Tumor_Sample_Barcode'
else:
    # Try to find BA#### pattern
    for col in df_mut.columns:
        if df_mut[col].astype(str).str.contains(r'BA\d{4}', na=False).any():
            sample_col = col
            break

print(f"Sample ID column: {sample_col}")

# Identify gene column - prefer 'symbol' over 'gene' for gene names
if 'symbol' in df_mut.columns:
    gene_col = 'symbol'
elif 'Hugo_Symbol' in df_mut.columns:
    gene_col = 'Hugo_Symbol'
else:
    gene_cols = [col for col in df_mut.columns if 'gene' in col.lower()]
    gene_col = gene_cols[0] if gene_cols else df_mut.columns[0]

print(f"Gene column: {gene_col}")
print()

# Count unique samples
unique_samples = df_mut[sample_col].nunique()
sample_list = df_mut[sample_col].unique()

print(f"SAMPLE COVERAGE:")
print(f"  Unique samples with mutations: {unique_samples:,}")
print(f"  Total mutation calls: {len(df_mut):,}")
print(f"  Mean mutations per sample: {len(df_mut)/unique_samples:.1f}")
print()

# Count unique genes
unique_genes = df_mut[gene_col].nunique()
gene_list = df_mut[gene_col].unique()

print(f"GENE COVERAGE:")
print(f"  Unique genes with mutations: {unique_genes:,}")
print()

# Determine panel type
if unique_genes < 100:
    panel_type = "Targeted panel (small)"
elif unique_genes < 500:
    panel_type = "Targeted panel (medium)"
elif unique_genes < 1000:
    panel_type = "Large targeted panel or gene panel"
else:
    panel_type = "Whole Exome Sequencing (WES) or comprehensive panel"

print(f"  Inferred sequencing type: {panel_type}")
print()

# ============================================================================
# 3. MUTATION TYPE ANALYSIS
# ============================================================================

print("3. MUTATION TYPE ANALYSIS")
print("-" * 80)

# Find mutation type column - prefer variant_classification
if 'variant_classification' in df_mut.columns:
    mut_type_col = 'variant_classification'
elif 'Variant_Classification' in df_mut.columns:
    mut_type_col = 'Variant_Classification'
else:
    type_cols = [col for col in df_mut.columns if 'type' in col.lower() or 'class' in col.lower()]
    mut_type_col = type_cols[0] if type_cols else None

if mut_type_col:
    mut_type_col = mut_type_col
    print(f"Mutation type column: {mut_type_col}")
    print()

    # Count mutation types
    mut_types = df_mut[mut_type_col].value_counts()
    print(f"MUTATION TYPES (n={len(mut_types)}):")
    for mut_type, count in mut_types.items():
        pct = count / len(df_mut) * 100
        print(f"  {mut_type}: {count:,} ({pct:.1f}%)")
    print()

    # Categorize into SNV, Indel, etc.
    snv_keywords = ['SNP', 'SNV', 'Missense', 'Nonsense', 'Silent', 'substitution']
    indel_keywords = ['Indel', 'Insertion', 'Deletion', 'Frame_Shift', 'In_Frame']

    mut_type_str = df_mut[mut_type_col].astype(str)

    n_snv = mut_type_str.str.contains('|'.join(snv_keywords), case=False, na=False).sum()
    n_indel = mut_type_str.str.contains('|'.join(indel_keywords), case=False, na=False).sum()
    n_other = len(df_mut) - n_snv - n_indel

    print("CATEGORIZED MUTATION TYPES:")
    print(f"  SNVs (single nucleotide variants): {n_snv:,} ({n_snv/len(df_mut)*100:.1f}%)")
    print(f"  Indels (insertions/deletions): {n_indel:,} ({n_indel/len(df_mut)*100:.1f}%)")
    print(f"  Other/Structural: {n_other:,} ({n_other/len(df_mut)*100:.1f}%)")
    print()
else:
    print("No mutation type column found")
    mut_type_col = None
    print()

# ============================================================================
# 4. VARIANT ALLELE FREQUENCY (VAF) ANALYSIS
# ============================================================================

print("4. VARIANT ALLELE FREQUENCY (VAF) ANALYSIS")
print("-" * 80)

# Find VAF column - prefer t_vaf (tumor VAF)
if 't_vaf' in df_mut.columns:
    vaf_col = 't_vaf'
elif 'vaf' in df_mut.columns:
    vaf_col = 'vaf'
else:
    vaf_cols = [col for col in df_mut.columns if 'vaf' in col.lower() or 'freq' in col.lower()]
    vaf_col = vaf_cols[0] if vaf_cols else None

if vaf_col:
    vaf_col = vaf_col
    print(f"VAF column found: {vaf_col}")

    # Convert to numeric
    vaf_values = pd.to_numeric(df_mut[vaf_col], errors='coerce')
    vaf_clean = vaf_values.dropna()

    if len(vaf_clean) > 0:
        print(f"\nVAF STATISTICS:")
        print(f"  N values: {len(vaf_clean):,}")
        print(f"  Mean VAF: {vaf_clean.mean():.3f}")
        print(f"  Median VAF: {vaf_clean.median():.3f}")
        print(f"  SD: {vaf_clean.std():.3f}")
        print(f"  Range: [{vaf_clean.min():.3f}, {vaf_clean.max():.3f}]")
        print()

        # VAF distribution categories
        low_vaf = (vaf_clean < 0.1).sum()
        mid_vaf = ((vaf_clean >= 0.1) & (vaf_clean < 0.4)).sum()
        high_vaf = (vaf_clean >= 0.4).sum()

        print(f"VAF DISTRIBUTION:")
        print(f"  Low VAF (<10%): {low_vaf:,} ({low_vaf/len(vaf_clean)*100:.1f}%)")
        print(f"  Mid VAF (10-40%): {mid_vaf:,} ({mid_vaf/len(vaf_clean)*100:.1f}%)")
        print(f"  High VAF (≥40%): {high_vaf:,} ({high_vaf/len(vaf_clean)*100:.1f}%)")
        print()
    else:
        vaf_col = None
        print("  VAF column found but no valid numeric values")
        print()
else:
    vaf_col = None
    print("No VAF column found")
    print()

# ============================================================================
# 5. AML DRIVER MUTATION ANALYSIS
# ============================================================================

print("5. AML DRIVER MUTATION ANALYSIS")
print("-" * 80)

# Define AML driver genes
driver_genes = {
    'FLT3': ['ITD', 'TKD'],  # FLT3 has specific subtypes
    'NPM1': [],
    'DNMT3A': [],
    'IDH1': [],
    'IDH2': [],
    'TET2': [],
    'RUNX1': [],
    'TP53': [],
    'NRAS': [],
    'KRAS': [],
    'CEBPA': [],
    'WT1': [],
    'ASXL1': [],
    'SRSF2': [],
    'STAG2': [],
    'KIT': [],
    'PTPN11': [],
    'CBL': [],
    'PHF6': [],
    'BCOR': []
}

driver_results = []

for gene, subtypes in driver_genes.items():
    # Get mutations for this gene
    gene_muts = df_mut[df_mut[gene_col].str.upper() == gene.upper()]

    if len(gene_muts) > 0:
        n_samples_mutated = gene_muts[sample_col].nunique()
        freq = n_samples_mutated / unique_samples * 100
        n_mutations = len(gene_muts)

        # Get mutation types for this gene
        if mut_type_col:
            mut_types_gene = gene_muts[mut_type_col].value_counts().to_dict()
            top_mut_type = gene_muts[mut_type_col].value_counts().index[0] if len(gene_muts) > 0 else 'N/A'
        else:
            mut_types_gene = {}
            top_mut_type = 'N/A'

        # Get VAF for this gene
        if vaf_col:
            gene_vaf = pd.to_numeric(gene_muts[vaf_col], errors='coerce').dropna()
            mean_vaf = gene_vaf.mean() if len(gene_vaf) > 0 else np.nan
            median_vaf = gene_vaf.median() if len(gene_vaf) > 0 else np.nan
        else:
            mean_vaf = np.nan
            median_vaf = np.nan

        # Special handling for FLT3 subtypes
        if gene == 'FLT3':
            # Look for ITD and TKD in mutation description columns
            desc_cols = [col for col in gene_muts.columns if 'desc' in col.lower() or
                        'protein' in col.lower() or 'hgvs' in col.lower()]

            itd_count = 0
            tkd_count = 0

            for col in desc_cols:
                itd_count += gene_muts[col].astype(str).str.contains('ITD', case=False, na=False).sum()
                tkd_count += gene_muts[col].astype(str).str.contains('TKD|D835', case=False, na=False).sum()

            if itd_count > 0 or tkd_count > 0:
                print(f"  {gene}: {n_samples_mutated} samples ({freq:.1f}%)")
                if itd_count > 0:
                    print(f"    - ITD: {itd_count} mutations")
                if tkd_count > 0:
                    print(f"    - TKD: {tkd_count} mutations")
            else:
                print(f"  {gene}: {n_samples_mutated} samples ({freq:.1f}%) - subtype not specified")
        else:
            print(f"  {gene}: {n_samples_mutated} samples ({freq:.1f}%), {n_mutations} mutations")

        if mut_type_col:
            print(f"    Most common type: {top_mut_type}")
        if vaf_col and not np.isnan(mean_vaf):
            print(f"    Mean VAF: {mean_vaf:.3f}, Median VAF: {median_vaf:.3f}")

        driver_results.append({
            'gene': gene,
            'n_samples_mutated': n_samples_mutated,
            'frequency_pct': freq,
            'n_mutations': n_mutations,
            'top_mutation_type': top_mut_type,
            'mean_vaf': mean_vaf,
            'median_vaf': median_vaf
        })
    else:
        print(f"  {gene}: 0 samples (0.0%)")
        driver_results.append({
            'gene': gene,
            'n_samples_mutated': 0,
            'frequency_pct': 0.0,
            'n_mutations': 0,
            'top_mutation_type': 'N/A',
            'mean_vaf': np.nan,
            'median_vaf': np.nan
        })

print()

# ============================================================================
# 6. TOP MUTATED GENES
# ============================================================================

print("6. TOP 30 MOST FREQUENTLY MUTATED GENES")
print("-" * 80)

# Count samples per gene
gene_sample_counts = df_mut.groupby(gene_col)[sample_col].nunique().sort_values(ascending=False)
top_30_genes = gene_sample_counts.head(30)

print(f"{'Rank':<6} {'Gene':<15} {'Samples':<10} {'Frequency':<12} {'Total Mutations':<15}")
print("-" * 65)

top_genes_data = []
for rank, (gene, n_samples) in enumerate(top_30_genes.items(), 1):
    freq = n_samples / unique_samples * 100
    n_muts = len(df_mut[df_mut[gene_col] == gene])
    print(f"{rank:<6} {gene:<15} {n_samples:<10} {freq:>6.1f}%{'':<5} {n_muts:<15}")

    top_genes_data.append({
        'rank': rank,
        'gene': gene,
        'n_samples_mutated': n_samples,
        'frequency_pct': freq,
        'total_mutations': n_muts
    })

print()

# ============================================================================
# 7. MUTATION BURDEN ANALYSIS
# ============================================================================

print("7. MUTATION BURDEN PER SAMPLE")
print("-" * 80)

# Count mutations per sample
mutations_per_sample = df_mut[sample_col].value_counts()

print(f"MUTATION BURDEN STATISTICS:")
print(f"  Mean mutations/sample: {mutations_per_sample.mean():.1f}")
print(f"  Median mutations/sample: {mutations_per_sample.median():.1f}")
print(f"  SD: {mutations_per_sample.std():.1f}")
print(f"  Range: [{mutations_per_sample.min()}, {mutations_per_sample.max()}]")
print()

# Distribution
low_burden = (mutations_per_sample <= 2).sum()
mid_burden = ((mutations_per_sample > 2) & (mutations_per_sample <= 10)).sum()
high_burden = (mutations_per_sample > 10).sum()

print(f"MUTATION BURDEN DISTRIBUTION:")
print(f"  Low (≤2 mutations): {low_burden} samples ({low_burden/unique_samples*100:.1f}%)")
print(f"  Medium (3-10 mutations): {mid_burden} samples ({mid_burden/unique_samples*100:.1f}%)")
print(f"  High (>10 mutations): {high_burden} samples ({high_burden/unique_samples*100:.1f}%)")
print()

# Identify samples with zero mutations (not in the mutation file)
# These would be samples with only germline or QC issues
print(f"NOTE: Samples in mutation file: {unique_samples}")
print(f"      These represent samples with at least one somatic mutation detected")
print()

# ============================================================================
# 8. CO-OCCURRENCE ANALYSIS
# ============================================================================

print("8. MUTATION CO-OCCURRENCE ANALYSIS")
print("-" * 80)

# Create binary mutation matrix for top 20 genes
top_20_genes = gene_sample_counts.head(20).index.tolist()

# Create sample x gene matrix
mutation_matrix = pd.DataFrame(0,
                              index=sample_list,
                              columns=top_20_genes)

for _, row in df_mut.iterrows():
    sample = row[sample_col]
    gene = row[gene_col]
    if gene in top_20_genes:
        mutation_matrix.loc[sample, gene] = 1

# Calculate co-occurrence matrix
cooccurrence_matrix = mutation_matrix.T.dot(mutation_matrix)

# Convert to frequency (percentage of samples)
cooccurrence_freq = cooccurrence_matrix / unique_samples * 100

print(f"Created co-occurrence matrix for top 20 genes")
print(f"Matrix dimensions: {cooccurrence_matrix.shape}")
print()

# Find top co-occurring pairs
cooccurrence_pairs = []
for i, gene1 in enumerate(top_20_genes):
    for j, gene2 in enumerate(top_20_genes):
        if i < j:  # Upper triangle only
            count = cooccurrence_matrix.loc[gene1, gene2]
            freq = cooccurrence_freq.loc[gene1, gene2]
            if count > 0:
                cooccurrence_pairs.append({
                    'gene1': gene1,
                    'gene2': gene2,
                    'n_samples': count,
                    'frequency_pct': freq
                })

# Sort by frequency
cooccurrence_pairs_df = pd.DataFrame(cooccurrence_pairs).sort_values('frequency_pct',
                                                                      ascending=False)

print("TOP 10 CO-OCCURRING MUTATION PAIRS:")
print(f"{'Gene 1':<15} {'Gene 2':<15} {'Samples':<10} {'Frequency':<12}")
print("-" * 55)
for _, row in cooccurrence_pairs_df.head(10).iterrows():
    print(f"{row['gene1']:<15} {row['gene2']:<15} {row['n_samples']:<10.0f} {row['frequency_pct']:>6.1f}%")
print()

# ============================================================================
# 9. SAVE RESULTS
# ============================================================================

print("9. SAVING RESULTS...")
print("-" * 80)

# 1. Mutation summary
mutation_summary = {
    'metric': [
        'Total mutation calls',
        'Unique samples with mutations',
        'Unique genes with mutations',
        'Mean mutations per sample',
        'Median mutations per sample',
        'Sequencing type inferred',
        'File format',
        'SNVs count',
        'Indels count',
        'Mean VAF',
        'Median VAF'
    ],
    'value': [
        f"{len(df_mut):,}",
        f"{unique_samples:,}",
        f"{unique_genes:,}",
        f"{mutations_per_sample.mean():.1f}",
        f"{mutations_per_sample.median():.1f}",
        panel_type,
        file_format,
        f"{n_snv:,}" if mut_type_col else 'N/A',
        f"{n_indel:,}" if mut_type_col else 'N/A',
        f"{vaf_clean.mean():.3f}" if vaf_col and len(vaf_clean) > 0 else 'N/A',
        f"{vaf_clean.median():.3f}" if vaf_col and len(vaf_clean) > 0 else 'N/A'
    ]
}

summary_df = pd.DataFrame(mutation_summary)
output_file = OUTPUT_DIR / "mutation_summary.csv"
summary_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 2. Top mutated genes
top_genes_df = pd.DataFrame(top_genes_data)
output_file = OUTPUT_DIR / "top_mutated_genes.csv"
top_genes_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 3. Driver mutation frequencies
driver_df = pd.DataFrame(driver_results).sort_values('frequency_pct', ascending=False)
output_file = OUTPUT_DIR / "driver_mutation_frequencies.csv"
driver_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 4. Mutation burden per sample
burden_df = pd.DataFrame({
    'sample_id': mutations_per_sample.index,
    'n_mutations': mutations_per_sample.values
}).sort_values('n_mutations', ascending=False)
output_file = OUTPUT_DIR / "mutation_burden_per_sample.csv"
burden_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 5. Co-occurrence matrix
output_file = OUTPUT_DIR / "mutation_cooccurrence_matrix.csv"
cooccurrence_freq.to_csv(output_file)
print(f"✓ Saved: {output_file.name}")

# 6. Co-occurrence pairs
output_file = OUTPUT_DIR / "mutation_cooccurrence_pairs.csv"
cooccurrence_pairs_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

print()

# ============================================================================
# 10. VISUALIZATIONS
# ============================================================================

print("10. GENERATING VISUALIZATIONS...")
print("-" * 80)

# Figure 1: Top 20 mutated genes
plt.figure(figsize=(12, 8))
top_20_data = gene_sample_counts.head(20)
plt.barh(range(len(top_20_data)),
         top_20_data.values / unique_samples * 100,
         color='steelblue')
plt.yticks(range(len(top_20_data)), top_20_data.index)
plt.xlabel('Frequency (%)')
plt.title(f'Top 20 Most Frequently Mutated Genes\n({unique_samples:,} samples)')
plt.gca().invert_yaxis()
plt.grid(axis='x', alpha=0.3)
plt.tight_layout()
output_file = FIGURE_DIR / "top_mutated_genes.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 2: Mutation burden distribution
plt.figure(figsize=(10, 6))
plt.hist(mutations_per_sample.values, bins=50, edgecolor='black', alpha=0.7, color='coral')
plt.xlabel('Number of Mutations per Sample')
plt.ylabel('Number of Samples')
plt.title('Mutation Burden Distribution')
plt.axvline(mutations_per_sample.mean(), color='red', linestyle='--',
            label=f'Mean: {mutations_per_sample.mean():.1f}')
plt.axvline(mutations_per_sample.median(), color='orange', linestyle='--',
            label=f'Median: {mutations_per_sample.median():.1f}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
output_file = FIGURE_DIR / "mutation_burden_distribution.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 3: VAF distribution (if available)
if vaf_col and len(vaf_clean) > 0:
    plt.figure(figsize=(10, 6))
    plt.hist(vaf_clean.values, bins=50, edgecolor='black', alpha=0.7, color='green')
    plt.xlabel('Variant Allele Frequency (VAF)')
    plt.ylabel('Number of Mutations')
    plt.title('VAF Distribution')
    plt.axvline(vaf_clean.mean(), color='red', linestyle='--',
                label=f'Mean: {vaf_clean.mean():.3f}')
    plt.axvline(vaf_clean.median(), color='orange', linestyle='--',
                label=f'Median: {vaf_clean.median():.3f}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    output_file = FIGURE_DIR / "vaf_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file.name}")
    plt.close()

# Figure 4: Co-occurrence heatmap
plt.figure(figsize=(14, 12))
sns.heatmap(cooccurrence_freq, annot=True, fmt='.1f', cmap='YlOrRd',
            square=True, cbar_kws={'label': 'Co-occurrence Frequency (%)'})
plt.title('Mutation Co-occurrence Matrix (Top 20 Genes)\n% of samples with both mutations')
plt.tight_layout()
output_file = FIGURE_DIR / "mutation_cooccurrence_heatmap.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 5: Driver mutations bar plot
plt.figure(figsize=(12, 8))
driver_df_plot = driver_df[driver_df['frequency_pct'] > 0].sort_values('frequency_pct')
plt.barh(range(len(driver_df_plot)), driver_df_plot['frequency_pct'], color='purple', alpha=0.7)
plt.yticks(range(len(driver_df_plot)), driver_df_plot['gene'])
plt.xlabel('Frequency (%)')
plt.title('AML Driver Mutation Frequencies')
plt.grid(axis='x', alpha=0.3)
plt.tight_layout()
output_file = FIGURE_DIR / "driver_mutation_frequencies.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

print()

# ============================================================================
# 11. FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("TASK 2.4 COMPLETE: MUTATION DATA ANALYSIS")
print("=" * 80)
print()
print("SUMMARY:")
print(f"  ✓ Samples with mutations: {unique_samples:,}")
print(f"  ✓ Genes with mutations: {unique_genes:,}")
print(f"  ✓ Total mutation calls: {len(df_mut):,}")
print(f"  ✓ Sequencing type: {panel_type}")
print(f"  ✓ Mean mutations/sample: {mutations_per_sample.mean():.1f}")
print(f"  ✓ Top mutated gene: {top_30_genes.index[0]} ({top_30_genes.iloc[0]/unique_samples*100:.1f}%)")
print()
print("OUTPUTS SAVED:")
print(f"  • mutation_summary.csv")
print(f"  • top_mutated_genes.csv (top 30)")
print(f"  • driver_mutation_frequencies.csv (20 AML drivers)")
print(f"  • mutation_burden_per_sample.csv")
print(f"  • mutation_cooccurrence_matrix.csv (top 20 genes)")
print(f"  • mutation_cooccurrence_pairs.csv")
print(f"  • top_mutated_genes.png")
print(f"  • mutation_burden_distribution.png")
if vaf_col and len(vaf_clean) > 0:
    print(f"  • vaf_distribution.png")
print(f"  • mutation_cooccurrence_heatmap.png")
print(f"  • driver_mutation_frequencies.png")
print()
print("=" * 80)
