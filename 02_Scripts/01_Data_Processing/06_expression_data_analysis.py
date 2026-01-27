"""
Task 2.3: Expression Data Analysis
Efficient analysis of the large gene expression file (269MB)

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
expr_file = DATA_DIR / "beataml_expression.txt"

print("=" * 80)
print("TASK 2.3: EXPRESSION DATA ANALYSIS")
print("=" * 80)
print(f"\nAnalyzing: {expr_file.name}")
print(f"File size: {expr_file.stat().st_size / (1024**2):.2f} MB")
print()

# ============================================================================
# 1. EFFICIENT DATA LOADING
# ============================================================================

print("1. LOADING DATA (this may take a moment for 269MB file)...")
print("-" * 80)

# Load with optimized settings
df_expr = pd.read_csv(expr_file, sep='\t', low_memory=False)

print(f"✓ Data loaded successfully")
print(f"  Shape: {df_expr.shape[0]:,} rows × {df_expr.shape[1]:,} columns")
print(f"  Memory usage: {df_expr.memory_usage(deep=True).sum() / (1024**2):.2f} MB")
print()

# ============================================================================
# 2. DATA STRUCTURE ANALYSIS
# ============================================================================

print("2. DATA STRUCTURE")
print("-" * 80)

# Identify metadata vs expression columns
metadata_cols = []
sample_cols = []

for col in df_expr.columns:
    # Sample columns typically start with 'BA' for BeatAML
    if col.startswith('BA'):
        sample_cols.append(col)
    else:
        metadata_cols.append(col)

print(f"Metadata columns ({len(metadata_cols)}): {', '.join(metadata_cols[:10])}")
if len(metadata_cols) > 10:
    print(f"  ... and {len(metadata_cols) - 10} more")
print()

print(f"Sample columns ({len(sample_cols)}): First 10: {', '.join(sample_cols[:10])}")
print()

# Determine data orientation
num_genes = len(df_expr)
num_samples = len(sample_cols)

print(f"DATA ORIENTATION:")
print(f"  Genes (rows): {num_genes:,}")
print(f"  Samples (columns): {num_samples:,}")
print()

# ============================================================================
# 3. NORMALIZATION & DATA TYPE DETECTION
# ============================================================================

print("3. NORMALIZATION & DATA TYPE")
print("-" * 80)

# Sample some expression values for analysis
expr_values = df_expr[sample_cols].values.flatten()
expr_values_clean = expr_values[~np.isnan(expr_values)]

# Statistical summary
print(f"Expression value statistics:")
print(f"  Min: {np.min(expr_values_clean):.4f}")
print(f"  Max: {np.max(expr_values_clean):.4f}")
print(f"  Mean: {np.mean(expr_values_clean):.4f}")
print(f"  Median: {np.median(expr_values_clean):.4f}")
print(f"  SD: {np.std(expr_values_clean):.4f}")
print()

# Check for characteristics
has_negative = np.any(expr_values_clean < 0)
has_zeros = np.any(expr_values_clean == 0)
max_value = np.max(expr_values_clean)

print(f"Data characteristics:")
print(f"  Negative values: {'Yes' if has_negative else 'No'}")
print(f"  Zero values: {'Yes' if has_zeros else 'No'}")
print(f"  Likely log-transformed: {'Yes (has negatives)' if has_negative else 'Possibly'}")
print()

# Infer normalization method
if max_value < 25 and has_negative:
    inferred_type = "Log2-transformed (FPKM/TPM)"
elif max_value > 100000:
    inferred_type = "Raw counts"
elif max_value < 1000:
    inferred_type = "FPKM or TPM (normalized)"
else:
    inferred_type = "Unknown - needs metadata check"

print(f"Inferred data type: {inferred_type}")
print()

# ============================================================================
# 4. GENE DETECTION ANALYSIS
# ============================================================================

print("4. GENE DETECTION ANALYSIS")
print("-" * 80)

# Define detection threshold (typical: > 0 for log-transformed, > 1 for non-log)
if has_negative:
    detection_threshold = -10  # Very low for log-transformed (essentially any expression)
else:
    detection_threshold = 1

print(f"Using detection threshold: > {detection_threshold}")
print()

# Calculate genes detected per sample
genes_detected_per_sample = []
sample_names = []

for sample in sample_cols:
    n_detected = (df_expr[sample] > detection_threshold).sum()
    genes_detected_per_sample.append(n_detected)
    sample_names.append(sample)

gene_detection_df = pd.DataFrame({
    'sample_id': sample_names,
    'genes_detected': genes_detected_per_sample,
    'detection_rate': [x/num_genes*100 for x in genes_detected_per_sample]
})

print(f"Gene detection per sample:")
print(f"  Mean genes detected: {np.mean(genes_detected_per_sample):,.0f}")
print(f"  Median genes detected: {np.median(genes_detected_per_sample):,.0f}")
print(f"  Min genes detected: {np.min(genes_detected_per_sample):,}")
print(f"  Max genes detected: {np.max(genes_detected_per_sample):,}")
print(f"  SD: {np.std(genes_detected_per_sample):,.0f}")
print()

# Identify potential QC failures (samples with unusually low detection)
q1 = np.percentile(genes_detected_per_sample, 25)
q3 = np.percentile(genes_detected_per_sample, 75)
iqr = q3 - q1
low_detection_threshold = q1 - 1.5 * iqr

low_detection_samples = gene_detection_df[
    gene_detection_df['genes_detected'] < low_detection_threshold
]

print(f"Potential QC failures (genes detected < {low_detection_threshold:.0f}):")
print(f"  Number of samples: {len(low_detection_samples)}")
if len(low_detection_samples) > 0:
    print(f"  Sample IDs: {', '.join(low_detection_samples['sample_id'].head(10).tolist())}")
    if len(low_detection_samples) > 10:
        print(f"    ... and {len(low_detection_samples) - 10} more")
print()

# Calculate detection rate per gene
print("Calculating gene detection across samples...")
genes_detected_across_samples = (df_expr[sample_cols] > detection_threshold).sum(axis=1)
gene_detection_rate = genes_detected_across_samples / num_samples * 100

print(f"Gene detection across samples:")
print(f"  Genes detected in >90% samples: {(gene_detection_rate > 90).sum():,}")
print(f"  Genes detected in >50% samples: {(gene_detection_rate > 50).sum():,}")
print(f"  Genes detected in >10% samples: {(gene_detection_rate > 10).sum():,}")
print(f"  Genes detected in >1% samples: {(gene_detection_rate > 1).sum():,}")
print(f"  Genes never detected: {(gene_detection_rate == 0).sum():,}")
print()

# ============================================================================
# 5. MISSING DATA ANALYSIS
# ============================================================================

print("5. MISSING DATA ANALYSIS")
print("-" * 80)

# Check for NaN values
total_values = df_expr[sample_cols].size
missing_values = df_expr[sample_cols].isna().sum().sum()
missing_pct = missing_values / total_values * 100

print(f"Missing data:")
print(f"  Total expression values: {total_values:,}")
print(f"  Missing (NaN) values: {missing_values:,}")
print(f"  Missing percentage: {missing_pct:.4f}%")
print()

# Check if any samples have high missingness
sample_missing = df_expr[sample_cols].isna().sum()
if sample_missing.max() > 0:
    samples_with_missing = sample_missing[sample_missing > 0]
    print(f"Samples with missing values: {len(samples_with_missing)}")
    print(f"  Max missing per sample: {samples_with_missing.max():,} ({samples_with_missing.max()/num_genes*100:.2f}%)")
else:
    print("No missing values detected - excellent data quality! ✓")
print()

# ============================================================================
# 6. SAMPLE QC METRICS
# ============================================================================

print("6. SAMPLE QC METRICS")
print("-" * 80)

# Calculate QC metrics for each sample
qc_metrics = []

for sample in sample_cols:
    sample_expr = df_expr[sample].dropna()

    qc_metrics.append({
        'sample_id': sample,
        'n_genes': len(sample_expr),
        'n_detected': (sample_expr > detection_threshold).sum(),
        'detection_rate': (sample_expr > detection_threshold).sum() / len(sample_expr) * 100,
        'mean_expr': sample_expr.mean(),
        'median_expr': sample_expr.median(),
        'sd_expr': sample_expr.std(),
        'min_expr': sample_expr.min(),
        'max_expr': sample_expr.max(),
        'missing_values': df_expr[sample].isna().sum()
    })

qc_df = pd.DataFrame(qc_metrics)

print(f"Sample QC summary:")
print(f"  Mean expression per sample:")
print(f"    Mean: {qc_df['mean_expr'].mean():.4f}")
print(f"    Range: [{qc_df['mean_expr'].min():.4f}, {qc_df['mean_expr'].max():.4f}]")
print()

# ============================================================================
# 7. DATA DISTRIBUTION VISUALIZATION
# ============================================================================

print("7. GENERATING VISUALIZATIONS...")
print("-" * 80)

# Figure 1: Overall expression distribution
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
# Sample for histogram (all values would be too many)
sample_for_hist = np.random.choice(expr_values_clean, size=min(100000, len(expr_values_clean)), replace=False)
plt.hist(sample_for_hist, bins=100, edgecolor='black', alpha=0.7)
plt.xlabel('Expression Value')
plt.ylabel('Frequency')
plt.title(f'Expression Distribution\n({num_genes:,} genes × {num_samples} samples)')
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.boxplot(sample_for_hist, vert=True)
plt.ylabel('Expression Value')
plt.title('Expression Value Distribution (Boxplot)')
plt.grid(True, alpha=0.3)

plt.tight_layout()
output_file = FIGURE_DIR / "expression_distribution.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 2: Gene detection per sample
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.hist(genes_detected_per_sample, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
plt.xlabel('Number of Genes Detected')
plt.ylabel('Number of Samples')
plt.title(f'Gene Detection per Sample\n(threshold: > {detection_threshold})')
plt.axvline(np.mean(genes_detected_per_sample), color='red', linestyle='--',
            label=f'Mean: {np.mean(genes_detected_per_sample):,.0f}')
plt.axvline(np.median(genes_detected_per_sample), color='orange', linestyle='--',
            label=f'Median: {np.median(genes_detected_per_sample):,.0f}')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.hist(gene_detection_rate, bins=50, edgecolor='black', alpha=0.7, color='coral')
plt.xlabel('Detection Rate Across Samples (%)')
plt.ylabel('Number of Genes')
plt.title('Gene Detection Rate Across All Samples')
plt.grid(True, alpha=0.3)

plt.tight_layout()
output_file = FIGURE_DIR / "gene_detection_patterns.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 3: Sample QC heatmap
plt.figure(figsize=(10, 6))
qc_metrics_subset = qc_df[['mean_expr', 'detection_rate']].head(50)  # First 50 samples
sns.heatmap(qc_metrics_subset.T, cmap='viridis', cbar_kws={'label': 'Value'})
plt.xlabel('Sample Index')
plt.ylabel('QC Metric')
plt.title('Sample QC Metrics (First 50 Samples)')
plt.tight_layout()
output_file = FIGURE_DIR / "sample_qc_heatmap.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

print()

# ============================================================================
# 8. SAVE RESULTS
# ============================================================================

print("8. SAVING RESULTS...")
print("-" * 80)

# 1. Expression data summary
summary_data = {
    'metric': [
        'Number of genes',
        'Number of samples',
        'Total data points',
        'Metadata columns',
        'Min expression value',
        'Max expression value',
        'Mean expression value',
        'Median expression value',
        'SD expression value',
        'Has negative values',
        'Has zero values',
        'Inferred data type',
        'Missing values (N)',
        'Missing percentage',
        'Mean genes detected per sample',
        'Median genes detected per sample',
        'Samples with low detection',
        'Genes detected in >90% samples',
        'Genes detected in >50% samples',
        'Genes never detected'
    ],
    'value': [
        f"{num_genes:,}",
        f"{num_samples}",
        f"{total_values:,}",
        f"{len(metadata_cols)}",
        f"{np.min(expr_values_clean):.4f}",
        f"{np.max(expr_values_clean):.4f}",
        f"{np.mean(expr_values_clean):.4f}",
        f"{np.median(expr_values_clean):.4f}",
        f"{np.std(expr_values_clean):.4f}",
        'Yes' if has_negative else 'No',
        'Yes' if has_zeros else 'No',
        inferred_type,
        f"{missing_values:,}",
        f"{missing_pct:.4f}%",
        f"{np.mean(genes_detected_per_sample):,.0f}",
        f"{np.median(genes_detected_per_sample):,.0f}",
        f"{len(low_detection_samples)}",
        f"{(gene_detection_rate > 90).sum():,}",
        f"{(gene_detection_rate > 50).sum():,}",
        f"{(gene_detection_rate == 0).sum():,}"
    ]
}

summary_df = pd.DataFrame(summary_data)
output_file = OUTPUT_DIR / "expression_data_summary.csv"
summary_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 2. Gene detection per sample
output_file = OUTPUT_DIR / "gene_detection_per_sample.csv"
gene_detection_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 3. Sample QC metrics
output_file = OUTPUT_DIR / "expression_sample_qc.csv"
qc_df.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")

# 4. Gene-level summary (add gene info if available in metadata)
if 'gene_symbol' in df_expr.columns or 'gene' in df_expr.columns:
    gene_col = 'gene_symbol' if 'gene_symbol' in df_expr.columns else 'gene'
    gene_summary = pd.DataFrame({
        'gene': df_expr[gene_col] if gene_col in df_expr.columns else df_expr.index,
        'detection_rate': gene_detection_rate,
        'n_samples_detected': genes_detected_across_samples,
        'mean_expression': df_expr[sample_cols].mean(axis=1),
        'median_expression': df_expr[sample_cols].median(axis=1),
        'sd_expression': df_expr[sample_cols].std(axis=1)
    })
    output_file = OUTPUT_DIR / "gene_level_summary.csv"
    gene_summary.to_csv(output_file, index=False)
    print(f"✓ Saved: {output_file.name}")

print()

# ============================================================================
# 9. FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("TASK 2.3 COMPLETE: EXPRESSION DATA ANALYSIS")
print("=" * 80)
print()
print("SUMMARY:")
print(f"  ✓ Data structure: {num_genes:,} genes × {num_samples} samples")
print(f"  ✓ Data type: {inferred_type}")
print(f"  ✓ Data quality: {100-missing_pct:.4f}% complete")
print(f"  ✓ Gene detection: Mean {np.mean(genes_detected_per_sample):,.0f} per sample")
print(f"  ✓ QC flags: {len(low_detection_samples)} samples with low detection")
print()
print("OUTPUTS SAVED:")
print(f"  • expression_data_summary.csv")
print(f"  • gene_detection_per_sample.csv")
print(f"  • expression_sample_qc.csv")
if 'gene_symbol' in df_expr.columns or 'gene' in df_expr.columns:
    print(f"  • gene_level_summary.csv")
print(f"  • expression_distribution.png")
print(f"  • gene_detection_patterns.png")
print(f"  • sample_qc_heatmap.png")
print()
print("=" * 80)
