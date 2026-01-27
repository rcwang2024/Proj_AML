"""
Phase 4: Data Quality Assessment
Task 4.2: Expression Data Quality Metrics

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
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
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/02_QC_Reports")
FIGURE_DIR = Path("D:/Projects/Project_AML/04_Figures/01_QC_Figures")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TASK 4.2: EXPRESSION DATA QUALITY METRICS")
print("=" * 80)
print()

# ============================================================================
# STEP 1: LOAD EXPRESSION DATA
# ============================================================================

print("STEP 1: LOADING EXPRESSION DATA")
print("-" * 80)

expr_file = DATA_DIR / "beataml_expression.txt"
df_expr = pd.read_csv(expr_file, sep='\t', low_memory=False)

sample_cols = [col for col in df_expr.columns if col.startswith('BA')]
expr_matrix = df_expr[sample_cols].values.T  # samples × genes

print(f"✓ Expression data loaded")
print(f"  Samples: {len(sample_cols)}")
print(f"  Genes: {expr_matrix.shape[1]:,}")
print()

# ============================================================================
# STEP 2: PCA-BASED OUTLIER DETECTION
# ============================================================================

print("STEP 2: PCA-BASED OUTLIER DETECTION")
print("-" * 80)

# Standardize and perform PCA
scaler = StandardScaler()
expr_scaled = scaler.fit_transform(expr_matrix)

pca = PCA(n_components=10)
pca_coords = pca.fit_transform(expr_scaled)

# Detect outliers using Mahalanobis-like distance in PC space
from scipy.spatial.distance import mahalanobis

# Calculate center and covariance
pc_center = pca_coords.mean(axis=0)
pc_cov = np.cov(pca_coords.T)

# Calculate distances
try:
    pc_cov_inv = np.linalg.inv(pc_cov)
    distances = [mahalanobis(sample, pc_center, pc_cov_inv) for sample in pca_coords]
except:
    # Fallback to Euclidean if singular
    distances = np.sqrt(((pca_coords - pc_center) ** 2).sum(axis=1))

# Identify outliers (> 3 SD from mean)
dist_mean = np.mean(distances)
dist_std = np.std(distances)
threshold = dist_mean + 3 * dist_std

pca_outliers = np.array(sample_cols)[np.array(distances) > threshold].tolist()

print(f"PCA outlier detection:")
print(f"  Mean distance: {dist_mean:.2f}")
print(f"  SD: {dist_std:.2f}")
print(f"  Threshold (mean + 3*SD): {threshold:.2f}")
print(f"  Outliers detected: {len(pca_outliers)}")
if pca_outliers:
    print(f"  Outlier samples: {', '.join(pca_outliers[:10])}")
    if len(pca_outliers) > 10:
        print(f"    ... and {len(pca_outliers)-10} more")
print()

# ============================================================================
# STEP 3: SAMPLE CORRELATION ANALYSIS
# ============================================================================

print("STEP 3: SAMPLE CORRELATION ANALYSIS")
print("-" * 80)

print("Calculating sample-sample correlation matrix...")
# Calculate pairwise correlation (this may take a moment)
corr_matrix = np.corrcoef(expr_matrix)

# Calculate median correlation for each sample
median_corr = np.median(corr_matrix, axis=1)
mean_corr = np.mean(corr_matrix, axis=1)

print(f"✓ Correlation matrix calculated")
print(f"  Mean of median correlations: {np.mean(median_corr):.3f}")
print(f"  Range: [{np.min(median_corr):.3f}, {np.max(median_corr):.3f}]")

# Identify low correlation samples (median < 0.8)
low_corr_mask = median_corr < 0.8
low_corr_samples = np.array(sample_cols)[low_corr_mask].tolist()

print(f"\nLow correlation samples (median < 0.8): {len(low_corr_samples)}")
if low_corr_samples:
    print(f"  Samples: {', '.join(low_corr_samples[:10])}")
    if len(low_corr_samples) > 10:
        print(f"    ... and {len(low_corr_samples)-10} more")
print()

# ============================================================================
# STEP 4: GENE DETECTION OUTLIERS
# ============================================================================

print("STEP 4: GENE DETECTION OUTLIERS")
print("-" * 80)

# Load gene detection data from previous analysis
gene_det_file = Path("D:/Projects/Project_AML/03_Results/01_Processed_Data/gene_detection_per_sample.csv")

if gene_det_file.exists():
    df_gene_det = pd.read_csv(gene_det_file)

    # Identify outliers (< mean - 2*SD)
    mean_det = df_gene_det['genes_detected'].mean()
    std_det = df_gene_det['genes_detected'].std()
    threshold_low = mean_det - 2 * std_det

    gene_det_outliers = df_gene_det[df_gene_det['genes_detected'] < threshold_low]['sample_id'].tolist()

    print(f"Gene detection outliers:")
    print(f"  Mean genes detected: {mean_det:.0f}")
    print(f"  SD: {std_det:.0f}")
    print(f"  Threshold (mean - 2*SD): {threshold_low:.0f}")
    print(f"  Outliers: {len(gene_det_outliers)}")
else:
    print("⚠ Gene detection file not found - skipping this check")
    gene_det_outliers = []

print()

# ============================================================================
# STEP 5: HIERARCHICAL CLUSTERING OUTLIERS
# ============================================================================

print("STEP 5: HIERARCHICAL CLUSTERING FOR OUTLIERS")
print("-" * 80)

print("Performing hierarchical clustering...")
# Use correlation distance
correlation_dist = 1 - corr_matrix
np.fill_diagonal(correlation_dist, 0)

# Hierarchical clustering
linkage_matrix = linkage(pdist(expr_matrix, metric='correlation'), method='average')

# Identify outliers based on linkage height
# Samples with very high linkage distances may be outliers
from scipy.cluster.hierarchy import fcluster
max_d = 0.5  # threshold for outlier detection
clusters = fcluster(linkage_matrix, max_d, criterion='distance')

# Identify singletons or small clusters as potential outliers
cluster_counts = pd.Series(clusters).value_counts()
small_clusters = cluster_counts[cluster_counts <= 3].index.tolist()

hclust_outliers = []
for cluster_id in small_clusters:
    cluster_samples = np.array(sample_cols)[clusters == cluster_id].tolist()
    hclust_outliers.extend(cluster_samples)

print(f"✓ Hierarchical clustering complete")
print(f"  Clusters identified: {len(np.unique(clusters))}")
print(f"  Small clusters (≤3 samples): {len(small_clusters)}")
print(f"  Potential outliers: {len(hclust_outliers)}")
print()

# ============================================================================
# STEP 6: COMBINE OUTLIER RESULTS
# ============================================================================

print("STEP 6: COMBINING OUTLIER DETECTION RESULTS")
print("-" * 80)

# Combine all outlier lists
all_outliers = set(pca_outliers + low_corr_samples + gene_det_outliers + hclust_outliers)

# Create outlier summary
outlier_summary = []
for sample in sample_cols:
    is_pca_outlier = sample in pca_outliers
    is_low_corr = sample in low_corr_samples
    is_gene_det_outlier = sample in gene_det_outliers
    is_hclust_outlier = sample in hclust_outliers

    n_flags = sum([is_pca_outlier, is_low_corr, is_gene_det_outlier, is_hclust_outlier])

    if n_flags > 0:
        idx = sample_cols.index(sample)
        outlier_summary.append({
            'sample_id': sample,
            'pca_outlier': is_pca_outlier,
            'low_correlation': is_low_corr,
            'median_correlation': median_corr[idx],
            'gene_detection_outlier': is_gene_det_outlier,
            'clustering_outlier': is_hclust_outlier,
            'n_outlier_flags': n_flags,
            'pca_distance': distances[idx] if idx < len(distances) else np.nan
        })

df_outliers = pd.DataFrame(outlier_summary).sort_values('n_outlier_flags', ascending=False)

print(f"Total unique outliers: {len(all_outliers)}")
print(f"  PCA outliers: {len(pca_outliers)}")
print(f"  Low correlation: {len(low_corr_samples)}")
print(f"  Gene detection: {len(gene_det_outliers)}")
print(f"  Clustering: {len(hclust_outliers)}")
print()

print(f"Samples with multiple flags:")
for n_flags in range(4, 0, -1):
    n_samples = (df_outliers['n_outlier_flags'] == n_flags).sum()
    if n_samples > 0:
        print(f"  {n_flags} flags: {n_samples} samples")

print()

# Save outlier list
output_file = OUTPUT_DIR / "expression_outliers.csv"
df_outliers.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")
print()

# ============================================================================
# STEP 7: CREATE VISUALIZATIONS
# ============================================================================

print("STEP 7: CREATING VISUALIZATIONS")
print("-" * 80)

# Figure 1: PCA outliers plot
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# PC1 vs PC2 with outliers highlighted
ax = axes[0]
outlier_mask = np.array([s in pca_outliers for s in sample_cols])
ax.scatter(pca_coords[~outlier_mask, 0], pca_coords[~outlier_mask, 1],
          alpha=0.5, s=20, c='steelblue', label='Normal')
ax.scatter(pca_coords[outlier_mask, 0], pca_coords[outlier_mask, 1],
          alpha=0.8, s=50, c='red', marker='X', label='PCA Outliers')
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('PCA with Outliers Highlighted', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Distance plot
ax = axes[1]
ax.hist(distances, bins=50, alpha=0.7, edgecolor='black')
ax.axvline(threshold, color='red', linestyle='--', linewidth=2,
          label=f'Threshold: {threshold:.1f}')
ax.set_xlabel('Mahalanobis Distance')
ax.set_ylabel('Number of Samples')
ax.set_title('Distribution of PCA Distances', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = FIGURE_DIR / "pca_outliers.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 2: Sample correlation heatmap (subset)
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Heatmap of median correlations
ax = axes[0]
sorted_idx = np.argsort(median_corr)
im = ax.imshow([median_corr[sorted_idx]], aspect='auto', cmap='RdYlGn', vmin=0.5, vmax=1.0)
ax.set_yticks([])
ax.set_xlabel('Samples (sorted by median correlation)')
ax.set_title('Median Sample Correlation', fontweight='bold')
plt.colorbar(im, ax=ax, label='Median Correlation')

# Histogram of median correlations
ax = axes[1]
ax.hist(median_corr, bins=50, alpha=0.7, edgecolor='black', color='steelblue')
ax.axvline(0.8, color='red', linestyle='--', linewidth=2, label='QC Threshold (0.8)')
ax.axvline(np.median(median_corr), color='orange', linestyle='--', linewidth=2,
          label=f'Median: {np.median(median_corr):.3f}')
ax.set_xlabel('Median Correlation')
ax.set_ylabel('Number of Samples')
ax.set_title('Distribution of Median Sample Correlations', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = FIGURE_DIR / "sample_correlation_metrics.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 3: Full correlation heatmap (downsampled for visualization)
print("Creating correlation heatmap (this may take a moment)...")
# Downsample for visualization if too many samples
if len(sample_cols) > 100:
    step = len(sample_cols) // 100
    idx_subset = np.arange(0, len(sample_cols), step)
    corr_subset = corr_matrix[idx_subset][:, idx_subset]
    sample_subset = [sample_cols[i] for i in idx_subset]
else:
    corr_subset = corr_matrix
    sample_subset = sample_cols

plt.figure(figsize=(12, 10))
sns.heatmap(corr_subset, cmap='RdYlGn', center=0.9, vmin=0.5, vmax=1.0,
           xticklabels=False, yticklabels=False, cbar_kws={'label': 'Correlation'})
plt.title(f'Sample-Sample Correlation Heatmap\n(showing {len(sample_subset)} of {len(sample_cols)} samples)',
         fontweight='bold', fontsize=14)
plt.tight_layout()
output_file = FIGURE_DIR / "sample_correlation_heatmap.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("TASK 4.2 COMPLETE: EXPRESSION DATA QUALITY METRICS")
print("=" * 80)
print()
print("SUMMARY:")
print(f"  ✓ Analyzed {len(sample_cols)} samples")
print(f"  ✓ Multiple outlier detection methods applied")
print()
print("OUTLIER DETECTION RESULTS:")
print(f"  • PCA-based outliers: {len(pca_outliers)}")
print(f"  • Low correlation (< 0.8): {len(low_corr_samples)}")
print(f"  • Gene detection outliers: {len(gene_det_outliers)}")
print(f"  • Clustering outliers: {len(hclust_outliers)}")
print(f"  • Total unique outliers: {len(all_outliers)}")
print()
print("QUALITY METRICS:")
print(f"  • Mean median correlation: {np.mean(median_corr):.3f}")
print(f"  • Samples < 0.8 correlation: {(median_corr < 0.8).sum()}")
print()
print("OUTPUTS SAVED:")
print("  • expression_outliers.csv")
print("  • pca_outliers.png")
print("  • sample_correlation_metrics.png")
print("  • sample_correlation_heatmap.png")
print()
print("=" * 80)
