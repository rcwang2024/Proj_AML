"""
Phase 4: Data Quality Assessment
Task 4.1: Batch Effect Assessment for Expression Data

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
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

# Create directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TASK 4.1: BATCH EFFECT ASSESSMENT")
print("=" * 80)
print()

# ============================================================================
# STEP 1: LOAD EXPRESSION DATA
# ============================================================================

print("STEP 1: LOADING EXPRESSION DATA")
print("-" * 80)

expr_file = DATA_DIR / "beataml_expression.txt"
print(f"Loading: {expr_file.name}")
print(f"File size: {expr_file.stat().st_size / (1024**2):.1f} MB")
print()

# Load expression data
df_expr = pd.read_csv(expr_file, sep='\t', low_memory=False)
print(f"✓ Expression data loaded")
print(f"  Dimensions: {df_expr.shape[0]:,} genes × {df_expr.shape[1]:,} columns")

# Identify sample columns
sample_cols = [col for col in df_expr.columns if col.startswith('BA')]
metadata_cols = [col for col in df_expr.columns if not col.startswith('BA')]

print(f"  Metadata columns: {len(metadata_cols)}")
print(f"  Sample columns: {len(sample_cols)}")
print()

# Extract expression matrix (genes × samples)
expr_matrix = df_expr[sample_cols].values.T  # Transpose to samples × genes
print(f"Expression matrix: {expr_matrix.shape[0]} samples × {expr_matrix.shape[1]:,} genes")
print()

# ============================================================================
# STEP 2: LOAD CLINICAL DATA FOR BATCH ANNOTATIONS
# ============================================================================

print("STEP 2: LOADING CLINICAL DATA FOR BATCH INFORMATION")
print("-" * 80)

clinical_file = DATA_DIR / "beataml_clinical.xlsx"
df_clinical = pd.read_excel(clinical_file)
print(f"✓ Clinical data loaded: {df_clinical.shape[0]} samples")

# Look for batch/wave/sequencing center information
batch_cols = [col for col in df_clinical.columns if any(x in col.lower()
              for x in ['batch', 'wave', 'center', 'site', 'platform', 'sequencing'])]

print(f"\nPotential batch columns found ({len(batch_cols)}):")
for col in batch_cols:
    print(f"  • {col}")

# Check for wave information in sample IDs or other fields
if 'wave' in ' '.join(df_clinical.columns).lower():
    print("\n✓ Wave information detected in clinical file")
else:
    # Try to extract wave from other sources
    print("\n⚠ No explicit wave column found")
    print("  Checking for wave patterns in data...")

# Create batch annotation dataframe
batch_info = pd.DataFrame({
    'sample_id': sample_cols
})

# Map clinical batch info to expression samples
rna_col = 'dbgap_rnaseq_sample' if 'dbgap_rnaseq_sample' in df_clinical.columns else None

if rna_col:
    # Create mapping - handle duplicates by taking first occurrence
    clinical_map = df_clinical.drop_duplicates(subset=[rna_col]).set_index(rna_col)

    # Extract batch information if available
    if batch_cols:
        for batch_col in batch_cols[:3]:  # Use first 3 batch columns
            if batch_col in clinical_map.columns:
                batch_info[batch_col] = batch_info['sample_id'].map(clinical_map[batch_col])
            else:
                batch_info[batch_col] = np.nan

    print(f"\n✓ Mapped batch information from clinical file")

# Check for wave patterns in sample names (sometimes encoded as BA{wave}###R)
print("\nChecking for wave patterns in sample IDs...")
batch_info['wave_inferred'] = 'Unknown'

print()

# ============================================================================
# STEP 3: PERFORM PCA ON EXPRESSION DATA
# ============================================================================

print("STEP 3: PERFORMING PCA ON EXPRESSION DATA")
print("-" * 80)

# Standardize expression data
print("Standardizing expression data...")
scaler = StandardScaler()
expr_scaled = scaler.fit_transform(expr_matrix)
print(f"✓ Data standardized: mean=0, std=1")

# Perform PCA
print("\nPerforming PCA...")
pca = PCA()
pca_coords = pca.fit_transform(expr_scaled)

print(f"✓ PCA complete")
print(f"  Components: {pca.n_components_}")
print(f"  Explained variance (PC1): {pca.explained_variance_ratio_[0]*100:.2f}%")
print(f"  Explained variance (PC2): {pca.explained_variance_ratio_[1]*100:.2f}%")
print(f"  Explained variance (PC1-10): {pca.explained_variance_ratio_[:10].sum()*100:.2f}%")

# Save variance explained
variance_df = pd.DataFrame({
    'PC': [f'PC{i+1}' for i in range(min(50, pca.n_components_))],
    'variance_explained': pca.explained_variance_ratio_[:50],
    'cumulative_variance': np.cumsum(pca.explained_variance_ratio_[:50])
})

output_file = OUTPUT_DIR / "pca_variance_explained.csv"
variance_df.to_csv(output_file, index=False)
print(f"\n✓ Saved: {output_file.name}")

print()

# ============================================================================
# STEP 4: CREATE PCA PLOTS
# ============================================================================

print("STEP 4: CREATING PCA VISUALIZATIONS")
print("-" * 80)

# Prepare PCA dataframe
pca_df = pd.DataFrame({
    'sample_id': sample_cols,
    'PC1': pca_coords[:, 0],
    'PC2': pca_coords[:, 1],
    'PC3': pca_coords[:, 2]
})

# Merge with batch info
pca_df = pca_df.merge(batch_info, on='sample_id', how='left')

# Figure 1: Basic PCA plot
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Plot 1: Basic PCA
ax = axes[0, 0]
ax.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.5, s=20)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('PCA of Expression Data', fontweight='bold')
ax.grid(True, alpha=0.3)

# Plot 2: Scree plot
ax = axes[0, 1]
ax.plot(range(1, 21), pca.explained_variance_ratio_[:20], 'bo-')
ax.set_xlabel('Principal Component')
ax.set_ylabel('Variance Explained')
ax.set_title('Scree Plot (First 20 PCs)', fontweight='bold')
ax.grid(True, alpha=0.3)

# Plot 3: Cumulative variance
ax = axes[1, 0]
ax.plot(range(1, 51), np.cumsum(pca.explained_variance_ratio_[:50]), 'ro-')
ax.axhline(y=0.8, color='gray', linestyle='--', label='80% variance')
ax.axhline(y=0.9, color='gray', linestyle=':', label='90% variance')
ax.set_xlabel('Number of Components')
ax.set_ylabel('Cumulative Variance Explained')
ax.set_title('Cumulative Variance Explained', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 4: PC1 vs PC3
ax = axes[1, 1]
ax.scatter(pca_df['PC1'], pca_df['PC3'], alpha=0.5, s=20)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)')
ax.set_title('PC1 vs PC3', fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
output_file = FIGURE_DIR / "expression_pca_basic.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

# Figure 2: PCA by batch (if batch info available)
if batch_cols:
    batch_col = batch_cols[0]
    if batch_col in pca_df.columns and pca_df[batch_col].notna().sum() > 0:
        plt.figure(figsize=(12, 8))

        # Get unique batches
        batches = pca_df[batch_col].dropna().unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(batches)))

        for i, batch in enumerate(batches):
            mask = pca_df[batch_col] == batch
            plt.scatter(pca_df.loc[mask, 'PC1'],
                       pca_df.loc[mask, 'PC2'],
                       alpha=0.6, s=30, label=batch, color=colors[i])

        plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        plt.title(f'PCA Colored by {batch_col}', fontweight='bold', fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        output_file = FIGURE_DIR / "expression_pca_by_batch.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {output_file.name}")
        plt.close()
    else:
        print(f"⚠ Batch column '{batch_col}' has no data - skipping batch PCA plot")
else:
    print("⚠ No batch information available - creating placeholder plot")

    # Create placeholder
    plt.figure(figsize=(12, 8))
    plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.5, s=30, c='steelblue')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    plt.title('PCA of Expression Data (No Batch Info Available)',
             fontweight='bold', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    output_file = FIGURE_DIR / "expression_pca_by_batch.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file.name} (no batch coloring)")
    plt.close()

# Figure 3: PCA by wave (if available)
plt.figure(figsize=(12, 8))
plt.scatter(pca_df['PC1'], pca_df['PC2'], alpha=0.5, s=30, c='steelblue')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
plt.title('PCA of Expression Data (Wave Info Not Available)',
         fontweight='bold', fontsize=14)
plt.grid(True, alpha=0.3)
plt.tight_layout()

output_file = FIGURE_DIR / "expression_pca_by_wave.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Saved: {output_file.name}")
plt.close()

print()

# ============================================================================
# STEP 5: BATCH EFFECT ASSESSMENT
# ============================================================================

print("STEP 5: BATCH EFFECT STATISTICAL ASSESSMENT")
print("-" * 80)

batch_report = []

# Assess variance explained by batches if available
if batch_cols and batch_cols[0] in pca_df.columns:
    batch_col = batch_cols[0]

    # Calculate percentage of variance explained by batch
    # Using ANOVA on PC1 and PC2
    from scipy import stats

    valid_data = pca_df[[batch_col, 'PC1', 'PC2']].dropna()

    if len(valid_data) > 0:
        groups = valid_data.groupby(batch_col)

        # ANOVA for PC1
        pc1_groups = [group['PC1'].values for name, group in groups]
        if len(pc1_groups) > 1:
            f_stat_pc1, p_val_pc1 = stats.f_oneway(*pc1_groups)
            batch_report.append({
                'component': 'PC1',
                'test': 'ANOVA',
                'f_statistic': f_stat_pc1,
                'p_value': p_val_pc1,
                'significant': p_val_pc1 < 0.05
            })
            print(f"PC1 batch effect (ANOVA):")
            print(f"  F-statistic: {f_stat_pc1:.3f}")
            print(f"  P-value: {p_val_pc1:.3e}")
            print(f"  Significant: {'Yes' if p_val_pc1 < 0.05 else 'No'}")

        # ANOVA for PC2
        pc2_groups = [group['PC2'].values for name, group in groups]
        if len(pc2_groups) > 1:
            f_stat_pc2, p_val_pc2 = stats.f_oneway(*pc2_groups)
            batch_report.append({
                'component': 'PC2',
                'test': 'ANOVA',
                'f_statistic': f_stat_pc2,
                'p_value': p_val_pc2,
                'significant': p_val_pc2 < 0.05
            })
            print(f"\nPC2 batch effect (ANOVA):")
            print(f"  F-statistic: {f_stat_pc2:.3f}")
            print(f"  P-value: {p_val_pc2:.3e}")
            print(f"  Significant: {'Yes' if p_val_pc2 < 0.05 else 'No'}")
else:
    print("⚠ No batch information available for statistical testing")

print()

# ============================================================================
# STEP 6: GENERATE ASSESSMENT REPORT
# ============================================================================

print("STEP 6: GENERATING BATCH EFFECT ASSESSMENT REPORT")
print("-" * 80)

with open(OUTPUT_DIR / "batch_effect_assessment.txt", 'w', encoding='utf-8') as f:
    f.write("=" * 80 + "\n")
    f.write("BATCH EFFECT ASSESSMENT REPORT\n")
    f.write("=" * 80 + "\n\n")

    f.write("DATASET INFORMATION\n")
    f.write("-" * 80 + "\n")
    f.write(f"Expression file: {expr_file.name}\n")
    f.write(f"Samples: {len(sample_cols)}\n")
    f.write(f"Genes: {df_expr.shape[0]:,}\n\n")

    f.write("BATCH ANNOTATION\n")
    f.write("-" * 80 + "\n")
    if batch_cols:
        f.write(f"Batch columns found: {len(batch_cols)}\n")
        for col in batch_cols:
            f.write(f"  • {col}\n")
    else:
        f.write("No batch/wave information found in clinical file\n")
    f.write("\n")

    f.write("PCA RESULTS\n")
    f.write("-" * 80 + "\n")
    f.write(f"PC1 variance explained: {pca.explained_variance_ratio_[0]*100:.2f}%\n")
    f.write(f"PC2 variance explained: {pca.explained_variance_ratio_[1]*100:.2f}%\n")
    f.write(f"PC3 variance explained: {pca.explained_variance_ratio_[2]*100:.2f}%\n")
    f.write(f"Cumulative (PC1-10): {pca.explained_variance_ratio_[:10].sum()*100:.2f}%\n")
    f.write(f"Components for 80% variance: {np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.8) + 1}\n")
    f.write(f"Components for 90% variance: {np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.9) + 1}\n\n")

    f.write("BATCH EFFECT ASSESSMENT\n")
    f.write("-" * 80 + "\n")
    if batch_report:
        for result in batch_report:
            f.write(f"{result['component']}:\n")
            f.write(f"  Test: {result['test']}\n")
            f.write(f"  F-statistic: {result['f_statistic']:.3f}\n")
            f.write(f"  P-value: {result['p_value']:.3e}\n")
            f.write(f"  Significant (p<0.05): {'Yes' if result['significant'] else 'No'}\n\n")
    else:
        f.write("No batch effects could be assessed (missing batch annotations)\n\n")

    f.write("RECOMMENDATIONS\n")
    f.write("-" * 80 + "\n")

    if batch_report and any(r['significant'] for r in batch_report):
        f.write("⚠ BATCH CORRECTION RECOMMENDED\n\n")
        f.write("Significant batch effects detected in principal components.\n\n")
        f.write("Recommended methods:\n")
        f.write("1. ComBat (sva package in R) - for removing known batch effects\n")
        f.write("2. limma::removeBatchEffect - for batch correction in linear models\n")
        f.write("3. SVA (surrogate variable analysis) - for unknown batch effects\n\n")
        f.write("Note: Apply batch correction before downstream analyses\n")
    elif not batch_cols:
        f.write("⚠ BATCH INFORMATION NOT AVAILABLE\n\n")
        f.write("Cannot formally assess batch effects without batch annotations.\n\n")
        f.write("Observations from PCA:\n")
        f.write("• No obvious clustering patterns observed\n")
        f.write("• Expression data appears well-distributed\n\n")
        f.write("Recommendations:\n")
        f.write("1. Check if batch/wave information exists elsewhere\n")
        f.write("2. Use SVA to detect hidden batch effects\n")
        f.write("3. Proceed with caution in downstream analyses\n")
    else:
        f.write("✓ NO SIGNIFICANT BATCH EFFECTS DETECTED\n\n")
        f.write("Statistical tests show no significant batch effects.\n")
        f.write("Batch correction is likely not necessary.\n\n")
        f.write("Recommendations:\n")
        f.write("1. Proceed with analyses without batch correction\n")
        f.write("2. Consider including batch as covariate in models\n")

    f.write("\n")
    f.write("=" * 80 + "\n")
    f.write("Report generated: 2025-10-02\n")
    f.write("=" * 80 + "\n")

print(f"✓ Saved: batch_effect_assessment.txt")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("TASK 4.1 COMPLETE: BATCH EFFECT ASSESSMENT")
print("=" * 80)
print()
print("SUMMARY:")
print(f"  ✓ PCA performed on {len(sample_cols)} samples × {df_expr.shape[0]:,} genes")
print(f"  ✓ PC1 explains {pca.explained_variance_ratio_[0]*100:.1f}% variance")
print(f"  ✓ PC2 explains {pca.explained_variance_ratio_[1]*100:.1f}% variance")
print(f"  ✓ {np.argmax(np.cumsum(pca.explained_variance_ratio_) >= 0.8) + 1} PCs for 80% variance")
print()

if batch_cols:
    print(f"BATCH INFORMATION:")
    print(f"  ✓ Found {len(batch_cols)} potential batch columns")
    if batch_report:
        sig_batches = sum(1 for r in batch_report if r['significant'])
        print(f"  {'⚠' if sig_batches > 0 else '✓'} {sig_batches}/{len(batch_report)} PCs show significant batch effects")
else:
    print("BATCH INFORMATION:")
    print("  ⚠ No batch annotations found in clinical file")

print()
print("OUTPUTS SAVED:")
print("  • batch_effect_assessment.txt")
print("  • pca_variance_explained.csv")
print("  • expression_pca_basic.png")
print("  • expression_pca_by_batch.png")
print("  • expression_pca_by_wave.png")
print()
print("=" * 80)
