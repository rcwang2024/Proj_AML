"""
Task 2.1: Drug Response Data Structure Analysis
================================================
Comprehensive analysis of BeatAML drug response data including:
- Data structure and dimensions
- Missing data patterns
- Drug and sample statistics
- AUC distribution analysis
- Quality metrics
- Comparison with pilot dataset

Author: Data Analysis Pipeline
Date: 2025-10-02
Project: AML Multi-Omics Integration - Phase 2
"""

import os
import sys
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def load_drug_data(filepath):
    """Load drug response data."""
    print(f"Loading drug response data from: {filepath}")
    df = pd.read_csv(filepath, sep='\t', low_memory=False)
    print(f"✓ Loaded: {df.shape[0]:,} rows × {df.shape[1]:,} columns\n")
    return df

def analyze_data_structure(df):
    """Analyze basic data structure."""
    print("=" * 80)
    print("DATA STRUCTURE ANALYSIS")
    print("=" * 80)

    print(f"\nDimensions: {df.shape[0]:,} rows × {df.shape[1]:,} columns")

    # Display columns
    print("\nColumn names:")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")

    # Data types
    print("\nData types:")
    print(df.dtypes.value_counts())

    # Check format (wide vs long)
    print("\nData Format Assessment:")

    # Check if it's long format (typical columns: sample_id, drug, value)
    likely_long = any(col.lower() in ['inhibitor', 'drug', 'compound'] for col in df.columns)

    if likely_long:
        print("  Format: LONG format (each row is a sample-drug combination)")
    else:
        print("  Format: WIDE format (each column is a drug, each row is a sample)")

    return likely_long

def analyze_samples(df, is_long_format):
    """Analyze unique samples."""
    print("\n" + "=" * 80)
    print("SAMPLE ANALYSIS")
    print("=" * 80)

    # Identify sample ID column
    sample_cols = [col for col in df.columns if 'sample' in col.lower()]

    if sample_cols:
        sample_col = sample_cols[0]
        if 'rnaseq' in sample_col.lower():
            sample_col = [c for c in sample_cols if 'rnaseq' in c.lower()][0]

        unique_samples = df[sample_col].dropna().unique()
        n_unique = len(unique_samples)

        print(f"\nSample ID column: {sample_col}")
        print(f"Unique samples: {n_unique:,}")
        print(f"\nSample ID format examples:")
        for i, sample in enumerate(unique_samples[:5], 1):
            print(f"  {i}. {sample}")

        return sample_col, unique_samples
    else:
        print("\n⚠ Warning: Could not identify sample ID column")
        return None, None

def analyze_drugs(df, is_long_format):
    """Analyze drugs/inhibitors."""
    print("\n" + "=" * 80)
    print("DRUG/INHIBITOR ANALYSIS")
    print("=" * 80)

    # Identify drug column
    drug_cols = [col for col in df.columns if any(x in col.lower() for x in ['inhibitor', 'drug', 'compound'])]

    if drug_cols:
        drug_col = drug_cols[0]
        unique_drugs = df[drug_col].dropna().unique()
        n_unique = len(unique_drugs)

        print(f"\nDrug column: {drug_col}")
        print(f"Unique drugs/inhibitors: {n_unique:,}")

        print(f"\nTop 10 drugs (by frequency):")
        drug_counts = df[drug_col].value_counts().head(10)
        for i, (drug, count) in enumerate(drug_counts.items(), 1):
            print(f"  {i}. {drug}: {count:,} measurements")

        return drug_col, unique_drugs
    else:
        # If wide format, drugs are column names (exclude ID columns)
        id_cols = [col for col in df.columns if any(x in col.lower() for x in ['sample', 'subject', 'id', 'age'])]
        drug_columns = [col for col in df.columns if col not in id_cols]

        print(f"\nFormat: Wide (drugs as columns)")
        print(f"Number of drug columns: {len(drug_columns)}")
        print(f"\nFirst 10 drug columns:")
        for i, drug in enumerate(drug_columns[:10], 1):
            print(f"  {i}. {drug}")

        return None, drug_columns

def analyze_metrics(df):
    """Analyze available metrics."""
    print("\n" + "=" * 80)
    print("AVAILABLE METRICS")
    print("=" * 80)

    # Look for metric columns
    metric_keywords = ['auc', 'ic50', 'ec50', 'viability', 'response', 'sensitivity']
    metric_cols = [col for col in df.columns if any(kw in col.lower() for kw in metric_keywords)]

    if metric_cols:
        print(f"\nFound {len(metric_cols)} metric columns:")
        for col in metric_cols:
            n_values = df[col].notna().sum()
            print(f"  • {col}: {n_values:,} values")
    else:
        print("\n⚠ No obvious metric columns found")

    return metric_cols

def analyze_missing_data(df, drug_col, sample_col):
    """Analyze missing data patterns."""
    print("\n" + "=" * 80)
    print("MISSING DATA ANALYSIS")
    print("=" * 80)

    # Overall missing
    total_cells = df.shape[0] * df.shape[1]
    missing_cells = df.isna().sum().sum()
    pct_missing = (missing_cells / total_cells) * 100

    print(f"\nOverall missing data: {missing_cells:,} / {total_cells:,} ({pct_missing:.2f}%)")

    # Missing by column
    print("\nMissing data by column:")
    missing_by_col = df.isna().sum().sort_values(ascending=False)
    missing_by_col_pct = (missing_by_col / len(df)) * 100

    for col, n_missing in missing_by_col.head(10).items():
        pct = missing_by_col_pct[col]
        print(f"  {col}: {n_missing:,} ({pct:.1f}%)")

    # Create missing data summary
    missing_summary = pd.DataFrame({
        'Column': df.columns,
        'Missing_Count': df.isna().sum().values,
        'Missing_Percent': (df.isna().sum() / len(df) * 100).values
    }).sort_values('Missing_Percent', ascending=False)

    return missing_summary

def analyze_auc_distribution(df):
    """Analyze AUC value distribution."""
    print("\n" + "=" * 80)
    print("AUC DISTRIBUTION ANALYSIS")
    print("=" * 80)

    # Find AUC column
    auc_cols = [col for col in df.columns if 'auc' in col.lower()]

    if not auc_cols:
        print("\n⚠ No AUC column found")
        return None

    auc_col = auc_cols[0]
    auc_values = df[auc_col].dropna()

    print(f"\nAUC column: {auc_col}")
    print(f"Number of AUC values: {len(auc_values):,}")

    # Summary statistics
    print("\nSummary Statistics:")
    print(f"  Mean: {auc_values.mean():.4f}")
    print(f"  Median: {auc_values.median():.4f}")
    print(f"  Std Dev: {auc_values.std():.4f}")
    print(f"  Min: {auc_values.min():.4f}")
    print(f"  Max: {auc_values.max():.4f}")
    print(f"  Q1 (25%): {auc_values.quantile(0.25):.4f}")
    print(f"  Q3 (75%): {auc_values.quantile(0.75):.4f}")

    # Detect outliers (using IQR method)
    Q1 = auc_values.quantile(0.25)
    Q3 = auc_values.quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    outliers = auc_values[(auc_values < lower_bound) | (auc_values > upper_bound)]
    print(f"\nOutliers (IQR method):")
    print(f"  Lower bound: {lower_bound:.4f}")
    print(f"  Upper bound: {upper_bound:.4f}")
    print(f"  Number of outliers: {len(outliers):,} ({len(outliers)/len(auc_values)*100:.2f}%)")

    return auc_col, auc_values

def analyze_sample_drug_counts(df, drug_col, sample_col):
    """Analyze how many drugs tested per sample."""
    print("\n" + "=" * 80)
    print("SAMPLE DRUG TESTING COMPLETENESS")
    print("=" * 80)

    if drug_col and sample_col:
        # Count drugs per sample
        drugs_per_sample = df.groupby(sample_col)[drug_col].nunique()

        print(f"\nDrugs tested per sample:")
        print(f"  Mean: {drugs_per_sample.mean():.1f}")
        print(f"  Median: {drugs_per_sample.median():.0f}")
        print(f"  Min: {drugs_per_sample.min()}")
        print(f"  Max: {drugs_per_sample.max()}")

        print(f"\nDistribution of drug testing:")
        print(drugs_per_sample.describe())

        return drugs_per_sample
    else:
        print("\n⚠ Cannot calculate without drug and sample columns")
        return None

def create_visualizations(df, auc_col, auc_values, output_dir):
    """Create visualization plots."""
    print("\n" + "=" * 80)
    print("CREATING VISUALIZATIONS")
    print("=" * 80)

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True, parents=True)

    # 1. AUC histogram
    if auc_values is not None:
        plt.figure(figsize=(10, 6))
        plt.hist(auc_values, bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('AUC Value')
        plt.ylabel('Frequency')
        plt.title('Distribution of Drug Response AUC Values')
        plt.axvline(auc_values.mean(), color='red', linestyle='--', label=f'Mean: {auc_values.mean():.2f}')
        plt.axvline(auc_values.median(), color='green', linestyle='--', label=f'Median: {auc_values.median():.2f}')
        plt.legend()
        plt.tight_layout()

        outfile = fig_dir / "auc_distribution.png"
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved: {outfile}")

        # 2. AUC boxplot
        plt.figure(figsize=(8, 6))
        plt.boxplot(auc_values, vert=True)
        plt.ylabel('AUC Value')
        plt.title('AUC Value Distribution (Boxplot)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        outfile = fig_dir / "auc_boxplot.png"
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved: {outfile}")

def save_outputs(df, missing_summary, drugs_per_sample, drug_col, sample_col, output_dir):
    """Save analysis outputs."""
    print("\n" + "=" * 80)
    print("SAVING OUTPUTS")
    print("=" * 80)

    # 1. Drug summary
    if drug_col:
        drug_summary = df.groupby(drug_col).size().reset_index(name='n_measurements')
        drug_summary = drug_summary.sort_values('n_measurements', ascending=False)

        outfile = output_dir / "drug_response_summary.csv"
        drug_summary.to_csv(outfile, index=False)
        print(f"✓ Saved: {outfile}")

    # 2. Sample drug counts
    if drugs_per_sample is not None:
        sample_drug_df = drugs_per_sample.reset_index()
        sample_drug_df.columns = ['sample_id', 'n_drugs_tested']

        outfile = output_dir / "samples_drug_counts.csv"
        sample_drug_df.to_csv(outfile, index=False)
        print(f"✓ Saved: {outfile}")

    # 3. Missing data summary
    outfile = output_dir / "drug_missing_data_summary.csv"
    missing_summary.to_csv(outfile, index=False)
    print(f"✓ Saved: {outfile}")

def main():
    """Main analysis function."""

    # Set paths
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "01_Data" / "BeatAML_Downloaded_Data"
    output_dir = project_root / "03_Results" / "01_Processed_Data"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("DRUG RESPONSE DATA STRUCTURE ANALYSIS")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Data Directory: {data_dir}")
    print(f"Output Directory: {output_dir}")
    print("=" * 80)
    print()

    # Load data
    drug_file = data_dir / "beataml_drug_auc.txt"
    df = load_drug_data(drug_file)

    # Analyze structure
    is_long_format = analyze_data_structure(df)

    # Analyze samples
    sample_col, unique_samples = analyze_samples(df, is_long_format)

    # Analyze drugs
    drug_col, unique_drugs = analyze_drugs(df, is_long_format)

    # Analyze metrics
    metric_cols = analyze_metrics(df)

    # Analyze missing data
    missing_summary = analyze_missing_data(df, drug_col, sample_col)

    # Analyze AUC distribution
    auc_col, auc_values = analyze_auc_distribution(df)

    # Analyze sample drug counts
    drugs_per_sample = analyze_sample_drug_counts(df, drug_col, sample_col)

    # Create visualizations
    if auc_values is not None:
        create_visualizations(df, auc_col, auc_values, output_dir)

    # Save outputs
    save_outputs(df, missing_summary, drugs_per_sample, drug_col, sample_col, output_dir)

    print("\n" + "=" * 80)
    print("✓ DRUG RESPONSE ANALYSIS COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
