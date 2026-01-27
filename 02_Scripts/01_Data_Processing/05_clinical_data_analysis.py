"""
Task 2.2: Clinical Data Structure Analysis
===========================================
Comprehensive analysis of BeatAML clinical data including:
- Demographics (age, sex)
- Disease characteristics
- Survival data
- Mutation data
- Risk stratification
- Treatment information
- Missing data patterns

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

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def load_clinical_data(filepath):
    """Load clinical data."""
    print(f"Loading clinical data from: {filepath}")
    df = pd.read_excel(filepath)
    print(f"✓ Loaded: {df.shape[0]:,} rows × {df.shape[1]:,} columns\n")
    return df

def list_all_variables(df):
    """List all clinical variables."""
    print("=" * 80)
    print("CLINICAL VARIABLES")
    print("=" * 80)

    print(f"\nTotal variables: {df.shape[1]}")
    print(f"\nAll column names:")

    for i, col in enumerate(df.columns, 1):
        dtype = df[col].dtype
        n_unique = df[col].nunique()
        n_missing = df[col].isna().sum()
        pct_missing = (n_missing / len(df)) * 100

        print(f"  {i:2d}. {col:50s} | Type: {str(dtype):10s} | Unique: {n_unique:5d} | Missing: {pct_missing:5.1f}%")

def analyze_demographics(df):
    """Analyze demographics: age and sex."""
    print("\n" + "=" * 80)
    print("DEMOGRAPHICS ANALYSIS")
    print("=" * 80)

    # Age analysis
    age_cols = [col for col in df.columns if 'age' in col.lower()]

    if age_cols:
        print(f"\nAge columns found: {age_cols}")

        for age_col in age_cols:
            if df[age_col].notna().sum() > 0:
                # Convert to numeric, coercing errors
                age_data = pd.to_numeric(df[age_col], errors='coerce').dropna()

                if len(age_data) == 0:
                    continue

                print(f"\n{age_col}:")
                print(f"  N: {len(age_data):,}")
                print(f"  Mean: {age_data.mean():.1f} years")
                print(f"  Median: {age_data.median():.1f} years")
                print(f"  Std Dev: {age_data.std():.1f} years")
                print(f"  Range: [{age_data.min():.0f}, {age_data.max():.0f}] years")
                print(f"  Q1: {age_data.quantile(0.25):.1f} years")
                print(f"  Q3: {age_data.quantile(0.75):.1f} years")

    # Sex/Gender analysis
    sex_cols = [col for col in df.columns if any(x in col.lower() for x in ['sex', 'gender'])]

    if sex_cols:
        print(f"\nSex/Gender columns found: {sex_cols}")

        for sex_col in sex_cols:
            if df[sex_col].notna().sum() > 0:
                sex_counts = df[sex_col].value_counts()

                print(f"\n{sex_col}:")
                for val, count in sex_counts.items():
                    pct = (count / df[sex_col].notna().sum()) * 100
                    print(f"  {val}: {count:,} ({pct:.1f}%)")

def analyze_disease_characteristics(df):
    """Analyze disease characteristics."""
    print("\n" + "=" * 80)
    print("DISEASE CHARACTERISTICS")
    print("=" * 80)

    # De novo vs secondary
    denovo_cols = [col for col in df.columns if 'denovo' in col.lower()]
    if denovo_cols:
        print(f"\nDe novo AML:")
        for col in denovo_cols:
            counts = df[col].value_counts()
            print(f"  {col}:")
            for val, count in counts.items():
                pct = (count / df[col].notna().sum()) * 100
                print(f"    {val}: {count:,} ({pct:.1f}%)")

    # Relapse status
    relapse_cols = [col for col in df.columns if 'relapse' in col.lower()]
    if relapse_cols:
        print(f"\nRelapse status:")
        for col in relapse_cols:
            counts = df[col].value_counts()
            print(f"  {col}:")
            for val, count in counts.items():
                pct = (count / df[col].notna().sum()) * 100
                print(f"    {val}: {count:,} ({pct:.1f}%)")

    # Prior malignancy
    prior_mal_cols = [col for col in df.columns if 'priorMalignancy' in col]
    if prior_mal_cols:
        print(f"\nPrior malignancy:")
        for col in prior_mal_cols[:3]:  # Show first 3
            if df[col].notna().sum() > 0:
                counts = df[col].value_counts()
                print(f"  {col}:")
                for val, count in counts.head(5).items():
                    pct = (count / df[col].notna().sum()) * 100
                    print(f"    {val}: {count:,} ({pct:.1f}%)")

    # Disease type/diagnosis
    dx_cols = [col for col in df.columns if col.startswith('dx') or col.startswith('specific')]
    if dx_cols:
        print(f"\nDiagnosis columns: {len(dx_cols)} found")
        print(f"  Examples: {dx_cols[:3]}")

def analyze_survival(df):
    """Analyze survival data."""
    print("\n" + "=" * 80)
    print("SURVIVAL DATA ANALYSIS")
    print("=" * 80)

    # Vital status
    vital_cols = [col for col in df.columns if 'vital' in col.lower()]
    if vital_cols:
        print(f"\nVital Status:")
        for col in vital_cols:
            counts = df[col].value_counts()
            total = df[col].notna().sum()

            print(f"  {col}:")
            for val, count in counts.items():
                pct = (count / total) * 100
                print(f"    {val}: {count:,} ({pct:.1f}%)")

    # Overall survival
    os_cols = [col for col in df.columns if 'overall' in col.lower() and 'survival' in col.lower()]
    if os_cols:
        print(f"\nOverall Survival:")
        for col in os_cols:
            os_data = df[col].dropna()

            if len(os_data) > 0:
                print(f"  {col}:")
                print(f"    N: {len(os_data):,}")
                print(f"    Mean: {os_data.mean():.1f}")
                print(f"    Median: {os_data.median():.1f}")
                print(f"    Std Dev: {os_data.std():.1f}")
                print(f"    Range: [{os_data.min():.0f}, {os_data.max():.0f}]")

    # Cause of death
    cod_cols = [col for col in df.columns if 'cause' in col.lower() and 'death' in col.lower()]
    if cod_cols:
        print(f"\nCause of Death:")
        for col in cod_cols:
            counts = df[col].value_counts()
            print(f"  {col}: {df[col].notna().sum():,} values")
            print(f"  Top 5 causes:")
            for val, count in counts.head(5).items():
                print(f"    {val}: {count:,}")

def analyze_mutations(df):
    """Analyze mutation data in clinical file."""
    print("\n" + "=" * 80)
    print("MUTATION DATA IN CLINICAL FILE")
    print("=" * 80)

    # Common AML mutations
    mutation_genes = ['FLT3', 'NPM1', 'RUNX1', 'ASXL1', 'TP53', 'DNMT3A', 'IDH1', 'IDH2',
                      'TET2', 'NRAS', 'KRAS', 'CEBPA']

    mutation_cols = [col for col in df.columns if any(gene in col for gene in mutation_genes)]

    if mutation_cols:
        print(f"\nMutation columns found: {len(mutation_cols)}")

        for col in mutation_cols:
            if df[col].notna().sum() > 0:
                counts = df[col].value_counts()
                total = df[col].notna().sum()

                print(f"\n{col}:")
                print(f"  Total assessed: {total:,}")
                for val, count in counts.head(5).items():
                    pct = (count / total) * 100
                    print(f"    {val}: {count:,} ({pct:.1f}%)")
    else:
        print("\n⚠ No mutation columns found in clinical file")

def analyze_risk_stratification(df):
    """Analyze risk stratification."""
    print("\n" + "=" * 80)
    print("RISK STRATIFICATION")
    print("=" * 80)

    # ELN risk
    eln_cols = [col for col in df.columns if 'eln' in col.lower()]
    if eln_cols:
        print(f"\nELN Risk Classification:")
        for col in eln_cols:
            counts = df[col].value_counts()
            total = df[col].notna().sum()

            print(f"  {col}:")
            print(f"    Total classified: {total:,}")
            for val, count in counts.items():
                pct = (count / total) * 100
                print(f"      {val}: {count:,} ({pct:.1f}%)")

    # Cytogenetic risk (from karyotype)
    karyo_cols = [col for col in df.columns if 'karyotype' in col.lower() or 'cytogenetic' in col.lower()]
    if karyo_cols:
        print(f"\nCytogenetic data available:")
        for col in karyo_cols:
            n_available = df[col].notna().sum()
            print(f"  {col}: {n_available:,} samples")

def analyze_treatment(df):
    """Analyze treatment information."""
    print("\n" + "=" * 80)
    print("TREATMENT INFORMATION")
    print("=" * 80)

    # Treatment columns
    tx_cols = [col for col in df.columns if any(x in col.lower() for x in ['treatment', 'therapy', 'regimen'])]

    if tx_cols:
        print(f"\nTreatment columns found: {len(tx_cols)}")

        # Induction therapy response
        response_cols = [col for col in tx_cols if 'response' in col.lower() and 'induction' in col.lower()]
        if response_cols:
            print(f"\nInduction Therapy Response:")
            for col in response_cols:
                counts = df[col].value_counts()
                total = df[col].notna().sum()

                print(f"  {col}: {total:,} patients")
                for val, count in counts.head(5).items():
                    pct = (count / total) * 100
                    print(f"    {val}: {count:,} ({pct:.1f}%)")

        # Treatment types
        type_cols = [col for col in tx_cols if 'type' in col.lower()][:3]
        if type_cols:
            print(f"\nTreatment types available: {len(type_cols)} columns")

def analyze_missing_data(df):
    """Comprehensive missing data analysis."""
    print("\n" + "=" * 80)
    print("MISSING DATA ANALYSIS")
    print("=" * 80)

    # Overall missing
    total_cells = df.shape[0] * df.shape[1]
    missing_cells = df.isna().sum().sum()
    pct_missing = (missing_cells / total_cells) * 100

    print(f"\nOverall: {missing_cells:,} / {total_cells:,} ({pct_missing:.2f}%)")

    # Missing by variable
    missing_by_var = df.isna().sum().sort_values(ascending=False)
    missing_pct = (missing_by_var / len(df)) * 100

    print(f"\nVariables with >50% missing data:")
    high_missing = missing_pct[missing_pct > 50]
    for var, pct in high_missing.items():
        print(f"  {var}: {pct:.1f}%")

    print(f"\nVariables with 0% missing data: {len(missing_pct[missing_pct == 0])}")

    # Create missing data summary
    missing_summary = pd.DataFrame({
        'Variable': df.columns,
        'Missing_Count': df.isna().sum().values,
        'Missing_Percent': (df.isna().sum() / len(df) * 100).values,
        'Data_Type': df.dtypes.values
    }).sort_values('Missing_Percent', ascending=False)

    return missing_summary

def create_visualizations(df, output_dir):
    """Create visualization plots."""
    print("\n" + "=" * 80)
    print("CREATING VISUALIZATIONS")
    print("=" * 80)

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True, parents=True)

    # 1. Age distribution
    age_cols = [col for col in df.columns if 'ageAtDiagnosis' in col]
    if age_cols:
        age_col = age_cols[0]
        age_data = df[age_col].dropna()

        if len(age_data) > 0:
            plt.figure(figsize=(10, 6))
            plt.hist(age_data, bins=30, edgecolor='black', alpha=0.7)
            plt.xlabel('Age at Diagnosis (years)')
            plt.ylabel('Frequency')
            plt.title('Age Distribution at AML Diagnosis')
            plt.axvline(age_data.mean(), color='red', linestyle='--', label=f'Mean: {age_data.mean():.1f}')
            plt.axvline(age_data.median(), color='green', linestyle='--', label=f'Median: {age_data.median():.1f}')
            plt.legend()
            plt.tight_layout()

            outfile = fig_dir / "age_distribution.png"
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"✓ Saved: {outfile}")

    # 2. Missing data heatmap (top 30 variables with most missing)
    missing_data = df.isna().sum().sort_values(ascending=False).head(30)
    if len(missing_data) > 0:
        plt.figure(figsize=(10, 12))
        missing_pct = (missing_data / len(df)) * 100
        plt.barh(range(len(missing_pct)), missing_pct.values)
        plt.yticks(range(len(missing_pct)), missing_pct.index, fontsize=8)
        plt.xlabel('Missing Data (%)')
        plt.title('Top 30 Variables by Missing Data Percentage')
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()

        outfile = fig_dir / "clinical_missing_data.png"
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved: {outfile}")

def save_outputs(df, missing_summary, output_dir):
    """Save analysis outputs."""
    print("\n" + "=" * 80)
    print("SAVING OUTPUTS")
    print("=" * 80)

    # 1. Clinical data summary
    summary_data = []
    for col in df.columns:
        summary_data.append({
            'Variable': col,
            'Data_Type': str(df[col].dtype),
            'N_Total': len(df),
            'N_Missing': df[col].isna().sum(),
            'Pct_Missing': (df[col].isna().sum() / len(df)) * 100,
            'N_Unique': df[col].nunique()
        })

    summary_df = pd.DataFrame(summary_data)

    outfile = output_dir / "clinical_data_summary.csv"
    summary_df.to_csv(outfile, index=False)
    print(f"✓ Saved: {outfile}")

    # 2. Variable completeness
    outfile = output_dir.parent / "02_QC_Reports" / "clinical_completeness.csv"
    outfile.parent.mkdir(exist_ok=True, parents=True)
    missing_summary.to_csv(outfile, index=False)
    print(f"✓ Saved: {outfile}")

    # 3. Demographics table
    demo_data = []

    # Age
    age_cols = [col for col in df.columns if 'ageAtDiagnosis' in col]
    if age_cols:
        age_data = df[age_cols[0]].dropna()
        demo_data.append({
            'Variable': 'Age at Diagnosis',
            'N': len(age_data),
            'Mean': age_data.mean(),
            'Median': age_data.median(),
            'SD': age_data.std(),
            'Min': age_data.min(),
            'Max': age_data.max()
        })

    # Sex
    sex_cols = [col for col in df.columns if 'consensus_sex' in col]
    if sex_cols:
        sex_counts = df[sex_cols[0]].value_counts()
        for sex, count in sex_counts.items():
            demo_data.append({
                'Variable': f'Sex: {sex}',
                'N': count,
                'Mean': np.nan,
                'Median': np.nan,
                'SD': np.nan,
                'Min': np.nan,
                'Max': np.nan
            })

    if demo_data:
        demo_df = pd.DataFrame(demo_data)
        outfile = output_dir / "demographics_table.csv"
        demo_df.to_csv(outfile, index=False)
        print(f"✓ Saved: {outfile}")

def main():
    """Main analysis function."""

    # Set paths
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "01_Data" / "BeatAML_Downloaded_Data"
    output_dir = project_root / "03_Results" / "01_Processed_Data"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("CLINICAL DATA STRUCTURE ANALYSIS")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Data Directory: {data_dir}")
    print(f"Output Directory: {output_dir}")
    print("=" * 80)
    print()

    # Load data
    clinical_file = data_dir / "beataml_clinical.xlsx"
    df = load_clinical_data(clinical_file)

    # List all variables
    list_all_variables(df)

    # Analyze demographics
    analyze_demographics(df)

    # Analyze disease characteristics
    analyze_disease_characteristics(df)

    # Analyze survival
    analyze_survival(df)

    # Analyze mutations
    analyze_mutations(df)

    # Analyze risk stratification
    analyze_risk_stratification(df)

    # Analyze treatment
    analyze_treatment(df)

    # Analyze missing data
    missing_summary = analyze_missing_data(df)

    # Create visualizations
    create_visualizations(df, output_dir)

    # Save outputs
    save_outputs(df, missing_summary, output_dir)

    print("\n" + "=" * 80)
    print("✓ CLINICAL DATA ANALYSIS COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
