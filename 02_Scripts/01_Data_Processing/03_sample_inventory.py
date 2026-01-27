"""
Task 1.3: Create Sample Inventory and Overlap Analysis
=======================================================
This script creates a comprehensive inventory of samples across all datasets
and analyzes overlap for multi-omics integration.

Outputs:
- Sample inventory table
- Overlap analysis (Venn-style)
- Integration-ready sample list

Author: Data Analysis Pipeline
Date: 2025-10-02
Project: AML Multi-Omics Integration
"""

import os
import sys
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

def extract_sample_ids(df, dataset_name):
    """Extract sample IDs from a dataframe based on dataset type."""

    if dataset_name == "expression":
        # Sample IDs are column names (skip first 4 annotation columns)
        sample_ids = [col for col in df.columns[4:] if col.startswith('BA')]
        id_column = "column_names"
        id_format = "BA####R"

    elif dataset_name == "drug_auc":
        # Sample IDs in dbgap_rnaseq_sample column
        sample_ids = df['dbgap_rnaseq_sample'].dropna().unique().tolist()
        id_column = "dbgap_rnaseq_sample"
        id_format = "BA####R"

    elif dataset_name == "clinical":
        # Sample IDs might be in dbgap_rnaseq_sample or similar
        if 'dbgap_rnaseq_sample' in df.columns:
            sample_ids = df['dbgap_rnaseq_sample'].dropna().unique().tolist()
            id_column = "dbgap_rnaseq_sample"
            id_format = "BA####R"
        elif 'dbgap_subject_id' in df.columns:
            sample_ids = df['dbgap_subject_id'].dropna().unique().tolist()
            id_column = "dbgap_subject_id"
            id_format = "numeric"
        else:
            sample_ids = []
            id_column = "unknown"
            id_format = "unknown"

    elif dataset_name == "mutations":
        # Sample IDs in dbgap_sample_id column (DNA samples)
        if 'dbgap_sample_id' in df.columns:
            sample_ids = df['dbgap_sample_id'].dropna().unique().tolist()
            id_column = "dbgap_sample_id"
            id_format = "BA####D"
        else:
            sample_ids = []
            id_column = "unknown"
            id_format = "unknown"

    elif dataset_name == "raw_inhibitor":
        # Sample IDs in dbgap_rnaseq_sample column
        sample_ids = df['dbgap_rnaseq_sample'].dropna().unique().tolist()
        id_column = "dbgap_rnaseq_sample"
        id_format = "BA####R"

    else:
        sample_ids = []
        id_column = "unknown"
        id_format = "unknown"

    return sample_ids, id_column, id_format

def load_datasets(data_dir):
    """Load all datasets and extract sample information."""

    print("Loading datasets...")
    print("=" * 80)

    datasets = {}

    # 1. Expression data
    print("\n1. Loading expression data...")
    try:
        expr = pd.read_csv(data_dir / "beataml_expression.txt", sep='\t', low_memory=False)
        datasets['expression'] = expr
        print(f"   ✓ Loaded: {expr.shape}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        datasets['expression'] = None

    # 2. Drug AUC data
    print("\n2. Loading drug AUC data...")
    try:
        drug_auc = pd.read_csv(data_dir / "beataml_drug_auc.txt", sep='\t', low_memory=False)
        datasets['drug_auc'] = drug_auc
        print(f"   ✓ Loaded: {drug_auc.shape}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        datasets['drug_auc'] = None

    # 3. Clinical data
    print("\n3. Loading clinical data...")
    try:
        clinical = pd.read_excel(data_dir / "beataml_clinical.xlsx")
        datasets['clinical'] = clinical
        print(f"   ✓ Loaded: {clinical.shape}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        datasets['clinical'] = None

    # 4. Mutations data
    print("\n4. Loading mutations data...")
    try:
        mutations = pd.read_csv(data_dir / "beataml_mutations.txt", sep='\t', low_memory=False)
        datasets['mutations'] = mutations
        print(f"   ✓ Loaded: {mutations.shape}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        datasets['mutations'] = None

    # 5. Raw inhibitor data
    print("\n5. Loading raw inhibitor data...")
    try:
        raw_inhib = pd.read_csv(data_dir / "beataml_raw_inhibitor.txt", sep='\t', low_memory=False)
        datasets['raw_inhibitor'] = raw_inhib
        print(f"   ✓ Loaded: {raw_inhib.shape}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        datasets['raw_inhibitor'] = None

    print("\n" + "=" * 80)
    return datasets

def create_sample_inventory(datasets):
    """Create comprehensive sample inventory."""

    print("\nCreating sample inventory...")
    print("=" * 80)

    inventory = []
    sample_sets = {}

    for dataset_name, df in datasets.items():
        if df is None:
            print(f"\n✗ Skipping {dataset_name} (not loaded)")
            continue

        print(f"\nProcessing {dataset_name}...")

        # Extract sample IDs
        sample_ids, id_column, id_format = extract_sample_ids(df, dataset_name)

        n_samples = len(sample_ids)

        print(f"  ID Column: {id_column}")
        print(f"  ID Format: {id_format}")
        print(f"  N Samples: {n_samples:,}")

        # Store sample set for overlap analysis
        sample_sets[dataset_name] = set(sample_ids)

        # Add to inventory
        inventory.append({
            'data_type': dataset_name,
            'n_samples': n_samples,
            'sample_id_column': id_column,
            'id_format': id_format,
            'date_analyzed': datetime.now().strftime('%Y-%m-%d')
        })

    print("\n" + "=" * 80)

    return pd.DataFrame(inventory), sample_sets

def analyze_overlap(sample_sets):
    """Analyze overlap between datasets."""

    print("\nAnalyzing sample overlap...")
    print("=" * 80)

    # Get RNA-seq sample sets (expression, drug_auc, raw_inhibitor, clinical)
    rna_samples = {}
    dna_samples = {}

    for name, samples in sample_sets.items():
        # Check sample format
        if samples and len(samples) > 0:
            sample_list = list(samples)
            if isinstance(sample_list[0], str) and 'R' in sample_list[0]:
                rna_samples[name] = samples
            elif isinstance(sample_list[0], str) and 'D' in sample_list[0]:
                dna_samples[name] = samples
            else:
                # Numeric IDs - try to map
                pass

    print("\nRNA-seq based datasets:")
    for name, samples in rna_samples.items():
        print(f"  {name}: {len(samples)} samples")

    print("\nDNA-seq based datasets:")
    for name, samples in dna_samples.items():
        print(f"  {name}: {len(samples)} samples")

    # Calculate overlaps for RNA samples
    if len(rna_samples) >= 2:
        print("\n" + "-" * 80)
        print("RNA Sample Overlaps:")
        print("-" * 80)

        names = list(rna_samples.keys())
        for i, name1 in enumerate(names):
            for name2 in names[i+1:]:
                overlap = rna_samples[name1] & rna_samples[name2]
                only1 = rna_samples[name1] - rna_samples[name2]
                only2 = rna_samples[name2] - rna_samples[name1]

                print(f"\n{name1} ∩ {name2}:")
                print(f"  Overlap: {len(overlap)} samples")
                print(f"  Only in {name1}: {len(only1)} samples")
                print(f"  Only in {name2}: {len(only2)} samples")

    # Find core multi-omics samples (expression + drug + clinical)
    if 'expression' in rna_samples and 'drug_auc' in rna_samples:
        core_samples = rna_samples['expression'] & rna_samples['drug_auc']

        if 'clinical' in rna_samples:
            core_samples = core_samples & rna_samples['clinical']

        print("\n" + "=" * 80)
        print("CORE MULTI-OMICS SAMPLES")
        print("=" * 80)
        print(f"Samples with Expression + Drug Response: {len(core_samples)}")

        return core_samples

    return set()

def create_integration_table(datasets, core_samples):
    """Create integration-ready sample table."""

    print("\nCreating integration table...")
    print("=" * 80)

    integration_data = []

    for sample in sorted(core_samples):
        row = {'sample_id': sample}

        # Check presence in each dataset
        row['has_expression'] = sample in extract_sample_ids(datasets['expression'], 'expression')[0]
        row['has_drug_response'] = sample in extract_sample_ids(datasets['drug_auc'], 'drug_auc')[0]

        if datasets['clinical'] is not None:
            clinical_samples = extract_sample_ids(datasets['clinical'], 'clinical')[0]
            row['has_clinical'] = sample in clinical_samples
        else:
            row['has_clinical'] = False

        # Check for mutations (need to map RNA -> DNA ID)
        # BA####R -> BA####D
        if datasets['mutations'] is not None:
            dna_id = sample.replace('R', 'D')
            mutation_samples = extract_sample_ids(datasets['mutations'], 'mutations')[0]
            row['has_mutations'] = dna_id in mutation_samples
            if row['has_mutations']:
                row['dna_sample_id'] = dna_id
        else:
            row['has_mutations'] = False

        integration_data.append(row)

    df = pd.DataFrame(integration_data)

    print(f"\nIntegration table created: {len(df)} samples")
    print("\nData availability:")
    for col in ['has_expression', 'has_drug_response', 'has_clinical', 'has_mutations']:
        if col in df.columns:
            n = df[col].sum()
            pct = (n / len(df)) * 100
            print(f"  {col.replace('has_', '')}: {n} ({pct:.1f}%)")

    # Count complete cases
    complete_multi_omics = df[df[['has_expression', 'has_drug_response', 'has_mutations']].all(axis=1)]
    print(f"\nComplete multi-omics (Expression + Drug + Mutations): {len(complete_multi_omics)} samples")

    return df

def main():
    """Main function."""

    # Set paths
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "01_Data" / "BeatAML_Downloaded_Data"
    output_dir = project_root / "03_Results" / "01_Processed_Data"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("BeatAML Sample Inventory and Overlap Analysis")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Data Directory: {data_dir}")
    print(f"Output Directory: {output_dir}")
    print("=" * 80)

    # Load datasets
    datasets = load_datasets(data_dir)

    # Create sample inventory
    inventory_df, sample_sets = create_sample_inventory(datasets)

    # Display inventory
    print("\n" + "=" * 80)
    print("SAMPLE INVENTORY")
    print("=" * 80)
    print(inventory_df.to_string(index=False))

    # Save inventory
    inventory_file = output_dir / "sample_inventory.csv"
    inventory_df.to_csv(inventory_file, index=False)
    print(f"\n✓ Saved inventory to: {inventory_file}")

    # Analyze overlap
    core_samples = analyze_overlap(sample_sets)

    # Create integration table
    if len(core_samples) > 0:
        integration_df = create_integration_table(datasets, core_samples)

        # Save integration table
        integration_file = output_dir / "sample_integration_table.csv"
        integration_df.to_csv(integration_file, index=False)
        print(f"\n✓ Saved integration table to: {integration_file}")

        # Save core sample list
        core_file = output_dir / "core_multi_omics_samples.txt"
        with open(core_file, 'w') as f:
            f.write("# Core Multi-Omics Samples (Expression + Drug Response)\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Total samples: {len(core_samples)}\n\n")
            for sample in sorted(core_samples):
                f.write(f"{sample}\n")
        print(f"✓ Saved core sample list to: {core_file}")

    print("\n" + "=" * 80)
    print("✓ Sample inventory analysis complete!")
    print("=" * 80)

if __name__ == "__main__":
    main()
