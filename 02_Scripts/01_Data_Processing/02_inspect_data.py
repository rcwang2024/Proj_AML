"""
Task 1.2: Inspect BeatAML Data Files in Detail
==============================================
This script provides comprehensive inspection of each downloaded BeatAML file:
- File format and structure
- Dimensions (rows x columns)
- Column names and data types
- First 10 rows preview
- Basic statistics
- Missing data analysis

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

def format_bytes(size):
    """Convert bytes to human readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} TB"

def get_basic_stats(df, col):
    """Get basic statistics for a column."""
    stats = {}
    try:
        if pd.api.types.is_numeric_dtype(df[col]):
            stats['min'] = df[col].min()
            stats['max'] = df[col].max()
            stats['mean'] = df[col].mean()
            stats['median'] = df[col].median()
        else:
            stats['unique'] = df[col].nunique()
            top_val = df[col].value_counts().head(1)
            if len(top_val) > 0:
                stats['top_value'] = top_val.index[0]
                stats['top_count'] = top_val.values[0]
    except:
        pass
    return stats

def inspect_file(filepath, filename, output_file):
    """Inspect a single data file."""

    print("\n" + "=" * 80)
    print(f"INSPECTING: {filename}")
    print("=" * 80)

    # Write to output
    output_file.write("\n" + "=" * 80 + "\n")
    output_file.write(f"File: {filename}\n")
    output_file.write("=" * 80 + "\n\n")

    # File size
    file_size = os.path.getsize(filepath)
    print(f"File Size: {format_bytes(file_size)}")
    output_file.write(f"File Size: {format_bytes(file_size)}\n")

    # Determine file type and read
    try:
        if filename.endswith('.xlsx'):
            print("File Format: Excel (.xlsx)")
            output_file.write("File Format: Excel (.xlsx)\n")
            df = pd.read_excel(filepath)
            file_format = "Excel"
        elif filename.endswith('.txt'):
            # Try to detect delimiter
            with open(filepath, 'r') as f:
                first_line = f.readline()
                if '\t' in first_line:
                    delimiter = '\t'
                    format_name = "Tab-delimited text"
                elif ',' in first_line:
                    delimiter = ','
                    format_name = "Comma-separated text"
                else:
                    delimiter = None
                    format_name = "Text file"

            print(f"File Format: {format_name}")
            output_file.write(f"File Format: {format_name}\n")

            if delimiter:
                df = pd.read_csv(filepath, sep=delimiter, low_memory=False)
            else:
                df = pd.read_csv(filepath, low_memory=False)
            file_format = "Text"
        else:
            print("File Format: Unknown")
            output_file.write("File Format: Unknown\n")
            return

    except Exception as e:
        print(f"✗ ERROR reading file: {e}")
        output_file.write(f"ERROR reading file: {e}\n")
        return

    # Dimensions
    n_rows, n_cols = df.shape
    print(f"Dimensions: {n_rows:,} rows × {n_cols:,} columns")
    output_file.write(f"Dimensions: {n_rows:,} rows × {n_cols:,} columns\n\n")

    # Column information
    print("\nColumn Information:")
    print("-" * 80)
    output_file.write("Column Information:\n")
    output_file.write("-" * 80 + "\n")

    col_info = []
    for col in df.columns:
        dtype = str(df[col].dtype)
        n_missing = df[col].isna().sum()
        pct_missing = (n_missing / len(df)) * 100

        col_info.append({
            'Column': col,
            'Type': dtype,
            'Missing': n_missing,
            'Missing_%': f"{pct_missing:.1f}%"
        })

        print(f"  {col}")
        print(f"    Type: {dtype}")
        print(f"    Missing: {n_missing:,} ({pct_missing:.1f}%)")

        # Basic stats
        stats = get_basic_stats(df, col)
        if 'min' in stats:
            print(f"    Range: [{stats['min']:.2f}, {stats['max']:.2f}]")
            print(f"    Mean: {stats['mean']:.2f}, Median: {stats['median']:.2f}")
        elif 'unique' in stats:
            print(f"    Unique values: {stats['unique']:,}")
            if 'top_value' in stats:
                print(f"    Most common: '{stats['top_value']}' ({stats['top_count']:,} times)")

    # Save column info as dataframe
    col_df = pd.DataFrame(col_info)
    output_file.write(col_df.to_string(index=False) + "\n\n")

    # First 10 rows preview
    print("\nFirst 10 Rows Preview:")
    print("-" * 80)
    preview = df.head(10)
    print(preview.to_string())

    output_file.write("First 10 Rows Preview:\n")
    output_file.write("-" * 80 + "\n")
    output_file.write(preview.to_string() + "\n\n")

    # Overall missing data
    total_cells = n_rows * n_cols
    total_missing = df.isna().sum().sum()
    pct_missing = (total_missing / total_cells) * 100

    print(f"\nOverall Missing Data: {total_missing:,} / {total_cells:,} ({pct_missing:.2f}%)")
    output_file.write(f"Overall Missing Data: {total_missing:,} / {total_cells:,} ({pct_missing:.2f}%)\n")

    # Data quality issues
    print("\nData Quality Check:")
    output_file.write("\nData Quality Check:\n")

    issues = []

    # Check for duplicate rows
    n_duplicates = df.duplicated().sum()
    if n_duplicates > 0:
        issues.append(f"⚠ {n_duplicates:,} duplicate rows found")
    else:
        print("  ✓ No duplicate rows")
        output_file.write("  No duplicate rows\n")

    # Check for columns with all missing data
    all_missing_cols = [col for col in df.columns if df[col].isna().all()]
    if all_missing_cols:
        issues.append(f"⚠ {len(all_missing_cols)} columns with all missing data: {', '.join(all_missing_cols)}")
    else:
        print("  ✓ No columns with all missing data")
        output_file.write("  No columns with all missing data\n")

    # Check for columns with >50% missing
    high_missing_cols = [col for col in df.columns if (df[col].isna().sum() / len(df)) > 0.5]
    if high_missing_cols:
        issues.append(f"⚠ {len(high_missing_cols)} columns with >50% missing data")
        print(f"  ⚠ {len(high_missing_cols)} columns with >50% missing data:")
        output_file.write(f"  Columns with >50% missing data:\n")
        for col in high_missing_cols:
            pct = (df[col].isna().sum() / len(df)) * 100
            print(f"      - {col} ({pct:.1f}% missing)")
            output_file.write(f"    - {col} ({pct:.1f}% missing)\n")
    else:
        print("  ✓ No columns with >50% missing data")
        output_file.write("  No columns with >50% missing data\n")

    if issues:
        for issue in issues:
            print(f"  {issue}")
            output_file.write(f"  {issue}\n")

    print("\n" + "=" * 80)
    output_file.write("\n")

def main():
    """Main inspection function."""

    # Set paths
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "01_Data" / "BeatAML_Downloaded_Data"
    output_dir = project_root / "03_Results" / "02_QC_Reports"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file_path = output_dir / "data_inspection_summary.txt"

    # Files to inspect
    files_to_inspect = [
        'beataml_expression.txt',
        'beataml_drug_auc.txt',
        'beataml_clinical.xlsx',
        'beataml_mutations.txt',
        'beataml_raw_inhibitor.txt',
        'beataml_drug_families.xlsx'
    ]

    print("=" * 80)
    print("BeatAML Data Files - Detailed Inspection")
    print("=" * 80)
    print(f"Inspection Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Output File: {output_file_path}")
    print("=" * 80)

    with open(output_file_path, 'w', encoding='utf-8') as output_file:
        output_file.write("BeatAML Data Files - Detailed Inspection\n")
        output_file.write("=" * 80 + "\n")
        output_file.write(f"Inspection Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        output_file.write("=" * 80 + "\n")

        for filename in files_to_inspect:
            filepath = data_dir / filename
            if filepath.exists():
                inspect_file(filepath, filename, output_file)
            else:
                print(f"\n✗ File not found: {filename}")
                output_file.write(f"\nFile not found: {filename}\n")

    print("\n" + "=" * 80)
    print(f"✓ Inspection complete! Report saved to:")
    print(f"  {output_file_path}")
    print("=" * 80)

if __name__ == "__main__":
    main()
