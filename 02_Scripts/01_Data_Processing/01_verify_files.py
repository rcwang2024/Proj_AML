"""
Task 1.1: Verify Downloaded BeatAML Data Files
==============================================
This script verifies that all expected BeatAML data files are present,
checks their integrity, and reports file sizes and basic properties.

Author: Data Analysis Pipeline
Date: 2025-10-02
Project: AML Multi-Omics Integration
"""

import os
import sys
from pathlib import Path
from datetime import datetime
import hashlib

# Set UTF-8 encoding for Windows console
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Expected files with their approximate sizes (in MB)
EXPECTED_FILES = {
    'beataml_expression.txt': 269,
    'beataml_drug_auc.txt': 19,
    'beataml_clinical.xlsx': 0.5,
    'beataml_mutations.txt': 3.5,
    'beataml_raw_inhibitor.txt': 48,
    'beataml_drug_families.xlsx': 0.1
}

def get_file_size_mb(filepath):
    """Get file size in megabytes."""
    try:
        size_bytes = os.path.getsize(filepath)
        size_mb = size_bytes / (1024 * 1024)
        return size_mb
    except OSError as e:
        return None

def check_file_readable(filepath):
    """Check if file can be opened and read."""
    try:
        with open(filepath, 'rb') as f:
            # Try to read first 1KB
            f.read(1024)
        return True
    except Exception as e:
        return False

def calculate_md5(filepath, chunk_size=8192):
    """Calculate MD5 checksum of file (for integrity check)."""
    try:
        md5 = hashlib.md5()
        with open(filepath, 'rb') as f:
            while chunk := f.read(chunk_size):
                md5.update(chunk)
        return md5.hexdigest()
    except Exception as e:
        return None

def verify_files(data_dir, output_log):
    """Main verification function."""

    print("=" * 80)
    print("BeatAML Data Files Verification")
    print("=" * 80)
    print(f"Data Directory: {data_dir}")
    print(f"Verification Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    print()

    results = []
    all_passed = True

    for filename, expected_size_mb in EXPECTED_FILES.items():
        filepath = os.path.join(data_dir, filename)

        print(f"Checking: {filename}")
        print("-" * 60)

        # Check existence
        exists = os.path.exists(filepath)
        print(f"  Exists: {'✓ YES' if exists else '✗ NO'}")

        if not exists:
            results.append({
                'filename': filename,
                'exists': False,
                'size_mb': None,
                'readable': False,
                'status': 'MISSING'
            })
            all_passed = False
            print(f"  Status: MISSING")
            print()
            continue

        # Check size
        actual_size_mb = get_file_size_mb(filepath)
        print(f"  Expected Size: ~{expected_size_mb} MB")
        print(f"  Actual Size: {actual_size_mb:.2f} MB")

        # Size tolerance check (±20%)
        size_ok = (actual_size_mb >= expected_size_mb * 0.8 and
                   actual_size_mb <= expected_size_mb * 1.2)
        print(f"  Size Check: {'✓ PASS' if size_ok else '⚠ WARNING'}")

        # Check readability
        readable = check_file_readable(filepath)
        print(f"  Readable: {'✓ YES' if readable else '✗ NO'}")

        # Calculate checksum (for smaller files only)
        checksum = None
        if actual_size_mb < 50:  # Only for files < 50MB
            checksum = calculate_md5(filepath)
            print(f"  MD5 Checksum: {checksum[:16]}..." if checksum else "  MD5 Checksum: Failed")

        # Overall status
        if readable and size_ok:
            status = "OK"
            print(f"  Status: ✓ OK")
        elif readable:
            status = "WARNING"
            print(f"  Status: ⚠ WARNING (size mismatch)")
            all_passed = False
        else:
            status = "ERROR"
            print(f"  Status: ✗ ERROR (not readable)")
            all_passed = False

        results.append({
            'filename': filename,
            'exists': exists,
            'size_mb': actual_size_mb,
            'readable': readable,
            'checksum': checksum,
            'status': status
        })

        print()

    # Summary
    print("=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)

    ok_count = sum(1 for r in results if r['status'] == 'OK')
    warning_count = sum(1 for r in results if r['status'] == 'WARNING')
    error_count = sum(1 for r in results if r['status'] in ['ERROR', 'MISSING'])

    print(f"Total Files Checked: {len(EXPECTED_FILES)}")
    print(f"  ✓ OK: {ok_count}")
    print(f"  ⚠ WARNING: {warning_count}")
    print(f"  ✗ ERROR/MISSING: {error_count}")
    print()

    if all_passed:
        print("Overall Status: ✓ ALL CHECKS PASSED")
    else:
        print("Overall Status: ⚠ SOME ISSUES DETECTED")

    print("=" * 80)

    # Write log file
    with open(output_log, 'w') as f:
        f.write("BeatAML Data Files Verification Log\n")
        f.write("=" * 80 + "\n")
        f.write(f"Verification Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Data Directory: {data_dir}\n")
        f.write("=" * 80 + "\n\n")

        for result in results:
            f.write(f"File: {result['filename']}\n")
            f.write(f"  Exists: {result['exists']}\n")
            size_str = f"{result['size_mb']:.2f}" if result['size_mb'] is not None else 'N/A'
            f.write(f"  Size (MB): {size_str}\n")
            f.write(f"  Readable: {result['readable']}\n")
            if result.get('checksum'):
                f.write(f"  MD5: {result['checksum']}\n")
            f.write(f"  Status: {result['status']}\n")
            f.write("\n")

        f.write("\nSummary:\n")
        f.write(f"  Total Files: {len(EXPECTED_FILES)}\n")
        f.write(f"  OK: {ok_count}\n")
        f.write(f"  WARNING: {warning_count}\n")
        f.write(f"  ERROR/MISSING: {error_count}\n")
        f.write(f"  Overall: {'PASSED' if all_passed else 'ISSUES DETECTED'}\n")

    print(f"\nLog file saved to: {output_log}")

    return all_passed

if __name__ == "__main__":
    # Set paths
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "01_Data" / "BeatAML_Downloaded_Data"
    output_log = project_root / "06_Documentation" / "Data_Analysis_Log.txt"

    # Create output directory if needed
    output_log.parent.mkdir(parents=True, exist_ok=True)

    # Run verification
    success = verify_files(str(data_dir), str(output_log))

    # Exit with appropriate code
    sys.exit(0 if success else 1)
