#!/usr/bin/env python3
"""
Extract UUID to barcode mapping for TARGET-AML samples
"""

import requests
import json
import pandas as pd
import time

# GDC API endpoint
API_URL = "https://api.gdc.cancer.gov/files"

# Query for TARGET-AML RNA-seq files
filters = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "cases.project.project_id",
                "value": ["TARGET-AML"]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "files.data_type",
                "value": ["Gene Expression Quantification"]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "files.analysis.workflow_type",
                "value": ["STAR - Counts"]
            }
        }
    ]
}

# Fields to retrieve
fields = [
    "file_id",
    "file_name",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.sample_id",
    "cases.samples.submitter_id",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
    "cases.samples.portions.analytes.aliquots.submitter_id"
]

# Build query
params = {
    "filters": json.dumps(filters),
    "fields": ",".join(fields),
    "format": "JSON",
    "size": 10000  # Get all results
}

print("Querying GDC API for TARGET-AML file metadata...")
response = requests.get(API_URL, params=params)

if response.status_code != 200:
    print(f"Error: API returned status code {response.status_code}")
    print(response.text)
    exit(1)

data = response.json()
hits = data['data']['hits']

print(f"Found {len(hits)} files")

# Extract mappings
mappings = []

for hit in hits:
    file_id = hit.get('file_id', '')
    file_name = hit.get('file_name', '')

    # Get case information
    cases = hit.get('cases', [])
    if cases:
        case = cases[0]  # Usually one case per file
        case_id = case.get('case_id', '')
        case_submitter_id = case.get('submitter_id', '')

        # Get sample information
        samples = case.get('samples', [])
        if samples:
            sample = samples[0]  # Usually one sample
            sample_id = sample.get('sample_id', '')
            sample_submitter_id = sample.get('submitter_id', '')

            # Get aliquot information (most specific level)
            portions = sample.get('portions', [])
            aliquot_id = ''
            aliquot_submitter_id = ''

            if portions:
                for portion in portions:
                    analytes = portion.get('analytes', [])
                    for analyte in analytes:
                        aliquots = analyte.get('aliquots', [])
                        if aliquots:
                            aliquot = aliquots[0]
                            aliquot_id = aliquot.get('aliquot_id', '')
                            aliquot_submitter_id = aliquot.get('submitter_id', '')
                            break
                    if aliquot_id:
                        break

            mappings.append({
                'file_id': file_id,
                'file_name': file_name,
                'case_id': case_id,
                'case_submitter_id': case_submitter_id,
                'sample_id': sample_id,
                'sample_submitter_id': sample_submitter_id,
                'aliquot_id': aliquot_id,
                'aliquot_submitter_id': aliquot_submitter_id
            })

# Create DataFrame
df = pd.DataFrame(mappings)

print(f"\nCreated mapping for {len(df)} files")
print("\nFirst few mappings:")
print(df.head())

# Save mapping
output_file = "03_Results/18_TARGET_Validation/uuid_barcode_mapping.csv"
df.to_csv(output_file, index=False)
print(f"\nSaved mapping to: {output_file}")

# Check overlap with clinical data
clinical = pd.read_csv("03_Results/18_TARGET_Validation/target_aml_clinical.csv")
clinical_barcodes = set(clinical['sample_barcode'].str.strip())

print(f"\nClinical barcodes: {len(clinical_barcodes)}")

# Try matching at different levels
case_matches = df[df['case_submitter_id'].isin(clinical_barcodes)]
sample_matches = df[df['sample_submitter_id'].isin(clinical_barcodes)]
aliquot_matches = df[df['aliquot_submitter_id'].isin(clinical_barcodes)]

print(f"Matches by case_submitter_id: {len(case_matches)}")
print(f"Matches by sample_submitter_id: {len(sample_matches)}")
print(f"Matches by aliquot_submitter_id: {len(aliquot_matches)}")

# Check if clinical barcodes are case IDs or sample IDs
print("\nSample clinical barcodes:")
print(list(clinical_barcodes)[:5])

print("\nSample GDC case_submitter_ids:")
print(df['case_submitter_id'].head().tolist())

print("\nSample GDC sample_submitter_ids:")
print(df['sample_submitter_id'].head().tolist())

print("\n✓ Mapping extraction complete")
