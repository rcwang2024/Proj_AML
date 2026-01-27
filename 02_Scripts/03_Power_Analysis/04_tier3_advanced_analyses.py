"""
Phase 5: Statistical Power and Feasibility Assessment
Task 5.2 Extended: TIER 3 Advanced Integrative Analyses

Author: AML Multi-Omics Project
Date: 2025-10-02
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Paths
DATA_DIR = Path("D:/Projects/Project_AML/01_Data/BeatAML_Downloaded_Data")
PROCESSED_DIR = Path("D:/Projects/Project_AML/03_Results/01_Processed_Data")
OUTPUT_DIR = Path("D:/Projects/Project_AML/03_Results/03_Power_Analysis")
REPORT_DIR = Path("D:/Projects/Project_AML/05_Reports")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TIER 3: ADVANCED INTEGRATIVE ANALYSES (EXPLORATORY)")
print("=" * 80)
print()

# ============================================================================
# LOAD DATA AND CALCULATE SAMPLE SIZES
# ============================================================================

print("LOADING DATA FOR TIER 3 ANALYSES")
print("-" * 80)

# Master mapping
master_map = pd.read_csv(PROCESSED_DIR / "master_sample_id_mapping.csv")

# Calculate cohort sizes
n_gold_standard = ((master_map['has_expression']) &
                   (master_map['has_mutations']) &
                   (master_map['has_drug_response']) &
                   (master_map['has_clinical'])).sum()

n_expr_drug = ((master_map['has_expression']) & (master_map['has_drug_response'])).sum()
n_mut_drug = ((master_map['has_mutations']) & (master_map['has_drug_response'])).sum()
n_expr_mut_drug = ((master_map['has_expression']) &
                   (master_map['has_mutations']) &
                   (master_map['has_drug_response'])).sum()

n_expression = master_map['has_expression'].sum()
n_mutations = master_map['has_mutations'].sum()
n_drug_response = master_map['has_drug_response'].sum()

print(f"✓ Sample sizes calculated")
print(f"  Gold standard (all 4 types): n = {n_gold_standard}")
print(f"  Expression + Drug: n = {n_expr_drug}")
print(f"  Mutations + Drug: n = {n_mut_drug}")
print(f"  Expression + Mutations + Drug: n = {n_expr_mut_drug}")
print()

# ============================================================================
# ANALYSIS 8: MULTI-OMICS NETWORK ANALYSIS
# ============================================================================

print("=" * 80)
print("ANALYSIS 3.4: MULTI-OMICS NETWORK ANALYSIS")
print("=" * 80)
print()

analysis_8 = {
    'tier': 3,
    'analysis_id': '3.4',
    'analysis_name': 'Multi-Omics Network Analysis',
    'goal': 'Integrate all 4 data types into unified biological network',
    'methods': 'Network analysis, pathway integration, graph-based methods',
    'required_n': 50,
    'available_n': n_gold_standard,
    'data_requirements': 'Complete quad-omics (Expression + Mutations + Drug + Clinical)',
    'timeline_weeks': '4-6 weeks',
    'scripts_location': '02_Scripts/06_Integration/',
}

print(f"Goal: {analysis_8['goal']}")
print(f"Method: {analysis_8['methods']}")
print()

print(f"Sample Size:")
print(f"  Required: n ≥ {analysis_8['required_n']} with complete quad-omics")
print(f"  Available: n = {analysis_8['available_n']}")
print()

# Network analysis approaches
print("NETWORK ANALYSIS APPROACHES:")
print()

network_approaches = []

# 1. Gene-gene co-expression networks
print("  1. Gene-Gene Co-expression Networks")
print(f"     Samples needed: ≥50 (have {n_gold_standard})")
power_coexpr = 0.85 if n_gold_standard >= 100 else 0.70
feasible_coexpr = n_gold_standard >= 50
print(f"     Power: {power_coexpr:.2f} {'✓' if feasible_coexpr else '⚠'}")
print(f"     Methods: WGCNA, correlation networks, module detection")
network_approaches.append(('Co-expression Networks', power_coexpr, feasible_coexpr))
print()

# 2. Mutation-expression regulatory networks
print("  2. Mutation-Expression Regulatory Networks")
print(f"     Samples needed: ≥50 with mutations + expression")
n_mut_expr = ((master_map['has_mutations']) & (master_map['has_expression'])).sum()
print(f"     Available: n = {n_mut_expr}")
power_mutreg = 0.80 if n_mut_expr >= 100 else 0.65
feasible_mutreg = n_mut_expr >= 50
print(f"     Power: {power_mutreg:.2f} {'✓' if feasible_mutreg else '⚠'}")
print(f"     Methods: Directed graphs, causal inference")
network_approaches.append(('Mutation-Expression Networks', power_mutreg, feasible_mutreg))
print()

# 3. Drug-target-pathway networks
print("  3. Drug-Target-Pathway Networks")
print(f"     Samples needed: ≥50 with drug + expression")
print(f"     Available: n = {n_expr_drug}")
power_drugtarget = 0.80 if n_expr_drug >= 100 else 0.65
feasible_drugtarget = n_expr_drug >= 50
print(f"     Power: {power_drugtarget:.2f} {'✓' if feasible_drugtarget else '⚠'}")
print(f"     Methods: Pathway enrichment, drug target databases (DGIdb, DrugBank)")
network_approaches.append(('Drug-Target Networks', power_drugtarget, feasible_drugtarget))
print()

# 4. Integrated multi-layer networks
print("  4. Integrated Multi-Layer Networks")
print(f"     Samples needed: ≥50 with all 4 data types")
print(f"     Available: n = {n_gold_standard}")
power_multilayer = 0.75 if n_gold_standard >= 100 else 0.60
feasible_multilayer = n_gold_standard >= 50
print(f"     Power: {power_multilayer:.2f} {'✓' if feasible_multilayer else '⚠'}")
print(f"     Methods: Similarity Network Fusion (SNF), multi-omics factor analysis")
network_approaches.append(('Multi-Layer Networks', power_multilayer, feasible_multilayer))
print()

# Overall feasibility
mean_power_network = np.mean([a[1] for a in network_approaches])
feasible_network_count = sum([1 for a in network_approaches if a[2]])

analysis_8['statistical_power'] = mean_power_network
analysis_8['feasibility'] = 'YES' if feasible_network_count >= 3 else 'LIMITED'
analysis_8['justification'] = f"n={n_gold_standard} quad-omics, {feasible_network_count}/{len(network_approaches)} approaches feasible"
analysis_8['network_approaches'] = '; '.join([a[0] for a in network_approaches if a[2]])

deliverables_8 = [
    "Gene-gene co-expression networks (WGCNA modules)",
    "Mutation-expression regulatory networks",
    "Drug-target-pathway integration maps",
    "Multi-layer network visualizations (Cytoscape)",
    "Hub gene/node identification",
    "Network-based biomarker discovery",
    "Pathway enrichment analysis",
    "Network communities/modules linked to clinical outcomes"
]
analysis_8['deliverables'] = '; '.join(deliverables_8)

print("SUMMARY:")
print(f"  Statistical power: {analysis_8['statistical_power']:.2f}")
print(f"  Feasibility: {analysis_8['feasibility']}")
print(f"  Justification: {analysis_8['justification']}")
print(f"  Feasible approaches: {analysis_8['network_approaches']}")
print(f"  Timeline: {analysis_8['timeline_weeks']}")
print()

# ============================================================================
# ANALYSIS 9: DRUG MECHANISM DISCOVERY
# ============================================================================

print("=" * 80)
print("ANALYSIS 3.5: DRUG MECHANISM DISCOVERY")
print("=" * 80)
print()

analysis_9 = {
    'tier': 3,
    'analysis_id': '3.5',
    'analysis_name': 'Drug Mechanism Discovery',
    'goal': 'Understand molecular changes associated with drug response',
    'methods': 'Correlation analysis, differential expression, pathway enrichment',
    'required_n': 100,  # For robust correlation
    'available_n': n_expr_drug,
    'data_requirements': 'Expression + Drug Response',
    'timeline_weeks': '3-4 weeks',
    'scripts_location': '02_Scripts/05_Drug_Response/',
}

print(f"Goal: {analysis_9['goal']}")
print(f"Method: {analysis_9['methods']}")
print()

print(f"Sample Size:")
print(f"  Required: n ≥ {analysis_9['required_n']} with expression + drug")
print(f"  Available: n = {analysis_9['available_n']}")
print()

# Analysis components
print("ANALYSIS COMPONENTS:")
print()

mechanism_components = []

# 1. Gene-drug correlation analysis
print("  1. Gene-Drug Correlation Analysis")
print(f"     Goal: Identify genes whose expression predicts drug sensitivity")
print(f"     Samples: n = {n_expr_drug}")
print(f"     Method: Spearman/Pearson correlation, FDR correction")
# For correlation, need ~100 samples for power ~0.8 to detect r=0.3
power_corr = 0.85 if n_expr_drug >= 100 else 0.70
feasible_corr = n_expr_drug >= 80
print(f"     Power: {power_corr:.2f} {'✓' if feasible_corr else '⚠'}")
mechanism_components.append(('Gene-Drug Correlation', power_corr, feasible_corr))
print()

# 2. Pathway-level drug associations
print("  2. Pathway-Level Drug Associations")
print(f"     Goal: Identify pathways enriched in drug-sensitive vs resistant")
print(f"     Samples: n = {n_expr_drug}")
print(f"     Method: GSEA, pathway enrichment (KEGG, Reactome)")
power_pathway = 0.80 if n_expr_drug >= 100 else 0.65
feasible_pathway = n_expr_drug >= 80
print(f"     Power: {power_pathway:.2f} {'✓' if feasible_pathway else '⚠'}")
mechanism_components.append(('Pathway Associations', power_pathway, feasible_pathway))
print()

# 3. Drug response signatures
print("  3. Drug Response Signatures")
print(f"     Goal: Define gene signatures predictive of drug response")
print(f"     Samples: n = {n_expr_drug}")
print(f"     Method: Split into sensitive/resistant, identify DEGs")
# Split into tertiles (top/bottom 33%) for sensitive vs resistant
n_per_group = int(n_expr_drug * 0.33)
print(f"     Per group: ~{n_per_group} (top/bottom tertiles)")
power_signature = 0.85 if n_per_group >= 50 else 0.70
feasible_signature = n_per_group >= 30
print(f"     Power: {power_signature:.2f} {'✓' if feasible_signature else '⚠'}")
mechanism_components.append(('Drug Response Signatures', power_signature, feasible_signature))
print()

# 4. Mutation-drug mechanism integration
print("  4. Mutation-Drug Mechanism Integration")
print(f"     Goal: Link mutations to drug mechanism of action")
print(f"     Samples: n = {n_expr_mut_drug} (with mutations)")
print(f"     Method: Stratify by mutation, analyze drug response patterns")
power_mutmech = 0.75 if n_expr_mut_drug >= 100 else 0.60
feasible_mutmech = n_expr_mut_drug >= 50
print(f"     Power: {power_mutmech:.2f} {'✓' if feasible_mutmech else '⚠'}")
mechanism_components.append(('Mutation-Drug Mechanisms', power_mutmech, feasible_mutmech))
print()

# Overall feasibility
mean_power_mech = np.mean([c[1] for c in mechanism_components])
feasible_mech_count = sum([1 for c in mechanism_components if c[2]])

analysis_9['statistical_power'] = mean_power_mech
analysis_9['feasibility'] = 'YES' if feasible_mech_count >= 3 else 'LIMITED'
analysis_9['justification'] = f"n={n_expr_drug} expr+drug, {feasible_mech_count}/{len(mechanism_components)} components feasible"
analysis_9['analysis_components'] = '; '.join([c[0] for c in mechanism_components if c[2]])

deliverables_9 = [
    "Gene-drug correlation matrix (top 1000 genes × all drugs)",
    "Pathway enrichment results for each drug class",
    "Drug response signatures (sensitive vs resistant)",
    "Volcano plots for key drugs",
    "Drug mechanism hypotheses based on molecular profiles",
    "Target validation (literature comparison)",
    "Drug repurposing opportunities"
]
analysis_9['deliverables'] = '; '.join(deliverables_9)

print("SUMMARY:")
print(f"  Statistical power: {analysis_9['statistical_power']:.2f}")
print(f"  Feasibility: {analysis_9['feasibility']}")
print(f"  Justification: {analysis_9['justification']}")
print(f"  Feasible components: {analysis_9['analysis_components']}")
print(f"  Timeline: {analysis_9['timeline_weeks']}")
print()

# ============================================================================
# ANALYSIS 10: PERSONALIZED TREATMENT RECOMMENDATION FRAMEWORK
# ============================================================================

print("=" * 80)
print("ANALYSIS 3.6: PERSONALIZED TREATMENT RECOMMENDATION FRAMEWORK")
print("=" * 80)
print()

analysis_10 = {
    'tier': 3,
    'analysis_id': '3.6',
    'analysis_name': 'Personalized Treatment Recommendation Framework',
    'goal': 'Develop framework for patient-specific drug selection based on molecular profile',
    'methods': 'Integrate molecular profile → drug prediction → clinical recommendation',
    'required_n': 100,
    'available_n': n_gold_standard,
    'data_requirements': 'Complete quad-omics (for training); Expression + Mutations minimum (for deployment)',
    'timeline_weeks': '4-6 weeks',
    'scripts_location': '02_Scripts/08_Personalized_Medicine/',
}

print(f"Goal: {analysis_10['goal']}")
print(f"Method: {analysis_10['methods']}")
print()

print(f"Sample Size:")
print(f"  Training cohort (complete quad-omics): n = {n_gold_standard}")
print(f"  Deployment cohort (expression + mutations): n = {n_expr_mut_drug}")
print()

# Framework components
print("FRAMEWORK COMPONENTS:")
print()

framework_components = []

# 1. Patient molecular profile integration
print("  1. Patient Molecular Profile Integration")
print(f"     Goal: Aggregate all molecular data into patient profile")
print(f"     Training samples: n = {n_gold_standard}")
print(f"     Components:")
print(f"       - Molecular subtype assignment")
print(f"       - Key mutation status (FLT3, NPM1, TP53, etc.)")
print(f"       - Gene expression features (top variance genes)")
print(f"       - Clinical features (age, ELN risk, etc.)")
power_profile = 0.85 if n_gold_standard >= 100 else 0.70
feasible_profile = n_gold_standard >= 50
print(f"     Power: {power_profile:.2f} {'✓' if feasible_profile else '⚠'}")
framework_components.append(('Profile Integration', power_profile, feasible_profile))
print()

# 2. Drug prediction model (from Analysis 1.4)
print("  2. Drug Response Prediction Models")
print(f"     Goal: Predict drug sensitivity from molecular profile")
print(f"     Training samples: n = {n_gold_standard}")
print(f"     Models: Random Forest, Elastic Net, XGBoost")
print(f"     Based on Analysis 1.4 results")
power_prediction = 0.85 if n_gold_standard >= 150 else 0.70
feasible_prediction = n_gold_standard >= 100
print(f"     Power: {power_prediction:.2f} {'✓' if feasible_prediction else '⚠'}")
framework_components.append(('Drug Prediction Models', power_prediction, feasible_prediction))
print()

# 3. Decision support algorithm
print("  3. Decision Support Algorithm")
print(f"     Goal: Rank drugs for each patient based on predicted sensitivity")
print(f"     Method:")
print(f"       - Predict AUC for all drugs")
print(f"       - Rank by sensitivity (lower AUC = more sensitive)")
print(f"       - Filter by FDA approval status")
print(f"       - Consider mutation-drug associations")
print(f"       - Incorporate clinical contraindications")
power_decision = 0.80
feasible_decision = n_gold_standard >= 100
print(f"     Power: {power_decision:.2f} {'✓' if feasible_decision else '⚠'}")
framework_components.append(('Decision Support', power_decision, feasible_decision))
print()

# 4. Clinical validation framework
print("  4. Clinical Validation Framework")
print(f"     Goal: Validate predictions against clinical outcomes")
print(f"     Validation samples: Use held-out test set")
print(f"     Metrics:")
print(f"       - Prediction accuracy")
print(f"       - Clinical concordance")
print(f"       - Survival benefit (if treatment data available)")
power_validation = 0.75 if n_gold_standard >= 100 else 0.60
feasible_validation = n_gold_standard >= 100
print(f"     Power: {power_validation:.2f} {'✓' if feasible_validation else '⚠'}")
framework_components.append(('Clinical Validation', power_validation, feasible_validation))
print()

# Overall feasibility
mean_power_framework = np.mean([c[1] for c in framework_components])
feasible_framework_count = sum([1 for c in framework_components if c[2]])

analysis_10['statistical_power'] = mean_power_framework
analysis_10['feasibility'] = 'YES' if feasible_framework_count >= 3 else 'LIMITED'
analysis_10['justification'] = f"n={n_gold_standard} for training, {feasible_framework_count}/{len(framework_components)} components feasible"
analysis_10['framework_components'] = '; '.join([c[0] for c in framework_components if c[2]])

deliverables_10 = [
    "Patient molecular profile dashboard (web interface or R Shiny)",
    "Drug ranking algorithm with confidence scores",
    "Validation report (accuracy, concordance)",
    "Clinical decision support tool",
    "User manual and deployment guide",
    "Case studies (example patients with recommendations)",
    "Comparison with standard-of-care treatment selection"
]
analysis_10['deliverables'] = '; '.join(deliverables_10)

print("SUMMARY:")
print(f"  Statistical power: {analysis_10['statistical_power']:.2f}")
print(f"  Feasibility: {analysis_10['feasibility']}")
print(f"  Justification: {analysis_10['justification']}")
print(f"  Feasible components: {analysis_10['framework_components']}")
print(f"  Timeline: {analysis_10['timeline_weeks']}")
print()

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("=" * 80)
print("SAVING TIER 3 ANALYSES")
print("=" * 80)
print()

# Create DataFrame
tier3_analyses = [analysis_8, analysis_9, analysis_10]
df_tier3 = pd.DataFrame(tier3_analyses)

# Save to CSV
output_file = OUTPUT_DIR / "tier3_advanced_analyses.csv"
df_tier3.to_csv(output_file, index=False)
print(f"✓ Saved: {output_file.name}")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 80)
print("TIER 3 ANALYSES COMPLETE")
print("=" * 80)
print()
print(f"TIER 3: {len(tier3_analyses)} advanced analyses")
print()
for analysis in tier3_analyses:
    feasible_icon = '✓' if analysis['feasibility'] == 'YES' else '⚠'
    print(f"  {feasible_icon} {analysis['analysis_id']}: {analysis['analysis_name']}")
    print(f"     n={analysis['available_n']}, Power={analysis['statistical_power']:.2f}, {analysis['feasibility']}")
print()

feasible_tier3 = sum([1 for a in tier3_analyses if a['feasibility'] == 'YES'])
print(f"FEASIBILITY: {feasible_tier3}/{len(tier3_analyses)} fully feasible")
print()
print("OUTPUTS SAVED:")
print(f"  {output_file.name}")
print()
print("=" * 80)
