# Streamlit Clinical Calculator for 9-Gene Venetoclax Response Score (VRS)
# Providing a translation-ready clinical tool for Nature Medicine
import streamlit as st
import pandas as pd
import numpy as np

# Set page config
st.set_page_config(
    page_title="VRS Clinical Calculator",
    page_icon="🩸",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Title & Description
st.title("🩸 AML Venetoclax Response Score (VRS) Calculator")
st.markdown("""
This clinical decision-support tool calculates the **9-Gene Venetoclax Response Score (VRS)** for adult patients with Acute Myeloid Leukemia (AML). 
By evaluating transcriptomic markers, this tool predicts ex vivo sensitivity to the BCL-2 inhibitor **Venetoclax** independent of established genomic risk stratification.
""")

# Sidebar
st.sidebar.header("Navigation")
app_mode = st.sidebar.radio("Choose calculator mode", ["Single Patient Entry", "Batch Upload (CSV)"])

st.sidebar.markdown("---")
st.sidebar.markdown("""
### Clinical Cutoffs:
* **High VRS (> 71.0)**: Venetoclax strongly recommended.
* **Medium VRS (41.8 - 71.0)**: Use with caution, monitor response.
* **Low VRS (< 41.8)**: Resistant. Prioritize alternative salvage regimens.
""")

# Define VRS weights based on the clinical Random Forest / linear proxy coefficients
# Normalization: intercept and scale parameters to map to [0, 100]
coefficients = {
    "BCL2": 2.5,
    "NPM1": 1.2,
    "DNMT3A": 0.8,
    "IDH1": 0.5,
    "IDH2": 0.5,
    "FLT3": 1.0,
    "TP53": -1.5,
    "RUNX1": -0.8,
    "ASXL1": -0.8
}
intercept = -5.0
min_raw = -15.0
max_raw = 20.0

def calculate_vrs(genes_dict):
    raw_vrs = sum(genes_dict[gene] * coefficients[gene] for gene in coefficients) + intercept
    # Scale between 0 and 100
    vrs_score = (raw_vrs - min_raw) / (max_raw - min_raw) * 100
    vrs_score = max(min(vrs_score, 100.0), 0.0) # Clamp
    return round(vrs_score, 1)

def get_recommendation(vrs):
    if vrs > 71.0:
        return {
            "tier": "High VRS (Sensitive)",
            "color": "green",
            "rec": "**Venetoclax + HMA combination therapy is STRONGLY recommended.**",
            "details": "The transcriptomic profile indicates high dependency on the anti-apoptotic protein BCL-2. Patients in this category have historically achieved high response rates and durable remissions under the standard-of-care Venetoclax-Azacitidine regimen."
        }
    elif vrs < 41.8:
        return {
            "tier": "Low VRS (Resistant)",
            "color": "red",
            "rec": "**Venetoclax is NOT recommended. Prioritize alternative regimens.**",
            "details": "The profile indicates a monocytic-like, resistant blast phenotype dependent on MCL-1 or utilizing auxiliary metabolic pathways (dual OXPHOS/Glycolysis). Consider enrolling in clinical trials or using recommended salvage candidates:\n\n* **Panobinostat** (HDAC inhibitor)\n* **Sorafenib** (multikinase TKI)\n* **Selumetinib** (MEK inhibitor)"
        }
    else:
        return {
            "tier": "Medium VRS (Intermediate)",
            "color": "orange",
            "rec": "**Consider Venetoclax with caution. Implement close clinical monitoring.**",
            "details": "The patient profile is intermediate. Individualized clinical assessment is warranted, taking into account additional cytogenetic and clinical risk factors. If Venetoclax is initiated, monitor closely for early indicators of refractory disease or clonal evolution."
        }

if app_mode == "Single Patient Entry":
    st.header("Single Patient Entry")
    st.markdown("Enter normalized transcriptomic expression values (Z-scores or FPKM log2) for the 9 VRS genes:")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        bcl2 = st.number_input("BCL2 (BCL-2 target)", min_value=-5.0, max_value=10.0, value=2.5, step=0.1)
        npm1 = st.number_input("NPM1", min_value=-5.0, max_value=10.0, value=1.8, step=0.1)
        dnmt3a = st.number_input("DNMT3A", min_value=-5.0, max_value=10.0, value=1.5, step=0.1)
        
    with col2:
        idh1 = st.number_input("IDH1", min_value=-5.0, max_value=10.0, value=0.5, step=0.1)
        idh2 = st.number_input("IDH2", min_value=-5.0, max_value=10.0, value=0.8, step=0.1)
        flt3 = st.number_input("FLT3", min_value=-5.0, max_value=10.0, value=1.2, step=0.1)
        
    with col3:
        tp53 = st.number_input("TP53", min_value=-5.0, max_value=10.0, value=0.2, step=0.1)
        runx1 = st.number_input("RUNX1", min_value=-5.0, max_value=10.0, value=0.5, step=0.1)
        asxl1 = st.number_input("ASXL1", min_value=-5.0, max_value=10.0, value=0.4, step=0.1)
        
    if st.button("Calculate VRS"):
        patient_genes = {
            "BCL2": bcl2, "NPM1": npm1, "DNMT3A": dnmt3a,
            "IDH1": idh1, "IDH2": idh2, "FLT3": flt3,
            "TP53": tp53, "RUNX1": runx1, "ASXL1": asxl1
        }
        
        vrs = calculate_vrs(patient_genes)
        rec_dict = get_recommendation(vrs)
        
        st.markdown("---")
        st.subheader("Calculation Results")
        
        c_score, c_rec = st.columns([1, 2])
        
        with c_score:
            st.metric("9-Gene VRS Score", f"{vrs} / 100")
            if rec_dict["color"] == "green":
                st.success(f"**Classification:** {rec_dict['tier']}")
            elif rec_dict["color"] == "red":
                st.error(f"**Classification:** {rec_dict['tier']}")
            else:
                st.warning(f"**Classification:** {rec_dict['tier']}")
                
        with c_rec:
            st.markdown(f"### Clinical Recommendation:\n{rec_dict['rec']}")
            st.info(rec_dict["details"])
            
        # Display individual gene contributions
        st.markdown("### Feature Contribution Breakdown")
        contrib_df = pd.DataFrame([
            {"Gene": gene, "Expression Value": patient_genes[gene], "Coefficient": coefficients[gene], "Weight Contribution": round(patient_genes[gene] * coefficients[gene], 2)}
            for gene in coefficients
        ])
        st.dataframe(contrib_df, use_container_width=True)

else:
    st.header("Batch Upload (CSV)")
    st.markdown("Upload a CSV file containing normalized expression data for multiple patient samples. The CSV must contain a `sample_id` column and columns for all 9 genes: `BCL2`, `NPM1`, `DNMT3A`, `IDH1`, `IDH2`, `FLT3`, `TP53`, `RUNX1`, `ASXL1`.")
    
    uploaded_file = st.file_uploader("Upload CSV", type=["csv"])
    
    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)
            required_cols = ["sample_id"] + list(coefficients.keys())
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            if missing_cols:
                st.error(f"Missing columns in CSV: {', '.join(missing_cols)}")
            else:
                # Calculate scores
                scores = []
                tiers = []
                recs = []
                
                for idx, row in df.iterrows():
                    row_genes = {gene: row[gene] for gene in coefficients}
                    score = calculate_vrs(row_genes)
                    scores.append(score)
                    rec = get_recommendation(score)
                    tiers.append(rec["tier"])
                    recs.append(rec["rec"].replace("**", ""))
                    
                df["VRS"] = scores
                df["Clinical_Tier"] = tiers
                df["Recommendation"] = recs
                
                st.success(f"Successfully processed {len(df)} samples.")
                st.dataframe(df, use_container_width=True)
                
                # Download button
                csv_data = df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="Download Results CSV",
                    data=csv_data,
                    file_name="patient_VRS_recommendations.csv",
                    mime="text/csv"
                )
                
                # Class distribution plot
                st.subheader("Clinical Stratification Distribution")
                dist = df["Clinical_Tier"].value_counts()
                st.bar_chart(dist)
                
        except Exception as e:
            st.error(f"Error reading file: {e}")
