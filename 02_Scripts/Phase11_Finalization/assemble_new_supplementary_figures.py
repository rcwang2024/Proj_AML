# ==============================================================================
# PHASE 11: FINALIZATION
# Script: assemble_new_supplementary_figures.py
# Purpose: Assemble consolidated high-fidelity Figure S7 (DepMap target and salvage)
#          and Figure S9 (proteomic validation).
# ==============================================================================
import fitz
import os

def assemble_s7_and_s9():
    print("=== STARTING NEW SUPPLEMENTARY FIGURES ASSEMBLY ===")
    
    src_dir = r"d:\Proj_AML\04_Figures\Supplementary_Figures"
    dest_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures"
    os.makedirs(dest_dir, exist_ok=True)
    
    # --------------------------------------------------------------------------
    # 1. Assemble Figure S7 (Panels A, B, C, D)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S7...")
    s7_src_path = os.path.join(src_dir, "Supplementary_Figure_S7.pdf")
    
    if os.path.exists(s7_src_path):
        doc_s7 = fitz.open(s7_src_path)
        page_s7 = doc_s7[0]
        
        dest_s7_pdf = os.path.join(dest_dir, "FigureS7.pdf")
        dest_s7_png = os.path.join(dest_dir, "FigureS7.png")
        
        doc_s7.save(dest_s7_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s7.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s7_png)
        doc_s7.close()
        print(f"  [SUCCESS] Finalized Figure S7 PDF & PNG in: {dest_dir}")
    else:
        print(f"  [ERROR] Source component for Figure S7 (Supplementary_Figure_S7.pdf) not found!")
        
    # --------------------------------------------------------------------------
    # 2. Assemble Figure S8 (Proteomics Cohort Validation)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S8...")
    s8_src_path = os.path.join(src_dir, "Supplementary_Figure_S9.pdf") # Formerly S9 proteomics, now repurposed as S8
    
    if os.path.exists(s8_src_path):
        # Repurposed: Proteomic bulk-level validation of pVRS score
        doc_s8 = fitz.open(s8_src_path)
        page_s8 = doc_s8[0]
        
        dest_s8_pdf = os.path.join(dest_dir, "FigureS8.pdf")
        dest_s8_png = os.path.join(dest_dir, "FigureS8.png")
        
        doc_s8.save(dest_s8_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s8.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s8_png)
        doc_s8.close()
        print(f"  [SUCCESS] Finalized Figure S8 PDF & PNG in: {dest_dir}")
    else:
        print(f"  [ERROR] Source component for Figure S8 (Supplementary_Figure_S9.pdf) not found!")

    # --------------------------------------------------------------------------
    # 2.5. Assemble Figure S9 (European Cohorts Validation Boxplot)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S9...")
    s9_src_path = r"d:\Proj_AML\04_Figures\Phase14_AMLSG_Validation\Figure5_VRS_validation_european.pdf"
    
    if os.path.exists(s9_src_path):
        doc_s9 = fitz.open(s9_src_path)
        page_s9 = doc_s9[0]
        
        dest_s9_pdf = os.path.join(dest_dir, "FigureS9.pdf")
        dest_s9_png = os.path.join(dest_dir, "FigureS9.png")
        
        doc_s9.save(dest_s9_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s9.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s9_png)
        doc_s9.close()
        print(f"  [SUCCESS] Finalized Figure S9 PDF & PNG in: {dest_dir}")
    else:
        print(f"  [ERROR] Source component for Figure S9 (Figure5_VRS_validation_european.pdf) not found!")

    # --------------------------------------------------------------------------
    # 3. Assemble Figure S11 (Historical trial validation)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S11...")
    s11_src_path = os.path.join(src_dir, "Supplementary_Figure_S11.pdf")
    
    if os.path.exists(s11_src_path):
        doc_s11 = fitz.open(s11_src_path)
        page_s11 = doc_s11[0]
        
        dest_s11_pdf = os.path.join(dest_dir, "FigureS11.pdf")
        dest_s11_png = os.path.join(dest_dir, "FigureS11.png")
        
        doc_s11.save(dest_s11_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s11.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s11_png)
        doc_s11.close()
        print(f"  [SUCCESS] Finalized Figure S11 PDF & PNG in: {dest_dir}")
    else:
        print(f"  [ERROR] Source component for Figure S11 (Supplementary_Figure_S11.pdf) not found!")
        
    # --------------------------------------------------------------------------
    # 4. Assemble Figure S12 (CITE-seq Single-Cell ADT Validation)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S12...")
    s12_src_path = os.path.join(src_dir, "Supplementary_Figure_S12.pdf")
    
    if os.path.exists(s12_src_path):
        doc_s12 = fitz.open(s12_src_path)
        page_s12 = doc_s12[0]
        
        dest_s12_pdf = os.path.join(dest_dir, "FigureS12.pdf")
        dest_s12_png = os.path.join(dest_dir, "FigureS12.png")
        
        doc_s12.save(dest_s12_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s12.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s12_png)
        doc_s12.close()
        print(f"  [SUCCESS] Finalized Figure S12 PDF & PNG in: {dest_dir}")
    # --------------------------------------------------------------------------
    # 5. Assemble Figure S13 (CITE-seq Single-Cell Longitudinal Dynamics)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S13...")
    s13_src_path = os.path.join(src_dir, "Supplementary_Figure_S13.pdf")
    
    if os.path.exists(s13_src_path):
        doc_s13 = fitz.open(s13_src_path)
        page_s13 = doc_s13[0]
        
        dest_s13_pdf = os.path.join(dest_dir, "FigureS13.pdf")
        dest_s13_png = os.path.join(dest_dir, "FigureS13.png")
        
        doc_s13.save(dest_s13_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s13.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s13_png)
        doc_s13.close()
        print(f"  [SUCCESS] Finalized Figure S13 PDF & PNG in: {dest_dir}")
    else:
        print(f"  [ERROR] Source component for Figure S13 (Supplementary_Figure_S13.pdf) not found!")
        
    # --------------------------------------------------------------------------
    # 6. Assemble Figure S14 (CITE-seq Single-Cell UMAP Projections)
    # --------------------------------------------------------------------------
    print("\nAssembling Figure S14...")
    s14_src_path = os.path.join(src_dir, "Supplementary_Figure_S14.pdf")
    
    if os.path.exists(s14_src_path):
        doc_s14 = fitz.open(s14_src_path)
        page_s14 = doc_s14[0]
        
        dest_s14_pdf = os.path.join(dest_dir, "FigureS14.pdf")
        dest_s14_png = os.path.join(dest_dir, "FigureS14.png")
        
        doc_s14.save(dest_s14_pdf)
        
        # Save high-res PNG (600 DPI equivalent)
        pix = page_s14.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
        pix.save(dest_s14_png)
        doc_s14.close()
        print(f"  [SUCCESS] Finalized Figure S14 PDF & PNG in: {dest_dir}")
    else:
        print(f"  [ERROR] Source component for Figure S14 (Supplementary_Figure_S14.pdf) not found!")
        
    print("\n=== NEW SUPPLEMENTARY FIGURES ASSEMBLY COMPLETED ===")

if __name__ == "__main__":
    assemble_s7_and_s9()
