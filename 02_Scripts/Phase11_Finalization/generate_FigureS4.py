import subprocess
import os
import shutil
import fitz

def main():
    print("=== STARTING FIGURE S4 GENERATION AND SYNCHRONIZATION ===")
    
    # 1. Run all three component generator scripts
    scripts = [
        (r"d:\Proj_AML\02_Scripts\Phase2_Validation\12_immune_cell_deconvolution_simplified.R", "MCP-counter Heatmap (Panel A)"),
        (r"d:\Proj_AML\02_Scripts\Phase11_Finalization\generate_S4A_header.R", "Panel A Header"),
        (r"d:\Proj_AML\02_Scripts\Phase11_Finalization\04_generate_FigureS4_HighFid.R", "Boxplots & Correlation (Panel B/C)")
    ]
    
    for path, desc in scripts:
        print(f"\nRunning generator: {desc} ({path})")
        subprocess.run(["Rscript", path], check=True)
        
    # 2. Canvas assembly settings
    heatmap_path = r"d:\Proj_AML\04_Figures\15_Immune_Deconvolution\immune_cell_heatmap.pdf"
    header_path = r"d:\Proj_AML\04_Figures\15_Immune_Deconvolution\FigureS4_HeaderA.pdf"
    boxplot_path = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts\s4_pB_upgraded.pdf"
    corr_path = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts\s4_pC.pdf"
    
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS4.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS4.png"
    
    print("\nAssembling high-fidelity Figure S4 canvas (2000pt x 2100pt)...")
    doc_out = fitz.open()
    # 3-panel layout: True Wide Heatmap top, Boxplot/Corr bottom
    # Stretched height to 2100 to fit enlarged Plot A
    page_out = doc_out.new_page(width=2000, height=2100)
    
    margin = 50
    gap = 50
    h_heatmap = 900 
    w_half = (2000 - 2*margin - gap) / 2
    h_bottom = 850
    
    # 2.1 Place Header
    rect_header = fitz.Rect(margin, margin, 2000-margin, margin + 120)
    src_header = fitz.open(header_path)
    page_out.show_pdf_page(rect_header, src_header, 0)
    src_header.close()
    print("  [OK] Placed Panel A Header")

    # 2.2 Place Heatmap (Stretched and high row text legibility)
    rect_a = fitz.Rect(margin, margin + 120, 2000-margin, margin + 120 + h_heatmap)
    src_a = fitz.open(heatmap_path)
    page_out.show_pdf_page(rect_a, src_a, 0)
    src_a.close()
    print("  [OK] Placed Panel A Heatmap")

    # 2.3 Place Boxplots
    rect_b = fitz.Rect(margin, margin + 120 + h_heatmap + 50, margin + w_half, margin + 120 + h_heatmap + 50 + h_bottom)
    src_b = fitz.open(boxplot_path)
    page_out.show_pdf_page(rect_b, src_b, 0)
    src_b.close()
    print("  [OK] Placed Panel B Boxplots")
    
    # 2.4 Place Correlation
    rect_c = fitz.Rect(margin + w_half + gap, margin + 120 + h_heatmap + 50, 2000-margin, margin + 120 + h_heatmap + 50 + h_bottom)
    src_c = fitz.open(corr_path)
    page_out.show_pdf_page(rect_c, src_c, 0)
    src_c.close()
    print("  [OK] Placed Panel C Correlation")
    
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    doc_out.save(dest_pdf_1)
    
    # Render PNG preview at ultra high res (300 DPI)
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png_1)
    doc_out.close()
    print(f"  [OK] Successfully compiled primary consolidated Figure S4: {dest_pdf_1}")
    
    # 3. Synchronize with target folders
    targets = [
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS4.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS4.png"
        ),
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS4.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS4.png"
        )
    ]
    
    def copy_file_safe(src, dst):
        try:
            if os.path.exists(dst) and os.path.samefile(src, dst):
                print(f"  [INFO] Junction/link detected: {dst} is identical to {src}. No copy needed.")
            else:
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(src, dst)
                print(f"  [OK] Copied to: {dst}")
        except Exception:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)
            print(f"  [OK] Copied to: {dst} (via backup handler)")

    print("\nSynchronizing outputs across workspace targets...")
    for pdf_tgt, png_tgt in targets:
        copy_file_safe(dest_pdf_1, pdf_tgt)
        copy_file_safe(dest_png_1, png_tgt)
        
    print("\n=== FIGURE S4 COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
