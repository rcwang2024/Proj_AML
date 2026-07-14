import subprocess
import os
import shutil
import fitz

def main():
    print("=== STARTING FIGURE S6 GENERATION AND SYNCHRONIZATION ===")
    
    # 1. Run R script to generate raw component PDFs
    r_script = r"d:\Proj_AML\02_Scripts\Phase11_Finalization\06_generate_FigureS6_HighFid.R"
    print(f"Running component generator: {r_script}")
    subprocess.run(["Rscript", r_script], check=True)
    
    # 2. Canvas assembly settings
    width, height = 2000, 1800
    print("\nAssembling high-fidelity Figure S6 canvas (2000pt x 1800pt)...")
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=width, height=height)
    
    margin = 50
    w_half = (2000 - 3*margin) / 2
    row_spacing = 100
    
    base_path = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts"
    
    # Row 1: A and B
    h_top = 900
    rect_a = fitz.Rect(margin, margin, margin + w_half, margin + h_top)
    rect_b = fitz.Rect(2000 - margin - w_half, margin, 2000 - margin, margin + h_top)
    
    # Row 2: C (Full Width)
    y2 = margin + h_top + row_spacing
    h_bottom = 600
    rect_c = fitz.Rect(margin, y2, 2000 - margin, y2 + h_bottom)
    
    panels = [
        ("s6_pA.pdf", rect_a),
        ("s6_pB.pdf", rect_b),
        ("s6_pC.pdf", rect_c)
    ]
    
    for filename, rect in panels:
        path = os.path.join(base_path, filename)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Component not found: {path}")
            
        src_doc = fitz.open(path)
        page_out.show_pdf_page(rect, src_doc, 0)
        src_doc.close()
        print(f"  [OK] Placed {filename}")
        
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS6.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS6.png"
    
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    doc_out.save(dest_pdf_1)
    
    # Render PNG preview at ultra high res (300 DPI)
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png_1)
    doc_out.close()
    print(f"  [OK] Successfully compiled primary consolidated Figure S6: {dest_pdf_1}")
    
    # 3. Synchronize with target folders
    targets = [
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS6.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS6.png"
        ),
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS6.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS6.png"
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
        
    print("\n=== FIGURE S6 COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
