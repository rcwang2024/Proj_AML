import subprocess
import os
import shutil
import fitz

def get_best_fit_rect(quad_rect, panel_w, panel_h):
    quad_w = quad_rect.x1 - quad_rect.x0
    quad_h = quad_rect.y1 - quad_rect.y0
    scale = min(quad_w / panel_w, quad_h / panel_h)
    new_w = panel_w * scale
    new_h = panel_h * scale
    dx = (quad_w - new_w) / 2
    dy = (quad_h - new_h) / 2
    return fitz.Rect(quad_rect.x0 + dx, quad_rect.y0 + dy, quad_rect.x0 + dx + new_w, quad_rect.y0 + dy + new_h)

def main():
    print("=== STARTING FIGURE S2 GENERATION AND SYNCHRONIZATION (PEDIATRIC TARGET-AML) ===")
    
    # 1. Run R script to generate raw component PDFs
    r_script = r"d:\Proj_AML\02_Scripts\Phase11_Finalization\generate_FigureS2_HighFid_Pediatric.R"
    print(f"Running component generator: {r_script}")
    subprocess.run(["Rscript", r_script], check=True)
    
    # 2. Copy compiled pediatric KM plot to primary destination
    src_pdf = r"d:\Proj_AML\04_Figures\18_TARGET_Validation\HighFid\FigureS2_Pediatric_KM.pdf"
    src_png = r"d:\Proj_AML\04_Figures\18_TARGET_Validation\HighFid\FigureS2_Pediatric_KM.png"
    
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS2.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS2.png"
    
    print("\nCopying high-fidelity Figure S2...")
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    shutil.copy(src_pdf, dest_pdf_1)
    shutil.copy(src_png, dest_png_1)
    print(f"  [OK] Copied S2 PDF to {dest_pdf_1}")
    print(f"  [OK] Copied S2 PNG to {dest_png_1}")
    
    # 3. Synchronize with target folders
    targets = [
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS2.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS2.png"
        ),
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS2.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS2.png"
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
        
    print("\n=== FIGURE S2 COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
