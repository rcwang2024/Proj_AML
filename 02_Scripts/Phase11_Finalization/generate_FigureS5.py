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
    print("=== STARTING FIGURE S5 GENERATION AND SYNCHRONIZATION ===")
    
    # 1. Run R script to generate raw component PDFs
    r_script = r"d:\Proj_AML\02_Scripts\Phase11_Finalization\05_generate_FigureS5_HighFid.R"
    print(f"Running component generator: {r_script}")
    subprocess.run(["Rscript", r_script], check=True)
    
    # 2. Canvas assembly settings
    path_in = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts"
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS5.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS5.png"
    
    panels = [
        "s5_pA.pdf",
        "s5_pB.pdf",
        "s5_pC.pdf",
        "s5_pD.pdf"
    ]
    
    print("\nAssembling high-fidelity Figure S5 canvas (2000pt x 1800pt)...")
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=2000, height=1800)
    
    margin = 50
    w = 925
    h = 825
    
    for i, p_name in enumerate(panels):
        p_path = os.path.join(path_in, p_name)
        if not os.path.exists(p_path):
            raise FileNotFoundError(f"Component not found: {p_path}")
            
        row = i // 2
        col = i % 2
        
        rect = fitz.Rect(
            margin + col * (w + margin),
            margin + row * (h + margin),
            margin + col * (w + margin) + w,
            margin + row * (h + margin) + h
        )
        
        src = fitz.open(p_path)
        actual_w = src[0].rect.width
        actual_h = src[0].rect.height
        fit_rect = get_best_fit_rect(rect, actual_w, actual_h)
        page_out.show_pdf_page(fit_rect, src, 0)
        src.close()
        print(f"  [OK] Placed {p_name}")
        
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    doc_out.save(dest_pdf_1)
    
    # Render PNG preview at ultra high res (300 DPI)
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png_1)
    doc_out.close()
    print(f"  [OK] Successfully compiled primary consolidated Figure S5: {dest_pdf_1}")
    
    # 3. Synchronize with target folders
    targets = [
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS5.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS5.png"
        ),
        (
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS5.pdf",
            r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS5.png"
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
        
    print("\n=== FIGURE S5 COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
