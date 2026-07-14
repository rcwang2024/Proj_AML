import fitz
import os
import shutil


def get_best_fit_rect(quad_rect, panel_w, panel_h):
    quad_w = quad_rect.x1 - quad_rect.x0
    quad_h = quad_rect.y1 - quad_rect.y0
    scale = min(quad_w / panel_w, quad_h / panel_h)
    new_w = panel_w * scale
    new_h = panel_h * scale
    dx = (quad_w - new_w) / 2
    dy = (quad_h - new_h) / 2
    return fitz.Rect(quad_rect.x0 + dx, quad_rect.y0 + dy, quad_rect.x0 + dx + new_w, quad_rect.y0 + dy + new_h)

def assemble_s3_high_fid():
    path_in = r"d:\Proj_AML\04_Figures\12_ELN_Comparison\HighFid"
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS3.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS3.png"
    dest_pdf_2 = r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS3.pdf"
    dest_png_2 = r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS3.png"
    
    # Compilation folder for LaTeX
    dest_pdf_3 = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS3.pdf"
    dest_png_3 = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS3.png"
    
    # 2x2 Layout
    panels = [
        "s3_pA.pdf", # Top-Left
        "s3_pB.pdf", # Top-Right
        "s3_pC.pdf", # Bottom-Left
        "s3_pD.pdf"  # Bottom-Right
    ]
    
    doc_out = fitz.open()
    # Square page for 2x2 (2000pt x 2000pt)
    page_out = doc_out.new_page(width=2000, height=2000)
    
    margin = 50
    w = 925
    h = 925
    
    for i, p_name in enumerate(panels):
        p_path = os.path.join(path_in, p_name)
        if not os.path.exists(p_path):
            print(f"Missing {p_path}")
            continue
            
        col = i % 2
        row = i // 2
        
        rect = fitz.Rect(
            margin + col*(w + margin),
            margin + row*(h + margin),
            margin + col*(w + margin) + w,
            margin + row*(h + margin) + h
        )
        
        src = fitz.open(p_path)
        
        actual_w = src[0].rect.width
        actual_h = src[0].rect.height
        fit_rect = get_best_fit_rect(rect, actual_w, actual_h)
        page_out.show_pdf_page(fit_rect, src, 0)
        src.close()
        
    # Ensure destination directories exist
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    os.makedirs(os.path.dirname(dest_pdf_2), exist_ok=True)
    os.makedirs(os.path.dirname(dest_pdf_3), exist_ok=True)
    
    # Save first output (Submission_Hub)
    doc_out.save(dest_pdf_1)
    
    # Render PNG
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png_1)
    print(f"Successfully finalized High-Fidelity Figure S3 in Submission_Hub: {dest_pdf_1}")
    
    # Close document
    doc_out.close()
    
    # Copy helper to handle potential samefile (symlinks/junctions)
    def copy_file_safe(src, dst, label):
        try:
            if os.path.samefile(src, dst):
                print(f"Junction/link detected: {dst} points to the same file as {src}. No copy needed.")
            else:
                shutil.copy(src, dst)
                print(f"Successfully copied {label} to: {dst}")
        except Exception:
            shutil.copy(src, dst)
            print(f"Successfully copied {label} to: {dst}")

    # Copy to Submission_Hub_2026-05-16/03_Supplementary_Figures/
    copy_file_safe(dest_pdf_1, dest_pdf_2, "PDF")
    copy_file_safe(dest_png_1, dest_png_2, "PNG")
    
    # Copy to 01_Manuscript_Source/Supplementary/
    copy_file_safe(dest_pdf_1, dest_pdf_3, "PDF (LaTeX Compilation)")
    copy_file_safe(dest_png_1, dest_png_3, "PNG (LaTeX Compilation)")

if __name__ == "__main__":
    assemble_s3_high_fid()
