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

def copy_file_safe(src, dst):
    try:
        if not os.path.exists(src):
            print(f"  [WARNING] Source not found: {src}")
            return False
        if os.path.exists(dst) and os.path.samefile(src, dst):
            return True
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception as e:
        print(f"  [WARNING] Copy failed from {src} to {dst}: {e}")
        return False

def assemble_s6():
    print("=== ASSEMBLING HIGH-FIDELITY FIGURE S6 (Metabolism) ===")
    
    width, height = 2000, 1800
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=width, height=height)
    
    margin = 50
    w_half = (2000 - 3*margin) / 2
    row_spacing = 100
    
    base_path = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts"
    
    h_top = 900
    rect_a = fitz.Rect(margin, margin, margin + w_half, margin + h_top)
    rect_b = fitz.Rect(2000 - margin - w_half, margin, 2000 - margin, margin + h_top)
    
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
        if os.path.exists(path):
            src_doc = fitz.open(path)
            actual_w = src_doc[0].rect.width
            actual_h = src_doc[0].rect.height
            fit_rect = get_best_fit_rect(rect, actual_w, actual_h)
            page_out.show_pdf_page(fit_rect, src_doc, 0)
            src_doc.close()
            print(f"  [OK] Inserted {filename}")
        else:
            print(f"  [WARN] Missing {filename}")
            
    output_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures"
    dest_pdf = os.path.join(output_dir, "FigureS6.pdf")
    dest_png = os.path.join(output_dir, "FigureS6.png")
    
    os.makedirs(output_dir, exist_ok=True)
    doc_out.save(dest_pdf)
    
    pix = page_out.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
    pix.save(dest_png)
    doc_out.close()
    
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary"
    copy_file_safe(dest_pdf, os.path.join(latex_dir, "FigureS6.pdf"))
    copy_file_safe(dest_png, os.path.join(latex_dir, "FigureS6.png"))
    
    print(f"Successfully finalized High-Fidelity Figure S6: {dest_pdf}")

if __name__ == "__main__":
    assemble_s6()
