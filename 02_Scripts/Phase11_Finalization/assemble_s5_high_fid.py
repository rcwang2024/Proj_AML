import fitz
import os


def get_best_fit_rect(quad_rect, panel_w, panel_h):
    quad_w = quad_rect.x1 - quad_rect.x0
    quad_h = quad_rect.y1 - quad_rect.y0
    scale = min(quad_w / panel_w, quad_h / panel_h)
    new_w = panel_w * scale
    new_h = panel_h * scale
    dx = (quad_w - new_w) / 2
    dy = (quad_h - new_h) / 2
    return fitz.Rect(quad_rect.x0 + dx, quad_rect.y0 + dy, quad_rect.x0 + dx + new_w, quad_rect.y0 + dy + new_h)

def assemble_s5():
    print("=== ASSEMBLING HIGH-FIDELITY FIGURE S5 (Benchmarking) ===")
    
    path_in = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts"
    dest_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS5.pdf"
    dest_png = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS5.png"
    
    # 2x2 Layout
    panels = [
        "s5_pA.pdf",
        "s5_pB.pdf",
        "s5_pC.pdf",
        "s5_pD.pdf"
    ]
    
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=2000, height=1800)
    
    margin = 50
    w = 925
    h = 825
    
    for i, p_name in enumerate(panels):
        p_path = os.path.join(path_in, p_name)
        if not os.path.exists(p_path):
            print(f"Missing {p_path}")
            continue
            
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
        print(f"  [OK] Inserted {p_name}")
        
    os.makedirs(os.path.dirname(dest_pdf), exist_ok=True)
    doc_out.save(dest_pdf)
    
    # Ultra high res PNG for review (600 DPI equivalent)
    pix = page_out.get_pixmap(matrix=fitz.Matrix(600/72, 600/72))
    pix.save(dest_png)
    print(f"Successfully finalized High-Fidelity Figure S5: {dest_pdf}")
    doc_out.close()

if __name__ == "__main__":
    assemble_s5()
