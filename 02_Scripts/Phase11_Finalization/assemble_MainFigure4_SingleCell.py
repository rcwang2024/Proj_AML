# ==============================================================================
# PHASE 11: FINALIZATION
# Script: assemble_MainFigure4_SingleCell.py
# Purpose: Assemble consolidated high-fidelity Main Figure 4 (Single-Cell Resolution
#          and Lineage Dynamics of Venetoclax Resistance)
# ==============================================================================
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

def assemble_main_figure4():
    print("=== ASSEMBLING MAIN FIGURE 4 (SINGLE-CELL RESOLUTION) ===")
    path_in = r"d:\Proj_AML\04_Figures\Phase11_Finalization"
    dest_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures"
    os.makedirs(dest_dir, exist_ok=True)
    
    dest_pdf = os.path.join(dest_dir, "Figure4_Consolidated.pdf")
    dest_png = os.path.join(dest_dir, "Figure4_Consolidated.png")
    
    panels = [
        "FigureS18_sc_UMAP_celltypes.pdf",
        "FigureS18_sc_UMAP_VRS.pdf",
        "FigureS18_sc_VRS_boxplot.pdf",
        "FigureS18_sc_dotplot.pdf",
        "FigureS18_sc_BCL2_MCL1_tradeoff.pdf"
    ]
    
    paths = [os.path.join(path_in, p) for p in panels]
    if any(not os.path.exists(p) for p in paths):
        print("  [ERROR] Main Figure 4 component files not found.")
        return
        
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=2000, height=2700)
    
    margin = 50
    w = 925
    h = 800
    gap = 50
    
    # Row 1: A and B
    for i in range(2):
        col = i % 2
        rect = fitz.Rect(margin + col*(w + gap), margin, margin + col*(w + gap) + w, margin + h)
        src = fitz.open(paths[i])
        page_out.show_pdf_page(get_best_fit_rect(rect, src[0].rect.width, src[0].rect.height), src, 0)
        src.close()
        
    # Row 2: C and D
    for i in range(2):
        col = i % 2
        rect = fitz.Rect(margin + col*(w + gap), margin + h + gap, margin + col*(w + gap) + w, margin + h + gap + h)
        src = fitz.open(paths[i+2])
        page_out.show_pdf_page(get_best_fit_rect(rect, src[0].rect.width, src[0].rect.height), src, 0)
        src.close()
        
    # Row 3: E (centered, width 1200)
    rect_e = fitz.Rect(400, margin + 2*(h + gap), 1600, margin + 2*(h + gap) + h)
    src_e = fitz.open(paths[4])
    page_out.show_pdf_page(get_best_fit_rect(rect_e, src_e[0].rect.width, src_e[0].rect.height), src_e, 0)
    src_e.close()
    
    doc_out.save(dest_pdf)
    
    # Save high-res PNG (300 DPI)
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png)
    doc_out.close()
    
    # Copy consolidated files to Manuscript Figures directory
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Figures"
    copy_file_safe(dest_pdf, os.path.join(latex_dir, "Figure4_Consolidated.pdf"))
    copy_file_safe(dest_png, os.path.join(latex_dir, "Figure4_Consolidated.png"))
    
    print(f"  [SUCCESS] Finalized Main Figure 4 PDF & PNG in {dest_dir} and copied to LaTeX directory.")

if __name__ == "__main__":
    assemble_main_figure4()
