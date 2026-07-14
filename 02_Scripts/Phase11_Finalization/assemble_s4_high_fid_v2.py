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

def assemble_s4_high_fid_v2():
    heatmap_path = r"d:\Proj_AML\04_Figures\15_Immune_Deconvolution\immune_cell_heatmap.pdf"
    boxplot_path = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts\s4_pB_upgraded.pdf"
    corr_path = r"d:\Proj_AML\05_Submission\Submission_Hub\05_Internal_Drafts\s4_pC.pdf"
    
    # Primary output paths
    dest_pdf_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS4.pdf"
    dest_png_1 = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS4.png"
    
    # Synchronization outputs
    dest_pdf_2 = r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS4.pdf"
    dest_png_2 = r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\03_Supplementary_Figures\FigureS4.png"
    
    # Compilation folder for LaTeX
    dest_pdf_3 = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS4.pdf"
    dest_png_3 = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary\FigureS4.png"
    
    doc_out = fitz.open()
    # 3-panel layout: True Wide Heatmap top, Boxplot/Corr bottom
    # Enlarged canvas height to 2100 to fit enlarged Plot A
    page_out = doc_out.new_page(width=2000, height=2100)
    
    margin = 50
    gap = 50
    # SIGNIFICANTLY ENLARGED Plot A (Heatmap) from 600 to 900 points (+50% increase)
    h_heatmap = 900 
    w_half = (2000 - 2*margin - gap) / 2
    h_bottom = 850
    
    # Panel A Header (R-Generated for 100% Identity)
    header_path = r"d:\Proj_AML\04_Figures\15_Immune_Deconvolution\FigureS4_HeaderA.pdf"
    rect_header = fitz.Rect(margin, margin, 2000-margin, margin + 120)
    src_header = fitz.open(header_path)
    page_out.show_pdf_page(rect_header, src_header, 0)
    src_header.close()

    # Panel A: Heatmap (True Wide & Vertically Expanded for high row text legibility)
    rect_a = fitz.Rect(margin, margin + 120, 2000-margin, margin + 120 + h_heatmap)
    src_a = fitz.open(heatmap_path)
    page_out.show_pdf_page(rect_a, src_a, 0)
    src_a.close()

    # Panel B: Boxplots (aligned below heatmap with 50pt gap)
    rect_b = fitz.Rect(margin, margin + 120 + h_heatmap + 50, margin + w_half, margin + 120 + h_heatmap + 50 + h_bottom)
    src_b = fitz.open(boxplot_path)
    page_out.show_pdf_page(rect_b, src_b, 0)
    src_b.close()
    
    # Panel C: Correlation
    rect_c = fitz.Rect(margin + w_half + gap, margin + 120 + h_heatmap + 50, 2000-margin, margin + 120 + h_heatmap + 50 + h_bottom)
    src_c = fitz.open(corr_path)
    page_out.show_pdf_page(rect_c, src_c, 0)
    src_c.close()
    
    # Ensure destination directories exist
    os.makedirs(os.path.dirname(dest_pdf_1), exist_ok=True)
    os.makedirs(os.path.dirname(dest_pdf_2), exist_ok=True)
    os.makedirs(os.path.dirname(dest_pdf_3), exist_ok=True)
    
    # Save primary PDF output
    doc_out.save(dest_pdf_1)
    
    # Render PNG
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png_1)
    print(f"Successfully finalized enlarged High-Fidelity Figure S4 in Submission_Hub: {dest_pdf_1}")
    
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
    assemble_s4_high_fid_v2()
