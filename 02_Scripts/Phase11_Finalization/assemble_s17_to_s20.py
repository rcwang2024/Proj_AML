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
            # Same physical file (junction/symlink)
            return True
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception as e:
        # Fallback in case of permissions or locks
        print(f"  [WARNING] Copy failed from {src} to {dst}: {e}")
        return False

def assemble_s17():
    print("Assembling Figure S17 (ELN Benchmarking)...")
    path_in = r"d:\Proj_AML\04_Figures\Phase11_Finalization"
    dest_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS17.pdf"
    dest_png = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS17.png"
    
    panel_a = os.path.join(path_in, "VRS_distribution_by_ELN.pdf")
    panel_b = os.path.join(path_in, "VRS_vs_AUC_stratified_by_ELN.pdf")
    
    if not os.path.exists(panel_a) or not os.path.exists(panel_b):
        print("  [ERROR] S17 component files not found.")
        return
        
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=1300, height=550)
    
    margin = 50
    # Place Panel A
    rect_a = fitz.Rect(margin, margin, margin + 450, 500)
    src_a = fitz.open(panel_a)
    page_out.show_pdf_page(get_best_fit_rect(rect_a, src_a[0].rect.width, src_a[0].rect.height), src_a, 0)
    src_a.close()
    
    # Place Panel B
    rect_b = fitz.Rect(margin + 450 + 50, margin, 1300 - margin, 500)
    src_b = fitz.open(panel_b)
    page_out.show_pdf_page(get_best_fit_rect(rect_b, src_b[0].rect.width, src_b[0].rect.height), src_b, 0)
    src_b.close()
    
    os.makedirs(os.path.dirname(dest_pdf), exist_ok=True)
    doc_out.save(dest_pdf)
    
    # Save high-res PNG
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png)
    doc_out.close()
    
    # Copy consolidated files to LaTeX directory
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary"
    copy_file_safe(dest_pdf, os.path.join(latex_dir, "FigureS17.pdf"))
    copy_file_safe(dest_png, os.path.join(latex_dir, "FigureS17.png"))
    
    print(f"  [SUCCESS] Finalized Figure S17 PDF & PNG and copied consolidated files.")

def assemble_s18():
    print("Assembling Figure S18 (Single-Cell Validation)...")
    path_in = r"d:\Proj_AML\04_Figures\Phase11_Finalization"
    dest_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS18.pdf"
    dest_png = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS18.png"
    
    panels = [
        "FigureS18_sc_UMAP_celltypes.pdf",
        "FigureS18_sc_UMAP_VRS.pdf",
        "FigureS18_sc_VRS_boxplot.pdf",
        "FigureS18_sc_dotplot.pdf",
        "FigureS18_sc_BCL2_MCL1_tradeoff.pdf"
    ]
    
    paths = [os.path.join(path_in, p) for p in panels]
    if any(not os.path.exists(p) for p in paths):
        print("  [ERROR] S18 component files not found.")
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
    
    os.makedirs(os.path.dirname(dest_pdf), exist_ok=True)
    doc_out.save(dest_pdf)
    
    # Save high-res PNG
    pix = page_out.get_pixmap(matrix=fitz.Matrix(300/72, 300/72))
    pix.save(dest_png)
    doc_out.close()
    
    # Copy consolidated files to LaTeX directory
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary"
    copy_file_safe(dest_pdf, os.path.join(latex_dir, "FigureS18.pdf"))
    copy_file_safe(dest_png, os.path.join(latex_dir, "FigureS18.png"))
            
    print(f"  [SUCCESS] Finalized Figure S18 PDF & PNG and copied consolidated files.")

def assemble_s19():
    print("Assembling Figure S19 (Noise Sensitivity)...")
    path_in = r"d:\Proj_AML\04_Figures\Phase11_Finalization"
    dest_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS19.pdf"
    dest_png = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS19.png"
    
    src_pdf = os.path.join(path_in, "VRS_noise_sensitivity.pdf")
    src_png = os.path.join(path_in, "VRS_noise_sensitivity.png")
    
    if os.path.exists(src_pdf):
        copy_file_safe(src_pdf, dest_pdf)
        if os.path.exists(src_png):
            copy_file_safe(src_png, dest_png)
            
        # Copy consolidated files to LaTeX directory
        latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary"
        copy_file_safe(dest_pdf, os.path.join(latex_dir, "FigureS19.pdf"))
        copy_file_safe(dest_png, os.path.join(latex_dir, "FigureS19.png"))
        print("  [SUCCESS] Finalized Figure S19 PDF & PNG and copied consolidated files.")
    else:
        print("  [ERROR] S19 source file not found.")

def assemble_s20():
    print("Assembling Figure S20 (VRS Web Calculator)...")
    src_png = r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\01_Manuscript_Source\Supplementary\FigureS20.png"
    dest_png = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS20.png"
    dest_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub\03_Supplementary_Figures\FigureS20.pdf"
    
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Supplementary"
    
    if os.path.exists(src_png):
        # Copy PNG files first
        copy_file_safe(src_png, dest_png)
        copy_file_safe(src_png, os.path.join(latex_dir, "FigureS20.png"))
        
        # Create PDF wrapper
        try:
            doc = fitz.open()
            img = fitz.open(src_png)
            rect = img[0].rect
            page = doc.new_page(width=rect.width, height=rect.height)
            page.insert_image(rect, filename=src_png)
            doc.save(dest_pdf)
            doc.close()
            img.close()
            
            # Copy PDF wrapper to LaTeX directory
            copy_file_safe(dest_pdf, os.path.join(latex_dir, "FigureS20.pdf"))
            print("  [SUCCESS] Finalized Figure S20 PDF & PNG and copied panels.")
        except Exception as e:
            print(f"  [WARNING] PyMuPDF PDF wrapper creation failed: {e}")
            # If pdf wrapper fails, just copy from fallback if exists
            fallback_pdf = r"d:\Proj_AML\05_Submission\Submission_Hub_2026-05-16\01_Manuscript_Source\Supplementary\FigureS20.pdf"
            copy_file_safe(fallback_pdf, dest_pdf)
            copy_file_safe(fallback_pdf, os.path.join(latex_dir, "FigureS20.pdf"))
            print("  [SUCCESS] Copied fallback Figure S20 PDF.")
    else:
        print("  [ERROR] S20 source file not found.")

def main():
    print("=== STARTING ASSEMBLY OF FIGURES S17 TO S20 ===")
    assemble_s17()
    assemble_s18()
    assemble_s19()
    assemble_s20()
    print("=== FIGURES S17 TO S20 ASSEMBLY COMPLETE ===")

if __name__ == "__main__":
    main()
