import fitz # PyMuPDF
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
    
    return fitz.Rect(
        quad_rect.x0 + dx, 
        quad_rect.y0 + dy, 
        quad_rect.x0 + dx + new_w, 
        quad_rect.y0 + dy + new_h
    )

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

def main():
    print("=== STARTING MAIN FIGURE 5 ASSEMBLY (VRS BIOMARKER) ===")
    print("Assembling high-fidelity 5-panel portrait Main Figure 5...")
    
    doc_out = fitz.open()
    page_out = doc_out.new_page(width=720, height=1150)
    
    panels = [
        {
            "id": "A",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\Figure5A_VRS_Histogram.pdf",
            "quad": fitz.Rect(25, 50, 345, 350),
            "label_pos": fitz.Point(25, 35)
        },
        {
            "id": "B",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\VRS_validation_beataml2_scatter.pdf",
            "quad": fitz.Rect(375, 50, 695, 350),
            "label_pos": fitz.Point(375, 35)
        },
        {
            "id": "C",
            "path": r"d:\Proj_AML\04_Figures\06_Drug_Response\drug_sensitivity_heatmap.pdf",
            "quad": fitz.Rect(25, 420, 400, 780),
            "label_pos": fitz.Point(25, 395)
        },
        {
            "id": "E",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\Figure5E_VIALE_Validation.pdf",
            "quad": fitz.Rect(410, 420, 695, 780),
            "label_pos": fitz.Point(410, 395)
        },
        {
            "id": "D",
            "path": r"d:\Proj_AML\04_Figures\29_ExternalValidation\Figure5D_Salvage_Boxplots.pdf",
            "quad": fitz.Rect(25, 830, 695, 1130),
            "label_pos": fitz.Point(25, 805)
        }
    ]
    
    for panel in panels:
        path = panel["path"]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Source panel for {panel['id']} not found: {path}")
            
        print(f"  Placing Panel {panel['id']} ({os.path.basename(path)})...")
        doc_in = fitz.open(path)
        actual_w = doc_in[0].rect.width
        actual_h = doc_in[0].rect.height
        fit_rect = get_best_fit_rect(panel["quad"], actual_w, actual_h)
        page_out.show_pdf_page(fit_rect, doc_in, 0)
        doc_in.close()
        
    print("  Drawing bold labels for Panel C and Panel E...")
    for panel in panels:
        if panel["id"] == "C":
            page_out.insert_text(
                panel["label_pos"], 
                "C. Global Ex Vivo Resistance Profile", 
                fontsize=9.6, 
                fontname="hebo",
                color=(0, 0, 0.54)
            )
            subtitle_pos = fitz.Point(panel["label_pos"].x, panel["label_pos"].y + 11)
            page_out.insert_text(
                subtitle_pos, 
                "Targeted drug sensitivity across patient cohort", 
                fontsize=8.4, 
                fontname="helv", 
                color=(0, 0, 0.54)
            )
        elif panel["id"] == "E":
            page_out.insert_text(
                panel["label_pos"], 
                "E. Clinical Trial Validation (VIALE-A)", 
                fontsize=9.6, 
                fontname="hebo",
                color=(0, 0, 0.54)
            )
            subtitle_pos = fitz.Point(panel["label_pos"].x, panel["label_pos"].y + 11)
            page_out.insert_text(
                subtitle_pos, 
                "VRS tracks the 24-gene resistance signature", 
                fontsize=8.4, 
                fontname="helv", 
                color=(0, 0, 0.54)
            )
        
    pdf_out = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures\Figure5_Consolidated.pdf"
    png_out = r"d:\Proj_AML\05_Submission\Submission_Hub\02_Main_Figures\Figure5_Consolidated.png"
    
    os.makedirs(os.path.dirname(pdf_out), exist_ok=True)
    doc_out.save(pdf_out)
    doc_out.close()
    
    doc_png = fitz.open(pdf_out)
    page_png = doc_png[0]
    pix = page_png.get_pixmap(dpi=300)
    pix.save(png_out)
    doc_png.close()
    
    latex_dir = r"d:\Proj_AML\05_Submission\Submission_Hub\01_Manuscript_Source\Figures"
    copy_file_safe(pdf_out, os.path.join(latex_dir, "Figure5_Consolidated.pdf"))
    copy_file_safe(png_out, os.path.join(latex_dir, "Figure5_Consolidated.png"))
    
    print("=== MAIN FIGURE 5 ASSEMBLY COMPLETED SUCCESSFULLY ===")

if __name__ == "__main__":
    main()
